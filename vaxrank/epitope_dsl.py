# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Topiary DSL node construction for epitope filtering and scoring.

Two entry points:

- :func:`build_filter_node` — returns a boolean ``DSLNode`` (or ``None``) for
  the configured ``EpitopeConfig``; used with ``topiary.ranking.apply_filter``
  to drop rows wholesale before scoring.

- :func:`build_score_node` — returns a numeric ``DSLNode`` whose ``.eval(ctx)``
  produces a ``pd.Series`` indexed by ``(source_sequence_name, peptide,
  peptide_offset, allele)`` group tuples.

When ``filter_expr`` / ``score_expr`` are set on the config, each is parsed
via :func:`topiary.ranking.parse`. When absent the nodes are built directly
from scalar fields (``binding_affinity_cutoff``, ``logistic_epitope_score_*``,
``scoring_mode``, ``percentile_rank_cutoff``). The default affinity-mode
score uses topiary 5.1's :class:`LogisticNormalizedExpr` so the output is in
``[0, 1]`` and matches the legacy
:meth:`vaxrank.epitope_prediction.EpitopePrediction.logistic_epitope_score`
byte-for-byte.

Multi-model inputs (e.g. LENS with both MHCflurry and netMHCpan producing
``pMHC_affinity``) are disambiguated via topiary 5.10's
``EvalContext(default_methods=...)`` and
``apply_filter(default_methods=...)``. Vaxrank computes the resolved
``default_methods`` dict by merging ``cfg.default_methods`` with an
auto-picked canonical (``mhcflurry`` > ``netmhcpan`` > alphabetical) for
Kinds the user hasn't specified, and forwards it to topiary.
"""

from __future__ import annotations

import functools
import logging

import pandas as pd


logger = logging.getLogger(__name__)


# Method name -> topiary Kind for external-input predictions (LENS / pVACseq
# import). Listed methods come from ``load_lens`` / ``load_pvacseq``;
# anything not in this map is assumed to be a pMHC_affinity predictor.
_METHOD_KIND_MAP = {
    "mhcflurry": "pMHC_affinity",
    "netmhcpan": "pMHC_affinity",
    "netmhcstabpan": "pMHC_stability",
    # pVACseq aggregates across predictors without naming one; treat it as
    # a generic affinity predictor so default DSL nodes light up.
    "pvacseq": "pMHC_affinity",
}


@functools.lru_cache(maxsize=64)
def _parse(expr):
    """Parse a DSL string into a topiary node. Cached so the same config's
    ``filter_expr`` / ``score_expr`` aren't re-parsed multiple times within
    one :func:`write_neoepitope_report` call (validator + build_filter_node
    + build_score_node all need the AST). Parsed nodes are immutable, so
    sharing is safe.
    """
    from topiary.ranking import parse
    return parse(expr)


def _kind_for_method(method_name):
    """Topiary Kind string for a prediction_method_name.

    Unknown methods default to pMHC_affinity so that a plain
    ``affinity.logistic(...)`` expression keeps working on external
    inputs even when the method name isn't in our registry.
    """
    if not method_name:
        return "pMHC_affinity"
    return _METHOD_KIND_MAP.get(str(method_name).lower(), "pMHC_affinity")


def predictions_to_topiary_df(predictions):
    """Convert a list of :class:`EpitopePrediction` into the topiary
    long-format DataFrame consumed by :class:`topiary.ranking.EvalContext`
    and :func:`topiary.ranking.apply_filter`.

    Each ``EpitopePrediction`` becomes one row; when multiple predictions
    share a (peptide, allele) pair (e.g. MHCflurry and netMHCpan rows
    from a LENS file), DSL expressions can select between them via
    ``affinity['mhcflurry']`` / ``affinity['netmhcpan']``, and when the
    predictions carry ``predictor_version`` they can also disambiguate
    by version via ``affinity['mhcflurry', '2.1.1']``. Unqualified refs
    (``affinity.value``) resolve to the Kind's default method — see
    :func:`resolve_default_methods`.
    """
    rows = []
    for p in predictions:
        tool = str(p.prediction_method_name or "")
        rows.append({
            "sample_name": "",
            "source_sequence_name": p.source_sequence or p.peptide_sequence,
            "peptide": p.peptide_sequence,
            "peptide_offset": p.offset,
            "peptide_length": len(p.peptide_sequence),
            "allele": p.allele,
            "n_flank": "",
            "c_flank": "",
            "prediction_method_name": tool,
            "predictor_version": p.predictor_version or "",
            "kind": _kind_for_method(tool),
            "value": p.ic50,
            "affinity": p.ic50,
            "percentile_rank": p.percentile_rank,
            "score": 0.0,
        })
    return pd.DataFrame(rows)


def epitopes_to_topiary_df(epitopes):
    """Convert a list of :class:`vaxrank.peptide_context.Epitope` into the
    topiary long-format DataFrame.

    Replacement for :func:`predictions_to_topiary_df` as part of the
    EpitopePrediction → Epitope migration (vaxrank#284). One row per
    leaf ``mhctools.Prediction`` in each Epitope's mutant context —
    the new types carry multi-axis predictions natively, so a single
    Epitope can emit N rows (one per allele × predictor × kind).

    Schema matches ``predictions_to_topiary_df`` so consumers
    (``apply_filter``, ``EvalContext``, ``score_predictions``) don't
    care which producer built the frame. Notable differences from
    the legacy producer:

    * ``score`` column carries the real ``mhctools.Prediction.score``
      value instead of the hardcoded ``0.0``. The legacy frame's
      ``score`` was never read by the default DSL (which references
      ``Affinity.value``), so populating it can only help.
    * ``n_flank`` / ``c_flank`` come from the parent ``PeptideContext``
      instead of being empty strings.
    * Comparator predictions (``epitope.wt``, ``nearest_self``, …) are
      NOT emitted as rows. The DSL frame is mutant-only by convention
      — comparator data stays on the Epitope for display. This matches
      legacy behavior (``predictions_to_topiary_df`` never emitted
      ``wt_*`` columns either).
    """
    rows = []
    for e in epitopes:
        ctx = e.mutant
        source_name = ctx.source_sequence or ctx.peptide_sequence
        for p in ctx.predictions_flat():
            rows.append({
                "sample_name": "",
                "source_sequence_name": source_name,
                "peptide": p.peptide,
                "peptide_offset": ctx.offset,
                "peptide_length": len(p.peptide),
                "allele": p.allele,
                "n_flank": ctx.n_flank,
                "c_flank": ctx.c_flank,
                "prediction_method_name": p.predictor_name,
                "predictor_version": p.predictor_version or "",
                "kind": p.kind,
                "value": p.value,
                "affinity": p.value,
                "percentile_rank": p.percentile_rank,
                "score": p.score,
            })
    return pd.DataFrame(rows)


def _legacy_predictions_to_epitopes(predictions):
    """Adapter: convert ``list[EpitopePrediction]`` to ``list[Epitope]``.

    Transitional bridge for the Phase 2 migration (vaxrank#284). Used
    by producers that still emit ``EpitopePrediction`` internally —
    they pipe through this adapter to expose the new shape to migrated
    consumers. Drops out of the codebase once every producer emits
    ``Epitope`` natively.

    Groups input by ``(peptide_sequence, source_sequence, offset)`` —
    each group becomes one ``Epitope``. Per-allele / per-predictor
    rows in the group populate the mutant ``PeptideContext.predictions``;
    when ``wt_peptide_sequence`` differs from the mutant peptide and
    ``wt_ic50`` is set, a parallel ``wt`` comparator context is built.

    ``mhctools.Prediction.score`` is set to ``0.0`` because
    ``EpitopePrediction`` doesn't carry a normalized score — matches
    the legacy frame's hardcoded value. Downstream code that calls
    ``best()`` (which ranks by score) will see ties and pick first;
    consumers needing a meaningful score should use ``predictions_for``
    + custom ranking. Real scores arrive when producers stop going
    through this adapter.
    """
    from mhctools.pred import Prediction
    from vaxrank.peptide_context import COMPARATOR_WT, Epitope, PeptideContext

    groups: dict = {}
    for p in predictions:
        key = (p.peptide_sequence, p.source_sequence, p.offset)
        groups.setdefault(key, []).append(p)

    epitopes = []
    for (peptide, source, offset), members in groups.items():
        first = members[0]
        mutant_preds = []
        # WT comparator is built whenever ANY group member carries a
        # non-None wt_ic50. The WT peptide sequence may be empty
        # (LENS doesn't emit one; some pVACseq exports omit the
        # ``WT Epitope Seq`` column) — we keep the empty string as
        # the codebase's "unknown WT peptide" sentinel rather than
        # dropping the WT IC50 signal entirely.
        wt_peptide = first.wt_peptide_sequence
        wt_preds = []
        has_wt_signal = any(ep.wt_ic50 is not None for ep in members)
        for ep in members:
            tool = str(ep.prediction_method_name or "")
            kind = _kind_for_method(tool)
            mutant_preds.append(Prediction(
                kind=kind,
                predictor_name=tool,
                predictor_version=ep.predictor_version or "",
                allele=ep.allele,
                peptide=ep.peptide_sequence,
                value=ep.ic50,
                score=0.0,
                percentile_rank=ep.percentile_rank,
            ))
            if has_wt_signal and ep.wt_ic50 is not None:
                wt_preds.append(Prediction(
                    kind=kind,
                    predictor_name=tool,
                    predictor_version=ep.predictor_version or "",
                    allele=ep.allele,
                    peptide=wt_peptide,
                    value=ep.wt_ic50,
                    score=0.0,
                    percentile_rank=None,
                ))

        mutant_ctx = PeptideContext(
            peptide_sequence=peptide,
            source_sequence=source,
            offset=offset,
            predictions=tuple(mutant_preds),
        )
        comparators = {}
        if wt_preds:
            comparators[COMPARATOR_WT] = PeptideContext(
                peptide_sequence=wt_peptide,
                source_sequence=source,
                offset=offset,
                predictions=tuple(wt_preds),
            )
        epitopes.append(Epitope(
            mutant=mutant_ctx,
            comparators=comparators,
            overlaps_mutation=first.overlaps_mutation,
            occurs_in_reference=first.occurs_in_reference,
        ))
    return epitopes


# Priority order for auto-picking a canonical method when the user hasn't
# set ``default_methods[kind]`` and the data exposes multiple models for
# that Kind. Listed by decreasing preference; any method not listed falls
# back to alphabetical.
_CANONICAL_METHOD_PREFERENCE = (
    "mhcflurry",
    "netmhcpan",
    "netmhcstabpan",
    "netmhcpan_el",
    "netmhcpan_ba",
)


def _auto_pick_canonical_method(methods):
    """Pick one method name from a non-empty set by priority order."""
    for pref in _CANONICAL_METHOD_PREFERENCE:
        if pref in methods:
            return pref
    return sorted(methods)[0]


def resolve_default_methods(cfg, topiary_df):
    """Return the ``default_methods`` dict to pass to topiary's
    :class:`~topiary.ranking.EvalContext` / :func:`apply_filter`.

    For every Kind with multiple models in ``topiary_df``:

    - if ``cfg.default_methods[kind]`` is set, use it;
    - otherwise auto-pick the canonical method (``mhcflurry`` >
      ``netmhcpan`` > ``netmhcstabpan`` > alphabetical) and log
      INFO so the user can make it explicit.

    Kinds with a single model are omitted (topiary doesn't need a
    default for them). Assumes ``validate_default_methods`` has
    already caught any typos in ``cfg.default_methods``.
    """
    if topiary_df.empty:
        return {}
    user_defaults = dict(cfg.default_methods or {})
    resolved = {}
    for kind, group in topiary_df.groupby("kind", sort=False):
        methods = set(group["prediction_method_name"].dropna().unique())
        if len(methods) <= 1:
            continue
        if kind in user_defaults:
            resolved[kind] = user_defaults[kind]
        else:
            picked = _auto_pick_canonical_method(methods)
            logger.info(
                "Kind %s has multiple methods %s and no default_methods "
                "entry; auto-picking %r. Set default_methods[%r] in the "
                "epitope config to override.",
                kind, sorted(methods), picked, kind)
            resolved[kind] = picked
    return resolved


def validate_default_methods(cfg, topiary_df):
    """Error if ``cfg.default_methods`` names a method that isn't in the data.

    Runs even for Kinds with a single model (where the default isn't
    actually consulted by topiary) so typos are caught eagerly.
    """
    if not cfg.default_methods:
        return
    if topiary_df.empty:
        return
    methods_by_kind = (
        topiary_df.groupby("kind")["prediction_method_name"]
        .apply(lambda s: set(s.dropna().unique()))
        .to_dict())
    for kind, method in cfg.default_methods.items():
        available = methods_by_kind.get(kind, set())
        if method not in available:
            raise ValueError(
                f"default_methods[{kind!r}]={method!r} but that method "
                f"is not present in the loaded predictions "
                f"(available for {kind}: {sorted(available)})")


def score_predictions(epitopes, cfg, *, topiary_df=None):
    """Score external-input epitopes using the configured Topiary DSL.

    Returns a ``pandas.Series`` of per-(peptide, allele) scores indexed
    by the group-key MultiIndex topiary uses internally
    (``source_sequence_name, peptide, peptide_offset, allele``). Callers
    merge this back onto their report DataFrame.

    Multi-model data (e.g. both MHCflurry and netMHCpan for
    ``pMHC_affinity``) is handled via topiary's
    ``EvalContext(default_methods=...)`` and
    ``apply_filter(default_methods=...)``: unqualified DSL references
    resolve to the per-Kind default (user-specified via
    ``cfg.default_methods`` or auto-picked canonical), while bracketed
    references like ``Affinity['netmhcpan']`` still pick up their
    specific rows.

    When ``cfg.score_expr`` / ``cfg.filter_expr`` are unset this applies
    the same default node that the main pipeline uses, which matches the
    legacy per-prediction logistic score byte-for-byte.

    Pass ``topiary_df`` to reuse an already-built frame (see
    :func:`epitopes_to_topiary_df`) rather than rebuilding from
    ``epitopes``.
    """
    from topiary.ranking import EvalContext, apply_filter

    df = (epitopes_to_topiary_df(epitopes)
          if topiary_df is None else topiary_df)
    if df.empty:
        return pd.Series(dtype=float)

    validate_default_methods(cfg, df)
    resolved = resolve_default_methods(cfg, df)

    filter_node = build_filter_node(cfg)
    if filter_node is not None:
        df = apply_filter(df, filter_node, default_methods=resolved)
        if df.empty:
            return pd.Series(dtype=float)

    score_node = build_score_node(cfg)
    ctx = EvalContext(df, default_methods=resolved)
    return (
        score_node.eval(ctx)
        .reindex(ctx.group_index)
        .fillna(0.0)
    )


def collect_dsl_references(node):
    """Walk a parsed DSL node and collect its external-data references.

    Returns a dict with two keys:

      ``columns`` - set of column names referenced via ``Column(...)``
      ``kinds``   - set of (kind, method, version) tuples from ``Field`` nodes

    ``method`` and ``version`` may be ``None`` when the formula used a
    kind without qualification (``affinity.value``). Column refs are
    already validated by topiary's ``apply_filter``; vaxrank only uses
    this to drive :func:`validate_dsl_against_predictions`.
    """
    from topiary.ranking.nodes import Column, Field

    columns = set()
    kinds = set()
    stack = [node]
    while stack:
        n = stack.pop()
        if n is None:
            continue
        if isinstance(n, Column):
            columns.add(n.col_name)
        elif isinstance(n, Field):
            kinds.add((n.kind, n.method, n.version))
        child_iter = getattr(n, "child_nodes", None)
        if callable(child_iter):
            stack.extend(child_iter())
    return {"columns": columns, "kinds": kinds}


def validate_dsl_against_predictions(cfg, epitopes, *, topiary_df=None):
    """Error early when ``filter_expr`` / ``score_expr`` references a
    predictor the loaded epitopes don't expose.

    Topiary's own ``apply_filter`` already validates ``Column(...)``
    references, but when a ``Field``'s (kind, method, version) isn't
    present in the data topiary can silently evaluate to NaN for that
    sub-expression. This function catches that up front and points the
    user at ``default_methods`` / ``--input-lens`` when appropriate.

    Pass ``topiary_df`` to reuse an already-built frame rather than
    rebuilding it from ``epitopes``.
    """
    nodes = []
    if cfg.filter_expr is not None:
        nodes.append(_parse(cfg.filter_expr))
    if cfg.score_expr is not None:
        nodes.append(_parse(cfg.score_expr))
    if not nodes:
        return

    df = (epitopes_to_topiary_df(epitopes)
          if topiary_df is None else topiary_df)
    available_kinds = set(df["kind"].unique()) if not df.empty else set()
    available_methods = (
        set(df["prediction_method_name"].dropna().unique())
        if not df.empty else set())
    # Drop the "no version known" sentinel so a formula pinning a specific
    # version against predictions that don't carry any version errors out,
    # and so the error message's version list matches what we actually
    # checked against.
    available_versions = set()
    if not df.empty:
        available_versions = {
            v for v in df["predictor_version"].dropna().unique() if v
        }

    for node in nodes:
        refs = collect_dsl_references(node)
        for kind, method, version in refs["kinds"]:
            if kind not in available_kinds:
                raise ValueError(
                    f"DSL expression references kind {kind!r} but loaded "
                    f"predictions only expose kinds {sorted(available_kinds)}. "
                    f"Check filter_expr / score_expr against the LENS "
                    f"file columns or input source."
                )
            if method is not None and method not in available_methods:
                raise ValueError(
                    f"DSL expression references predictor {method!r} but "
                    f"loaded predictions only expose methods "
                    f"{sorted(available_methods)}. Remove the bracket "
                    f"reference or load a data source that exposes this "
                    f"method."
                )
            if version is not None and version not in available_versions:
                raise ValueError(
                    f"DSL expression references predictor version "
                    f"{version!r} but loaded predictions only expose "
                    f"versions {sorted(available_versions)}. "
                    f"Drop the version from the bracket expression or "
                    f"update the LENS file / loader."
                )


def _default_score_node(cfg):
    from topiary.ranking import Affinity, Column, Const

    if cfg.scoring_mode == "percentile_rank":
        cutoff = cfg.percentile_rank_cutoff
        rank = Column("percentile_rank")
        # Match `max(0, 1 - rank / cutoff)` with the `rank >= cutoff → 0` gate.
        # The mask is strict `<`, so the DSL drops NaN (→ False) and ≥-cutoff
        # rows to 0 exactly as the legacy scorer does.
        linear = Const(1.0) - rank / cutoff
        return (rank < cutoff) * linear.clip(lo=0.0, hi=None)

    midpoint = cfg.logistic_epitope_score_midpoint
    width = cfg.logistic_epitope_score_width
    cutoff_mask = Affinity < cfg.binding_affinity_cutoff
    return cutoff_mask * Affinity.logistic_normalized(midpoint, width)


def build_filter_node(cfg):
    """Return the ``DSLNode`` used to filter predictions, or ``None`` for no-op.

    The default is no filter — legacy vaxrank never hard-dropped rows; it
    scored above-cutoff rows as 0 and then applied ``min_epitope_score`` as
    a gate. The default score node preserves that by multiplying by the
    cutoff mask. Users who want a hard pre-filter must set ``filter_expr``.
    """
    if cfg.filter_expr is not None:
        return _parse(cfg.filter_expr)
    return None


def build_score_node(cfg):
    """Return the ``DSLNode`` that computes the per-epitope score."""
    if cfg.score_expr is not None:
        return _parse(cfg.score_expr)
    return _default_score_node(cfg)
