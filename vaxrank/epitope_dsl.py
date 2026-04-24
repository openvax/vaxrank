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
    one :func:`write_neoepitope_report` call (validator + subsetter +
    filter-node + score-node all need the AST). Parsed nodes are
    immutable, so sharing is safe.
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
    from a LENS file loaded with ``--lens-predictor=all``), DSL
    expressions can select between them via
    ``affinity['mhcflurry']`` / ``affinity['netmhcpan']``, and when the
    predictions carry ``predictor_version`` they can also disambiguate
    by version via ``affinity['mhcflurry', '2.1.1']``.
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


def score_predictions(predictions, cfg, *, topiary_df=None):
    """Score external-input predictions using the configured Topiary DSL.

    Returns a ``pandas.Series`` of per-(peptide, allele) scores indexed by
    the group-key MultiIndex topiary uses internally
    (``source_sequence_name, peptide, peptide_offset, allele``). Callers
    merge this back onto their report DataFrame.

    When ``cfg.score_expr`` / ``cfg.filter_expr`` are unset this applies
    the same default node that the main pipeline uses, which matches the
    legacy :meth:`EpitopePrediction.logistic_epitope_score` byte-for-byte.

    Pass ``topiary_df`` to reuse an already-built frame (see
    :func:`predictions_to_topiary_df`) rather than rebuilding from
    ``predictions``; callers that validate + score back-to-back should
    build once and share.
    """
    from topiary.ranking import EvalContext, apply_filter

    df = (predictions_to_topiary_df(predictions)
          if topiary_df is None else topiary_df)
    if df.empty:
        return pd.Series(dtype=float)

    # Resolve Kinds with multiple models down to a single default method
    # per Kind (except where the DSL explicitly brackets other methods).
    df = subset_topiary_df_for_eval(df, cfg)

    filter_node = build_filter_node(cfg)
    if filter_node is not None:
        df = apply_filter(df, filter_node)
        if df.empty:
            return pd.Series(dtype=float)

    score_node = build_score_node(cfg)
    ctx = EvalContext(df)
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
    kind without qualification (``affinity.value``). The Column set is
    informational only: topiary's own ``apply_filter`` already validates
    Column references at eval time. The kinds set is what this module's
    :func:`validate_dsl_against_predictions` uses to catch missing
    predictor refs *before* we produce silent NaN scores.
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


def predictors_required_by_cfg(cfg):
    """Set of predictor method names referenced by ``filter_expr`` /
    ``score_expr`` (e.g. ``{"mhcflurry", "netmhcpan"}`` for a formula using
    ``affinity['mhcflurry']`` and ``affinity['netmhcpan']``).

    Returns an empty set when no formulas are set or when every
    ``Field`` reference is unqualified (``affinity.value`` with no
    bracket selection).
    """
    methods = set()
    for kind_to_methods in qualified_methods_by_kind(cfg).values():
        methods |= kind_to_methods
    return methods


def qualified_methods_by_kind(cfg):
    """Methods referenced via bracket selection, grouped by canonical Kind.

    Example: ``{"pMHC_affinity": {"mhcflurry", "netmhcpan"},
    "pMHC_stability": {"netmhcstabpan"}}`` for a formula that pins each.

    Unqualified references (``affinity.value`` with no bracket) don't
    contribute to any set but do imply the Kind needs a default-method
    resolution — see :func:`kinds_needing_method_disambiguation`.
    """
    out = {}
    for expr in (cfg.filter_expr, cfg.score_expr):
        if expr is None:
            continue
        refs = collect_dsl_references(_parse(expr))
        for kind, method, _version in refs["kinds"]:
            if method is not None:
                out.setdefault(kind, set()).add(method)
    return out


def kinds_needing_method_disambiguation(cfg):
    """Kinds whose DSL evaluation depends on seeing a single method per
    ``(source_sequence_name, peptide, peptide_offset, allele)`` group.

    A Kind needs disambiguation when:

    - the DSL references it unqualified (``affinity.value`` with no
      bracket selection), OR
    - ``filter_expr`` / ``score_expr`` are unset, because the synthesized
      default node implicitly references ``pMHC_affinity`` via
      :class:`~topiary.ranking.Affinity` (in affinity mode) or reads the
      ``percentile_rank`` column (in percentile_rank mode) — both of
      which break if the frame has multiple rows per group.

    Note that :func:`subset_topiary_df_for_eval` *also* falls back to
    the default method for any Kind whose bracket-referenced set is
    empty — see the ``not allowed`` branch there. So this function is
    not strictly "every Kind that gets a default"; it's "every Kind
    for which we must force resolution even when bracket refs exist".
    """
    out = set()
    for expr in (cfg.filter_expr, cfg.score_expr):
        if expr is None:
            continue
        refs = collect_dsl_references(_parse(expr))
        for kind, method, _version in refs["kinds"]:
            if method is None:
                out.add(kind)
    if cfg.filter_expr is None and cfg.score_expr is None:
        # The default node implicitly needs single-method rows for
        # pMHC_affinity (affinity mode uses unqualified Affinity;
        # percentile_rank mode reads percentile_rank which differs per
        # method even though it's a Column not a Field).
        out.add("pMHC_affinity")
    return out


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


def _resolve_default_for_kind(kind, methods, user_defaults):
    """Return the default method for a single Kind, user-overridable.

    Raises if the user specified a default that isn't in the data.
    Logs when we auto-pick a canonical.
    """
    if kind in user_defaults:
        choice = user_defaults[kind]
        if choice not in methods:
            raise ValueError(
                f"default_methods[{kind!r}]={choice!r} but that method "
                f"is not present in the loaded predictions "
                f"(available: {sorted(methods)})")
        return choice
    picked = _auto_pick_canonical_method(methods)
    logger.info(
        "Kind %s has multiple methods %s and no default_methods entry; "
        "auto-picking %r. Set default_methods[%r] in the epitope config "
        "to override.", kind, sorted(methods), picked, kind)
    return picked


def validate_default_methods(cfg, topiary_df):
    """Error if ``cfg.default_methods`` names a method that isn't in the data.

    This runs even for Kinds with a single model (where the default
    wouldn't actually be consulted during subsetting) so typos are
    caught rather than silently ignored.
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


def subset_topiary_df_for_eval(topiary_df, cfg):
    """Drop rows that would cause "Ambiguous Kind" errors at DSL eval.

    For each Kind with multiple models in the data, keep only rows
    whose method is (a) explicitly referenced by bracket in any DSL
    expression, or (b) the resolved default method (user-specified via
    ``default_methods[kind]`` or auto-picked). Kinds with a single
    model pass through unchanged.

    Lets users write unqualified ``Affinity.logistic(...)`` and have it
    resolve to the canonical affinity predictor while bracketed
    references like ``Affinity['netmhcpan']`` still pick up their
    specific rows. Raises if ``default_methods[kind]`` names a method
    that isn't in the data.
    """
    if topiary_df.empty:
        return topiary_df

    qualified = qualified_methods_by_kind(cfg)
    needs_disambiguation = kinds_needing_method_disambiguation(cfg)
    user_defaults = dict(cfg.default_methods or {})

    chunks = []
    for kind, group in topiary_df.groupby("kind", sort=False):
        methods = set(group["prediction_method_name"].dropna().unique())
        if len(methods) <= 1:
            chunks.append(group)
            continue

        allowed = set(qualified.get(kind, set())) & methods
        # Add the resolved default when the Kind is referenced unqualified
        # OR when no bracketed refs exist — either way we need a single
        # method to resolve to.
        if kind in needs_disambiguation or not allowed:
            allowed.add(_resolve_default_for_kind(kind, methods, user_defaults))
        chunks.append(group[group["prediction_method_name"].isin(allowed)])

    return pd.concat(chunks).reset_index(drop=True)


def validate_dsl_against_predictions(cfg, predictions, *, topiary_df=None):
    """Error early when ``filter_expr`` / ``score_expr`` reference a
    predictor the loaded predictions don't expose.

    Topiary's own ``apply_filter`` already validates ``Column(...)``
    references, but when a Field's (kind, method) isn't present in the
    data it silently evaluates to NaN — so users get "every score is
    zero" rather than a clear error. This function catches that case
    up front. The parsed node's ``Field.kind`` is already the
    canonical string (``"pMHC_affinity"`` etc.) so we match it
    directly against the values in the topiary DataFrame.

    Pass ``topiary_df`` to reuse an already-built frame rather than
    rebuilding it from ``predictions``.
    """
    nodes = []
    if cfg.filter_expr is not None:
        nodes.append(_parse(cfg.filter_expr))
    if cfg.score_expr is not None:
        nodes.append(_parse(cfg.score_expr))
    if not nodes:
        return

    df = (predictions_to_topiary_df(predictions)
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
                    f"Check filter_expr / score_expr against --lens-predictor "
                    f"and the LENS file columns."
                )
            if method is not None and method not in available_methods:
                raise ValueError(
                    f"DSL expression references predictor {method!r} but "
                    f"loaded predictions only expose methods "
                    f"{sorted(available_methods)}. Either pass "
                    f"--lens-predictor=all to emit every detected predictor "
                    f"or remove the reference from filter_expr / score_expr."
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
