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

import pandas as pd


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


def _parse(expr):
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
    resolution — see :func:`kinds_referenced_unqualified`.
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


def kinds_referenced_unqualified(cfg):
    """Set of Kinds the DSL references without bracket selection.

    Each such Kind needs a default method when the data exposes
    multiple models for it, otherwise topiary raises "Ambiguous" at
    eval time. The default score / filter node (when the user hasn't
    set ``score_expr`` / ``filter_expr``) also implicitly references
    ``pMHC_affinity`` unqualified.
    """
    out = set()
    for expr in (cfg.filter_expr, cfg.score_expr):
        if expr is None:
            continue
        refs = collect_dsl_references(_parse(expr))
        for kind, method, _version in refs["kinds"]:
            if method is None:
                out.add(kind)
    # The synthesized default node always touches pMHC_affinity unqualified
    # (``Affinity.logistic_normalized(...)`` in affinity mode), so even when
    # the user provides no formulas we still need a default for that Kind.
    if cfg.filter_expr is None and cfg.score_expr is None:
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


def resolve_default_methods(cfg, topiary_df):
    """Return the effective ``default_methods`` dict for each Kind with
    multiple models present in ``topiary_df``.

    Merges ``cfg.default_methods`` (explicit user choice) with an
    auto-picked canonical for Kinds the user didn't specify. Kinds with
    only one method in the data are omitted (no disambiguation needed).
    Only covers Kinds actually present in the DataFrame.
    """
    import logging
    logger = logging.getLogger(__name__)

    if topiary_df.empty:
        return dict(cfg.default_methods or {})

    user_defaults = dict(cfg.default_methods or {})
    resolved = {}
    for kind, group in topiary_df.groupby("kind"):
        methods = set(group["prediction_method_name"].dropna().unique())
        if len(methods) <= 1:
            continue
        if kind in user_defaults:
            if user_defaults[kind] not in methods:
                raise ValueError(
                    f"default_methods[{kind!r}]={user_defaults[kind]!r} "
                    f"but that method is not present in the loaded "
                    f"predictions (available: {sorted(methods)})")
            resolved[kind] = user_defaults[kind]
        else:
            picked = _auto_pick_canonical_method(methods)
            logger.info(
                "Kind %s has multiple methods %s and no default_methods "
                "entry; auto-picking %r as the default for unqualified "
                "references. Set default_methods[%r] in your epitope "
                "config to override.",
                kind, sorted(methods), picked, kind)
            resolved[kind] = picked
    return resolved


def subset_topiary_df_for_eval(topiary_df, cfg):
    """Drop rows that would cause "Ambiguous Kind" errors at DSL eval.

    For each Kind with multiple models in the data, keep only rows
    whose method is (a) explicitly referenced by bracket in any DSL
    expression, or (b) the resolved default method (user-specified or
    auto-picked). Kinds with a single model pass through unchanged.

    This lets users write unqualified ``Affinity.logistic(...)`` and
    have it resolve to the canonical affinity predictor while
    bracketed references like ``Affinity['netmhcpan']`` still pick up
    their specific rows.
    """
    if topiary_df.empty:
        return topiary_df

    qualified = qualified_methods_by_kind(cfg)
    unqualified_kinds = kinds_referenced_unqualified(cfg)
    resolved_defaults = resolve_default_methods(cfg, topiary_df)

    chunks = []
    for kind, group in topiary_df.groupby("kind", sort=False):
        methods = set(group["prediction_method_name"].dropna().unique())
        if len(methods) <= 1:
            chunks.append(group)
            continue

        allowed = set(qualified.get(kind, set())) & methods
        # Unqualified refs (including the default node's implicit
        # ``Affinity`` when no formulas are set) need a default method.
        if kind in unqualified_kinds or not allowed:
            default = resolved_defaults.get(kind)
            if default is not None:
                allowed.add(default)
        chunks.append(group[group["prediction_method_name"].isin(allowed)])

    import pandas as pd
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
