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


def score_predictions(predictions, cfg):
    """Score external-input predictions using the configured Topiary DSL.

    Returns a ``pandas.Series`` of per-(peptide, allele) scores indexed by
    the group-key MultiIndex topiary uses internally
    (``source_sequence_name, peptide, peptide_offset, allele``). Callers
    merge this back onto their report DataFrame.

    When ``cfg.score_expr`` / ``cfg.filter_expr`` are unset this applies
    the same default node that the main pipeline uses, which matches the
    legacy :meth:`EpitopePrediction.logistic_epitope_score` byte-for-byte.
    """
    from topiary.ranking import EvalContext, apply_filter

    df = predictions_to_topiary_df(predictions)
    if df.empty:
        return pd.Series(dtype=float)

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


def validate_dsl_against_predictions(cfg, predictions):
    """Error early when ``filter_expr`` / ``score_expr`` reference a
    predictor the loaded predictions don't expose.

    Topiary's own ``apply_filter`` already validates ``Column(...)``
    references, but when a Field's (kind, method) isn't present in the
    data it silently evaluates to NaN — so users get "every score is
    zero" rather than a clear error. This function catches that case
    up front. The parsed node's ``Field.kind`` is already the
    canonical string (``"pMHC_affinity"`` etc.) so we match it
    directly against the values in the topiary DataFrame.
    """
    nodes = []
    if cfg.filter_expr is not None:
        nodes.append(_parse(cfg.filter_expr))
    if cfg.score_expr is not None:
        nodes.append(_parse(cfg.score_expr))
    if not nodes:
        return

    df = predictions_to_topiary_df(predictions)
    available_kinds = set(df["kind"].unique()) if not df.empty else set()
    available_methods = (
        set(df["prediction_method_name"].dropna().unique())
        if not df.empty else set())
    available_versions = (
        set(df["predictor_version"].dropna().unique())
        if not df.empty else set())

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
                    f"versions {sorted(v for v in available_versions if v)}. "
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
