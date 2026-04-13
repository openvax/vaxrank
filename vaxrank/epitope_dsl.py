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


def _parse(expr):
    from topiary.ranking import parse
    return parse(expr)


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
