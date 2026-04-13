# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Configuration for epitope scoring and filtering.

This module defines the EpitopeConfig class which controls how epitopes
are scored based on MHC binding affinity and filtered for vaccine peptide
selection.

The default epitope scoring uses a logistic function to transform IC50
binding affinity values into a normalized score between 0 and 1:

    rescaled = (ic50 - midpoint) / width
    score = (1 / (1 + exp(rescaled))) / normalizer

Where normalizer ensures score ≈ 1.0 when ic50 ≈ 0:

    normalizer = 1 / (1 + exp(-midpoint / width))

Parameters:
- midpoint: IC50 value (nM) around which scores transition (default: 350 nM)
- width: Controls steepness of the scoring curve (default: 150)

Lower IC50 values indicate stronger binding and result in higher scores.

Users can override filtering and scoring with topiary DSL expressions via
``filter_expr`` and ``score_expr``. When either is omitted the scalar
fields above define the default via :mod:`vaxrank.epitope_dsl`.
"""

from typing import Optional

import msgspec

from .config.defaults import (
    DEFAULT_BINDING_AFFINITY_CUTOFF,
    DEFAULT_LOGISTIC_EPITOPE_SCORE_MIDPOINT,
    DEFAULT_LOGISTIC_EPITOPE_SCORE_WIDTH,
    DEFAULT_MIN_EPITOPE_SCORE,
    DEFAULT_PERCENTILE_RANK_CUTOFF,
)


class EpitopeConfig(msgspec.Struct, frozen=True, forbid_unknown_fields=True):
    """
    Configuration parameters for epitope scoring and filtering.

    Attributes
    ----------
    logistic_epitope_score_midpoint : float
        IC50 value (in nM) at which the epitope score equals 0.5.
        Default: 350.0 nM

    logistic_epitope_score_width : float
        Width parameter controlling the steepness of the logistic scoring
        function.
        Default: 150.0

    min_epitope_score : float
        Minimum normalized score threshold. Epitopes below this are excluded.
        Default: 0.00001

    binding_affinity_cutoff : float
        Maximum IC50 value (in nM) to consider.
        Default: 5000.0 nM

    scoring_mode : str
        How to score epitopes. Options:
        - "affinity": Use IC50 binding affinity with logistic scoring (default)
        - "percentile_rank": Use percentile rank directly (lower = better,
          inverted to produce higher scores for better binders)

    percentile_rank_cutoff : float
        The percentile rank at or above which an epitope is treated as a
        non-binder in percentile-rank scoring mode.

    filter_expr : str, optional
        Topiary DSL string used to drop epitope rows wholesale, e.g.
        ``"affinity <= 500 & affinity.rank <= 2.0"``. When omitted the
        filter is synthesized from ``binding_affinity_cutoff`` /
        ``percentile_rank_cutoff`` so default behavior is unchanged.

    score_expr : str, optional
        Topiary DSL string producing a per-peptide-allele score, e.g.
        ``"affinity.logistic(350, 150)"``. When omitted the score node
        is synthesized from the logistic / percentile-rank scalar fields
        (preserving the ``1/(1+exp(-mid/width))`` normalizer that keeps
        default scores byte-identical with pre-5.0 topiary behavior).
    """

    logistic_epitope_score_midpoint: float = DEFAULT_LOGISTIC_EPITOPE_SCORE_MIDPOINT
    logistic_epitope_score_width: float = DEFAULT_LOGISTIC_EPITOPE_SCORE_WIDTH
    min_epitope_score: float = DEFAULT_MIN_EPITOPE_SCORE
    binding_affinity_cutoff: float = DEFAULT_BINDING_AFFINITY_CUTOFF
    scoring_mode: str = "affinity"
    percentile_rank_cutoff: float = DEFAULT_PERCENTILE_RANK_CUTOFF
    filter_expr: Optional[str] = None
    score_expr: Optional[str] = None

    def __post_init__(self):
        if self.logistic_epitope_score_midpoint <= 0:
            raise ValueError(
                f"logistic_epitope_score_midpoint must be positive, "
                f"got {self.logistic_epitope_score_midpoint}"
            )
        if self.logistic_epitope_score_width <= 0:
            raise ValueError(
                f"logistic_epitope_score_width must be positive, "
                f"got {self.logistic_epitope_score_width}"
            )
        if self.min_epitope_score < 0:
            raise ValueError(
                f"min_epitope_score must be non-negative, "
                f"got {self.min_epitope_score}"
            )
        if self.binding_affinity_cutoff <= 0:
            raise ValueError(
                f"binding_affinity_cutoff must be positive, "
                f"got {self.binding_affinity_cutoff}"
            )
        if self.scoring_mode not in ("affinity", "percentile_rank"):
            raise ValueError(
                f"scoring_mode must be 'affinity' or 'percentile_rank', "
                f"got '{self.scoring_mode}'"
            )
        if self.percentile_rank_cutoff <= 0 or self.percentile_rank_cutoff > 100:
            raise ValueError(
                f"percentile_rank_cutoff must be in (0, 100], "
                f"got {self.percentile_rank_cutoff}"
            )
        if self.filter_expr is not None or self.score_expr is not None:
            # Parse eagerly so malformed DSL strings fail at config load
            # rather than mid-pipeline when predict_epitopes runs.
            from topiary.ranking import parse as _topiary_parse
            for name, expr in (
                ("filter_expr", self.filter_expr),
                ("score_expr", self.score_expr),
            ):
                if expr is None:
                    continue
                try:
                    _topiary_parse(expr)
                except Exception as exc:
                    raise ValueError(
                        f"Invalid {name} {expr!r}: {exc}"
                    ) from exc
