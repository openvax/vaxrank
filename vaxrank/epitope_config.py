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

The epitope scoring uses a logistic function to transform IC50 binding
affinity values into a normalized score between 0 and 1:

    rescaled = (ic50 - midpoint) / width
    score = (1 / (1 + exp(rescaled))) / normalizer

Where normalizer ensures score ≈ 1.0 when ic50 ≈ 0:

    normalizer = 1 / (1 + exp(-midpoint / width))

Parameters:
- midpoint: IC50 value (nM) around which scores transition (default: 350 nM)
- width: Controls steepness of the scoring curve (default: 150)

Lower IC50 values indicate stronger binding and result in higher scores.
"""

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
    """

    logistic_epitope_score_midpoint: float = DEFAULT_LOGISTIC_EPITOPE_SCORE_MIDPOINT
    logistic_epitope_score_width: float = DEFAULT_LOGISTIC_EPITOPE_SCORE_WIDTH
    min_epitope_score: float = DEFAULT_MIN_EPITOPE_SCORE
    binding_affinity_cutoff: float = DEFAULT_BINDING_AFFINITY_CUTOFF
    scoring_mode: str = "affinity"
    percentile_rank_cutoff: float = DEFAULT_PERCENTILE_RANK_CUTOFF
