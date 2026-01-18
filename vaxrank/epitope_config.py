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

    score = 1 / (1 + (ic50 / midpoint) ^ (4 * ln(3) / width))

Where:
- midpoint: IC50 value at which score equals 0.5 (default: 350 nM)
- width: Controls steepness of the scoring curve (default: 150)

Lower IC50 values indicate stronger binding and result in higher scores.
"""

import msgspec

# Default minimum epitope score threshold
# Epitopes with scores below this are filtered out
DEFAULT_MIN_EPITOPE_SCORE = 0.00001

# Default maximum IC50 value to consider
# Epitopes with IC50 above this are considered non-binders
DEFAULT_BINDING_AFFINITY_CUTOFF = 5000.0


class EpitopeConfig(msgspec.Struct, frozen=True):
    """
    Configuration parameters for epitope scoring and filtering.

    This immutable struct contains all parameters needed to score epitopes
    based on their MHC binding affinity predictions and filter them for
    vaccine peptide selection.

    Attributes
    ----------
    logistic_epitope_score_midpoint : float
        IC50 value (in nM) at which the epitope score equals 0.5.
        Lower values make the scoring more stringent.
        Default: 350.0 nM

    logistic_epitope_score_width : float
        Width parameter controlling the steepness of the logistic scoring
        function. Smaller values create a sharper transition between
        high and low scores.
        Default: 150.0

    min_epitope_score : float
        Minimum normalized score threshold. Epitopes with scores below
        this value are filtered out and not considered for vaccine
        peptide selection.
        Default: 0.00001

    binding_affinity_cutoff : float
        Maximum IC50 value (in nM) to consider. Epitopes with predicted
        binding affinity above this threshold are considered non-binders
        and excluded.
        Default: 5000.0 nM

    Examples
    --------
    >>> config = EpitopeConfig()
    >>> config.logistic_epitope_score_midpoint
    350.0

    >>> strict_config = EpitopeConfig(min_epitope_score=0.01)
    >>> strict_config.min_epitope_score
    0.01
    """

    logistic_epitope_score_midpoint: float = 350.0
    logistic_epitope_score_width: float = 150.0
    min_epitope_score: float = DEFAULT_MIN_EPITOPE_SCORE
    binding_affinity_cutoff: float = DEFAULT_BINDING_AFFINITY_CUTOFF
