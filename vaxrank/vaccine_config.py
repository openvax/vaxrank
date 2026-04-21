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
Configuration for vaccine peptide assembly.

This module defines the VaccineConfig class which controls how vaccine
peptides are assembled from epitope predictions. Vaccine peptides are
longer sequences (typically 25 amino acids) that contain one or more
predicted MHC-binding epitopes spanning the mutated residue(s).

The vaccine peptide selection considers:
1. Peptide length requirements for synthesis and delivery
2. Position of the mutation within the peptide
3. Number of strong MHC-binding epitopes contained
4. Total predicted immunogenicity score
"""

from typing import Optional

import msgspec

from .config.defaults import (
    DEFAULT_MANUFACTURABILITY_MAX_C_TERMINAL_HYDROPATHY,
    DEFAULT_MANUFACTURABILITY_MAX_KMER_HYDROPATHY_HIGH_PRIORITY,
    DEFAULT_MANUFACTURABILITY_MAX_KMER_HYDROPATHY_LOW_PRIORITY,
    DEFAULT_MANUFACTURABILITY_MIN_KMER_HYDROPATHY,
    DEFAULT_MAX_PEPTIDE_LENGTH,
    DEFAULT_MAX_VACCINE_PEPTIDES_PER_VARIANT,
    DEFAULT_MIN_PEPTIDE_LENGTH,
    DEFAULT_NUM_MUTANT_EPITOPES_TO_KEEP,
    DEFAULT_PADDING_AROUND_MUTATION,
    DEFAULT_PREFERRED_PEPTIDE_LENGTH,
    DEFAULT_VACCINE_PEPTIDE_SCORE_FRACTION_OF_BEST,
)
from .manufacturability import MANUFACTURABILITY_RULE_REGISTRY
from .ranking import RANKING_RULE_REGISTRY


COMBINED_SCORE_MODES = (
    "sqrt_reads_times_epitope",
    "reads_times_epitope",
    "epitope_only",
)
DEFAULT_COMBINED_SCORE_MODE = "sqrt_reads_times_epitope"


class VaccineConfig(msgspec.Struct, frozen=True, forbid_unknown_fields=True):
    """
    Configuration parameters for vaccine peptide assembly.

    This immutable struct contains all parameters needed to assemble
    epitope predictions into vaccine peptides suitable for synthesis
    and therapeutic use.

    Attributes
    ----------
    preferred_peptide_length : int
        Preferred length of vaccine peptides in amino acids.
        Default: 25

    min_peptide_length : int
        Minimum acceptable vaccine peptide length in amino acids.
        Default: 25

    max_peptide_length : int
        Maximum acceptable vaccine peptide length in amino acids.
        Default: 25

    padding_around_mutation : int
        Number of off-center window positions to consider when selecting
        vaccine peptides. A value of 5 means considering windows where
        the mutation is anywhere from position 5 to position 20 in a
        25-mer peptide.
        Default: 5

    max_vaccine_peptides_per_variant : int
        Maximum number of vaccine peptides to select for each variant.
        Multiple peptides may be selected if they have sufficiently
        different epitope content.
        Default: 1

    num_mutant_epitopes_to_keep : int
        Maximum number of epitope predictions to retain per variant
        during processing. Higher values increase computational cost
        but may capture more potential vaccine targets. Set to 0 to
        keep all predicted mutant epitopes.
        Default: 1000

    score_fraction_of_best : float
        Keep vaccine peptide candidates whose score is at least this fraction
        of the best candidate for the variant before lexicographic tie-breaking.
        Default: 0.99

    max_c_terminal_hydropathy : float
        Maximum acceptable mean hydropathy (GRAVY) score for the 7
        C-terminal residues.  Peptides exceeding this are penalised
        during manufacturability ranking.
        Default: 1.5

    min_kmer_hydropathy : float
        Minimum acceptable max-7mer GRAVY score.  Peptides with all
        windows below this floor are penalised (too hydrophilic).
        Default: 0.0

    max_kmer_hydropathy_low_priority : float
        Low-priority upper bound on any 7-mer GRAVY window.  Applied
        as a tie-breaker after higher-priority constraints.
        Default: 1.5

    max_kmer_hydropathy_high_priority : float
        High-priority upper bound on any 7-mer GRAVY window.  Exceeding
        this indicates a severely hydrophobic region that will be
        difficult to synthesise.
        Default: 2.5

    Examples
    --------
    >>> config = VaccineConfig()
    >>> config.preferred_peptide_length
    25

    >>> long_peptide_config = VaccineConfig(preferred_peptide_length=30)
    >>> long_peptide_config.preferred_peptide_length
    30

    >>> multi_peptide_config = VaccineConfig(max_vaccine_peptides_per_variant=3)
    >>> multi_peptide_config.max_vaccine_peptides_per_variant
    3
    """

    preferred_peptide_length: int = DEFAULT_PREFERRED_PEPTIDE_LENGTH
    min_peptide_length: int = DEFAULT_MIN_PEPTIDE_LENGTH
    max_peptide_length: int = DEFAULT_MAX_PEPTIDE_LENGTH
    padding_around_mutation: int = DEFAULT_PADDING_AROUND_MUTATION
    max_vaccine_peptides_per_variant: int = DEFAULT_MAX_VACCINE_PEPTIDES_PER_VARIANT
    num_mutant_epitopes_to_keep: int = DEFAULT_NUM_MUTANT_EPITOPES_TO_KEEP
    score_fraction_of_best: float = DEFAULT_VACCINE_PEPTIDE_SCORE_FRACTION_OF_BEST
    max_c_terminal_hydropathy: float = DEFAULT_MANUFACTURABILITY_MAX_C_TERMINAL_HYDROPATHY
    min_kmer_hydropathy: float = DEFAULT_MANUFACTURABILITY_MIN_KMER_HYDROPATHY
    max_kmer_hydropathy_low_priority: float = DEFAULT_MANUFACTURABILITY_MAX_KMER_HYDROPATHY_LOW_PRIORITY
    max_kmer_hydropathy_high_priority: float = DEFAULT_MANUFACTURABILITY_MAX_KMER_HYDROPATHY_HIGH_PRIORITY
    manufacturability_rules: Optional[tuple[str, ...]] = None
    combined_score_mode: str = DEFAULT_COMBINED_SCORE_MODE
    ranking_rules: Optional[tuple[str, ...]] = None

    def __post_init__(self):
        if self.preferred_peptide_length < 1:
            raise ValueError(
                f"preferred_peptide_length must be at least 1, "
                f"got {self.preferred_peptide_length}"
            )
        if self.min_peptide_length < 1:
            raise ValueError(
                f"min_peptide_length must be at least 1, "
                f"got {self.min_peptide_length}"
            )
        if self.max_peptide_length < 1:
            raise ValueError(
                f"max_peptide_length must be at least 1, "
                f"got {self.max_peptide_length}"
            )
        if self.min_peptide_length > self.max_peptide_length:
            raise ValueError(
                f"min_peptide_length ({self.min_peptide_length}) must be "
                f"<= max_peptide_length ({self.max_peptide_length})"
            )
        if not (self.min_peptide_length
                <= self.preferred_peptide_length
                <= self.max_peptide_length):
            raise ValueError(
                f"preferred_peptide_length ({self.preferred_peptide_length}) must be "
                f"between min_peptide_length ({self.min_peptide_length}) and "
                f"max_peptide_length ({self.max_peptide_length})"
            )
        if self.padding_around_mutation < 0:
            raise ValueError(
                f"padding_around_mutation must be non-negative, "
                f"got {self.padding_around_mutation}"
            )
        if self.max_vaccine_peptides_per_variant < 0:
            raise ValueError(
                f"max_vaccine_peptides_per_variant must be non-negative, "
                f"got {self.max_vaccine_peptides_per_variant}"
            )
        if self.num_mutant_epitopes_to_keep < 0:
            raise ValueError(
                f"num_mutant_epitopes_to_keep must be non-negative, "
                f"got {self.num_mutant_epitopes_to_keep}"
            )
        if not (0 < self.score_fraction_of_best <= 1):
            raise ValueError(
                f"score_fraction_of_best must be in (0, 1], "
                f"got {self.score_fraction_of_best}"
            )
        if self.manufacturability_rules is not None:
            for rule_name in self.manufacturability_rules:
                if rule_name not in MANUFACTURABILITY_RULE_REGISTRY:
                    raise ValueError(
                        f"Unknown manufacturability rule '{rule_name}'. "
                        f"Available: {sorted(MANUFACTURABILITY_RULE_REGISTRY)}"
                    )
        if self.combined_score_mode not in COMBINED_SCORE_MODES:
            raise ValueError(
                f"combined_score_mode must be one of {COMBINED_SCORE_MODES}, "
                f"got '{self.combined_score_mode}'"
            )
        if self.ranking_rules is not None:
            for rule_name in self.ranking_rules:
                if rule_name not in RANKING_RULE_REGISTRY:
                    raise ValueError(
                        f"Unknown ranking rule '{rule_name}'. "
                        f"Available: {sorted(RANKING_RULE_REGISTRY)}"
                    )
