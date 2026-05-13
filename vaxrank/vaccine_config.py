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

from .combined_score_dsl import parse_combined_score_expr
from .config.defaults import (
    DEFAULT_MAX_PEPTIDE_LENGTH,
    DEFAULT_MAX_VACCINE_PEPTIDES_PER_VARIANT,
    DEFAULT_MIN_PEPTIDE_LENGTH,
    DEFAULT_NUM_MUTANT_EPITOPES_TO_KEEP,
    DEFAULT_PADDING_AROUND_MUTATION,
    DEFAULT_PREFERRED_PEPTIDE_LENGTH,
    DEFAULT_VACCINE_PEPTIDE_SCORE_FRACTION_OF_BEST,
)
from .ranking import RANKING_RULE_REGISTRY


# Single source of truth for the default vaccine-peptide combined score.
# Everything that computes ``combined_score`` reads this constant (or a
# user-supplied override) via the DSL — no parallel hardcoded branch.
DEFAULT_COMBINED_SCORE_EXPR = "sqrt(n_alt_reads) * target_epitope_score"


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

    num_target_epitopes_to_keep : int
        Maximum number of epitope predictions to retain per variant
        during processing. Higher values increase computational cost
        but may capture more potential vaccine targets. Set to 0 to
        keep all predicted mutant epitopes.
        Default: 1000

    score_fraction_of_best : float
        Keep vaccine peptide candidates whose score is at least this fraction
        of the best candidate for the variant before lexicographic tie-breaking.
        Default: 0.99

    Manufacturability config (synthesis-difficulty thresholds + rules)
    moved off VaccineConfig in 2.20 — it's peptide-only, so it lives
    on :class:`vaxrank.manufacturability_config.ManufacturabilityConfig`
    and is threaded into the ranker as a separate parameter alongside
    VaccineConfig. See :func:`vaxrank.core_logic.run_vaxrank` and
    :class:`vaxrank.peptide.PeptideConstructConfig`.

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
    num_target_epitopes_to_keep: int = DEFAULT_NUM_MUTANT_EPITOPES_TO_KEEP
    score_fraction_of_best: float = DEFAULT_VACCINE_PEPTIDE_SCORE_FRACTION_OF_BEST
    # DSL expression evaluated to produce each vaccine peptide's
    # ``combined_score``. See :mod:`vaxrank.combined_score_dsl` for
    # grammar + bindings (``target_epitope_score``, ``n_alt_reads``,
    # ``expression_score``, etc.). The default is
    # :data:`DEFAULT_COMBINED_SCORE_EXPR` — the legacy
    # ``sqrt(n_alt_reads) * target_epitope_score`` formula. There is
    # no separate mode enum: the expression IS the mechanism.
    combined_score_expr: str = DEFAULT_COMBINED_SCORE_EXPR
    ranking_rules: Optional[tuple[str, ...]] = None
    # When True (default, legacy behavior), variants whose vaccine peptides
    # contain no mutant-overlapping epitopes are dropped from the report
    # entirely. Set False to keep every variant for which any vaccine
    # peptide could be assembled, regardless of epitope coverage —
    # useful when the report itself is downstream input for an
    # independent epitope predictor.
    require_target_epitopes_in_variant: bool = True

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
        if self.num_target_epitopes_to_keep < 0:
            raise ValueError(
                f"num_target_epitopes_to_keep must be non-negative, "
                f"got {self.num_target_epitopes_to_keep}"
            )
        if not (0 < self.score_fraction_of_best <= 1):
            raise ValueError(
                f"score_fraction_of_best must be in (0, 1], "
                f"got {self.score_fraction_of_best}"
            )
        # Validate the expression at config construction so a malformed
        # YAML / CLI override fails fast, not on first vaccine peptide.
        # ``parse_combined_score_expr`` performs the same AST-allowlist
        # check used at evaluation time.
        parse_combined_score_expr(self.combined_score_expr)
        if self.ranking_rules is not None:
            for rule_name in self.ranking_rules:
                if rule_name not in RANKING_RULE_REGISTRY:
                    raise ValueError(
                        f"Unknown ranking rule '{rule_name}'. "
                        f"Available: {sorted(RANKING_RULE_REGISTRY)}"
                    )
