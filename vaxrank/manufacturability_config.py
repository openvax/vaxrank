# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Peptide-synthesis difficulty rules + thresholds, as a standalone
config struct.

Pre-2.20 these knobs lived as flat fields on
:class:`vaxrank.vaccine_config.VaccineConfig`. Post-2.20 they live
on their own struct, attached to
:class:`vaxrank.peptide.PeptideConstructConfig` (peptide-only —
solid-phase synthesis rules don't apply to in-vivo-translated mRNA
constructs). The ranker pulls them in by accepting a
``manufacturability_config`` parameter alongside ``vaccine_config``;
when ``--vaccine-type`` doesn't include peptide, the parameter is
``None`` and the manufacturability sentinel inside
``ranking_rules`` is a no-op tie-break.
"""

from __future__ import annotations

from typing import Optional

import msgspec

from .config.defaults import (
    DEFAULT_MANUFACTURABILITY_MAX_C_TERMINAL_HYDROPATHY,
    DEFAULT_MANUFACTURABILITY_MAX_KMER_HYDROPATHY_HIGH_PRIORITY,
    DEFAULT_MANUFACTURABILITY_MAX_KMER_HYDROPATHY_LOW_PRIORITY,
    DEFAULT_MANUFACTURABILITY_MIN_KMER_HYDROPATHY,
)
from .manufacturability import MANUFACTURABILITY_RULE_REGISTRY


class ManufacturabilityConfig(
    msgspec.Struct, frozen=True, forbid_unknown_fields=True
):
    """Configuration for peptide-synthesis difficulty scoring.

    Attributes
    ----------
    max_c_terminal_hydropathy : float
        Maximum acceptable mean GRAVY score for the 7 C-terminal
        residues. Peptides exceeding this are penalised during
        manufacturability ranking.
    min_kmer_hydropathy : float
        Minimum acceptable max-7mer GRAVY score. Peptides with all
        windows below this floor are penalised (too hydrophilic).
    max_kmer_hydropathy_low_priority : float
        Low-priority upper bound on any 7-mer GRAVY window; applied
        as a tie-breaker after higher-priority constraints.
    max_kmer_hydropathy_high_priority : float
        High-priority upper bound on any 7-mer GRAVY window;
        exceeding this indicates a severely hydrophobic region
        that's difficult to synthesise.
    rules : Optional[tuple[str, ...]]
        Ordered list of synthesis-difficulty rules applied when the
        ``manufacturability`` sentinel appears in
        ``vaccine_peptides.ranking_rules``. ``None`` means use the
        registry's default rule set
        (:data:`vaxrank.manufacturability.DEFAULT_MANUFACTURABILITY_RULES`).
    """

    max_c_terminal_hydropathy: float = (
        DEFAULT_MANUFACTURABILITY_MAX_C_TERMINAL_HYDROPATHY)
    min_kmer_hydropathy: float = DEFAULT_MANUFACTURABILITY_MIN_KMER_HYDROPATHY
    max_kmer_hydropathy_low_priority: float = (
        DEFAULT_MANUFACTURABILITY_MAX_KMER_HYDROPATHY_LOW_PRIORITY)
    max_kmer_hydropathy_high_priority: float = (
        DEFAULT_MANUFACTURABILITY_MAX_KMER_HYDROPATHY_HIGH_PRIORITY)
    rules: Optional[tuple[str, ...]] = None

    def __post_init__(self):
        if self.rules is not None:
            for rule_name in self.rules:
                if rule_name not in MANUFACTURABILITY_RULE_REGISTRY:
                    raise ValueError(
                        f"Unknown manufacturability rule "
                        f"'{rule_name}'. Available: "
                        f"{sorted(MANUFACTURABILITY_RULE_REGISTRY)}")

    def thresholds_dict(self) -> dict[str, float]:
        """Convenience: the four threshold knobs as a dict, in the
        shape :class:`vaxrank.vaccine_peptide.VaccinePeptide`
        expects (``manufacturability_thresholds`` kwarg)."""
        return {
            "max_c_terminal_hydropathy": self.max_c_terminal_hydropathy,
            "min_kmer_hydropathy": self.min_kmer_hydropathy,
            "max_kmer_hydropathy_low_priority":
                self.max_kmer_hydropathy_low_priority,
            "max_kmer_hydropathy_high_priority":
                self.max_kmer_hydropathy_high_priority,
        }
