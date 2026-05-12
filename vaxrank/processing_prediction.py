# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""``ProcessingPrediction`` — a per-(peptide, source_sequence)
proteasomal-cleavage prediction.

Lives on its own axis from MHC binding because the two prediction
kinds are semantically different:

- :class:`mhctools.Prediction` (carried inside ``vaxrank.CandidateEpitope``)
  records **(peptide, allele) MHC-binding** scores (output of an
  ``mhctools.BindingPredictor`` — pMHC affinity / presentation /
  stability).
- ``ProcessingPrediction`` is a **(peptide, source_sequence)
  proteasomal-cleavage** score (output of an
  ``mhctools.ProcessingPredictor`` — no allele axis, depends on
  the peptide's flanking context within its source protein).

Pre-2.22 vaxrank annotated EpitopePrediction objects in place by
adding ``pepsickle_*`` attributes — that conflated the two
prediction kinds. ``ProcessingPrediction`` (this module) is the
post-2.22 canonical record; consumers join in by
``(peptide_sequence, source_sequence, predictor_name)`` at
render time.

Issue: openvax/vaxrank#272.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class ProcessingPrediction:
    """One ``ProcessingPredictor`` score for a (peptide,
    source_sequence) pair.

    The composite ``processing_score`` is the geometric mean of
    ``c_term_cleavage_prob`` and ``(1 - max_internal_cut_prob)`` —
    bounded [0, 1], 1 = ideal release, 0 = no clean release OR
    near-certain internal destruction. Geometric mean rather than
    raw product so a balanced ``(0.6, 0.6)`` row scores ~0.6
    instead of 0.36.

    Attributes
    ----------
    peptide_sequence : str
        The MHC-ligand peptide whose proteasomal cleavage was scored.
    source_sequence : str
        The protein context the peptide was scored within.
    predictor_name : str
        Lowercase predictor identifier (e.g. ``'pepsickle'``). Future
        per-position cleavage predictors (NetChop, PAProC, …) plug in
        with their own name; report writers join across predictors by
        ``(peptide_sequence, source_sequence, predictor_name)``.
    predictor_version : Optional[str]
        Predictor version string when the predictor exposes one,
        else ``None``.
    c_term_cleavage_prob : float
        Probability the proteasome cuts at the ligand's C-terminus
        (clean release). Range [0, 1].
    max_internal_cut_prob : float
        Peak cleavage probability strictly inside the ligand
        (high → ligand is destroyed before reaching MHC). Range
        [0, 1].
    processing_score : float
        Composite ``sqrt(c_term_cleavage_prob *
        (1 - max_internal_cut_prob))``.
    """

    peptide_sequence: str
    source_sequence: str
    predictor_name: str
    predictor_version: Optional[str] = None
    c_term_cleavage_prob: float = 0.0
    max_internal_cut_prob: float = 0.0
    processing_score: float = 0.0

    def key(self) -> tuple:
        """Stable join key used by report writers to look up the
        ProcessingPrediction for a given mutant CandidateEpitope at render
        time. Includes the predictor name so a future second
        per-position cleavage predictor (NetChop, …) lands in the
        same map without colliding."""
        return (
            self.peptide_sequence, self.source_sequence,
            self.predictor_name)
