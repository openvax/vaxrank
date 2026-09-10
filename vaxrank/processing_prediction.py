# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""``ProcessingPrediction`` — a per-(peptide, source_sequence, offset)
proteasomal-cleavage prediction.

Lives on its own axis from MHC binding because the two prediction
kinds are semantically different:

- :class:`mhctools.Prediction` (carried inside ``vaxrank.CandidateEpitope``)
  records **(peptide, allele) MHC-binding** scores (output of an
  ``mhctools.BindingPredictor`` — pMHC affinity / presentation /
  stability).
- ``ProcessingPrediction`` is a **(peptide, source_sequence, offset)
  proteasomal-cleavage** score (output of an
  ``mhctools.ProcessingPredictor`` — no allele axis, depends on
  the peptide's flanking context within its source protein).

Pre-2.22 vaxrank annotated flat record objects in place by
adding ``pepsickle_*`` attributes — that conflated the two
prediction kinds. ``ProcessingPrediction`` (this module) is the
post-2.22 canonical record; consumers join in by
``(peptide_sequence, source_sequence, peptide_offset, predictor_name)`` at
render time.

Issue: openvax/vaxrank#272.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Real
from typing import Optional


@dataclass(frozen=True)
class ProcessingPrediction:
    """One ``ProcessingPredictor`` score for a (peptide,
    source_sequence, peptide_offset) occurrence.

    The composite ``processing_score`` is the geometric mean of
    ``c_term_cleavage_prob`` and ``(1 - max_internal_cut_prob)`` —
    bounded [0, 1] but not a calibrated probability of release,
    destruction or presentation. It is unavailable at sequence endpoints.
    Geometric mean rather than
    raw product so a balanced ``(0.6, 0.6)`` row scores ~0.6
    instead of 0.36.

    Attributes
    ----------
    peptide_sequence : str
        The MHC-ligand peptide whose proteasomal cleavage was scored.
    source_sequence : str
        The protein context the peptide was scored within.
    peptide_offset : int
        Resolved zero-based start of this occurrence in source_sequence.
        Repeated peptides retain separate local processing predictions.
    predictor_name : str
        Lowercase predictor identifier (e.g. ``'pepsickle'``). Future
        per-position cleavage predictors (NetChop, PAProC, …) plug in
        with their own name; report writers join across predictors by
        ``(peptide_sequence, source_sequence, peptide_offset, predictor_name)``.
    predictor_version : Optional[str]
        Predictor version string when the predictor exposes one,
        else ``None``.
    c_term_cleavage_prob : Optional[float]
        Probability the proteasome cuts at the ligand's C-terminus
        boundary. Range [0, 1]; None at a sequence endpoint. An endpoint
        may be a cropped window, not an exposed molecular terminus.
    max_internal_cut_prob : float
        Peak cleavage probability strictly inside the ligand
        (model evidence, not proof of ligand destruction). Range [0, 1].
    processing_score : Optional[float]
        Composite ``sqrt(c_term_cleavage_prob *
        (1 - max_internal_cut_prob))``.
    """

    peptide_sequence: str
    source_sequence: str
    predictor_name: str
    peptide_offset: int
    predictor_version: Optional[str] = None
    c_term_cleavage_prob: Optional[float] = None
    max_internal_cut_prob: float = 0.0
    processing_score: Optional[float] = None
    c_boundary_status: str = ""

    def __post_init__(self):
        end = self.peptide_offset + len(self.peptide_sequence)
        if (type(self.peptide_offset) is not int or self.peptide_offset < 0
                or not self.peptide_sequence
                or self.source_sequence[self.peptide_offset:end] != self.peptide_sequence):
            raise ValueError("Processing prediction must match its exact source occurrence")
        status = ("sequence_endpoint" if end == len(self.source_sequence) else
                  ("observed" if self.c_term_cleavage_prob is not None else "unassessed"))
        if self.c_boundary_status and self.c_boundary_status != status:
            raise ValueError("Processing boundary status disagrees with its source occurrence")
        object.__setattr__(self, "c_boundary_status", status)
        if status == "sequence_endpoint" and self.c_term_cleavage_prob is not None:
            raise ValueError("Sequence endpoint cannot claim a C-boundary cleavage probability")
        if self.c_term_cleavage_prob is None and self.processing_score is not None:
            raise ValueError("Missing C-boundary evidence cannot produce a composite score")
        for value in (self.c_term_cleavage_prob, self.max_internal_cut_prob, self.processing_score):
            if value is not None and (isinstance(value, bool) or not isinstance(value, Real)
                                      or not math.isfinite(value) or not 0 <= value <= 1):
                raise ValueError("Processing scores must be finite probabilities in [0, 1] or None")

    def key(self) -> tuple:
        """Stable join key used by report writers to look up the
        ProcessingPrediction for a given target-ligand occurrence at render
        time. Includes the position and predictor name so a future second
        per-position cleavage predictor (NetChop, …) lands in the
        same map without colliding."""
        return (
            self.peptide_sequence, self.source_sequence,
            self.peptide_offset, self.predictor_name)
