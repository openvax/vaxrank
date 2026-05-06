# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Per-peptide context records used by ``VaccinePeptide``.

Two layers:

  ``PeptideContext`` — one peptide sequence + flanks + multi-axis
  binding predictions. Generic shape used for the antigenic
  candidate AND for any reference comparator (WT pair, nearest_self,
  nearest_vital_self, nearest_nonCTA, nearest_oncovirus, …).

  ``CandidateEpitope`` — one sliding-window position from a
  VaccinePeptide. Holds the mutant ``PeptideContext`` plus a
  ``comparators`` dict keyed by name. Future safety / homology
  features (#254 / #257 / #258) populate comparators with their
  respective contexts; today only ``'wt'`` is canonical.

## Why this shape

Replaces the legacy flat ``EpitopePrediction`` (single-axis
affinity-only record). Predictions live in ``mhctools.Prediction``
form so multi-axis data (pMHC_affinity / pMHC_presentation /
pMHC_stability / immunogenicity / proteasome_cleavage /
antigen_processing) flows through unchanged. Multi-predictor
support (post-#261) is native to mhctools' record — vaxrank just
holds them.

## Design choice: cross-predictor disambiguation is explicit

mhctools' ``Prediction.score`` is normalized so higher = better
*within one predictor*. Cross-predictor ``max(score)`` is
meaningless — mhcflurry's score scale and netmhcpan's score scale
are independent transforms. So ``best_for_kind(kind)`` raises
when there's no unambiguous answer and the caller hasn't passed
``predictor=``. Callers that legitimately want a per-predictor
loop iterate ``predictors_for_kind(kind)`` and call
``best_for_kind(kind, predictor=…)`` per predictor.

Issue: openvax/vaxrank#282 (replaces).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


# Reserved comparator names. Open-ended — callers can add new ones —
# but these are the canonical entries vaxrank knows about today or
# will populate as the related issues land. Documented here so
# downstream readers know what to expect.
COMPARATOR_WT = 'wt'                                # same position, ref allele
COMPARATOR_NEAREST_SELF = 'nearest_self'            # closest self-peptide (#254)
COMPARATOR_NEAREST_VITAL_SELF = 'nearest_vital_self'  # closest in vital tissues (#254)
COMPARATOR_NEAREST_NONCTA = 'nearest_nonCTA'        # closest non-CTA self (#257)
COMPARATOR_NEAREST_ONCOVIRUS = 'nearest_oncovirus'  # closest oncoviral peptide (#258)


@dataclass(frozen=True)
class PeptideContext:
    """One peptide sequence + flanking residues + (optional) MHC
    binding predictions.

    Generic — the same shape is used for the antigenic candidate
    (the mutant) and for every reference comparator (WT pair,
    nearest_self, …). When a comparator carries predictions, it's
    answering "does the patient's MHC also bind this comparator?"
    — which is the actual safety question for autoimmunity scoring.

    Predictions follow mhctools' multi-kind / multi-predictor /
    multi-allele shape: one ``mhctools.Prediction`` per
    (kind, predictor, allele) tuple. Stored as a flat tuple for
    serializability; structured access is via the methods below.
    """

    peptide_sequence: str
    n_flank: str = ''
    c_flank: str = ''
    # Source provenance — for the mutant: the assembled mutant
    # protein fragment. For ``nearest_self``: the matching self
    # protein. For ``nearest_oncovirus``: the viral genome / ORF.
    source_sequence: str = ''
    source_name: str = ''         # gene / virus / "self_proteome" / …
    offset: int = 0
    # tuple[mhctools.Prediction, ...] — left untyped so this module
    # stays import-light (mhctools pulls in heavy deps).
    predictions: tuple = ()

    # ------------------------------------------------------------------
    # Structured views — read-only, computed on demand.
    # ------------------------------------------------------------------

    def predictions_by_kind_and_predictor(self) -> dict:
        """Three-level nested view:
        ``{kind: {predictor_name: {allele: Prediction}}}``.

        Built fresh each call (the storage tuple is the source of
        truth). Use this when you need full structure; the
        ``best_*`` methods cover the common case.
        """
        nested: dict = {}
        for p in self.predictions:
            (nested.setdefault(p.kind, {})
                   .setdefault(p.predictor_name, {})[p.allele]) = p
        return nested

    def kinds(self) -> tuple:
        """Sorted tuple of distinct ``kind`` values in this context.
        Empty when the context has no predictions (e.g. a
        ``nearest_self`` comparator that's just a sequence with no
        binding probe yet)."""
        return tuple(sorted({p.kind for p in self.predictions}))

    def predictors_for_kind(self, kind: str) -> tuple:
        """Predictors that emitted ``kind`` for this peptide.
        Drives multi-predictor disambiguation in
        ``best_for_kind``."""
        return tuple(sorted({
            p.predictor_name for p in self.predictions
            if p.kind == kind}))

    def alleles_for(self, kind: str, predictor: str) -> tuple:
        """Alleles a given (kind, predictor) covers. Used by
        coverage / report code that walks per-allele evidence
        without forcing a single 'best' allele."""
        return tuple(sorted({
            p.allele for p in self.predictions
            if p.kind == kind and p.predictor_name == predictor
            and p.allele}))

    # ------------------------------------------------------------------
    # Best-of accessors — kind-named, predictor-aware.
    # ------------------------------------------------------------------

    def best_for_kind(self, kind: str, *,
                       predictor: Optional[str] = None):
        """Best prediction for ``kind`` across alleles for one
        predictor. Returns the ``mhctools.Prediction`` (carries
        allele, score, value, percentile_rank, …) or ``None`` when
        no record matches.

        Disambiguation:
        - If only one predictor emitted ``kind`` for this peptide,
          ``predictor`` may be omitted.
        - If multiple predictors emitted ``kind``, ``predictor`` is
          required — cross-predictor ``max(score)`` is meaningless
          since each predictor's score scale is its own.

        Raises ``ValueError`` when the call is ambiguous.
        """
        candidates = [p for p in self.predictions if p.kind == kind]
        if not candidates:
            return None
        predictor_set = {p.predictor_name for p in candidates}
        if predictor is None:
            if len(predictor_set) > 1:
                raise ValueError(
                    "%s has predictions from multiple predictors "
                    "(%s); pass predictor= to disambiguate. "
                    "Each predictor's score scale is its own — "
                    "cross-predictor max() is meaningless."
                    % (kind, sorted(predictor_set)))
            # Single predictor → unambiguous.
        else:
            candidates = [
                p for p in candidates if p.predictor_name == predictor]
            if not candidates:
                return None
        # mhctools normalizes ``score`` so higher = better within a
        # predictor regardless of whether the underlying axis is
        # IC50 (lower = better) or a probability (higher = better).
        # Best-across-alleles is unambiguously max(score).
        return max(candidates, key=lambda p: p.score)

    # Kind-specific helpers. Each forwards to ``best_for_kind``
    # with the same disambiguation contract.
    def best_affinity(self, *, predictor: Optional[str] = None):
        return self.best_for_kind('pMHC_affinity', predictor=predictor)

    def best_presentation(self, *, predictor: Optional[str] = None):
        return self.best_for_kind('pMHC_presentation', predictor=predictor)

    def best_stability(self, *, predictor: Optional[str] = None):
        return self.best_for_kind('pMHC_stability', predictor=predictor)

    def best_immunogenicity(self, *, predictor: Optional[str] = None):
        return self.best_for_kind('immunogenicity', predictor=predictor)

    def best_cleavage(self, *, predictor: Optional[str] = None):
        return self.best_for_kind('proteasome_cleavage', predictor=predictor)

    def best_antigen_processing(self, *, predictor: Optional[str] = None):
        return self.best_for_kind('antigen_processing', predictor=predictor)


@dataclass(frozen=True)
class CandidateEpitope:
    """One sliding-window peptide position from a VaccinePeptide,
    with its mutation context + reference comparators for safety
    / immunogenicity scoring.

    Comparator names are open-ended; vaxrank populates these as
    their respective issues land:

      ``'wt'``                  same position, reference allele
                                (drives mutation-specificity +
                                today's wt_ic50 column)
      ``'nearest_self'``        closest match in the patient's
                                reference proteome (#254 —
                                autoimmunity risk)
      ``'nearest_vital_self'``  closest match in vital-tissue
                                proteins (#254)
      ``'nearest_nonCTA'``      closest self-peptide that isn't
                                a cancer-testis antigen (#257)
      ``'nearest_oncovirus'``   closest match in oncovirus
                                reference genomes (#258)

    Each value is a ``PeptideContext`` carrying its own peptide
    sequence + (optionally) its own MHC predictions. When a
    comparator carries predictions, scorers can ask "does the
    patient's MHC also bind this comparator?" — the actual safety
    question for autoimmunity / off-target presentation.
    """

    mutant: PeptideContext
    # ``field(default_factory=dict)`` would normally be required
    # for mutable defaults, but the class is frozen so callers
    # construct dicts explicitly. Keeps the dataclass simple.
    comparators: dict = field(default_factory=dict)

    # Mutation-specific context. Lives at this level (not on
    # ``PeptideContext``) because it's about the mutant-vs-source
    # relationship, not a property of the peptide itself.
    overlaps_mutation: bool = False
    occurs_in_reference: bool = False

    # Convenience: most call sites today read the WT alongside the
    # mutant. Common-case shortcut.
    @property
    def wt(self) -> Optional[PeptideContext]:
        return self.comparators.get(COMPARATOR_WT)

    def comparator(self, name: str) -> Optional[PeptideContext]:
        """Generic accessor for any comparator. Returns ``None``
        when the comparator isn't present for this candidate
        (typical when a safety scorer hasn't run yet)."""
        return self.comparators.get(name)

    # Convenience pass-throughs to the mutant context — most
    # callers read these on the candidate level.
    @property
    def peptide_sequence(self) -> str:
        return self.mutant.peptide_sequence

    @property
    def length(self) -> int:
        return len(self.mutant.peptide_sequence)
