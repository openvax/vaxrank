# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Peptide and CandidateEpitope types used by ``VaccinePeptide``.

  ``Peptide`` — one amino-acid sequence + flanks + multi-axis
  binding predictions. Base shape used for reference comparator
  peptides (WT pair, nearest_self, nearest_vital_self, nearest_nonCTA,
  nearest_oncovirus, …) and as the parent of ``CandidateEpitope``.

  ``CandidateEpitope(Peptide)`` — the vaccine-candidate peptide
  itself plus a ``comparators`` dict of reference peptides + a few
  source-context flags. Inherits Peptide's fields, so accessors
  like ``ce.sequence`` / ``ce.predictions_for(...)`` work directly
  on the candidate (no ``.mutant`` indirection). Source-agnostic:
  works for mutation-derived neoepitopes, oncoviral peptides,
  CTAs, etc. Future safety / homology features (#254 / #257 /
  #258) populate the comparators dict with their respective
  peptides; today only ``'wt'`` is canonical.

## Storage shape

``predictions`` is a *natively nested* dict:
``{kind: {predictor_name: {predictor_version: tuple[Prediction, ...]}}}``.

The leaf tuple holds whatever ``mhctools.Prediction`` records exist
for that ``(kind, predictor, version)`` triple. Allele is a *record
property*, not a structural key — that way:

- ``proteasome_cleavage`` / ``antigen_processing`` (no allele) live
  cleanly as one-element leaves.
- ``pMHC_affinity`` / ``pMHC_stability`` (one record per allele)
  live as N-element leaves.
- ``pMHC_presentation`` lives as N-element leaves under mhctools'
  current per-allele emission (mhcflurry-pres deconvolutes across
  alleles internally; mhctools calls it per-allele to capture each
  allele's score — see ``mhctools/mhcflurry.py`` line 246).

The constructor accepts either nested form *or* a flat
``list``/``tuple`` of ``Prediction`` records for producer
ergonomics — flat input is auto-grouped on construction.
``predictions_flat()`` flattens back out for serialization or
iteration.

## Disambiguation contract

mhctools' ``Prediction.score`` is normalized so higher = better
*within one predictor*. Cross-predictor ``max(score)`` is
meaningless — mhcflurry's score scale and netmhcpan's score scale
are independent transforms. So ``best(kind)`` and
``predictions_for(kind)`` raise when multiple predictors emitted
``kind`` and the caller hasn't passed ``predictor=``.

Versions, however, *do* auto-resolve: when one predictor has
multiple ``predictor_version`` strings on file (e.g. mhcflurry
2.1.0 and 2.1.1 both ran), the leaf accessor defaults to the most
recent version (PEP 440 ordering). Pass ``version=`` to pin.

## Parametric scoring

``best`` uses ``Prediction.score`` (mhctools' normalized
higher-is-better axis). For a different ranking axis, use
``predictions_for(...)`` to get the leaf tuple and rank yourself —
the records expose ``value`` (e.g. IC50) and ``percentile_rank``
alongside ``score``.

Higher-level *peptide* ranking across the report is driven by the
topiary DSL (``EpitopeConfig.score_expr`` / row-frame
``apply_filter``); ``best`` is the per-record primitive and stays
deliberately simple.

## Kind aliases

Defers to ``topiary.ranking.KIND_ALIASES`` for canonical resolution
(``ba`` / ``el`` / ``aff`` / ``ic50`` / ``presentation`` / etc. →
canonical ``pMHC_*``). vaxrank does not maintain its own alias
table. Requires topiary >= 5.12.0 (public name).

Issue: openvax/vaxrank#282 (replaces).
"""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional

from packaging.version import InvalidVersion, Version
from topiary.ranking import KIND_ALIASES as _KIND_ALIASES

if TYPE_CHECKING:
    from mhctools.pred import Prediction
    from .vaccine_antigen import SelfReferenceMatch


def _resolve_kind(kind: str) -> str:
    """Map an alias / case-variant onto the canonical
    ``mhctools.Prediction.kind`` string. Defers to topiary's
    kind-alias table so vaxrank doesn't fork the canonical list.
    Unknown inputs pass through unchanged."""
    return _KIND_ALIASES.get(kind.lower(), kind)


def _nan_safe(x) -> tuple[bool, float]:
    """Wrap a float so it sorts deterministically even when NaN or
    None. Python's ``sorted`` is undefined on NaN keys (NaN
    comparisons return False both ways), and ``None`` doesn't
    compare with floats at all. Returns ``(is_missing, x_or_zero)``
    — NaN / None entries sort *after* finite ones; within the
    missing bucket they all tie at ``(True, 0.0)`` and stable sort
    preserves their relative order from the preceding key
    components, which is the best we can do."""
    if x is None:
        return (True, 0.0)
    is_nan = math.isnan(x)
    return (is_nan, 0.0 if is_nan else x)


def _sort_versions(versions) -> list:
    """Sort ``predictor_version`` strings: legacy / unparseable
    strings come first (lex-sorted), valid PEP 440 versions follow
    oldest → newest. The last element is therefore the *most recent*
    valid version when any are valid, else the last legacy string
    lex-sorted — that's what the leaf accessors default to."""
    parsed, fallback = [], []
    for v in versions:
        try:
            parsed.append((Version(v), v))
        except (InvalidVersion, TypeError):
            fallback.append(v)
    return sorted(fallback) + [v for _, v in sorted(parsed)]


def _group_predictions(predictions) -> dict:
    """Flat ``Sequence[Prediction]`` → nested
    ``{kind: {predictor: {version: tuple}}}``. Producer-friendly:
    most call sites already have a flat list of ``Prediction``s
    coming out of mhctools / topiary."""
    nested: dict = {}
    for p in predictions:
        (nested.setdefault(p.kind, {})
               .setdefault(p.predictor_name, {})
               .setdefault(p.predictor_version, [])
               .append(p))
    return {
        k: {pn: {pv: tuple(records) for pv, records in by_version.items()}
            for pn, by_version in by_predictor.items()}
        for k, by_predictor in nested.items()
    }


# Reserved comparator names. Open-ended — callers can add new ones —
# but these are the canonical entries vaxrank knows about today or
# will populate as the related issues land. Documented here so
# downstream readers know what to expect.
COMPARATOR_WT = 'wt'                                # same position, ref allele
COMPARATOR_NEAREST_SELF = 'nearest_self'            # closest self-peptide (#254)
COMPARATOR_NEAREST_VITAL_SELF = 'nearest_vital_self'  # closest in vital tissues (#254)
COMPARATOR_NEAREST_NONCTA = 'nearest_nonCTA'        # closest non-CTA self (#257)
COMPARATOR_NEAREST_ONCOVIRUS = 'nearest_oncovirus'  # closest oncoviral peptide (#258)


# Source-class constants — what the candidate was *derived from*.
# Open-ended (callers can introduce new classes); these are the
# canonical entries vaxrank knows about today.
SOURCE_CLASS_MUTATION = 'mutation'  # somatic-mutation-derived neoepitope
SOURCE_CLASS_VIRUS = 'virus'        # oncoviral / pathogen-derived
SOURCE_CLASS_SELF = 'self'          # self-protein-derived (CTAs, overexpressed
                                    # WT targets, lineage antigens, etc.)


@dataclass(frozen=True)
class Peptide:
    """One peptide sequence + flanking residues + (optional) MHC
    binding predictions.

    Generic — the same shape is used for reference comparator
    peptides (WT pair, nearest_self, …) and is the base class for
    ``CandidateEpitope`` (the vaccine-candidate itself). When a
    comparator carries predictions, it's answering "does the
    patient's MHC also bind this comparator?" — which is the actual
    safety question for autoimmunity scoring.

    ``sequence`` may be empty for an *anonymous comparator* — pVACseq
    inputs carry ``IC50 WT`` without ``WT Epitope Seq`` for some
    rows, and the resulting WT ``Peptide`` keeps the predictions
    while leaving the sequence blank. Renderers fall back to a
    blank display; predictions are accessed via
    ``predictions_flat()``.

    Hashability: instances are hashable by their *position identity*
    — ``(sequence, source_sequence, offset)`` — not by the full
    field set. The ``predictions`` dict is mutable so it can't
    participate in a hash, but two peptides at the same position
    are the natural dedup key (and match the grouping in
    ``candidate_epitopes_from_rows``). Equality is still
    field-wise, so the hash may collide for peptides with the same
    position but different predictions / flanks — that's a legal
    hash collision, resolved by ``__eq__``.

    See module docstring for the storage shape and disambiguation
    contract.
    """

    sequence: str
    n_flank: str = ''
    c_flank: str = ''
    # Source provenance — for the candidate: the assembled source
    # protein fragment (mutant protein, viral ORF, CTA protein, etc.).
    # For ``nearest_self``: the matching self protein. For
    # ``nearest_oncovirus``: the viral genome / ORF.
    source_sequence: str = ''
    source_name: str = ''         # gene / virus / "self_proteome" / …
    offset: int = 0
    # Nested storage; flat ``list[Prediction]`` is auto-grouped on
    # construction. Shape documented in the module docstring.
    predictions: dict = field(default_factory=dict)

    def __post_init__(self):
        # ``object.__setattr__`` is the escape hatch for writing to
        # a frozen-dataclass field.
        if isinstance(self.predictions, (list, tuple)):
            object.__setattr__(
                self, 'predictions', _group_predictions(self.predictions))

    def __hash__(self) -> int:
        # Position identity. The dataclass-generated hash would
        # include ``predictions`` (a mutable dict), which raises
        # TypeError. ``(sequence, source_sequence, offset)`` is the
        # natural dedup key — matches the grouping used by
        # ``candidate_epitopes_from_rows`` and lets these instances
        # land in sets / dict keys without surprise.
        return hash((self.sequence, self.source_sequence, self.offset))

    def kinds(self) -> tuple[str, ...]:
        """Sorted tuple of distinct ``kind`` values in this context.
        Empty when the context has no predictions (e.g. a
        ``nearest_self`` comparator that's just a sequence with no
        binding probe yet)."""
        return tuple(sorted(self.predictions))

    def predictors_for(self, kind: str) -> tuple[str, ...]:
        """Predictors that emitted ``kind`` for this peptide.
        Drives multi-predictor disambiguation in ``best`` and
        ``predictions_for``. Accepts kind aliases (``'ba'`` /
        ``'el'`` / …)."""
        canonical = _resolve_kind(kind)
        return tuple(sorted(self.predictions.get(canonical, {})))

    def versions_for(self, kind: str, *, predictor: str) -> tuple[str, ...]:
        """All versions of ``predictor`` that emitted ``kind``.

        Ordering: legacy / unparseable strings (lex-sorted) come
        first, followed by valid PEP 440 versions oldest → newest.
        That matches ``best`` / ``predictions_for`` defaulting to
        the *last* entry — i.e. the most-recent valid version when
        any exist, with legacy strings falling below.

        ``predictor`` is keyword-only for consistency with
        ``alleles_for`` / ``predictions_for``.
        """
        canonical = _resolve_kind(kind)
        by_version = self.predictions.get(canonical, {}).get(predictor, {})
        return tuple(_sort_versions(by_version))

    def alleles_for(self, kind: str, *,
                    predictor: Optional[str] = None,
                    version: Optional[str] = None) -> tuple[str, ...]:
        """Alleles attested for ``kind``, optionally restricted by
        ``predictor`` and / or ``version``. Filters out empty
        alleles (processing kinds emit ``allele=""``).

        Predictor- and version-agnostic by default — alleles are a
        set property, not predictor-specific, so when multiple
        predictors emitted ``kind`` we union their alleles rather
        than raise. This is unlike ``best`` / ``predictions_for``,
        where score scales differ across predictors so ranking
        requires explicit disambiguation.

        Explicit ``predictor=`` / ``version=`` arguments are
        validated against what's actually present — a typo raises
        ``ValueError`` rather than silently returning ``()``.
        Missing kind itself still returns ``()`` (asking about a
        kind we don't have isn't a typo).
        """
        canonical = _resolve_kind(kind)
        by_predictor = self.predictions.get(canonical, {})
        if not by_predictor:
            return ()
        if predictor is not None and predictor not in by_predictor:
            raise ValueError(
                f"{canonical} has no predictions from predictor "
                f"{predictor!r}; available: {sorted(by_predictor)}.")
        leaves = [
            (ver, records)
            for pn, by_ver in by_predictor.items()
            if predictor is None or pn == predictor
            for ver, records in by_ver.items()
        ]
        available_versions = {ver for ver, _ in leaves}
        if version is not None and version not in available_versions:
            raise ValueError(
                f"{canonical} has no predictions at version "
                f"{version!r}; available: "
                f"{sorted(available_versions)}.")
        return tuple(sorted({
            p.allele
            for ver, records in leaves
            if version is None or ver == version
            for p in records
            if p.allele
        }))

    def predictions_for(self, kind: str, *,
                        predictor: Optional[str] = None,
                        version: Optional[str] = None
                        ) -> tuple["Prediction", ...]:
        """All ``Prediction`` records for ``(kind, predictor,
        version)``. Empty tuple when nothing matches.

        Predictor disambiguation:
        - One predictor emitted ``kind`` → ``predictor`` may be
          omitted.
        - Multiple predictors emitted ``kind`` → ``predictor`` is
          required; raises ``ValueError`` otherwise.

        Version disambiguation:
        - Multiple versions on file → defaults to the most recent
          (PEP 440 ordering). Pass ``version=`` to pin.

        Use this to rank predictions on a non-default axis (e.g.
        ``min(ctx.predictions_for('affinity'), key=lambda p:
        p.percentile_rank)``) — ``best`` does the common case
        (``max`` over normalized ``score``) but isn't parametric.
        """
        canonical = _resolve_kind(kind)
        by_predictor = self.predictions.get(canonical, {})
        if not by_predictor:
            return ()
        if predictor is None:
            if len(by_predictor) > 1:
                raise ValueError(
                    f"{canonical} has predictions from multiple "
                    f"predictors ({sorted(by_predictor)}); pass "
                    f"predictor= to disambiguate. Each predictor's "
                    f"score scale is its own — cross-predictor "
                    f"ranking is not meaningful.")
            predictor = next(iter(by_predictor))
        by_version = by_predictor.get(predictor, {})
        if not by_version:
            return ()
        if version is None:
            version = _sort_versions(by_version)[-1]
        return by_version.get(version, ())

    def best(self, kind: str, *,
             predictor: Optional[str] = None,
             version: Optional[str] = None
             ) -> Optional["Prediction"]:
        """Best ``Prediction`` for ``(kind, predictor, version)``
        across alleles, ranked by ``Prediction.score``. Returns
        ``None`` when no record matches.

        ``kind`` accepts topiary's aliases (``'ba'`` / ``'el'`` /
        ``'aff'`` / ``'ic50'`` / ``'presentation'`` / …).

        For a different ranking axis (e.g. lowest IC50, lowest
        %-rank), call ``predictions_for(...)`` and rank the leaf
        tuple yourself — ``best`` is deliberately the common-case
        primitive only.

        CONTRACT: relies on mhctools normalizing
        ``Prediction.score`` so higher = better *within one
        predictor*. mhctools enforces this on producers; vaxrank
        consumes ``score`` directly. Cross-predictor max is *not*
        meaningful — handled by the ``predictor`` disambiguation
        in ``predictions_for``.
        """
        candidates = self.predictions_for(
            kind, predictor=predictor, version=version)
        if not candidates:
            return None
        return max(candidates, key=lambda p: p.score)

    # Kind-named sugar for the common ``best(kind)`` calls. For
    # predictor / version pinning, use ``best`` directly.
    def best_affinity(self): return self.best('pMHC_affinity')
    def best_presentation(self): return self.best('pMHC_presentation')
    def best_stability(self): return self.best('pMHC_stability')
    def best_immunogenicity(self): return self.best('immunogenicity')
    def best_cleavage(self): return self.best('proteasome_cleavage')
    def best_antigen_processing(self): return self.best('antigen_processing')

    def predictions_flat(self) -> tuple["Prediction", ...]:
        """Flatten the nested store back to a tuple of
        ``Prediction`` records. Sorted by ``(kind, predictor_name,
        predictor_version, allele, score, value, percentile_rank)``
        so the output is fully deterministic — even if a producer
        ever emits duplicate ``(kind, predictor, version, allele)``
        records, the score / value / %-rank tail breaks ties
        without falling back on input order. Float keys are wrapped
        with ``_nan_safe`` because mhctools sometimes emits NaN
        for ``percentile_rank`` and ``sorted`` is undefined on NaN."""
        flat = (
            p
            for by_predictor in self.predictions.values()
            for by_version in by_predictor.values()
            for records in by_version.values()
            for p in records)
        return tuple(sorted(flat, key=lambda p: (
            p.kind, p.predictor_name, p.predictor_version, p.allele,
            _nan_safe(p.score), _nan_safe(p.value),
            _nan_safe(p.percentile_rank))))

    def _slice_fits(self, start_offset: int, end_offset: int) -> bool:
        """True when this peptide fits inside the requested window —
        i.e. ``sliced()`` would not return None. Shared between
        ``Peptide.sliced`` and ``CandidateEpitope.sliced``."""
        return (self.offset >= start_offset
                and self.offset + len(self.sequence) <= end_offset)

    def sliced(self, start_offset: int,
               end_offset: int) -> Optional["Peptide"]:
        """Return a new ``Peptide`` whose source window is narrowed
        to ``[start_offset, end_offset)``. Returns ``None`` if the
        peptide doesn't fit inside the window.

        The peptide sequence, flanks, source name, and predictions
        are preserved — only ``source_sequence`` is sliced and
        ``offset`` is rebased relative to the new source start.
        """
        if not self._slice_fits(start_offset, end_offset):
            return None
        # Deep copy of predictions: the nested
        # ``{kind: {predictor: {version: tuple}}}`` chain is mutable
        # at every level above the leaf tuple. Sharing it across two
        # frozen instances would couple their state — a mutation
        # anywhere in the chain would leak between slices. The
        # ``Prediction`` leaves themselves are NamedTuple-like and
        # immutable; deepcopy is fast on the dict spine.
        return Peptide(
            sequence=self.sequence,
            n_flank=self.n_flank,
            c_flank=self.c_flank,
            source_sequence=self.source_sequence[start_offset:end_offset],
            source_name=self.source_name,
            offset=self.offset - start_offset,
            predictions=copy.deepcopy(self.predictions))


@dataclass(frozen=True)
class CandidateEpitope(Peptide):
    """A vaccine-candidate epitope. *Is-a* ``Peptide`` — carries the
    candidate's sequence + predictions directly — plus comparator
    peptides for safety / immunogenicity scoring and a couple of
    optional source-context flags.

    Source-agnostic: works for mutation-derived neoepitopes, oncoviral
    peptides, and self-derived candidates (CTAs etc.). Each source
    class populates the comparators it knows how to compute; the
    wrapper itself doesn't branch on source.

    Comparator names (the ``comparators`` dict keys) are open-ended;
    vaxrank populates these as their respective issues land:

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

    Each value is a bare ``Peptide``. When a comparator carries
    predictions, scorers can ask "does the patient's MHC also bind
    this comparator?" — the actual safety question for autoimmunity
    / off-target presentation.
    """

    comparators: dict = field(default_factory=dict)

    # Source provenance — what the candidate was derived from. One of
    # ``SOURCE_CLASS_MUTATION`` / ``SOURCE_CLASS_VIRUS`` /
    # ``SOURCE_CLASS_SELF`` (or None for legacy / unspecified inputs).
    # Producers set this; consumers can branch on it for source-aware
    # rendering or filtering, but most ranking logic stays
    # source-agnostic via the safety flags below.
    source_class: Optional[str] = None

    # Optional mutation-specific provenance — populated only when
    # ``source_class == SOURCE_CLASS_MUTATION``. Used internally by
    # ``predict_epitopes`` to gate WT-comparator generation and by
    # mutation-flavored reports / audits. Viral / self leave at False.
    overlaps_mutation: bool = False

    # Source-agnostic targetability. ``None`` means the legacy producer did
    # not supply a mask result and remains eligible for backwards-compatible
    # target/self splitting. New producers always supply a boolean.
    overlaps_targetable: Optional[bool] = None

    # Raw safety signal: this peptide's exact sequence appears
    # somewhere in the patient's reference proteome. True means a
    # cross-reactivity risk (the peptide matches a self-protein).
    # Conservative — flags CTA peptides as unsafe even when targeting
    # them is intended.
    occurs_in_reference: bool = False

    # CTA-aware safety signal: same as ``occurs_in_reference`` but with
    # source genes in oncoref's full CTA candidate universe excluded from
    # the reference set. Shared sequences in non-CTA genes remain present.
    # Consumers opt into CTA-friendly policy by reading this flag instead
    # of the raw one.
    occurs_in_non_CTA_reference: bool = False

    # Antigen-aware exact-self result. New producers populate this and
    # consumers read ``occurs_in_self_reference`` below. ``None`` retains the
    # raw-reference behavior for legacy files and external-input loaders.
    self_reference_match: Optional["SelfReferenceMatch"] = None

    # Patient alleles for which this candidate has input evidence. This is
    # usually identical to the non-empty alleles on prediction leaves, but
    # must be retained separately for allele-independent evidence such as
    # antigen processing. External formats like LENS repeat that evidence on
    # one row per patient allele while the canonical Prediction remains
    # allele-less.
    patient_alleles: tuple[str, ...] = ()

    # Stable identity supplied by external prediction loaders. It encodes
    # variant, gene/transcript, peptide context, peptide, and offset; see
    # ``ExternalPredictionKey``. Empty for native pipeline predictions, whose
    # position identity remains ``(source_sequence, sequence, offset)``.
    prediction_id: str = ""

    # Per-allele score for this epitope as computed by the configured
    # :class:`~vaxrank.epitope_config.EpitopeConfig` ``score_expr``
    # (see :mod:`vaxrank.epitope_dsl`). Keys are allele names exactly
    # as they appear on the predictions; values are the DSL output
    # summed over predictor leaves with the same (peptide, allele).
    # Populated by :func:`~vaxrank.epitope_logic.predict_epitopes`
    # (and equivalent loaders); empty when an epitope is loaded
    # without scores yet attached. Downstream consumers MUST read
    # from here and never recompute — the DSL is the single source
    # of truth.
    per_allele_scores: dict = field(default_factory=dict)

    # Why each allele was credited with this candidate's allele-independent
    # evidence — reported by the input, carried in the patient's genotype, or
    # selected by ranking the peptide's own allele-scoped predictions (with
    # the axis, predictor, value and rank that produced the selection).
    # Recorded rather than recomputed so a finished run can be reconstructed:
    # see :mod:`vaxrank.allele_evidence`. Empty when a candidate carries no
    # allele-independent evidence, which is the common case.
    allele_attributions: tuple = ()

    def __post_init__(self):
        """Normalize predictions and retain every explicitly known allele."""
        super().__post_init__()
        alleles = {
            allele
            for allele in self.patient_alleles
            if allele
        }
        alleles.update(
            prediction.allele
            for prediction in self.predictions_flat()
            if prediction.allele
        )
        alleles.update(
            allele for allele in self.per_allele_scores if allele)
        object.__setattr__(self, 'patient_alleles', tuple(sorted(alleles)))

    @property
    def epitope_score(self) -> float:
        """Per-epitope total score = sum of per-allele scores.

        This is the binding fed to the combined-score DSL as
        ``target_epitope_score`` (aggregated across an epitope set)
        and the value used by the vaccine-peptide ranker. Defined
        here exactly once; no parallel formula exists in the
        codebase.
        """
        return float(sum(self.per_allele_scores.values()))

    @property
    def occurs_in_self_reference(self) -> bool:
        """Exact-self status under this epitope's antigen policy."""
        if self.self_reference_match is not None:
            return self.self_reference_match.occurs
        return self.occurs_in_reference

    @property
    def is_targetable(self) -> bool:
        """Whether this epitope overlaps explicitly targetable content."""
        if self.overlaps_targetable is not None:
            return self.overlaps_targetable
        return True

    def __hash__(self) -> int:
        # Inherited from Peptide in spirit, but ``@dataclass(frozen=True)``
        # regenerates ``__hash__`` on every subclass — and ours would
        # include the mutable ``comparators`` dict. External candidates add
        # their provenance ID to the normal peptide-position identity.
        position = (self.sequence, self.source_sequence, self.offset)
        if not self.prediction_id:
            return hash(position)
        return hash((self.prediction_id,) + position)

    @property
    def prediction_group_source(self) -> str:
        """First component of the Topiary score-group identity."""
        return self.prediction_id or self.source_sequence or self.sequence

    @property
    def prediction_group_key(self) -> tuple[str, str, int]:
        """Exact source/peptide/offset identity used for score joins."""
        return self.prediction_group_source, self.sequence, self.offset

    # Convenience: most call sites today read the WT alongside the
    # candidate. Common-case shortcut.
    @property
    def wt(self) -> Optional[Peptide]:
        return self.comparators.get(COMPARATOR_WT)

    def comparator(self, name: str) -> Optional[Peptide]:
        """Generic accessor for any comparator. Returns ``None``
        when the comparator isn't present for this candidate
        (typical when a safety scorer hasn't run yet)."""
        return self.comparators.get(name)

    @property
    def length(self) -> int:
        return len(self.sequence)

    def sliced(self, start_offset: int,
               end_offset: int) -> Optional["CandidateEpitope"]:
        """Return a new ``CandidateEpitope`` whose source window is
        narrowed to ``[start_offset, end_offset)``. Returns ``None``
        if the candidate peptide doesn't fit inside the window.
        Comparators (WT etc.) keep their own source windows
        unchanged — they're independent peptides, not slices of
        the candidate's source.
        """
        if not self._slice_fits(start_offset, end_offset):
            return None
        # See ``Peptide.sliced`` for why predictions get a deep copy.
        # Comparators are independent ``Peptide`` instances and pass
        # through by reference — slicing the candidate doesn't slice
        # them. The sliced instance shares the same ``comparators``
        # dict as the source; that's the documented contract.
        return CandidateEpitope(
            sequence=self.sequence,
            n_flank=self.n_flank,
            c_flank=self.c_flank,
            source_sequence=self.source_sequence[start_offset:end_offset],
            source_name=self.source_name,
            offset=self.offset - start_offset,
            predictions=copy.deepcopy(self.predictions),
            comparators=self.comparators,
            source_class=self.source_class,
            overlaps_mutation=self.overlaps_mutation,
            overlaps_targetable=self.overlaps_targetable,
            occurs_in_reference=self.occurs_in_reference,
            occurs_in_non_CTA_reference=self.occurs_in_non_CTA_reference,
            self_reference_match=self.self_reference_match,
            patient_alleles=self.patient_alleles,
            prediction_id=self.prediction_id,
            per_allele_scores=dict(self.per_allele_scores),
            allele_attributions=self.allele_attributions)

    @classmethod
    def from_peptide(cls, peptide: "Peptide",
                     **extras) -> "CandidateEpitope":
        """Alternate constructor: spread a ``Peptide``'s fields into a
        ``CandidateEpitope`` along with the candidate-specific
        ``extras`` (``comparators``, ``overlaps_mutation``,
        ``occurs_in_reference``).

        Convenient when a producer or test already has a fully-built
        ``Peptide`` in hand. Equivalent to passing the peptide fields
        as keyword args to the regular constructor.
        """
        return cls(
            sequence=peptide.sequence,
            n_flank=peptide.n_flank,
            c_flank=peptide.c_flank,
            source_sequence=peptide.source_sequence,
            source_name=peptide.source_name,
            offset=peptide.offset,
            predictions=dict(peptide.predictions),
            **extras,
        )


def candidate_epitopes_from_rows(rows) -> list["CandidateEpitope"]:
    """Group leaf predictions into ``CandidateEpitope`` objects.

    Producers (``predict_epitopes``, LENS / pVACseq loaders) iterate
    row-shaped inputs (topiary frame, TSV) and build a flat list of
    row dicts; this function groups them by ``(prediction_id, peptide,
    source, offset)`` and emits one ``CandidateEpitope`` per group, in
    first-seen order.

    Each row is a mapping with the following keys:

      ``peptide``    str             — candidate peptide sequence
      ``source``     str             — source protein / transcript window
      ``source_name`` str, optional  — human-readable name of that window
                                       (transcript, variant, contig …)
      ``n_flank`` / ``c_flank``      str, optional — residues flanking the
                                       peptide in its source protein; carried
                                       because processing predictors score
                                       them
      ``offset``     int             — peptide offset within ``source``
      ``prediction_id`` str, optional — complete external provenance key
      ``mutant``     Prediction      — leaf prediction for one (allele,
                                       predictor, version) cell
      ``wt``         Prediction|None — parallel WT prediction, optional
      ``source_class``               str|None, default None
      ``overlaps_mutation``          bool,     default False
      ``overlaps_targetable``        bool|None, default None for legacy rows
      ``occurs_in_reference``        bool,     default False
      ``occurs_in_non_CTA_reference`` bool,    default False — when
                                     producers don't yet integrate a
                                     CTA database, pass the same
                                     value as ``occurs_in_reference``.
      ``self_reference_match``       SelfReferenceMatch|None — explicit
                                     antigen-aware exact-self result.
      ``patient_alleles``            iterable[str], optional — alleles
                                     explicitly associated with this input
                                     row even when ``mutant.allele`` is empty.
      ``per_allele_scores``          mapping[str, float], optional — scores
                                     already evaluated by the configured DSL.

    Semantics:

    - A ``wt`` whose ``peptide`` equals the mutant peptide is treated
      as no WT signal (self-comparator is meaningless) and dropped.
    - A ``wt`` whose ``peptide`` is empty (anonymous WT — pVACseq with
      ``IC50 WT`` but no ``WT Epitope Seq``) is kept; the resulting
      WT comparator ``Peptide`` carries the predictions but an empty
      ``sequence``. Downstream renderers display blank for the
      sequence and use ``predictions_flat()`` for the IC50.
    - ``overlaps_mutation`` / ``occurs_in_reference`` OR across rows
      in the same group: any True row makes the group True. The flags
      are properties of the (peptide, source, offset) position, not
      of an individual leaf record — if a producer disagrees across
      rows, OR preserves the True signal rather than depending on
      iteration order.
    """
    groups: dict = {}
    for row in rows:
        peptide = row['peptide']
        prediction_id = row.get('prediction_id') or ''
        key = (prediction_id, peptide, row['source'], row['offset'])
        slot = groups.get(key)
        if slot is None:
            slot = {
                'peptide': peptide,
                'source': row['source'],
                'source_name': row.get('source_name') or '',
                'n_flank': row.get('n_flank') or '',
                'c_flank': row.get('c_flank') or '',
                'offset': row['offset'],
                'prediction_id': prediction_id,
                'mutant_preds': [],
                'wt_preds': [],
                'wt_peptide': None,
                'source_class': row.get('source_class'),
                'overlaps_mutation': False,
                'overlaps_targetable': None,
                'occurs_in_reference': False,
                'occurs_in_non_CTA_reference': False,
                'self_reference_match': None,
                'patient_alleles': set(),
                # Accumulator for per-allele DSL scores. Same allele
                # may appear in multiple rows (different predictors);
                # all rows for one allele share the same score by
                # construction of the DSL eval, so last-write-wins.
                'per_allele_scores': {},
                'allele_attributions': (),
            }
            groups[key] = slot
        slot['mutant_preds'].append(row['mutant'])
        if row['mutant'].allele:
            slot['patient_alleles'].add(row['mutant'].allele)
        slot['patient_alleles'].update(
            allele for allele in (row.get('patient_alleles') or ()) if allele)
        for allele, score in (row.get('per_allele_scores') or {}).items():
            score = float(score)
            existing_score = slot['per_allele_scores'].get(allele)
            if existing_score is not None and existing_score != score:
                raise ValueError(
                    "Conflicting per-allele scores for one candidate epitope")
            slot['per_allele_scores'][allele] = score
        attributions = row.get('allele_attributions') or ()
        if attributions:
            existing = slot['allele_attributions']
            if existing and existing != tuple(attributions):
                raise ValueError(
                    "Conflicting allele attributions for one candidate "
                    "epitope")
            slot['allele_attributions'] = tuple(attributions)
        slot['overlaps_mutation'] |= bool(row.get('overlaps_mutation', False))
        overlaps_targetable = row.get('overlaps_targetable')
        if overlaps_targetable is not None:
            existing_targetable = slot['overlaps_targetable']
            if (
                existing_targetable is not None
                and existing_targetable != bool(overlaps_targetable)
            ):
                raise ValueError(
                    "Conflicting targetable-mask results for one candidate epitope"
                )
            slot['overlaps_targetable'] = bool(overlaps_targetable)
        slot['occurs_in_reference'] |= bool(row.get('occurs_in_reference', False))
        slot['occurs_in_non_CTA_reference'] |= bool(
            row.get('occurs_in_non_CTA_reference', False))
        self_reference_match = row.get('self_reference_match')
        if self_reference_match is not None:
            existing_match = slot['self_reference_match']
            if existing_match is not None and existing_match != self_reference_match:
                raise ValueError(
                    "Conflicting self-reference matches for one candidate epitope"
                )
            slot['self_reference_match'] = self_reference_match
        if 'allele_score' in row and row['mutant'].allele:
            # Allele-free leaves (antigen_processing, proteasome_cleavage)
            # carry a score for the peptide, not for an allele. Writing it
            # under '' would put a non-allele entry in a map whose contract
            # is explicitly per patient allele, and epitope_score sums the
            # map — so the candidate would score its processing evidence as
            # if it were an extra allele. attach_per_allele_scores applies
            # the same rule on the loader path.
            allele = row['mutant'].allele
            score = float(row['allele_score'])
            existing = slot['per_allele_scores'].get(allele)
            if existing is not None and existing != score:
                raise ValueError(
                    "Conflicting per-allele scores for one candidate epitope")
            slot['per_allele_scores'][allele] = score
        wt = row.get('wt')
        if wt is not None and wt.peptide and wt.peptide == peptide:
            # Self-WT: comparing the candidate against itself is meaningless.
            wt = None
        if wt is not None:
            slot['wt_preds'].append(wt)
            if slot['wt_peptide'] is None:
                slot['wt_peptide'] = wt.peptide

    out = []
    for slot in groups.values():
        comparators = {}
        if slot['wt_preds']:
            comparators[COMPARATOR_WT] = Peptide(
                sequence=slot['wt_peptide'] or '',
                source_sequence=slot['source'],
                offset=slot['offset'],
                predictions=tuple(slot['wt_preds']),
            )
        out.append(CandidateEpitope(
            sequence=slot['peptide'],
            source_sequence=slot['source'],
            source_name=slot['source_name'],
            n_flank=slot['n_flank'],
            c_flank=slot['c_flank'],
            offset=slot['offset'],
            predictions=tuple(slot['mutant_preds']),
            comparators=comparators,
            source_class=slot['source_class'],
            overlaps_mutation=slot['overlaps_mutation'],
            overlaps_targetable=slot['overlaps_targetable'],
            occurs_in_reference=slot['occurs_in_reference'],
            occurs_in_non_CTA_reference=slot['occurs_in_non_CTA_reference'],
            self_reference_match=slot['self_reference_match'],
            patient_alleles=tuple(sorted(slot['patient_alleles'])),
            prediction_id=slot['prediction_id'],
            per_allele_scores=dict(slot['per_allele_scores']),
            allele_attributions=slot['allele_attributions'],
        ))
    return out
