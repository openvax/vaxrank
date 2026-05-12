# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Per-peptide context records used by ``VaccinePeptide``.

Two layers:

  ``Peptide`` — one amino-acid sequence + flanks + multi-axis
  binding predictions. Generic shape used for the antigenic
  candidate AND for any reference comparator (WT pair, nearest_self,
  nearest_vital_self, nearest_nonCTA, nearest_oncovirus, …).

  ``CandidateEpitope`` — one sliding-window position from a
  VaccinePeptide. Holds the mutant ``Peptide`` plus a
  ``comparators`` dict keyed by name. Future safety / homology
  features (#254 / #257 / #258) populate comparators with their
  respective peptides; today only ``'wt'`` is canonical.

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

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional

from packaging.version import InvalidVersion, Version
from topiary.ranking import KIND_ALIASES as _KIND_ALIASES

if TYPE_CHECKING:
    from mhctools.pred import Prediction


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


@dataclass(frozen=True)
class Peptide:
    """One peptide sequence + flanking residues + (optional) MHC
    binding predictions.

    Generic — the same shape is used for the antigenic candidate
    (the mutant) and for every reference comparator (WT pair,
    nearest_self, …). When a comparator carries predictions, it's
    answering "does the patient's MHC also bind this comparator?"
    — which is the actual safety question for autoimmunity scoring.

    ``sequence`` may be empty for an *anonymous comparator* — pVACseq
    inputs carry ``IC50 WT`` without ``WT Epitope Seq`` for some
    rows, and the resulting WT ``Peptide`` keeps the predictions
    while leaving the sequence blank. Renderers fall back to a
    blank display; predictions are accessed via
    ``predictions_flat()``.

    See module docstring for the storage shape and disambiguation
    contract.
    """

    sequence: str
    n_flank: str = ''
    c_flank: str = ''
    # Source provenance — for the mutant: the assembled mutant
    # protein fragment. For ``nearest_self``: the matching self
    # protein. For ``nearest_oncovirus``: the viral genome / ORF.
    source: str = ''
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

    def sliced(self, start_offset: int,
               end_offset: int) -> Optional["Peptide"]:
        """Return a new ``Peptide`` whose source window is narrowed
        to ``[start_offset, end_offset)``. Returns ``None`` if the
        peptide doesn't fit inside the window.

        The peptide sequence, flanks, source name, and predictions
        are preserved — only ``source`` is sliced and ``offset`` is
        rebased relative to the new source start.
        """
        if self.offset < start_offset:
            return None
        if self.offset + len(self.sequence) > end_offset:
            return None
        # Defensive copy: the nested predictions dict is mutable, so
        # sharing it across two frozen instances would couple their
        # equality / hash semantics in surprising ways. Cheap to copy.
        return Peptide(
            sequence=self.sequence,
            n_flank=self.n_flank,
            c_flank=self.c_flank,
            source=self.source[start_offset:end_offset],
            source_name=self.source_name,
            offset=self.offset - start_offset,
            predictions=dict(self.predictions))


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

    Each value is a ``Peptide`` carrying its own sequence +
    (optionally) its own MHC predictions. When a comparator carries
    predictions, scorers can ask "does the patient's MHC also bind
    this comparator?" — the actual safety question for autoimmunity
    / off-target presentation.
    """

    mutant: Peptide
    comparators: dict = field(default_factory=dict)

    # Mutation-specific context. Lives at this level (not on
    # ``Peptide``) because it's about the mutant-vs-source
    # relationship, not a property of the peptide itself.
    overlaps_mutation: bool = False
    occurs_in_reference: bool = False

    # Convenience: most call sites today read the WT alongside the
    # mutant. Common-case shortcut.
    @property
    def wt(self) -> Optional[Peptide]:
        return self.comparators.get(COMPARATOR_WT)

    def comparator(self, name: str) -> Optional[Peptide]:
        """Generic accessor for any comparator. Returns ``None``
        when the comparator isn't present for this candidate
        (typical when a safety scorer hasn't run yet)."""
        return self.comparators.get(name)

    # Convenience pass-throughs to the mutant peptide — most
    # callers read these on the candidate level.
    @property
    def sequence(self) -> str:
        return self.mutant.sequence

    @property
    def length(self) -> int:
        return len(self.mutant.sequence)

    def sliced(self, start_offset: int,
               end_offset: int) -> Optional["CandidateEpitope"]:
        """Return a new ``CandidateEpitope`` whose mutant source window is
        narrowed to ``[start_offset, end_offset)``. Returns ``None``
        if the mutant peptide doesn't fit inside the window.
        Comparators (WT etc.) keep their own source windows
        unchanged — they're independent peptides, not slices of
        the mutant's source."""
        sliced_mutant = self.mutant.sliced(start_offset, end_offset)
        if sliced_mutant is None:
            return None
        return CandidateEpitope(
            mutant=sliced_mutant,
            comparators=self.comparators,
            overlaps_mutation=self.overlaps_mutation,
            occurs_in_reference=self.occurs_in_reference)


def candidate_epitopes_from_rows(rows) -> list["CandidateEpitope"]:
    """Group leaf predictions into ``CandidateEpitope`` objects.

    Producers (``predict_epitopes``, LENS / pVACseq loaders) iterate
    row-shaped inputs (topiary frame, TSV) and build a flat list of
    row dicts; this function groups them by ``(peptide, source,
    offset)`` and emits one ``CandidateEpitope`` per group, in
    first-seen order.

    Each row is a mapping with the following keys:

      ``peptide``    str             — mutant peptide sequence
      ``source``     str             — source protein / transcript window
      ``offset``     int             — peptide offset within ``source``
      ``mutant``     Prediction      — leaf prediction for one (allele,
                                       predictor, version) cell
      ``wt``         Prediction|None — parallel WT prediction, optional
      ``overlaps_mutation``    bool, default False
      ``occurs_in_reference``  bool, default False

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
        key = (peptide, row['source'], row['offset'])
        slot = groups.get(key)
        if slot is None:
            slot = {
                'peptide': peptide,
                'source': row['source'],
                'offset': row['offset'],
                'mutant_preds': [],
                'wt_preds': [],
                'wt_peptide': None,
                'overlaps_mutation': False,
                'occurs_in_reference': False,
            }
            groups[key] = slot
        slot['mutant_preds'].append(row['mutant'])
        slot['overlaps_mutation'] |= bool(row.get('overlaps_mutation', False))
        slot['occurs_in_reference'] |= bool(row.get('occurs_in_reference', False))
        wt = row.get('wt')
        if wt is not None and wt.peptide and wt.peptide == peptide:
            # Self-WT: comparing the mutant against itself is meaningless.
            wt = None
        if wt is not None:
            slot['wt_preds'].append(wt)
            if slot['wt_peptide'] is None:
                slot['wt_peptide'] = wt.peptide

    out = []
    for slot in groups.values():
        mutant_pep = Peptide(
            sequence=slot['peptide'],
            source=slot['source'],
            offset=slot['offset'],
            predictions=tuple(slot['mutant_preds']),
        )
        comparators = {}
        if slot['wt_preds']:
            comparators[COMPARATOR_WT] = Peptide(
                sequence=slot['wt_peptide'] or '',
                source=slot['source'],
                offset=slot['offset'],
                predictions=tuple(slot['wt_preds']),
            )
        out.append(CandidateEpitope(
            mutant=mutant_pep,
            comparators=comparators,
            overlaps_mutation=slot['overlaps_mutation'],
            occurs_in_reference=slot['occurs_in_reference'],
        ))
    return out
