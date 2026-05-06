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
are independent transforms. So ``best(kind)`` raises when multiple
predictors emitted ``kind`` and the caller hasn't passed
``predictor=``. Callers that legitimately want a per-predictor
loop iterate ``predictors_for_kind(kind)`` and call
``best(kind, predictor=…)`` per predictor.

Versions, however, *do* auto-resolve: when one predictor has
multiple ``predictor_version`` strings on file (e.g. mhcflurry
2.1.0 and 2.1.1 both ran), ``best`` defaults to the most recent
version (PEP 440 ordering). Pass ``version=`` to pin to a specific
version explicitly.

## Kind aliases

mhcflurry / netmhcpan / pVACseq / LENS use shorthand for the same
underlying kinds (``BA`` / ``EL`` / etc). ``best`` and friends
accept those aliases (case-insensitive) so callers don't have to
remember which canonical string vaxrank stores. The canonical
``pMHC_*`` strings keep working unchanged.

## Parametric scoring

``best`` accepts ``score_key`` — a ``Prediction -> float`` callable
that names the ranking axis (higher = better). Default is
:func:`default_score_key` (uses ``Prediction.score``, mhctools'
normalized higher-is-better axis). Pass a custom key for a
different axis (``lambda p: -p.percentile_rank`` for "lowest %-rank
wins") so the score-normalization assumption isn't load-bearing in
silent ways. Higher-level pipelines drive ranking via the topiary
DSL on row-frames (``EpitopeConfig.score_expr``); ``best`` stays
the per-record primitive.

Issue: openvax/vaxrank#282 (replaces).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import TYPE_CHECKING, Callable, Optional

from packaging.version import InvalidVersion, Version

if TYPE_CHECKING:
    from mhctools.pred import Prediction


# Lowercase alias → canonical mhctools.Prediction.kind. The canonical
# strings also resolve through this map (since ``str.lower`` is
# applied first), so callers can pass either form. Wrapped in a
# read-only ``MappingProxyType`` so accidental monkeypatching from
# tests doesn't bleed into other tests.
_KIND_ALIASES: MappingProxyType = MappingProxyType({
    # pMHC_affinity — "BA" in mhcflurry / netmhcpan parlance.
    'pmhc_affinity': 'pMHC_affinity',
    'affinity': 'pMHC_affinity',
    'ba': 'pMHC_affinity',
    'binding': 'pMHC_affinity',
    # pMHC_presentation — "EL" (elution likelihood) in mhcflurry-pres.
    'pmhc_presentation': 'pMHC_presentation',
    'presentation': 'pMHC_presentation',
    'el': 'pMHC_presentation',
    'elution': 'pMHC_presentation',
    # pMHC_stability.
    'pmhc_stability': 'pMHC_stability',
    'stability': 'pMHC_stability',
    # Single-word kinds — keep canonical and short alias both routed.
    'immunogenicity': 'immunogenicity',
    'antigen_processing': 'antigen_processing',
    'processing': 'antigen_processing',
    'ap': 'antigen_processing',
    'proteasome_cleavage': 'proteasome_cleavage',
    'cleavage': 'proteasome_cleavage',
    'proteasome': 'proteasome_cleavage',
    'tap_transport': 'tap_transport',
    'tap': 'tap_transport',
    'erap_trimming': 'erap_trimming',
    'erap': 'erap_trimming',
})


def _resolve_kind(kind: str) -> str:
    """Map an alias / case-variant onto the canonical
    ``mhctools.Prediction.kind`` string. Unknown inputs pass through
    unchanged (so future predictor kinds work without a registry
    update — they just won't have short aliases until added)."""
    return _KIND_ALIASES.get(kind.lower(), kind)


def _split_versions(versions) -> tuple:
    """Split a set of ``predictor_version`` strings into PEP 440
    parseable + fallback. Returns ``(parsed_oldest_to_newest,
    fallback_lex_sorted)`` — both are lists of the original strings.
    Shared helper for the version-disambiguation paths."""
    parsed = []
    fallback = []
    for v in versions:
        try:
            parsed.append((Version(v), v))
        except (InvalidVersion, TypeError):
            fallback.append(v)
    return ([v for _, v in sorted(parsed, key=lambda x: x[0])],
            sorted(fallback))


def _most_recent_version(versions) -> str:
    """Pick the most recent ``predictor_version`` from a set.

    PEP 440 ordering via ``packaging.version.Version``. Strings that
    don't parse (empty / non-semver / pre-#261 unset) sort below any
    valid version — so when both styles coexist, valid wins; when
    none parse, lexicographic ``max`` is the fallback. Returning a
    single string lets callers filter ``candidates`` by exact match.
    """
    parsed, fallback = _split_versions(versions)
    if parsed:
        return parsed[-1]
    return fallback[-1]


def default_score_key(prediction) -> float:
    """Default ranking key for ``PeptideContext.best`` and friends.

    CONTRACT: relies on mhctools normalizing ``Prediction.score`` so
    that **higher = better within one predictor** — across IC50
    (lower-is-better-natively) and probability (higher-is-better-
    natively) underlying axes. The contract is mhctools' to enforce
    on producers; vaxrank just consumes ``score`` directly.

    Override by passing ``score_key=`` to ``best`` if you want a
    different axis (e.g. ``score_key=lambda p: -p.percentile_rank``
    for "lowest %-rank wins") or a custom DSL-derived ranking. Higher
    return value = better, regardless of underlying scale.

    Cross-predictor comparison is *still* not meaningful even with a
    custom ``score_key`` — that disambiguation is enforced by
    ``best`` raising ``ValueError`` when multiple predictors emit
    ``kind`` and the caller hasn't pinned ``predictor=``.
    """
    return prediction.score


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
    # ``predictions`` is a flat tuple for serializability. Annotated
    # under ``TYPE_CHECKING`` to give static analyzers the full type
    # without importing ``mhctools.Prediction`` at runtime (mhctools
    # pulls in heavy deps).
    predictions: tuple["Prediction", ...] = ()

    # ------------------------------------------------------------------
    # Structured views — read-only, computed on demand.
    # ------------------------------------------------------------------

    def predictions_by_kind_and_predictor(self) -> dict:
        """Four-level nested view:
        ``{kind: {predictor_name: {predictor_version: {allele:
        Prediction}}}}``.

        Built fresh each call (the storage tuple is the source of
        truth). Use this when you need full structure; the
        ``best_*`` methods cover the common case. The version axis
        prevents collision when the same (kind, predictor, allele)
        was scored at multiple ``predictor_version`` values — every
        record stays addressable.
        """
        nested: dict = {}
        for p in self.predictions:
            (nested.setdefault(p.kind, {})
                   .setdefault(p.predictor_name, {})
                   .setdefault(p.predictor_version, {})[p.allele]) = p
        return nested

    def kinds(self) -> tuple:
        """Sorted tuple of distinct ``kind`` values in this context.
        Empty when the context has no predictions (e.g. a
        ``nearest_self`` comparator that's just a sequence with no
        binding probe yet)."""
        return tuple(sorted({p.kind for p in self.predictions}))

    def predictors_for_kind(self, kind: str) -> tuple:
        """Predictors that emitted ``kind`` for this peptide.
        Drives multi-predictor disambiguation in ``best``. Accepts
        kind aliases (``'ba'`` / ``'el'`` / …)."""
        canonical = _resolve_kind(kind)
        return tuple(sorted({
            p.predictor_name for p in self.predictions
            if p.kind == canonical}))

    def versions_for(self, kind: str, predictor: str) -> tuple:
        """Versions of ``predictor`` that emitted ``kind``, sorted
        oldest → newest by PEP 440. When multiple coexist, ``best``
        picks the last entry by default; this accessor exposes the
        full set for callers that want to walk every version.

        Unparseable strings (legacy / pre-#261) sort before parseable
        ones — so ``versions[-1]`` still yields the most recent valid
        version when both styles coexist.
        """
        canonical = _resolve_kind(kind)
        versions = {
            p.predictor_version for p in self.predictions
            if p.kind == canonical and p.predictor_name == predictor}
        parsed, fallback = _split_versions(versions)
        return tuple(fallback + parsed)

    def alleles_for(self, kind: str, predictor: str) -> tuple:
        """Alleles a given (kind, predictor) covers. Used by
        coverage / report code that walks per-allele evidence
        without forcing a single 'best' allele."""
        canonical = _resolve_kind(kind)
        return tuple(sorted({
            p.allele for p in self.predictions
            if p.kind == canonical and p.predictor_name == predictor
            and p.allele}))

    # ------------------------------------------------------------------
    # Best-of accessors — kind-named, predictor-aware.
    # ------------------------------------------------------------------

    def best(self, kind: str, *,
             predictor: Optional[str] = None,
             version: Optional[str] = None,
             score_key: Optional[Callable] = None):
        """Best prediction for ``kind`` across alleles. Returns the
        ``mhctools.Prediction`` (carries allele, score, value,
        percentile_rank, …) or ``None`` when no record matches.

        ``kind`` accepts aliases (case-insensitive): ``'ba'`` /
        ``'affinity'`` / ``'binding'`` → ``pMHC_affinity``;
        ``'el'`` / ``'presentation'`` / ``'elution'`` →
        ``pMHC_presentation``; ``'stability'``;
        ``'immunogenicity'``; ``'cleavage'`` / ``'proteasome'`` →
        ``proteasome_cleavage``; ``'processing'`` / ``'ap'`` →
        ``antigen_processing``; ``'tap'``; ``'erap'``. Canonical
        ``pMHC_*`` strings also work.

        Predictor disambiguation:
        - One predictor emitted ``kind`` → ``predictor`` may be
          omitted.
        - Multiple predictors emitted ``kind`` → ``predictor`` is
          required (cross-predictor ranking is meaningless under
          *any* score key, since each predictor's score scale is its
          own). Raises ``ValueError`` otherwise.

        Version disambiguation:
        - One ``predictor_version`` on file → unambiguous.
        - Multiple versions on file → defaults to the most recent
          (PEP 440 ordering). Pass ``version=`` to pin explicitly.
          Versions auto-resolve (vs. predictor's hard-raise) because
          the score scale is comparable across versions of the same
          predictor — just freshness differs.

        Score selection:
        - ``score_key`` is a callable ``Prediction -> float`` (higher
          return value wins). Default is :func:`default_score_key`,
          which uses ``Prediction.score`` (mhctools' normalized,
          higher-is-better axis). Pass a custom key when you want
          per-call ranking on a different axis — e.g.
          ``score_key=lambda p: -p.percentile_rank`` for
          "lowest %-rank wins". Higher-level pipelines use the
          topiary DSL (``EpitopeConfig.score_expr``) on a row-frame
          of predictions; ``best()`` is the per-record primitive.
        """
        canonical = _resolve_kind(kind)
        candidates = [p for p in self.predictions if p.kind == canonical]
        if not candidates:
            return None
        predictor_set = {p.predictor_name for p in candidates}
        if predictor is None:
            if len(predictor_set) > 1:
                raise ValueError(
                    "%s has predictions from multiple predictors "
                    "(%s); pass predictor= to disambiguate. "
                    "Each predictor's score scale is its own — "
                    "cross-predictor ranking is meaningless under "
                    "any score_key."
                    % (canonical, sorted(predictor_set)))
            # Single predictor → unambiguous.
        else:
            candidates = [
                p for p in candidates if p.predictor_name == predictor]
            if not candidates:
                return None
        version_set = {p.predictor_version for p in candidates}
        if version is None:
            if len(version_set) > 1:
                latest = _most_recent_version(version_set)
                candidates = [
                    p for p in candidates if p.predictor_version == latest]
        else:
            candidates = [
                p for p in candidates if p.predictor_version == version]
            if not candidates:
                return None
        if score_key is None:
            score_key = default_score_key
        return max(candidates, key=score_key)

    # Kind-specific helpers. Each forwards to ``best`` with the same
    # disambiguation contract — predictor required when ambiguous,
    # version auto-resolves to most recent, ``score_key`` overridable
    # for per-call ranking customization.
    def best_affinity(self, *, predictor: Optional[str] = None,
                      version: Optional[str] = None,
                      score_key: Optional[Callable] = None):
        return self.best(
            'pMHC_affinity', predictor=predictor, version=version,
            score_key=score_key)

    def best_presentation(self, *, predictor: Optional[str] = None,
                          version: Optional[str] = None,
                          score_key: Optional[Callable] = None):
        return self.best(
            'pMHC_presentation', predictor=predictor, version=version,
            score_key=score_key)

    def best_stability(self, *, predictor: Optional[str] = None,
                       version: Optional[str] = None,
                       score_key: Optional[Callable] = None):
        return self.best(
            'pMHC_stability', predictor=predictor, version=version,
            score_key=score_key)

    def best_immunogenicity(self, *, predictor: Optional[str] = None,
                            version: Optional[str] = None,
                            score_key: Optional[Callable] = None):
        return self.best(
            'immunogenicity', predictor=predictor, version=version,
            score_key=score_key)

    def best_cleavage(self, *, predictor: Optional[str] = None,
                      version: Optional[str] = None,
                      score_key: Optional[Callable] = None):
        return self.best(
            'proteasome_cleavage', predictor=predictor, version=version,
            score_key=score_key)

    def best_antigen_processing(self, *, predictor: Optional[str] = None,
                                version: Optional[str] = None,
                                score_key: Optional[Callable] = None):
        return self.best(
            'antigen_processing', predictor=predictor, version=version,
            score_key=score_key)


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
