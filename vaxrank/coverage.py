# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""MHC allele coverage measurement + greedy selection.

A vaccine's job is to engage a set of MHC alleles. Today vaxrank picks
top-k antigens by ``combined_score`` only — that can leave alleles
uncovered while others are over-engaged. This module adds the missing
piece: a measurement of which alleles each antigen covers, and a
greedy selector that prefers antigens which add the most uncovered
alleles before falling back to ``combined_score`` order.

## Two-axis evidence

For each (peptide, allele) pair we look at two complementary axes:

  presentation_percentile  — mhcflurry-pres: probability the peptide
                             is generated + presented at meaningful
                             abundance.
  percentile_rank          — netmhcpan-affinity: binding strength once
                             the peptide is in the ER.

A peptide-allele pair is "covered" at a given tier when the **better**
of the two axes (the smaller %-rank we have data for) is at-or-below
the tier threshold.

## Tiers

  strong  : %-rank ≤ 0.5  — top-priority for coverage spreading
  medium  : %-rank ≤ 1.0  — the selection objective
  low     : %-rank ≤ 2.0  — bonus signal, not driving selection
  none    : %-rank > 2.0  — not considered covered

## Animal-agnostic

Driven entirely off the patient's allele set
(``patient_info.mhc_alleles``). HLA, mouse H-2, swine SLA, … all work
the same way; predictors that don't have data for a given allele
simply don't contribute. No species-specific code anywhere.

## Selector

Greedy, lexicographic key per candidate antigen:

  (strong_alleles_added, medium_alleles_added, low_alleles_added,
   combined_score)

Pick the next antigen that maximizes new strong-tier coverage; ties
broken by medium, then low, then today's ``combined_score``. Once
every target allele is covered at strong tier, the selector
degenerates to pure ``combined_score`` order — so when coverage is
already complete, the rebalancer is a no-op.

Issue: openvax/vaxrank#269.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


# Default tier thresholds (% rank). Stable, hardcoded for now —
# expose as CLI knobs later if a user actually needs to tune them.
TIER_STRONG = 0.5
TIER_MEDIUM = 1.0
TIER_LOW = 2.0

# Names ordered most → least stringent. Shared between the report
# writers and the selector.
TIER_NAMES = ('strong', 'medium', 'low')


def _tier_for_percentile(pct: Optional[float]) -> Optional[str]:
    """Bucket one peptide-allele %-rank into a tier name, or
    ``None`` if it doesn't reach any tier."""
    if pct is None:
        return None
    if pct <= TIER_STRONG:
        return 'strong'
    if pct <= TIER_MEDIUM:
        return 'medium'
    if pct <= TIER_LOW:
        return 'low'
    return None


def _better_tier(a: Optional[str], b: Optional[str]) -> Optional[str]:
    """Return the more stringent of two tier names. ``None`` is the
    weakest. Order is: strong > medium > low > None."""
    rank = {'strong': 3, 'medium': 2, 'low': 1, None: 0}
    return a if rank[a] >= rank[b] else b


def _peptide_tier(p) -> Optional[str]:
    """Best tier achieved by one EpitopePrediction, taking the
    smaller of presentation_percentile and percentile_rank when
    both are available. Predictors that don't expose a given axis
    don't contribute (None on that axis is ignored)."""
    pres = getattr(p, 'presentation_percentile', None)
    aff = getattr(p, 'percentile_rank', None)
    return _better_tier(_tier_for_percentile(pres), _tier_for_percentile(aff))


@dataclass(frozen=True)
class AlleleCoverage:
    """Per-allele binding-evidence summary at one selection scope.

    Attributes
    ----------
    allele : str
        Allele name as it appears in the input (mhcgnomes-normalized
        upstream — e.g. ``'HLA-A*02:01'``, ``'H-2-Kb'``,
        ``'SLA-1*04:01'``).
    tier : Optional[str]
        Best tier achieved across all evidence in this scope.
        ``'strong' | 'medium' | 'low' | None``. ``None`` means
        uncovered at the tier-3 (low) cutoff.
    n_peptides : int
        Number of distinct peptide sequences in this scope that
        cleared the tier-3 (low) cutoff for this allele.
    best_presentation_pct : Optional[float]
        Smallest presentation percentile_rank seen for this allele,
        or ``None`` when no presentation predictor has data.
    best_affinity_pct : Optional[float]
        Smallest affinity percentile_rank seen for this allele,
        or ``None``.
    contributing_peptides : tuple[str, ...]
        Peptide sequences that hit at least the ``low`` tier
        (capped at 5 entries — full list lives in DEBUG logs).
    """

    allele: str
    tier: Optional[str]
    n_peptides: int
    best_presentation_pct: Optional[float]
    best_affinity_pct: Optional[float]
    contributing_peptides: tuple


def _min(a: Optional[float], b: Optional[float]) -> Optional[float]:
    """min() that treats None as 'no data' rather than raising."""
    if a is None:
        return b
    if b is None:
        return a
    return a if a < b else b


def compute_coverage(
    epitope_predictions,
    target_alleles,
) -> dict[str, AlleleCoverage]:
    """Build {allele: AlleleCoverage} for ``target_alleles`` against
    the evidence in ``epitope_predictions``.

    Every allele in ``target_alleles`` gets a row, even alleles with
    no predictions (tier=None, n=0). Alleles outside the target set
    are skipped — this is a coverage view, not a directory of every
    allele the predictor saw.

    ``epitope_predictions`` is the flat list shape that vaxrank uses
    everywhere post-#261 — one record per (peptide, allele,
    predictor). When multi-predictor data is in play (e.g.
    mhcflurry-presentation + netmhcpan-affinity for the same
    (peptide, allele)), the predictions land here as separate
    objects and ``compute_coverage`` aggregates evidence across
    them.
    """
    targets = set(target_alleles or [])
    if not targets:
        return {}
    # Two passes: collect best %-ranks per (peptide, allele), then
    # bucket into tiers. The first pass merges multi-predictor
    # evidence — for one (peptide, allele) we may see two records,
    # one with presentation and one with affinity.
    by_pep_allele: dict[tuple, dict] = {}
    for p in epitope_predictions:
        allele = getattr(p, 'allele', None)
        if allele not in targets:
            continue
        peptide = getattr(p, 'peptide_sequence', None) or ''
        if not peptide:
            continue
        key = (peptide, allele)
        slot = by_pep_allele.setdefault(
            key, {'pres': None, 'aff': None})
        slot['pres'] = _min(
            slot['pres'], getattr(p, 'presentation_percentile', None))
        slot['aff'] = _min(
            slot['aff'], getattr(p, 'percentile_rank', None))

    # Aggregate to per-allele coverage.
    out: dict[str, AlleleCoverage] = {}
    for allele in targets:
        peptides = []
        best_pres = None
        best_aff = None
        best_tier: Optional[str] = None
        for (peptide, a), slot in by_pep_allele.items():
            if a != allele:
                continue
            best_pres = _min(best_pres, slot['pres'])
            best_aff = _min(best_aff, slot['aff'])
            tier = _better_tier(
                _tier_for_percentile(slot['pres']),
                _tier_for_percentile(slot['aff']))
            if tier is not None:
                peptides.append(peptide)
                best_tier = _better_tier(best_tier, tier)
        out[allele] = AlleleCoverage(
            allele=allele,
            tier=best_tier,
            n_peptides=len(peptides),
            best_presentation_pct=best_pres,
            best_affinity_pct=best_aff,
            contributing_peptides=tuple(sorted(set(peptides))[:5]),
        )
    return out


def _antigen_predictions(vaccine_peptide):
    """All EpitopePredictions on this VaccinePeptide. Defensive:
    tolerates duck-typed test fixtures that don't carry the field."""
    return list(getattr(vaccine_peptide, 'mutant_epitope_predictions', []) or [])


def antigen_tier_per_allele(
    vaccine_peptide, target_alleles,
) -> dict[str, Optional[str]]:
    """Return ``{allele: tier}`` for one antigen — the best tier
    each target allele reaches via this antigen's epitope
    predictions. ``None`` for alleles the antigen doesn't cover at
    any tier."""
    coverage = compute_coverage(
        _antigen_predictions(vaccine_peptide), target_alleles)
    return {a: coverage[a].tier for a in target_alleles}


def summarize_construction_decisions(
    ranked_variants_with_vaccine_peptides,
    *,
    cap,
    target_alleles,
    selected=None,
):
    """Modality-agnostic "what landed in the vaccine?" summary used
    by the report writers (#269 + #270).

    Parameters
    ----------
    ranked_variants_with_vaccine_peptides : list[tuple]
        ``[(Variant, [VaccinePeptide])]`` — the score-ranked input.
    cap : int
        Number of antigens this modality can hold
        (``antigens_per_construct × max_constructs`` for both
        peptide and mRNA today; future modalities supply their own).
    target_alleles : list[str]
        Patient MHC alleles for the coverage view. Empty disables
        the per-allele coverage block.
    selected : optional list[tuple]
        The antigen list actually picked by the assembler (post
        coverage rebalancing). Defaults to ``ranked[:cap]`` for
        callers that don't run the rebalancer.

    Returns
    -------
    dict
        ``{
            'cap': int,
            'total_ranked': int,
            'selected': [{...}, ...],     # rank, gene_name, description,
                                          # combined_score, allele_tiers
            'dropped': [...],             # same row shape
            'coverage': [AlleleCoverage, ...],  # one per target_allele
          }``
    """
    if selected is None:
        selected = list(ranked_variants_with_vaccine_peptides[:cap])

    selected_ids = {id(item) for item in selected}
    rank_by_id = {
        id(item): i + 1
        for i, item in enumerate(ranked_variants_with_vaccine_peptides)
    }

    def _row(item):
        variant, peptides = item
        vp = peptides[0] if peptides else None
        gene_name = (
            getattr(vp.mutant_protein_fragment, 'gene_name', '')
            if vp is not None else '')
        rank = rank_by_id.get(id(item), 0)
        # Per-allele tier for this antigen — drives the
        # ``covers A*02:01 (strong)`` annotation in the report.
        allele_tiers = (
            antigen_tier_per_allele(vp, target_alleles)
            if (vp is not None and target_alleles) else {})
        # Manufacturability tuple for the per-modality block (used
        # by the peptide construction view; mRNA renders ignore it).
        # Structured as a flat dict so the template can render
        # without knowing the namedtuple's field names.
        scores = getattr(vp, 'manufacturability_scores', None) if vp else None
        manufacturability = (
            dict(scores._asdict()) if scores is not None
            and hasattr(scores, '_asdict') else {})
        return {
            'rank': rank,
            'gene_name': gene_name or '',
            'description': '%s_%s' % (
                gene_name or '?',
                getattr(variant, 'short_description', str(variant))),
            'combined_score': (
                float(vp.combined_score) if vp is not None else 0.0),
            'allele_tiers': {
                a: tier for a, tier in allele_tiers.items()
                if tier is not None
            },
            'manufacturability': manufacturability,
        }

    selected_rows = [_row(item) for item in selected]
    dropped_rows = [
        _row(item) for item in ranked_variants_with_vaccine_peptides
        if id(item) not in selected_ids
    ]

    # Coverage of the SELECTED pool: aggregate evidence across
    # every selected antigen's mutant-epitope predictions.
    coverage_records: list = []
    if target_alleles:
        all_predictions = []
        for variant, peptides in selected:
            for vp in (peptides or [])[:1]:
                all_predictions.extend(
                    getattr(vp, 'mutant_epitope_predictions', None) or [])
        coverage = compute_coverage(all_predictions, target_alleles)
        # Sort: covered (strong > medium > low) before uncovered;
        # within same tier, alphabetical.
        tier_rank = {'strong': 0, 'medium': 1, 'low': 2, None: 3}
        coverage_records = sorted(
            coverage.values(),
            key=lambda c: (tier_rank.get(c.tier, 3), c.allele))

    return {
        'cap': cap,
        'total_ranked': len(ranked_variants_with_vaccine_peptides),
        'selected': selected_rows,
        'dropped': dropped_rows,
        'coverage': coverage_records,
    }


def select_antigens_for_coverage(
    candidates,
    target_alleles,
    n_to_select,
    *,
    score_attr: str = 'combined_score',
):
    """Greedy coverage-aware selector.

    Picks ``n_to_select`` antigens from ``candidates`` (the
    ranked-by-score list of ``(variant, [VaccinePeptide])`` tuples)
    using the lexicographic key:

      (strong_alleles_added, medium_alleles_added, low_alleles_added,
       combined_score)

    Stops either when ``n_to_select`` is hit or when no remaining
    candidate adds new coverage at any tier — at which point the
    fallback runs and the rest of the slots are filled in
    ``combined_score`` order.

    No-op for empty ``target_alleles`` (returns ``candidates[:n_to_select]``
    unchanged) — the function is safe to call unconditionally.

    Returns the same shape as the input: a list of
    ``(variant, [VaccinePeptide])`` tuples in selection order.
    """
    if n_to_select <= 0:
        return []
    if not candidates:
        return []
    if not target_alleles:
        return list(candidates[:n_to_select])
    # Targets we'd like covered. We track tier achieved for each.
    target = list(target_alleles)
    achieved: dict[str, Optional[str]] = {a: None for a in target}

    def _new_tier_counts(antigen_tiers: dict[str, Optional[str]]):
        """How many *new* alleles this antigen would cover at each
        tier, given current ``achieved`` state."""
        new_strong = new_medium = new_low = 0
        for a, tier in antigen_tiers.items():
            cur = achieved[a]
            if tier == 'strong' and cur != 'strong':
                new_strong += 1
            elif tier == 'medium' and cur not in ('strong', 'medium'):
                new_medium += 1
            elif tier == 'low' and cur is None:
                new_low += 1
        return new_strong, new_medium, new_low

    # Pre-compute per-candidate (antigen) tier maps once.
    pool = list(candidates)
    tier_maps = []
    for variant, peptides in pool:
        # Antigen-level tier = best tier achieved by any peptide on
        # this variant; we use the top vaccine peptide as the
        # representative since constructs ship one-per-variant.
        if not peptides:
            tier_maps.append({a: None for a in target})
            continue
        tier_maps.append(antigen_tier_per_allele(peptides[0], target))

    selected: list = []
    selected_idx = set()
    while len(selected) < n_to_select and len(selected_idx) < len(pool):
        best_i = None
        best_key = None
        # Score ties broken by original (combined_score) rank — pool
        # was already sorted by score, so iterating in index order
        # gives us that for free.
        for i, (variant, peptides) in enumerate(pool):
            if i in selected_idx:
                continue
            new_strong, new_medium, new_low = _new_tier_counts(tier_maps[i])
            score = (
                float(getattr(peptides[0], score_attr, 0.0))
                if peptides else 0.0)
            # Lexicographic: max strong, then medium, then low,
            # then score (larger wins).
            key = (new_strong, new_medium, new_low, score, -i)
            if best_key is None or key > best_key:
                best_key = key
                best_i = i
        # If the best candidate adds no new coverage, the greedy
        # stage is done — fill remaining slots in score order.
        new_strong, new_medium, new_low = best_key[:3]
        if new_strong == 0 and new_medium == 0 and new_low == 0:
            break
        selected.append(pool[best_i])
        selected_idx.add(best_i)
        # Bookkeeping: update achieved tiers from this antigen.
        for a, tier in tier_maps[best_i].items():
            achieved[a] = _better_tier(achieved[a], tier)
    # Score-order fallback for the remaining slots.
    if len(selected) < n_to_select:
        for i, item in enumerate(pool):
            if i in selected_idx:
                continue
            selected.append(item)
            selected_idx.add(i)
            if len(selected) >= n_to_select:
                break
    return selected
