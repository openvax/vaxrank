# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Tests for ``vaxrank.coverage`` (#269): MHC allele coverage
measurement + greedy coverage-aware antigen selection.

Pins:
- Tier rule: best-of-two-axes (presentation_percentile, percentile_rank)
  bucketed at 0.5 / 1.0 / 2.0 % cutoffs into strong / medium / low.
- Animal-agnostic: same code path for HLA, mouse H-2, swine SLA — only
  thing that varies is which predictor has data on which allele.
- Greedy selector: lexicographic (strong, medium, low, score). Falls
  back to pure score order when no candidate adds coverage.
- Multi-predictor evidence aggregates: presentation %-rank from one
  record + affinity %-rank from another for the same (peptide, allele)
  combine (smaller of the two wins).
"""

from types import SimpleNamespace

from mhctools.pred import Prediction

from vaxrank.coverage import (
    AlleleCoverage,
    antigen_tier_per_allele,
    compute_coverage,
    select_antigens_for_coverage,
    summarize_construction_decisions,
)
from vaxrank.peptide_context import CandidateEpitope, Peptide


def _ep(peptide, allele, *, presentation_pct=None, affinity_pct=None,
        method=''):
    """One ``CandidateEpitope`` carrying one or two leaf ``Prediction`` records —
    a presentation-rank one (``pMHC_presentation``) and/or an affinity-
    rank one (``pMHC_affinity``). Matches the multi-kind shape
    ``compute_coverage`` operates on."""
    preds = []
    if presentation_pct is not None:
        preds.append(Prediction(
            kind='pMHC_presentation', predictor_name=method or 'pres',
            predictor_version='', allele=allele, peptide=peptide,
            value=None, score=0.0, percentile_rank=presentation_pct))
    if affinity_pct is not None:
        preds.append(Prediction(
            kind='pMHC_affinity', predictor_name=method or 'aff',
            predictor_version='', allele=allele, peptide=peptide,
            value=None, score=0.0, percentile_rank=affinity_pct))
    return CandidateEpitope(
        mutant=Peptide(
            sequence=peptide, predictions=tuple(preds)),
        overlaps_mutation=True, occurs_in_reference=False)


def _vp(epitopes, *, combined_score=1.0, gene_name='GENE'):
    """Minimal VaccinePeptide-like record (CandidateEpitope shape)."""
    return SimpleNamespace(
        mutant_epitopes=epitopes,
        mutant_protein_fragment=SimpleNamespace(
            gene_name=gene_name, amino_acids='AAAAAAA'),
        combined_score=combined_score,
        manufacturability_scores=None)


# -- compute_coverage --------------------------------------------------


def test_compute_coverage_buckets_into_correct_tiers():
    """0.4% → strong; 0.7% → medium; 1.5% → low; 3.0% → uncovered."""
    preds = [
        _ep('SIINFEKL', 'HLA-A*02:01', presentation_pct=0.4),
        _ep('AAAAAAAAA', 'HLA-B*07:02', presentation_pct=0.7),
        _ep('CCCCCCCCC', 'HLA-C*03:04', presentation_pct=1.5),
        _ep('DDDDDDDDD', 'HLA-A*68:01', presentation_pct=3.0),
    ]
    targets = ['HLA-A*02:01', 'HLA-B*07:02', 'HLA-C*03:04', 'HLA-A*68:01']
    coverage = compute_coverage(preds, targets)
    assert coverage['HLA-A*02:01'].tier == 'strong'
    assert coverage['HLA-B*07:02'].tier == 'medium'
    assert coverage['HLA-C*03:04'].tier == 'low'
    assert coverage['HLA-A*68:01'].tier is None  # 3% > 2.0 = uncovered


def test_compute_coverage_aggregates_across_predictors():
    """When mhcflurry says presentation_pct=1.5 and netmhcpan says
    affinity_pct=0.3 for the same (peptide, allele), the better axis
    wins and the allele lands in tier-strong."""
    preds = [
        _ep('SIINFEKL', 'HLA-A*02:01', presentation_pct=1.5,
            method='mhcflurry'),
        _ep('SIINFEKL', 'HLA-A*02:01', affinity_pct=0.3,
            method='netmhcpan'),
    ]
    coverage = compute_coverage(preds, ['HLA-A*02:01'])
    assert coverage['HLA-A*02:01'].tier == 'strong'


def test_compute_coverage_uncovered_alleles_appear_as_none():
    """Every target allele gets a row, even those with no
    predictions. Drives the report's ``✗ uncovered`` line."""
    preds = [_ep('SIINFEKL', 'HLA-A*02:01', presentation_pct=0.4)]
    coverage = compute_coverage(preds, ['HLA-A*02:01', 'HLA-B*07:02'])
    assert coverage['HLA-A*02:01'].tier == 'strong'
    assert coverage['HLA-B*07:02'].tier is None
    assert coverage['HLA-B*07:02'].n_peptides == 0


def test_compute_coverage_animal_agnostic_mouse_alleles():
    """Mouse H-2 alleles flow through the same code. mhcflurry
    typically has no mouse coverage so only affinity % rank is
    available; the model uses what it's given."""
    preds = [
        _ep('SIINFEKL', 'H-2-Kb', affinity_pct=0.3, method='netmhcpan'),
        _ep('AAAAAAAA', 'H-2-Db', affinity_pct=1.5, method='netmhcpan'),
    ]
    coverage = compute_coverage(preds, ['H-2-Kb', 'H-2-Db'])
    assert coverage['H-2-Kb'].tier == 'strong'
    assert coverage['H-2-Db'].tier == 'low'


def test_compute_coverage_skips_alleles_outside_target_set():
    """Alleles outside the patient's set are ignored — coverage is
    a per-patient view, not a directory of every allele the
    predictor saw."""
    preds = [
        _ep('SIINFEKL', 'HLA-A*02:01', presentation_pct=0.4),
        _ep('SIINFEKL', 'HLA-X*99:99', presentation_pct=0.1),
    ]
    coverage = compute_coverage(preds, ['HLA-A*02:01'])
    assert set(coverage) == {'HLA-A*02:01'}


def test_compute_coverage_empty_targets_returns_empty():
    """No targets → no coverage rows. Lets callers gate the report
    section on ``if patient_info.mhc_alleles``."""
    coverage = compute_coverage(
        [_ep('SIINFEKL', 'HLA-A*02:01', presentation_pct=0.4)],
        target_alleles=[])
    assert coverage == {}


# -- antigen_tier_per_allele ------------------------------------------


def test_antigen_tier_per_allele_takes_best_peptide():
    """One antigen has many sliding-window peptides; the antigen's
    tier for an allele is the *best* tier any of its peptides
    achieves."""
    vp = _vp([
        _ep('AAAAAAAA', 'HLA-A*02:01', presentation_pct=1.8),  # low
        _ep('SIINFEKL', 'HLA-A*02:01', presentation_pct=0.3),  # strong
        _ep('DDDDDDDD', 'HLA-A*02:01', presentation_pct=0.9),  # medium
    ])
    tiers = antigen_tier_per_allele(vp, ['HLA-A*02:01'])
    assert tiers['HLA-A*02:01'] == 'strong'


# -- select_antigens_for_coverage -------------------------------------


def test_selector_prefers_coverage_over_score():
    """When the top-scoring antigen covers no targets and a
    lower-scoring one covers a strong-tier target, the selector
    picks the lower-scoring one first."""
    high_score_uncovered = _vp(
        [_ep('AAAAAAAA', 'HLA-A*02:01', presentation_pct=3.0)],
        combined_score=10.0, gene_name='HIGH')
    low_score_strong = _vp(
        [_ep('SIINFEKL', 'HLA-A*02:01', presentation_pct=0.3)],
        combined_score=1.0, gene_name='LOW')
    candidates = [
        ('var1', [high_score_uncovered]),
        ('var2', [low_score_strong]),
    ]
    selected = select_antigens_for_coverage(
        candidates, ['HLA-A*02:01'], n_to_select=2)
    # Coverage-driving antigen comes first.
    assert selected[0][0] == 'var2'
    assert selected[1][0] == 'var1'


def test_selector_falls_back_to_score_when_coverage_satisfied():
    """Once every target allele is covered at strong tier, the
    selector reverts to pure score order — so an already-covered
    pool gives the same selection as today's behavior."""
    a = _vp(
        [_ep('SIINFEKL', 'HLA-A*02:01', presentation_pct=0.3)],
        combined_score=10.0, gene_name='A')
    b = _vp(
        [_ep('AAAAAAAA', 'HLA-A*02:01', presentation_pct=0.4)],
        combined_score=8.0, gene_name='B')
    c = _vp(
        [_ep('DDDDDDDD', 'HLA-A*02:01', presentation_pct=0.5)],
        combined_score=12.0, gene_name='C')
    candidates = [('vA', [a]), ('vB', [b]), ('vC', [c])]
    selected = select_antigens_for_coverage(
        candidates, ['HLA-A*02:01'], n_to_select=3)
    # First pick: any of the three covers strong; selector picks
    # the highest-scoring one (C). Then coverage is satisfied —
    # remaining slots fill in pure score order: A (10), B (8).
    names = [s[0] for s in selected]
    assert names[0] == 'vC'
    # Ordering after coverage is satisfied should mirror
    # combined_score among the rest.
    assert names == ['vC', 'vA', 'vB']


def test_selector_no_op_with_empty_targets():
    """Empty target_alleles → return ranked[:n] unchanged."""
    a = _vp([_ep('SIINFEKL', 'HLA-A*02:01', presentation_pct=0.3)],
            combined_score=10.0)
    b = _vp([_ep('AAAAAAAA', 'HLA-B*07:02', presentation_pct=0.4)],
            combined_score=8.0)
    candidates = [('vA', [a]), ('vB', [b])]
    selected = select_antigens_for_coverage(
        candidates, target_alleles=[], n_to_select=2)
    assert [s[0] for s in selected] == ['vA', 'vB']


def test_selector_spreads_across_multiple_alleles():
    """With three target alleles and three antigens each covering a
    different one, the selector picks all three (one per allele)
    rather than three of the same allele."""
    a02 = _vp(
        [_ep('SIINFEKL', 'HLA-A*02:01', presentation_pct=0.2)],
        combined_score=5.0, gene_name='A02')
    b07 = _vp(
        [_ep('AAAAAAAA', 'HLA-B*07:02', presentation_pct=0.2)],
        combined_score=5.0, gene_name='B07')
    c03 = _vp(
        [_ep('CCCCCCCC', 'HLA-C*03:04', presentation_pct=0.2)],
        combined_score=5.0, gene_name='C03')
    candidates = [('vA', [a02]), ('vB', [b07]), ('vC', [c03])]
    selected = select_antigens_for_coverage(
        candidates, ['HLA-A*02:01', 'HLA-B*07:02', 'HLA-C*03:04'],
        n_to_select=3)
    # All three picked — full spread.
    assert {s[0] for s in selected} == {'vA', 'vB', 'vC'}


# -- summarize_construction_decisions ---------------------------------


def test_summarize_includes_coverage_for_selected_pool():
    """The summary's ``coverage`` field is computed over the
    *selected* antigens, not the full ranked input."""
    selected_a = _vp(
        [_ep('SIINFEKL', 'HLA-A*02:01', presentation_pct=0.3)],
        combined_score=10.0, gene_name='SELECTED')
    dropped_b = _vp(
        [_ep('DDDDDDDD', 'HLA-B*07:02', presentation_pct=0.3)],
        combined_score=1.0, gene_name='DROPPED')
    ranked = [
        ('var1', [selected_a]),
        ('var2', [dropped_b]),
    ]
    summary = summarize_construction_decisions(
        ranked, cap=1, target_alleles=['HLA-A*02:01', 'HLA-B*07:02'])
    assert summary['cap'] == 1
    assert summary['total_ranked'] == 2
    assert len(summary['selected']) == 1
    assert len(summary['dropped']) == 1
    coverage_by_allele = {c.allele: c for c in summary['coverage']}
    assert coverage_by_allele['HLA-A*02:01'].tier == 'strong'
    # B*07 is only available on the dropped antigen, so the
    # selected-pool coverage is uncovered for B*07.
    assert coverage_by_allele['HLA-B*07:02'].tier is None


def test_ascii_report_renders_coverage_block_in_construction_section():
    """End-to-end render: a populated ``vaccine_constructions``
    block produces the new coverage section in the ASCII output —
    section header, per-allele coverage lines, per-antigen
    ``covers A (tier)`` annotations."""
    import io
    from vaxrank.report import _make_report

    template_data = {
        'patient_info': {'Patient ID': 'TEST'},
        'package_versions': {},
        'reviewers': [], 'final_review': '',
        'input_json_file': None,
        'include_manufacturability': False,
        'include_wt_epitopes': False,
        'args': [], 'variants': [],
        'mrna_ranking': None,
        'vaccine_constructions': {
            'mrna': {
                'antigens_per_construct': 5, 'max_constructs': 1,
                'cap': 5, 'total_ranked': 5,
                'selected': [
                    {'rank': 1, 'gene_name': 'GENE_A',
                     'description': 'GENE_A_chr1_100',
                     'combined_score': 10.0,
                     'allele_tiers': {'HLA-A*02:01': 'strong'},
                     'manufacturability': {}},
                ],
                'dropped': [],
                'coverage': [
                    AlleleCoverage(
                        allele='HLA-A*02:01', tier='strong',
                        n_peptides=2,
                        best_presentation_pct=0.3,
                        best_affinity_pct=0.5,
                        contributing_peptides=('SIINFEKL',)),
                    AlleleCoverage(
                        allele='HLA-B*07:02', tier=None,
                        n_peptides=0,
                        best_presentation_pct=None,
                        best_affinity_pct=None,
                        contributing_peptides=()),
                ],
            },
        },
    }
    f = io.StringIO()
    _make_report(template_data, f, 'templates/template.txt')
    rendered = f.getvalue()
    # Section header.
    assert 'Vaccine construction — mrna' in rendered
    # Per-allele coverage lines: one covered, one uncovered.
    assert 'HLA-A*02:01' in rendered
    assert 'HLA-B*07:02' in rendered
    assert 'strong' in rendered
    assert 'none' in rendered
    # Per-antigen tier annotation.
    assert 'covers HLA-A*02:01 (strong)' in rendered
    # Best presentation %-rank surfaces in the coverage detail line.
    assert '0.30' in rendered


def test_summarize_per_row_has_allele_tiers_and_manufacturability():
    """Each row in ``selected``/``dropped`` carries the per-allele
    tier map (drives the ``covers A*02:01 (strong)`` annotation)
    plus the manufacturability dict (peptide construction view
    renders it)."""
    vp = _vp(
        [_ep('SIINFEKL', 'HLA-A*02:01', presentation_pct=0.3)],
        combined_score=5.0, gene_name='GENEX')
    # Stub the manufacturability_scores like ManufacturabilityScores
    # would (a namedtuple with _asdict).
    from collections import namedtuple
    Mfg = namedtuple('Mfg', ['cysteine_count', 'cterm_hydropathy'])
    vp.manufacturability_scores = Mfg(0, 0.5)

    summary = summarize_construction_decisions(
        [('var1', [vp])], cap=1, target_alleles=['HLA-A*02:01'])
    row = summary['selected'][0]
    assert row['gene_name'] == 'GENEX'
    assert row['allele_tiers'] == {'HLA-A*02:01': 'strong'}
    assert row['manufacturability'] == {
        'cysteine_count': 0, 'cterm_hydropathy': 0.5}
