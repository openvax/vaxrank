"""Scientific comparisons remain legible and attributable in the report."""

from dataclasses import replace

import pytest
from mhctools.pred import Prediction

from vaxrank.candidate_epitope import CandidateEpitope, Peptide
from vaxrank.processing_prediction import ProcessingPrediction
from vaxrank.report import TemplateDataCreator
from vaxrank.report_display import epitope_display, marked_sequence, model_inventory
from vaxrank.vaccine_antigen import AminoAcidInterval, TargetableMask
from vaxrank.vaccine_peptide import VaccinePeptide

from .test_assembly_antigens import candidate
from .test_product_reports import render, template


def prediction(model='netMHCpan', version='4.2c', kind='pMHC_affinity', value=78.82, score=0.99):
    return Prediction(kind=kind, predictor_name=model, predictor_version=version,
                      peptide='KRFHATISF', allele='HLA-B*27:05', value=value,
                      score=score, percentile_rank=0.004)


def ranked_pair(wt='KRFHATVSF', include_flurry=False):
    antigen, _ = candidate(0, sequence='KRFHATISF', prefix='HG')
    antigen = replace(antigen, targetable_mask=TargetableMask((AminoAcidInterval(8, 9),)))
    predictions = [prediction(), prediction(kind='pMHC_presentation', value=None)]
    if include_flurry:
        predictions += [prediction('mhcflurry', '2.3', value=50.),
                        prediction('mhcflurry', '2.3', 'pMHC_presentation', None, .8),
                        replace(prediction('mhcflurry', '2.3', 'antigen_processing', None, .4), allele='')]
    epitope = CandidateEpitope(sequence='KRFHATISF', source_sequence=antigen.amino_acids,
                              offset=2, predictions=predictions, occurs_in_reference=False,
                              overlaps_targetable=True, per_allele_scores={'HLA-B*27:05': .94},
                              comparators={'wt': Peptide(sequence=wt, predictions=[prediction(value=89.19)])})
    return [(antigen, [VaccinePeptide(antigen=antigen, epitopes=[epitope])])]


def test_full_report_leads_with_evidence_and_attributes_scores_and_comparator():
    data = template(ranked_pair(), candidate_details=True)
    html = render(data, 'html')
    assert html.index('id="candidate-details"') < html.index('id="products"') < html.index('id="run-details"')
    assert 'KRFHAT<span class="mutant">I</span>SF' in html
    assert 'KRFHAT<span class="reference-change">V</span>SF' in html
    assert 'NetMHCpan Score_EL: <strong>0.990</strong>' in html
    assert 'Presentation rank: 0.004%' in html
    assert 'Version 4.2c' in html
    assert 'Target IC50: <strong>78.82 nM</strong>' in html
    assert 'WT IC50: 89.19 nM' in html
    assert '1 candidate window(s)' not in html
    assert 'No results in the reported candidates.' in html
    assert 'WT does not establish an observed patient-normal sequence.' in html
    text = render(data, 'txt')
    assert text.index('=== Candidate window details') < text.index('Emitted products')
    assert 'NetMHCpan: Score_EL' in text


def test_flurry_processing_and_pepsickle_are_separate_recorded_evidence():
    ranked = ranked_pair(include_flurry=True)
    epitope = ranked[0][1][0].epitopes[0]
    pp = ProcessingPrediction(epitope.sequence, epitope.source_sequence, 'pepsickle', 2,
                              predictor_version='1.0', c_term_cleavage_prob=.8,
                              max_internal_cut_prob=.2, processing_score=.8)
    html = render(template(ranked, candidate_details=True,
                           processing_predictions_by_key={pp.key(): pp}), 'html')
    assert 'mhcflurry presentation: <strong>0.800</strong>' in html
    assert 'mhcflurry processing: 0.400' in html
    assert 'Pepsickle C-term: 0.80' in html
    assert 'Pepsickle max internal: 0.20' in html
    assert 'Cleavage composite: 0.80' in html
    assert 'Recorded cleavage annotations' in html


@pytest.mark.parametrize('wt,note', [('', 'Sequence unavailable'), ('KRFHATVS', 'Different lengths; no residue alignment inferred.')])
def test_missing_and_unequal_comparators_are_explicit(wt, note):
    html = render(template(ranked_pair(wt), candidate_details=True), 'html')
    assert note in html
    assert 'class="reference-change">' not in html.split('<table class="pure-table target-epitopes">', 1)[1].split('</table>', 1)[0]


def test_repeated_sequence_does_not_guess_a_target_location():
    ranked = ranked_pair()
    antigen, peptides = ranked[0]
    epitope = replace(peptides[0].epitopes[0], source_sequence='other source', offset=0)
    display = epitope_display(epitope, antigen, {}, prediction())
    assert not display['location_known']
    assert str(display['target']) == epitope.sequence


def test_disjoint_target_intervals_and_junction_are_preserved_and_escaped():
    assert str(marked_sequence('ACDEFG', [(1, 2), (4, 5)])) == (
        'A<span class="mutant">C</span>DE<span class="mutant">F</span>G')
    text = str(marked_sequence('ACDE', [(2, 2)]))
    assert 'AC<span class="junction" title="Target junction">|</span>DE' == text
    assert '<script>' not in str(marked_sequence('<script>'))


def test_disabled_and_unavailable_cleavage_are_distinct():
    assert model_inventory([], {}, False)[-1]['status'].startswith('Disabled for this run')
    assert model_inventory([], {}, None)[-1]['status'].startswith('No cleavage results recorded')


def test_wt_affinity_matches_predictor_version_without_fallback():
    creator = TemplateDataCreator.__new__(TemplateDataCreator)
    creator.processing_predictions_by_key = {}
    epitope = ranked_pair()[0][1][0].epitopes[0]
    epitope = replace(epitope, comparators={'wt': Peptide(sequence='KRFHATVSF', predictions=[
        prediction(version='4.1', value=200.), prediction(version='4.2c', value=89.19)])})
    assert creator.epitope_data(epitope, prediction())['WT IC50'] == '89.19 nM'
    assert creator.epitope_data(epitope, prediction(version='4.0'))['WT IC50'] == 'No prediction'


def test_product_highlighting_uses_consistent_mapped_masks_without_changing_products():
    from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
    from vaxrank.report_display import product_target_intervals

    ranked = ranked_pair()
    products = assemble_peptide_constructs(ranked, PeptideConstructConfig())
    before = [(p.sequence, p.antigens) for p in products]
    data = template(ranked, {'peptide': products})
    view = data['product_groups'][0]['products'][0]
    assert 'class="mutant">I</span>' in view['sequence_html']
    assert view['display_name'].startswith('GENE0')
    assert view['antigens'][0]['genes'] == ''
    assert [(p.sequence, p.antigens) for p in products] == before
    # Same source sequence with differing target masks cannot acquire a union highlight.
    antigen, peptides = ranked[0]
    alternative = replace(antigen, targetable_mask=TargetableMask((AminoAcidInterval(2, 3),)))
    ranked.append((alternative, [VaccinePeptide(antigen=alternative, epitopes=peptides[0].epitopes)]))
    masks = product_target_intervals(products[0], ranked)
    assert list(masks.values()) == [None]


@pytest.mark.parametrize('interval', [(1, 2), (11, 12)])
def test_adjacent_target_residue_does_not_become_a_deletion_junction(interval):
    antigen, peptides = ranked_pair()[0]
    antigen = replace(antigen, targetable_mask=TargetableMask((AminoAcidInterval(*interval),)))
    display = epitope_display(peptides[0].epitopes[0], antigen, {}, prediction())
    assert str(display['target']) == 'KRFHATISF'


def test_product_with_an_unavailable_source_mask_is_not_partially_highlighted():
    from copy import deepcopy
    from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
    from vaxrank.report_display import product_target_intervals

    ranked = ranked_pair()
    product, = assemble_peptide_constructs(ranked, PeptideConstructConfig())
    antigens = deepcopy(product.antigens)
    missing = dict(antigens[0]['sources'][0], source_identifier='unavailable-source')
    antigens[0]['sources'].append(missing)
    product = replace(product, antigens=antigens)
    assert list(product_target_intervals(product, ranked).values()) == [None]
