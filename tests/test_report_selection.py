"""Reports explain saved decisions without changing products or inventing history."""

from copy import deepcopy
from dataclasses import replace

import pytest

from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs
from vaxrank.native_serialization import to_native_json
from vaxrank.patient_info import PatientInfo
from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
from vaxrank.report import TemplateDataCreator
from vaxrank.report_selection import selection_policy_view, window_view
from vaxrank.window_selection import WindowSelection, select_windows

from .test_assembly_antigens import candidate
from .test_product_reports import render, template
from .test_report_grouping import Links
from .test_window_selection import vaccine


@pytest.mark.parametrize('extension', ['html', 'txt'])
def test_effective_cap_wins_over_invocation_defaults_and_reporting_is_read_only(extension):
    ranked = [candidate(0)]
    products = assemble_peptide_constructs(ranked, PeptideConstructConfig(max_constructs=10))
    selection = {'run_configuration': {'effective_constructs': {'peptide': {'max_constructs': 20}}},
                 'construct_assembly': {'peptide': {'configuration': {'max_constructs': 10}}}}
    before = to_native_json((selection, products, [p.antigen for _, ps in ranked for p in ps]))
    data = TemplateDataCreator(ranked, PatientInfo('saved-cap'), final_review='', reviewers='',
        args_for_report={'manufacturability': False, 'wt_epitopes': True, 'peptide_max_constructs': 20},
        input_json_file=None, selection=selection, products={'peptide': products}).compute_template_data()
    assert ('Maximum constructs', '10') in data['selection_policy']['modalities']['peptide']['rows']
    report = render(data, extension)
    assert 'Invocation arguments, including parser defaults' in report
    assert 'Saved assembly configuration' in report
    assert before == to_native_json((selection, products, [p.antigen for _, ps in ranked for p in ps]))


def test_topiary_policy_filter_and_unknown_handling_are_shown():
    view = selection_policy_view({'run_configuration': {'effective_epitopes': {
        'selection_policy': {'score_by': 'presentation.score', 'filter_by': 'binding.rank < 2',
                             'unknown': 'exclude', 'min_score': .1}}}})
    rows = view['sections'][0]['rows']
    assert ('Explicit row filter', 'binding.rank < 2') in rows
    assert ('Unknown criterion outcome', 'exclude') in rows
    assert ('Per-allele score expression', 'presentation.score') in rows


def test_policy_labels_are_rendered_as_text():
    data = template(selection={'run_configuration': {'name': '<script>test</script>'}})
    report = render(data, 'html')
    assert '<script>test</script>' not in report
    assert '&lt;script&gt;test&lt;/script&gt;' in report


@pytest.mark.parametrize('tie', ['preferred_length', 'input_order'])
def test_window_audit_explains_exact_ties_without_changing_choice(tie):
    first = vaccine('ASIINFEKLAAA', [(1, 9, 1., False)])
    second = (replace(first, window_selection_audit={}) if tie == 'input_order' else
              vaccine('SIINFEKLAAA', [(0, 8, 1., False)]))
    audit = []
    selected = select_windows([first, second], WindowSelection(), preferred_length=11, limit=1, audit_out=audit)
    assert selected[0] is (first if tie == 'input_order' else second)
    alternatives = audit[0]['alternatives']
    winner, loser = sorted(alternatives, key=lambda a: a['selection_rank'])
    assert winner['selected'] and winner['decision_reason'] == 'within_window_limit'
    assert not loser['selected'] and loser['decision_reason'] == 'window_limit'
    assert winner['ordering_key'][0] == loser['ordering_key'][0] == -1.
    if tie == 'preferred_length':
        assert winner['ordering_key'][1] == 0 < loser['ordering_key'][1]
    else:
        assert winner['ordering_key'][:-1] == loser['ordering_key'][:-1]
        assert winner['ordering_key'][-1] == 0 < loser['ordering_key'][-1]
    view = window_view(selected[0].window_selection_audit, selected[0].amino_acids)
    assert len(view['alternatives']) == 1  # identical sequences can be distinct occurrences
    assert view['alternatives'][0]['decision_reason'] == 'window_limit'


@pytest.mark.parametrize('modality', ['peptide', 'mrna'])
def test_deduplication_precedes_capacity_and_report_links_all_sources(modality):
    ranked = [candidate(0), candidate(1, prefix='GGG'), candidate(2, 'GILGFVFTL')]
    audit = {}
    common = dict(min_antigen_length_aa=9, max_antigen_length_aa=9, max_constructs=1,
                  antigens_per_construct=1)
    if modality == 'peptide':
        products = assemble_peptide_constructs(ranked, PeptideConstructConfig(**common), selection_audit=audit)
    else:
        products = assemble_mrna_constructs(ranked, RNAConstructConfig(
            **common, optimize_linkers=False, codon_method='use_best_codon'), selection_audit=audit)
    emitted, excluded = audit['antigens']
    assert len(emitted['sources']) == 2 and emitted['outcome'] == 'emitted'
    assert excluded['outcome'] == 'construct_capacity' and excluded['products'] == []
    selection = {'construct_assembly': {modality: {'source_graph': to_native_json(
        dict(inputs=[], products={}, assembly_decisions=audit))}}}
    before = to_native_json((selection, products))
    data = template(ranked, {modality: products}, selection=selection)
    rows = data['product_groups'][0]['assembly_decisions']['rows']
    assert rows[0]['sources'] == ['S1', 'S2'] and rows[1]['sources'] == ['S3']
    assert data['ranked_candidates'][2]['assembly_outcomes'] == [
        dict(modality=modality, outcome='construct_capacity', shared=False)]
    for extension in ('html', 'txt'):
        report = render(data, extension)
        assert 'construct capacity' in report and 'shared exact sequence' in report
    links = Links(render(data, 'html'))
    assert links.targets <= set(links.identities)
    assert before == to_native_json((selection, products))


def test_window_policy_status_coordinates_and_links():
    antigen, peptides = candidate(0)
    products = assemble_peptide_constructs([(antigen, peptides)], PeptideConstructConfig(
        min_antigen_length_aa=9, max_antigen_length_aa=11, window_selection={}))
    data = template([(antigen, peptides)], {'peptide': products})
    decision, = data['window_decisions']
    parent = decision['parent']
    assert parent['sequence'][parent['start']:parent['end']] == products[0].sequence
    assert decision['sources'] == ['S1'] and decision['products'] == [products[0].name]
    assert data['product_groups'][0]['products'][0]['window_decisions'] == ['W1']
    assert decision['serum_status'] == 'Skipped: serum weight is zero'
    for extension in ('html', 'txt'):
        report = render(data, extension)
        assert 'not a completed self assessment' in report
        assert 'Skipped: serum weight is zero' in report
    links = Links(render(data, 'html'))
    assert links.targets <= set(links.identities)


def test_historical_window_metrics_do_not_invent_tie_decisions():
    peptide = vaccine('SIINFEKLAAA', [(0, 8, 1., False)])
    select_windows([peptide], WindowSelection(), preferred_length=11, limit=1)
    audit = deepcopy(peptide.window_selection_audit)
    del audit['decision']
    view = window_view(audit, peptide.amino_acids)
    assert view['decision'] is None and view['parent'] is None
    data = template()
    data['window_decisions'] = [dict(view, identifier='W1', sources=[], products=[])]
    assert 'Historical selection order and tie-break decisions were not recorded' in render(data, 'html')


def test_junction_status_uses_saved_optimization_and_actual_linkers():
    ranked = [candidate(0), candidate(1, 'GILGFVFTL')]
    product, = assemble_mrna_constructs(ranked, RNAConstructConfig(
        optimize_linkers=False, codon_method='use_best_codon'))
    data = template(ranked, {'mrna': [product]})
    assert data['product_groups'][0]['products'][0]['junction']['status'] == 'Not assessed: junction prediction disabled'
    optimized = replace(product, elements={**product.elements, 'junction_swap': dict(
        enabled=True, chosen=['AAY'], burden_strong=0, burden_mild=1,
        default_burden_strong=2, default_burden_mild=3, prediction={'endpoint': 'binding'})})
    data = template(ranked, {'mrna': [optimized]})
    view = data['product_groups'][0]['products'][0]['junction']
    assert view['status'] == 'Optimization evaluated'
    assert ('Default strong burden', '2') in view['rows']
    assert view['linkers'] == product.elements['linkers_per_junction']
    assert 'Independent per-junction benefit was not recorded' in render(data, 'txt')


def test_current_reranking_does_not_invent_historical_product_decisions():
    ranked = [candidate(0)]
    antigen, peptides = ranked[0]
    historical = assemble_peptide_constructs(ranked, PeptideConstructConfig())
    graph = to_native_json(dict(inputs=[dict(candidates=[dict(
        antigen=antigen, amino_acids=peptides[0].amino_acids, window_selection={})])], products={}))
    select_windows(peptides, WindowSelection(), preferred_length=13, limit=1)
    assert peptides[0].window_selection_audit['decision']
    data = template(ranked, {'peptide': historical},
        selection={'construct_assembly': {'peptide': {'source_graph': graph}}})
    assert data['window_decisions'] == []
    assert data['ranked_candidates'][0]['window_decisions'] == []
    assert 'A window-selection rationale was not recorded for this product' in render(data, 'txt')


def test_failed_window_eligibility_is_distinct_from_capacity():
    vp = vaccine('ASIINFEKL', [(0, 9, 1., False)])
    ranked, audit = [(vp.antigen, [vp])], []
    products = assemble_peptide_constructs(ranked, PeptideConstructConfig(
        min_antigen_length_aa=9, max_antigen_length_aa=9, n_terminal_acetylation=True,
        window_selection=dict(serum_weight=.5, serum_models=['dpp4-qpisa'])), window_audit=audit)
    assert products == [] and audit[0]['selected'] == []
    data = template(ranked, {'peptide': products},
        selection={'construct_assembly': {'peptide': {'window_audit': audit}}})
    assert data['ranked_candidates'][0]['assembly_outcomes'][0]['outcome'] == 'window_eligibility'
    assert data['window_decisions'][0]['sources'] == ['S1']
    assert 'No windows met the recorded eligibility constraints' in render(data, 'html')
    assert 'construct capacity' not in render(data, 'txt')
