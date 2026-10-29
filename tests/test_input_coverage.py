"""Original denominators survive native export independently of surviving rows."""

from dataclasses import replace

from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.input_coverage import record_direct_coverage
from vaxrank.input_scope import InputProvenance, InputScope

from .test_variant_counts import _isovar_result, _results


def test_direct_counts_keep_rejected_variants_and_replay_without_predictions(tmp_path):
    first, a = _isovar_result(coding=True, rna=True)
    second, b = _isovar_result(coding=True, rna=False)
    third, c = _isovar_result(coding=False, rna=True)
    unannotatable, _ = _isovar_result(coding=False, rna=False)
    dataset = EpitopeDataset.from_predictions(())
    dataset.provenance = (InputProvenance('original', 'vcf_bam', 'input.vcf; input.bam', 'test', InputScope()),)
    record = record_direct_coverage(dataset, _results([a, b, c], [first]),
                                   [first, second, third, unannotatable])
    assert record['counts'] == dict(loaded=4, annotatable=3, coding_effect=2,
                                   coding_with_rna_protein=1, with_candidate_windows=1)
    assert len(record['annotation_unavailable']) == 1
    assert len(record['variant_outcomes']) == 3
    assert not dataset.epitopes  # Missing downstream rows must not erase the input universe.
    path = tmp_path / 'native.tsv'
    dataset.save(path)
    restored = EpitopeDataset.load(path)
    assert restored.selection['source_evidence'] == dataset.selection['source_evidence']


def test_distinct_observations_of_one_input_are_not_overwritten():
    first, result = _isovar_result(coding=True, rna=True)
    dataset = EpitopeDataset.from_predictions(())
    dataset.provenance = (InputProvenance('original', 'vcf_bam', 'input', 'test', InputScope()),)
    record_direct_coverage(dataset, _results([result], [first]), [first])
    record_direct_coverage(dataset, _results([result]), [first])
    observations = [v['direct_input_coverage'] for v in dataset.selection['source_evidence'].values()]
    assert len(observations) == 2
    assert {o['counts']['with_candidate_windows'] for o in observations} == {0, 1}
    assert {o['input_id'] for o in observations} == {'original'}
    # Repeating the same observation does not create a second denominator.
    record_direct_coverage(dataset, _results([result]), [first])
    assert len(dataset.selection['source_evidence']) == 2
    dataset.provenance = (replace(dataset.provenance[0], input_id='other-library'),)
    record_direct_coverage(dataset, _results([result]), [first])
    assert len(dataset.selection['source_evidence']) == 3


def _coverage_fixture():
    from types import SimpleNamespace
    from vaxrank.antigen_input import AntigenDataset
    from vaxrank.external_input import ExternalRankingResult, ExternalVariantEntry
    from vaxrank.input_coverage import ranking_coverage
    from .test_assembly_antigens import candidate

    ranked = [candidate(0), candidate(1, 'GILGFVFTL')]
    held, _ = candidate(2)
    held = replace(held, kind='fusion', source_identifier='PARTNER1--PARTNER2',
                   tumor_specificity=replace(held.tumor_specificity, status='held_out',
                                              rationale_code='insufficient_rna_fragments,ambiguous_coding_candidate'))
    missing, _ = candidate(3)
    catalog = AntigenDataset({a.source_identifier: a for a in [ranked[0][0], ranked[1][0], held, missing]})
    epitopes = [e for _, peptides in ranked for p in peptides for e in p.epitopes]
    native = EpitopeDataset.from_predictions(epitopes)
    original = InputProvenance('original', 'antigens', 'original.json', 'a' * 64,
        InputScope(patient_id='P', sample_id='S', library_id='L', timepoint='T',
                   rna_evidence_extent='regional_subset'))
    native.provenance = (original,)
    native.antigens = {p.epitopes[0].prediction_group_source: a for a, peptides in ranked for p in peptides}
    native.selection['antigen_inputs'] = {'original': catalog.to_dict()['antigens']}
    report = SimpleNamespace(dataset=native, input_provenance=replace(original, input_id='replay',
        path='replay.tsv', source_format='epitopes', scope=InputScope(patient_id='P')))
    outcome = ExternalRankingResult(entries=tuple(ExternalVariantEntry(source=a, vaccine_peptides=tuple(p))
                                                   for a, p in ranked))
    return native, report, outcome, ranked, ranking_coverage(native, [report], [outcome], ranked)


def test_coverage_reconciles_held_out_missing_predictions_and_actual_product_cap(tmp_path):
    from vaxrank.input_coverage import coverage_view, ranking_coverage
    from vaxrank.peptide import assemble_peptide_constructs, PeptideConstructConfig
    from .test_product_reports import template

    native, report, outcome, ranked, coverage = _coverage_fixture()
    products = {'peptide': assemble_peptide_constructs(ranked, PeptideConstructConfig(max_constructs=1))}
    ranking = template(ranked, products)['ranked_candidates']
    view = coverage_view(coverage, {}, ranking, products)
    rows = {r['identity']: r for r in view['sources']}
    assert rows['PARTNER1--PARTNER2']['stage'] == 'Construct admission'
    assert 'ambiguous_coding_candidate' in rows['PARTNER1--PARTNER2']['reasons'][0]
    assert rows['candidate-3']['retained_epitopes'] == 0
    assert 'No retained epitope predictions' in rows['candidate-3']['reasons']
    assert rows['candidate-0']['selected'] is True
    assert rows['candidate-1']['selected'] is False
    assert rows['candidate-1']['rank'] == 2
    assert not view['direct']  # Historical source records cannot recover a VCF denominator.
    families = {f['name']: f for f in view['families']}
    assert families['CTAs']['status'] == families['Viral sequences']['status'] == 'No supplied evidence recorded'
    assert families['Small variants']['records'] == 3
    assert families['Small variants']['ranked'] == 2
    assert families['Small variants']['represented'] == 1
    assert coverage['original_inputs'][0]['scope']['rna_evidence_extent'] == 'regional_subset'
    assert coverage['current_inputs'][0]['scope']['sample_id'] is None
    path = tmp_path / 'native.tsv'
    native.save(path)
    report.dataset = EpitopeDataset.load(path)
    restored = ranking_coverage(report.dataset, [report], [outcome], ranked)
    assert restored == coverage


def test_coverage_unknown_selection_is_not_reported_as_zero():
    from vaxrank.input_coverage import coverage_view
    from vaxrank.peptide import assemble_peptide_constructs, PeptideConstructConfig
    from .test_product_reports import template
    _, _, _, ranked, coverage = _coverage_fixture()
    product, = assemble_peptide_constructs(ranked, PeptideConstructConfig(max_constructs=1))
    historical = replace(product, antigens=[])
    for products in ({}, {'peptide': [historical]}):
        view = coverage_view(coverage, {}, template(ranked, products)['ranked_candidates'], products)
        assert all(f['represented'] is None for f in view['families'])
    empty = coverage_view(coverage, {}, template(ranked, {'peptide': []})['ranked_candidates'], {'peptide': []})
    assert all(f['represented'] == 0 for f in empty['families'])


def test_unresolved_producer_models_are_visible_without_any_antigen_or_prediction():
    from vaxrank.input_coverage import coverage_view
    evidence = {'evidence-1': {'family': 'variant_orfs', 'decisions': [], 'isovar_result': {
        'sample_id': 'sample', 'rna_candidates': [], 'models': [dict(model_id='canonical-CDS',
        reason='outside_reference_orf', state='unresolved')], 'variant': {'contig': '1', 'start': 22}}}}
    view = coverage_view(dict(sources=[], original_inputs=[], current_inputs=[]), evidence, [], {})
    producer, = view['producers']
    assert producer['models'] == 1 and producer['candidates'] == producer['rna_reconstructions'] == 0
    assert producer['unresolved'][0]['identity'] == 'canonical-CDS'
    assert producer['unresolved'][0]['reason'] == 'outside_reference_orf'
    assert next(f for f in view['families'] if f['name'] == 'Additional variant ORFs')['status'] == 'Supplied'


def test_lineage_and_exclusion_ledger_are_in_both_report_formats():
    from vaxrank.patient_info import PatientInfo
    from vaxrank.report import TemplateDataCreator
    from .test_product_reports import render
    _, _, _, ranked, coverage = _coverage_fixture()
    data = TemplateDataCreator(ranked, PatientInfo('P', input_coverage=coverage), '', '',
        dict(manufacturability=False, wt_epitopes=False), None).compute_template_data()
    for extension in ('html', 'txt'):
        text = render(data, extension)
        for phrase in ('PARTNER1--PARTNER2', 'No retained epitope predictions', 'regional_subset',
                       'original.json', 'replay.tsv', 'Original input receipts', 'Current report inputs'):
            assert phrase in text


def test_extent_is_explicit_and_does_not_change_historical_input_identities():
    import pytest
    from vaxrank.input_scope import scope_dict, scope_from_mapping
    assert 'rna_evidence_extent' not in scope_dict(InputScope())
    scope = scope_from_mapping({'rna_evidence_extent': 'regional_subset'}, 'manifest')
    assert scope.rna_evidence_extent == 'regional_subset'
    with pytest.raises(ValueError, match='rna_evidence_extent'):
        InputScope(rna_evidence_extent='assumed_full')


def test_deferred_window_selection_is_not_a_zero_candidate_observation():
    first, result = _isovar_result(coding=True, rna=True)
    result.filter_values = {'min_alt_fragments': False, 'max_mismatches': True}
    dataset = EpitopeDataset.from_predictions(())
    dataset.provenance = (InputProvenance('original', 'vcf_bam', 'input', 'test', InputScope()),)
    record = record_direct_coverage(dataset, _results([result]), [first], window_selection_performed=False)
    assert record['counts']['with_candidate_windows'] is None
    assert record['variant_outcomes'][0]['has_vaccine_peptide'] is None
    assert record['counts']['coding_with_rna_protein'] == 1

    assert record['variant_outcomes'][0]['isovar_filter_results'] == result.filter_values


def test_real_direct_pipeline_and_native_replay_keep_original_coverage(tmp_path, monkeypatch):
    from pathlib import Path
    from vaxrank.cli.entry_point import run_cli
    data = Path(__file__).parent / 'data/b16.f10'
    native = tmp_path / 'native.tsv'
    run_cli(['--vcf', str(data / 'b16.vcf'), '--bam', str(data / 'b16.combined.bam'),
        '--mhc-predictor', 'random', '--mhc-alleles', 'H2-Kb,H2-Db', '--mhc-epitope-lengths', '8',
        '--no-processing-aware-annotation', '--vaccine-peptide-length', '15',
        '--output-epitopes', str(native), '--output-ascii-report', str(tmp_path / 'fresh.txt')])
    dataset = EpitopeDataset.load(native)
    observation, = (r['direct_input_coverage'] for r in dataset.selection['source_evidence'].values())
    assert observation['counts']['loaded'] == 5
    assert len(dataset.direct_sources) == observation['counts']['annotatable'] == 5
    assert observation['counts']['with_candidate_windows'] > 0
    assert observation['isovar_filter_thresholds']

    def forbidden(*args, **kwargs):
        raise AssertionError('Coverage replay cannot acquire predictions')
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    run_cli(['--input-epitopes', str(native), '--no-processing-aware-annotation',
        '--output-epitopes', str(tmp_path / 'replay.tsv'),
        '--output-ascii-report', str(tmp_path / 'replay.txt')])
    assert EpitopeDataset.load(tmp_path / 'replay.tsv').selection['source_evidence'] == dataset.selection['source_evidence']
    fresh, replay = ((tmp_path / (name + '.txt')).read_text() for name in ('fresh', 'replay'))
    def coverage(text):
        return text.split('=== Input coverage ===', 1)[1].split('Emitted products', 1)[0]
    assert coverage(fresh) == coverage(replay)
    assert 'Recorded sources: 5' in coverage(fresh)
