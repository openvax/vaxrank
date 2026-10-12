"""Scored evidence -> all windows -> per-core exposure, without reranking."""

from dataclasses import replace
import json

import pytest
from mhctools import vaccine_route_steps

from vaxrank.core_exposure import enumerate_exposure_windows, assess_dataset_windows
from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.cli.window_exposure import main
from .test_window_selection import vaccine


def dataset():
    vp = vaccine('AASIINFEKLGILGFVFTLAA', [(2, 10, 1., False), (10, 19, 3., False)])
    # The prediction identity is source-local and survives window slicing.
    epitopes = tuple(replace(epitope, prediction_id='core-%s' % i) for i, epitope in enumerate(vp.epitopes))
    return EpitopeDataset.from_predictions(epitopes, antigens={e.prediction_group_source: vp.antigen for e in epitopes})


def profile():
    return {'scenario': 'Synthetic sensitivity test, not human IM kinetics', 'times_hours': [0, 1, 48]}


def test_all_windows_keep_multiple_scored_occurrences_and_source_coordinates():
    evidence = dataset()
    windows, gaps = enumerate_exposure_windows(evidence, 19, 20)
    assert len(windows) == 5  # Three 19-mers and two 20-mers; no one-best pruning.
    assert not gaps
    request, result = assess_dataset_windows(evidence, profile(), min_length=19, max_length=20)
    assert len(result['candidates']) == 5
    for window in request['candidates']:
        assert len(window['cores']) == 2
        assert [core['score_contribution'] for core in window['cores']] == [1., 3.]
        assert window['metadata']['gene_name'] == 'PRAME'
        assert window['peptide']['source_start'] == window['metadata']['source_start']
        for core in window['cores']:
            assert window['peptide']['sequence'][core['start']:core['end']] == core['sequence']
            assert core['prediction_id'] in ('core-0', 'core-1')
    assert all(w['status'] == 'unassessed' for w in result['candidates'])
    assert evidence.epitopes[0].per_allele_scores == {'HLA-A*02:01': 1.}


def test_conditional_delivery_uses_existing_scores_and_one_microgram():
    evidence = dataset()
    sequence = evidence.epitopes[0].source_sequence
    scenario = profile()
    rates = {s.name: {'value_per_hour': 0, 'basis': 'disabled', 'source': 'Synthetic excluded route'}
             for s in vaccine_route_steps('slp', 'I')}
    for name in ('lymph_entry', 'lymph_transit', 'node_antigen_uptake', 'interstitium_clearance'):
        rates[name] = {'value_per_hour': 1, 'basis': 'assumed', 'source': 'Analytic competition oracle'}
    cuts = [{'compartment': compartment, 'start': 0, 'end': len(sequence), 'sequence': sequence,
             'channels': [], 'basis': 'assumed', 'source': 'Synthetic no-cut oracle'}
            for compartment in ('interstitium', 'afferent_lymph', 'node_fluid', 'node_endosome')]
    scenario['kinetics'] = {'I': {'tap_max_length': len(sequence), 'rates': rates, 'cleavage_assessments': cuts}}
    _, result = assess_dataset_windows(evidence, scenario, min_length=len(sequence), max_length=len(sequence))
    window = result['candidates'][0]
    assert window['dose']['amount'] == 1
    assert window['dose']['unit'] == 'ug'
    assert window['aggregate']['score_weighted_mean_node_entry_per_input'] == pytest.approx(.5)
    assert all(core['exposure']['node_apc_entry']['fraction_of_dosed_peptide_copies'] == pytest.approx(.5)
               for core in window['cores'])


def test_repeated_sequence_and_multiple_alleles_survive():
    vp = vaccine('AASIINFEKLASIINFEKLAA', [(2, 10, 1., False), (11, 19, 2., False)])
    epitopes = tuple(replace(e, prediction_id='repeat-%s' % i,
                             per_allele_scores={'HLA-A*02:01': float(i+1), 'HLA-B*07:02': .5}) for i, e in enumerate(vp.epitopes))
    data = EpitopeDataset.from_predictions(epitopes, antigens={e.prediction_group_source: vp.antigen for e in epitopes})
    request, _ = assess_dataset_windows(data, profile(), min_length=21, max_length=21)
    cores = request['candidates'][0]['cores']
    assert len(cores) == 4
    assert {core['start'] for core in cores} == {2, 11}
    assert {core['allele'] for core in cores} == {'HLA-A*02:01', 'HLA-B*07:02'}


def test_missing_context_and_modifications_remain_visible():
    data = dataset()
    del data.antigens['core-0']
    request, result = assess_dataset_windows(data, profile(), min_length=20, max_length=20, c_term='amidated')
    assert result['source_coverage_gaps'][0]['reason'] == 'missing_antigen_context'
    assert request['candidates'][0]['peptide']['c_term'] == 'amidated'
    assert 'unsupported C-terminal chemistry' in result['candidates'][0]['cores'][0]['missing_inputs'][0]


def test_unknown_profile_window_does_not_silently_drop_rates():
    scenario = profile()
    scenario['candidates'] = {'wrong-window-id': {'kinetics': {}}}
    with pytest.raises(ValueError, match='unknown candidate window'):
        assess_dataset_windows(dataset(), scenario, min_length=20, max_length=20)


def test_saved_dataset_cli_outputs_all_windows_and_protects_previous_run(tmp_path):
    data = dataset()
    source, scenario, output = tmp_path / 'native.tsv', tmp_path / 'profile.json', tmp_path / 'exposure'
    data.save(source)
    scenario.write_text(json.dumps(profile()))
    args = ['--input-epitopes', str(source), '--profile', str(scenario), '--min-length', '19', '--max-length', '20',
            '--n-term', 'free', '--c-term', 'free', '--output-dir', str(output)]
    main(args)
    result = json.loads((output / 'exposure.json').read_text())
    assert len(result['candidates']) == 5
    assert len((output / 'cores.csv').read_text().splitlines()) == 11
    assert len(result['input_files']) == 2
    assert 'no shorter binding core' in (output / 'report.md').read_text()
    original = (output / 'exposure.json').read_bytes()
    with pytest.raises(SystemExit):
        main(args)
    assert (output / 'exposure.json').read_bytes() == original


def test_explicit_core_interval_keeps_its_original_ligand_score_and_source():
    scenario = profile()
    scenario['core_intervals'] = [{'prediction_id': 'core-0', 'allele': 'HLA-A*02:01',
                                  'start': 1, 'end': 7, 'source': 'Synthetic coordinate test annotation'}]
    request, _ = assess_dataset_windows(dataset(), scenario, min_length=19, max_length=20)
    for window in request['candidates']:
        core = next(c for c in window['cores'] if c['prediction_id'] == 'core-0')
        assert core['sequence'] == 'IINFEK'
        assert core['scored_ligand'] == 'SIINFEKL'
        assert core['score_contribution'] == 1.
        assert core['start'] == 3 - window['metadata']['source_start']
        assert core['core_source'] == 'Synthetic coordinate test annotation'
    scenario['core_intervals'][0]['end'] = 99
    with pytest.raises(ValueError, match='within its scored ligand'):
        assess_dataset_windows(dataset(), scenario, min_length=19, max_length=20)


def test_overlapping_ligands_sharing_a_documented_core_keep_all_score_components():
    vp = vaccine('AASIINFEKLAA', [(1, 11, 1., False), (2, 10, 2., False)])
    epitopes = tuple(replace(e, prediction_id='shared-%s' % i) for i, e in enumerate(vp.epitopes))
    evidence = EpitopeDataset.from_predictions(epitopes, antigens={e.prediction_group_source: vp.antigen for e in epitopes})
    scenario = profile()
    scenario['core_intervals'] = [
        {'prediction_id': 'shared-0', 'allele': 'HLA-A*02:01', 'start': 1, 'end': 9, 'source': 'Synthetic annotation 1'},
        {'prediction_id': 'shared-1', 'allele': 'HLA-A*02:01', 'start': 0, 'end': 8, 'source': 'Synthetic annotation 2'}]
    request, _ = assess_dataset_windows(evidence, scenario, min_length=12, max_length=12)
    cores = request['candidates'][0]['cores']
    assert len(cores) == 1
    assert cores[0]['sequence'] == 'SIINFEKL'
    assert cores[0]['score_contribution'] == 3.
    assert len(cores[0]['score_components']) == 2
    assert {c['prediction_id'] for c in cores[0]['score_components']} == {'shared-0', 'shared-1'}


def test_identical_ligands_at_different_positions_can_document_one_shared_core():
    vp = vaccine('AAAAAA', [(0, 4, 1., False), (1, 5, 2., False)])
    epitopes = tuple(replace(e, prediction_id='position-%s' % i) for i, e in enumerate(vp.epitopes))
    evidence = EpitopeDataset.from_predictions(epitopes, antigens={e.prediction_group_source: vp.antigen for e in epitopes})
    scenario = profile()
    scenario['core_intervals'] = [
        {'prediction_id': 'position-0', 'allele': 'HLA-A*02:01', 'start': 1, 'end': 3, 'source': 'Synthetic annotation 1'},
        {'prediction_id': 'position-1', 'allele': 'HLA-A*02:01', 'start': 0, 'end': 2, 'source': 'Synthetic annotation 2'}]
    request, _ = assess_dataset_windows(evidence, scenario, min_length=6, max_length=6)
    core, = request['candidates'][0]['cores']
    assert (core['start'], core['end'], core['score_contribution']) == (1, 3, 3.)
    assert {c['scored_ligand_start'] for c in core['score_components']} == {0, 1}
