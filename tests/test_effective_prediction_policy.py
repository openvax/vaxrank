"""Effective policies and contextual annotations survive composed CLI replay."""

import json
import math
import pytest

from vaxrank.cli.entry_point import run_cli
from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.processing_prediction import ProcessingPrediction
from vaxrank.prediction_policy import retained_processing

from .input_scope_helpers import write_input_manifest
from .test_junction_prediction_policy import INPUTS, forbidden, run_external
from .test_junction_prediction_policy import junction_factory  # noqa: F401 -- shared pytest fixture


@pytest.mark.usefixtures('junction_factory')
def test_explicit_junction_endpoint_is_reported_and_replayed_without_inference(tmp_path, monkeypatch):
    directory = run_external(tmp_path, 'first', [
        '--mrna-junction-predictor', 'junction-test',
        '--mrna-junction-rank-kind', 'pMHC_presentation',
        '--mrna-junction-rank-strong', '0.7', '--mrna-junction-rank-mild', '3',
        '--output-html-report', str(tmp_path / 'first.html'), '--report-candidate-details'])
    first = EpitopeDataset.load(directory / 'candidate_predictions.tsv')
    stages = first.selection['run_configuration']['effective_predictions']
    assert stages['candidate']['status'] == 'recorded evidence; no new inference'
    assert stages['candidate']['policy']['unqualified_methods']['pMHC_affinity'] == 'mhcflurry'
    assert stages['junction']['status'] == 'executed'
    run, = stages['junction']['runs']
    assert run['rank_kind'] == 'pMHC_presentation'
    assert run['rank_kind_policy'] == 'explicit'
    assert run['rank_strong'] == .7 and run['rank_mild'] == 3
    assert run['coverage'] == 'complete'
    inventory, = run['configured_models']
    assert inventory['prediction_method_name'] == 'junction-test'
    assert inventory['predictor_version'] == '2.0'
    assert set(inventory['alleles']) == set(run['alleles'])
    summary = (directory / 'run_summary.txt').read_text()
    assert 'junction-test [2.0]: pMHC_presentation' in summary
    assert 'strong ≤ 0.7%' in summary
    assert 'Effective prediction configuration' in (directory / 'vaccine_report.txt').read_text()
    html = (tmp_path / 'first.html').read_text()
    assert html.index('id="products"') < html.index('id="run-details"') < html.index('id="effective-predictions"')
    assert 'Disabled for this run' in html
    assert 'junction-test [2.0]: pMHC_presentation' in html
    monkeypatch.setattr('vaxrank.cli.entry_point.mhc_binding_predictor_from_args', forbidden)
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    replay = tmp_path / 'replay'
    run_cli(['--input-epitopes', str(directory / 'candidate_predictions.tsv'),
             '--output-dir', str(replay), '--pdf-backend', 'weasyprint',
             '--no-processing-aware-annotation'])
    restored = EpitopeDataset.load(replay / 'candidate_predictions.tsv')
    assert restored.selection['construct_assembly'] == first.selection['construct_assembly']
    assert restored.selection['run_configuration']['effective_predictions']['junction']['status'] == (
        'replayed evidence; no new inference')
    assert json.loads((replay / 'manifest.json').read_text()) == json.loads(
        (directory / 'manifest.json').read_text())


def test_processing_occurrence_evidence_is_retained_without_reinference(tmp_path, monkeypatch):
    annotations = {}

    def annotate(ranked, predictions, **kwargs):
        epitope = next(e for e in predictions if e.source_sequence)
        record = ProcessingPrediction(
            peptide_sequence=epitope.sequence, source_sequence=epitope.source_sequence,
            peptide_offset=epitope.offset, predictor_name='pepsickle', predictor_version='test-1',
            c_term_cleavage_prob=.4, max_internal_cut_prob=.3,
            processing_score=math.sqrt(.4 * .7))
        annotations[record.key()] = record
        return 1, annotations

    monkeypatch.setattr('vaxrank.cli.entry_point.annotate_predictions_with_processing', annotate)
    source = tmp_path / 'first.tsv'
    run_cli(['--input-manifest', write_input_manifest(tmp_path / 'inputs.json', INPUTS[:1]),
             '--output-epitopes', str(source), '--processing-aware-annotation',
             '--pepsickle-threshold', '0.7'])
    first = EpitopeDataset.load(source)
    assert retained_processing([first]) == annotations
    stage = first.selection['run_configuration']['effective_predictions']['processing']
    assert stage['configuration']['threshold'] == .7
    assert stage['coverage'] == dict(annotated_occurrences=1, assessed_c_boundaries=1)
    monkeypatch.setattr('vaxrank.cli.entry_point.annotate_predictions_with_processing', forbidden)
    output = tmp_path / 'replay.tsv'
    run_cli(['--input-epitopes', str(source), '--output-epitopes', str(output)])
    replay = EpitopeDataset.load(output)
    assert retained_processing([replay]) == annotations
    replay_stage = replay.selection['run_configuration']['effective_predictions']['processing']
    assert replay_stage['status'] == 'replayed evidence; no new inference'
    assert replay_stage['configuration'] == stage['configuration']
    assert replay.epitopes == first.epitopes
