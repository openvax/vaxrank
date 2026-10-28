"""Annotation follows rendered mutation evidence, not the input flag spelling."""

from dataclasses import replace
from pathlib import Path
import shutil
from types import SimpleNamespace

import pytest

from vaxrank.cli.entry_point import main, run_cli
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dataset import ANTIGEN_COLUMN, EpitopeDataset
from vaxrank.native_serialization import to_native_json
from vaxrank.vaccine_antigen import (
    AminoAcidInterval, TargetableMask, TumorSpecificityAttestation, VaccineAntigen,
)

from .input_scope_helpers import write_input_manifest
from .test_epitope_dataset import evidence, forbidden
from .test_isovar_support_contract import reference
from .test_native_references import native_dataset, delete_original_reference


@pytest.fixture
def native_antigens(tmp_path):
    sequences = ['SIINFEKL', 'GILGFVFTL']
    antigens = [VaccineAntigen(
        kind=kind, amino_acids=sequence, source_identifier=kind + '-source',
        targetable_mask=TargetableMask((AminoAcidInterval(0, len(sequence)),)),
        tumor_specificity=TumorSpecificityAttestation(
            status='admitted', evidence_kind='synthetic_fixture',
            evidence_source='test_report_annotation_requirements',
            patient_specific=True, rationale_code='test_only'))
        for kind, sequence in zip(('CTA', 'viral'), sequences)]
    dataset = EpitopeDataset.from_topiary(evidence(**{
        ANTIGEN_COLUMN: [to_native_json(a) for a in antigens]}), sample_name='patient')
    dataset.config = EpitopeConfig(score_expr='1.0')
    path = tmp_path / 'antigens.tsv'
    dataset.save(path)
    return path


def input_args(route, path, tmp_path, reference_assembly='GRCh38'):
    if route == 'single':
        return ['--input-epitopes', str(path)]
    if route == 'external':
        return ['--external-input', 'epitopes=' + str(path)]
    return ['--input-manifest', write_input_manifest(
        tmp_path / 'manifest.json', [('epitopes', path)], patient_id='patient',
        reference_assembly=reference_assembly)]


@pytest.mark.parametrize('route', ['single', 'external', 'manifest'])
def test_native_nonmutation_reports_and_design_need_no_ensembl(
        native_antigens, tmp_path, monkeypatch, route):
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    directory = tmp_path / 'design'
    html = tmp_path / 'report.html'
    run_cli(input_args(route, native_antigens, tmp_path) + [
        '--output-dir', str(directory), '--output-html-report', str(html),
        '--vaccine-type', 'peptide', '--no-processing-aware-annotation', '--report-candidate-details'])
    text = html.read_text()
    assert 'CTA-source' in text and 'viral-source' in text
    assert 'Source provenance' in text
    fasta = (directory / 'vaccine.fasta').read_text()
    assert 'SIINFEKL' in fasta and 'GILGFVFTL' in fasta


def test_pooled_native_inputs_render_real_reports_without_ensembl(native_antigens, tmp_path):
    dataset = EpitopeDataset.load(native_antigens)
    paths = []
    for index, epitope in enumerate(dataset.epitopes):
        part = replace(dataset, epitopes=(epitope,))
        path = tmp_path / ('part%d.tsv' % index)
        part.save(path)
        paths.append(('epitopes', path))
    manifest = write_input_manifest(tmp_path / 'pool.json', paths, patient_id='patient')
    report = tmp_path / 'pool.txt'
    run_cli(['--input-manifest', manifest, '--output-ascii-report', str(report),
             '--no-processing-aware-annotation', '--report-candidate-details'])
    text = report.read_text()
    assert 'CTA-source' in text and 'viral-source' in text


@pytest.mark.parametrize('mixed', [False, True])
def test_mutation_and_mixed_reports_still_require_annotation(native_antigens, tmp_path, mixed):
    lens = Path(__file__).parent / 'data' / 'epitope_fixtures' / 'lens_example.tsv'
    paths = [('lens', lens)] + ([('epitopes', native_antigens)] if mixed else [])
    manifest = write_input_manifest(tmp_path / 'mixed.json', paths, patient_id='patient')
    report = tmp_path / 'report.txt'
    with pytest.raises(ValueError, match='mutation antigen.*--ensembl-release'):
        run_cli(['--input-manifest', manifest, '--output-ascii-report', str(report),
                 '--config-text', 'epitopes.score_expr=1.0', '--no-processing-aware-annotation', '--report-candidate-details'])
    assert not report.exists()


def test_mutations_filtered_out_of_report_do_not_require_annotation(tmp_path):
    lens = Path(__file__).parent / 'data' / 'epitope_fixtures' / 'lens_example.tsv'
    report = tmp_path / 'empty.txt'
    run_cli(['--external-input', 'lens=' + str(lens), '--output-ascii-report', str(report),
             '--config-text', 'epitopes.filter_expr=affinity.value < 0',
             '--no-processing-aware-annotation', '--report-candidate-details'])
    assert report.is_file()


@pytest.fixture
def mutation_native(tmp_path):
    genome = reference.__wrapped__(tmp_path, SimpleNamespace(param='+'))
    return native_dataset(tmp_path, genome)


@pytest.mark.parametrize('route', ['single', 'external', 'manifest'])
def test_native_mutation_reports_use_retained_annotation(mutation_native, tmp_path, monkeypatch, route):
    _, path = mutation_native
    delete_original_reference(tmp_path)
    monkeypatch.setattr('pyensembl.download_cache.DownloadCache._fetch', forbidden)
    report = tmp_path / 'mutation.txt'
    run_cli(input_args(route, path, tmp_path, reference_assembly='support-test-+') + [
        '--index-native-references', '--output-ascii-report', str(report),
        '--vaccine-peptide-length', '10', '--no-processing-aware-annotation', '--report-candidate-details'])
    text = report.read_text()
    assert 'Effect type: Substitution' in text and 'Effect description: p.K6Q' in text


def test_native_manifest_rejects_conflicting_retained_reference(mutation_native, tmp_path):
    _, path = mutation_native
    report = tmp_path / 'wrong-reference.txt'
    with pytest.raises(ValueError, match='conflicting reference_assembly'):
        run_cli(input_args('manifest', path, tmp_path, reference_assembly='GRCh38') + [
            '--output-ascii-report', str(report), '--no-processing-aware-annotation', '--report-candidate-details'])
    assert not report.exists()


def test_native_mutation_without_transcripts_requires_annotation(mutation_native, tmp_path):
    dataset, path = mutation_native
    dataset.mutation_fragments = {key: replace(fragment, supporting_reference_transcripts=[])
                                  for key, fragment in dataset.mutation_fragments.items()}
    dataset.save(path)
    report = tmp_path / 'unannotated.txt'
    with pytest.raises(ValueError, match='mutation antigen.*--ensembl-release'):
        run_cli(['--external-input', 'epitopes=' + str(path), '--output-ascii-report', str(report),
                 '--vaccine-peptide-length', '10', '--no-processing-aware-annotation', '--report-candidate-details'])
    assert not report.exists()


def test_missing_native_annotation_is_a_clear_cli_error(mutation_native, tmp_path, monkeypatch, capsys):
    _, path = mutation_native
    delete_original_reference(tmp_path)
    shutil.rmtree(str(path) + '.references')
    monkeypatch.setattr('pyensembl.download_cache.DownloadCache._fetch', forbidden)
    report = tmp_path / 'missing.txt'
    with pytest.raises(SystemExit) as error:
        main(['--input-epitopes', str(path), '--output-ascii-report', str(report),
              '--vaccine-peptide-length', '10', '--no-processing-aware-annotation', '--report-candidate-details'])
    assert error.value.code == 1
    message = capsys.readouterr().err
    assert 'vaxrank: error:' in message and '.references' in message
    assert 'Traceback' not in message and not report.exists()


@pytest.mark.parametrize('route', ['single', 'external', 'manifest'])
def test_compact_mutation_report_needs_no_annotation_queries(mutation_native, tmp_path, monkeypatch, route):
    _, path = mutation_native
    delete_original_reference(tmp_path)
    shutil.rmtree(str(path) + '.references')
    monkeypatch.setattr('pyensembl.download_cache.DownloadCache._fetch', forbidden)
    monkeypatch.setattr('vaxrank.report.TemplateDataCreator.effect_data', forbidden)
    report = tmp_path / 'compact.txt'
    run_cli(input_args(route, path, tmp_path, reference_assembly='support-test-+') + [
        '--output-ascii-report', str(report), '--vaccine-peptide-length', '10',
        '--no-processing-aware-annotation'])
    text = report.read_text()
    assert 'Candidate ranking' in text and 'AAAAQGG' in text
    assert 'Candidate window details' not in text
