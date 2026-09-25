"""Scope compatibility gates real imports, before any scoring or prediction."""

from dataclasses import asdict
import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from vaxrank.cli.arg_parser import parse_vaxrank_args
from vaxrank.cli.entry_point import main, write_run_summary
from vaxrank.epitope_io import load_predictions, save_predictions
from vaxrank.external_input import load_external_ranked
from vaxrank.external_rescoring import prepare_reports, external_inputs
from vaxrank.input_scope import read_input_manifest
from .input_scope_helpers import FIXTURE_SCOPE, write_input_manifest


DATA = Path(__file__).parent / 'data/epitope_fixtures'
INPUTS = [('lens', str(DATA / 'lens_example.tsv')),
          ('pvacseq', str(DATA / 'pvacseq_example.tsv'))]


def manifest_args(tmp_path, inputs=INPUTS, **scope):
    return SimpleNamespace(input_manifest=write_input_manifest(
        tmp_path / 'inputs.json', inputs, **scope))


def forbidden(*args, **kwargs):
    raise AssertionError('Invalid scope reached scoring, output or live prediction')


@pytest.mark.parametrize('field', ['patient_id', 'reference_assembly', 'mhc_alleles'])
def test_missing_shared_scope_is_rejected_before_scoring_or_export(tmp_path, monkeypatch, field):
    import vaxrank.external_rescoring as module
    monkeypatch.setattr(module, 'attach_per_allele_scores', forbidden)
    monkeypatch.setattr(module, 'save_predictions', forbidden)
    scope = {**FIXTURE_SCOPE, field: None}
    with pytest.raises(ValueError, match=field):
        prepare_reports(INPUTS, scopes=[scope, scope], input_predictions_path=tmp_path / 'out.tsv')


def test_output_patient_label_cannot_authorize_pooling():
    args = SimpleNamespace(external_input=[f'{fmt}={path}' for fmt, path in INPUTS],
                           output_patient_id='a-label')
    with pytest.raises(ValueError, match='--input-manifest'):
        load_external_ranked(args)


def test_complete_producer_declarations_allow_pooling_without_a_manifest(tmp_path):
    inputs = []
    for fmt, filename in INPUTS:
        frame = pd.read_csv(filename, sep='\t')
        for key, value in FIXTURE_SCOPE.items():
            frame[key] = json.dumps(value) if isinstance(value, list) else value
        path = tmp_path / f'{fmt}.tsv'
        frame.to_csv(path, sep='\t', index=False)
        inputs.append(f'{fmt}={path}')
    ranked, _, _, patient, _ = load_external_ranked(SimpleNamespace(external_input=inputs))
    assert ranked
    assert patient.patient_id == 'fixture-patient'
    assert all(p.manifest_path is None for p in patient.input_provenance)
    assert all(p.report_declarations['patient_id'] == ['fixture-patient']
               for p in patient.input_provenance)


@pytest.mark.parametrize('field,values', [
    ('patient_id', ['patient-one', 'patient-two']),
    ('annotation', ['ensembl:109', 'ensembl:110']),
])
def test_per_input_conflicts_are_rejected_without_a_shared_default(tmp_path, field, values):
    args = manifest_args(tmp_path)
    path = Path(args.input_manifest)
    document = json.loads(path.read_text())
    document.pop(field, None)
    for entry, value in zip(document['inputs'], values):
        entry[field] = value
    path.write_text(json.dumps(document))
    with pytest.raises(ValueError, match=f'conflicting {field}'):
        load_external_ranked(args)


@pytest.mark.parametrize('field,other', [
    ('patient_id', 'different-patient'), ('reference_assembly', 'GRCh37'),
    ('annotation', 'ensembl:110'), ('mhc_alleles', ['HLA-A*01:01']),
])
def test_conflicting_input_declarations_fail_before_models_or_scoring(tmp_path, monkeypatch, field, other):
    import mhctools.cli
    import vaxrank.external_rescoring as module
    path = tmp_path / 'inputs.json'
    write_input_manifest(path, INPUTS, annotation='ensembl:109')
    document = json.loads(path.read_text())
    document['inputs'][1][field] = other
    path.write_text(json.dumps(document))
    args = parse_vaxrank_args([
        '--input-manifest', str(path), '--external-predictions', 'fresh',
        '--mhc-predictor', 'random', '--mhc-alleles', 'HLA-A*02:01'])
    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', forbidden)
    monkeypatch.setattr(module, 'attach_per_allele_scores', forbidden)
    with pytest.raises(ValueError, match=f'conflicting {field}'):
        load_external_ranked(args)


@pytest.mark.parametrize('field,value', [
    ('patient_id', 'another-patient'), ('reference_assembly', 'GRCh37'),
    ('annotation', 'ensembl:110'), ('sample_id', 'another-sample'),
    ('mhc_alleles', json.dumps(['HLA-A*01:01'])),
])
def test_manifest_cannot_override_producer_declarations(tmp_path, field, value):
    source = pd.read_csv(INPUTS[1][1], sep='\t')
    source[field] = value
    path = tmp_path / 'producer.tsv'
    source.to_csv(path, sep='\t', index=False)
    args = manifest_args(tmp_path, [('pvacseq', str(path))],
                         annotation='ensembl:109', sample_id='first-sample')
    with pytest.raises(ValueError, match=f'conflicting {field}'):
        load_external_ranked(args)


def test_all_producer_rows_are_checked_including_late_lens_build_marker(tmp_path):
    frame = pd.read_csv(INPUTS[0][1], sep='\t').iloc[[0]]
    frame['origin_descriptor'] = 'Hsap38.chr1.1.10.+'
    frame = pd.concat([frame] * 502, ignore_index=True)
    frame.loc[501, 'origin_descriptor'] = 'Hsap37.chr1.1.10.+'
    path = tmp_path / 'late-conflict.tsv'
    frame.to_csv(path, sep='\t', index=False)
    with pytest.raises(ValueError, match='conflicting reference_assembly'):
        prepare_reports([('lens', str(path))])


def test_partial_allele_coverage_is_not_a_genotype_conflict(tmp_path):
    inputs = []
    for index, (fmt, filename) in enumerate(INPUTS):
        frame = pd.read_csv(filename, sep='\t')
        frame = frame.iloc[[index]]
        path = tmp_path / f'{fmt}.tsv'
        frame.to_csv(path, sep='\t', index=False)
        inputs.append((fmt, str(path)))
    args = manifest_args(tmp_path, inputs, mhc_alleles=['A0201', 'B0702', 'C0702'])
    ranked, report, predictions, patient, _ = load_external_ranked(args)
    assert ranked and predictions
    assert patient.mhc_alleles == ['HLA-A*02:01', 'HLA-B*07:02', 'HLA-C*07:02']
    assert [p.observed_mhc_alleles for p in patient.input_provenance] == [
        ('HLA-A*02:01',), ('HLA-B*07:02',)]
    assert set(report['input_patient_id']) == {'fixture-patient'}
    assert {a for e in predictions for a in e.patient_alleles} == {
        'HLA-A*02:01', 'HLA-B*07:02'}


def test_observed_allele_outside_declared_genotype_is_rejected(tmp_path):
    args = manifest_args(tmp_path, mhc_alleles=['HLA-A*02:01'])
    with pytest.raises(ValueError, match='reported alleles outside declared'):
        load_external_ranked(args)


def test_fresh_prediction_cannot_ignore_declared_genotype(tmp_path, monkeypatch):
    import mhctools.cli
    args = manifest_args(tmp_path)
    args.external_predictions = 'fresh'
    args.mhc_predictor = ['random']
    monkeypatch.setattr(mhctools.cli, 'mhc_alleles_from_args', lambda _: ['HLA-A*01:01'])
    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', forbidden)
    with pytest.raises(ValueError, match='Fresh prediction alleles outside'):
        load_external_ranked(args)


@pytest.mark.parametrize('genome,field', [
    (SimpleNamespace(reference_name='GRCh37', release=75), 'reference_assembly'),
    (SimpleNamespace(reference_name='GRCh38', release=109), 'annotation'),
])
def test_configured_genome_must_match_input_declarations(tmp_path, genome, field):
    args = manifest_args(tmp_path, annotation='ensembl:110')
    args.genome = genome
    with pytest.raises(ValueError, match=f'conflicting {field}'):
        load_external_ranked(args)


@pytest.mark.parametrize('suffix', ['csv', 'tsv'])
def test_samples_and_rna_observations_survive_selection_slice_and_native_reload(tmp_path, suffix):
    inputs = []
    for sample, depth in [('first', 100), ('second', 300)]:
        frame = pd.read_csv(INPUTS[1][1], sep='\t').iloc[[0]].copy()
        frame['RNA Depth'] = depth
        frame['RNA VAF'] = .5
        frame['sample_id'] = sample
        frame['library_id'] = 'shared-library'
        frame['timepoint'] = 'before-treatment'
        path = tmp_path / f'{sample}.tsv'
        frame.to_csv(path, sep='\t', index=False)
        inputs.append(('pvacseq', str(path)))
    ranked, _, predictions, patient, _ = load_external_ranked(manifest_args(tmp_path, inputs))
    assert len(ranked) == 1
    assert patient.num_somatic_variants == 1
    assert len(predictions) == 2
    selected = ranked[0][1][0]
    assert selected.mutant_protein_fragment.n_rna_alt == 150
    assert {p.scope.sample_id for p in patient.input_provenance} == {'first', 'second'}
    original = [asdict(e.input_provenance) for e in predictions]
    evidence = [e.input_evidence for e in predictions]
    assert [e[0]['n_rna_alt'] for e in evidence] == [50, 150]
    assert all(e[0]['rna_evidence_method'] for e in evidence)
    # Scope follows the representative into actual peptide-window objects.
    assert selected.epitopes[0].input_provenance.scope.sample_id == 'second'
    path = tmp_path / f'predictions.{suffix}'
    save_predictions(predictions, path)
    recovered = load_predictions(path)
    assert [asdict(e.input_provenance) for e in recovered] == original
    assert [e.input_evidence for e in recovered] == evidence
    # The patient/run metadata has the same complete serializable record.
    assert patient.from_json(patient.to_json()).input_provenance == patient.input_provenance


def test_legacy_single_report_remains_usable_with_explicitly_unknown_scope():
    ranked, _, predictions, patient, _ = load_external_ranked(
        SimpleNamespace(input_pvacseq=INPUTS[1][1]))
    assert ranked and predictions
    assert patient.input_provenance[0].scope.patient_id is None
    assert patient.input_provenance[0].scope.mhc_alleles is None
    assert all(e.input_provenance.scope.sample_id is None for e in predictions)


@pytest.mark.parametrize('fmt,path', INPUTS)
def test_declared_reference_reaches_variant_identity_without_guessing_annotation(tmp_path, fmt, path):
    ranked, _, predictions, _, _ = load_external_ranked(
        manifest_args(tmp_path, [(fmt, path)], reference_assembly='GRCh37'))
    assert ranked
    assert all(source.reference_name == 'GRCh37' for source, _ in ranked)
    assert all(e.input_provenance.scope.annotation is None for e in predictions)


def test_same_report_in_different_samples_has_distinct_saved_occurrence_ids(tmp_path):
    one = load_external_ranked(manifest_args(tmp_path, INPUTS[:1], sample_id='first'))[2]
    two = load_external_ranked(manifest_args(tmp_path, INPUTS[:1], sample_id='second'))[2]
    assert {e.prediction_id for e in one}.isdisjoint(e.prediction_id for e in two)
    assert one[0].input_provenance.content_sha256 == two[0].input_provenance.content_sha256


def test_manifest_relative_paths_and_unchanged_direct_cli(tmp_path):
    path = tmp_path / 'input.yaml'
    path.write_text('schema: vaxrank.input_manifest.v1\ninputs:\n'
                    '  - format: lens\n    path: report.tsv\n')
    args = parse_vaxrank_args(['--input-manifest=' + str(path)])
    assert external_inputs(args) == [('lens', str(tmp_path / 'report.tsv'))]
    direct = parse_vaxrank_args(['--vcf', 'calls.vcf', '--bam', 'rna.bam',
                                 '--mhc-predictor', 'random', '--mhc-alleles', 'H2-Kb'])
    assert direct.vcf == ['calls.vcf'] and direct.bam == 'rna.bam'
    assert not getattr(direct, 'input_manifest', None)


def test_manifest_cannot_hide_extra_cli_inputs(tmp_path):
    args = manifest_args(tmp_path)
    args.input_lens = INPUTS[0][1]
    with pytest.raises(ValueError, match='cannot be combined'):
        load_external_ranked(args)


@pytest.mark.parametrize('fmt,path', INPUTS)
def test_compressed_reports_preserve_scope_and_evidence(tmp_path, fmt, path):
    compressed = tmp_path / f'{fmt}.tsv.gz'
    pd.read_csv(path, sep='\t').to_csv(compressed, sep='\t', index=False)
    _, _, original, _, _ = load_external_ranked(manifest_args(tmp_path, [(fmt, path)]))
    ranked, _, recovered, patient, _ = load_external_ranked(
        manifest_args(tmp_path, [(fmt, str(compressed))]))
    assert ranked
    assert patient.input_provenance[0].scope.patient_id == 'fixture-patient'
    assert [e.predictions for e in recovered] == [e.predictions for e in original]
    assert [e.input_evidence for e in recovered] == [e.input_evidence for e in original]


def test_invalid_genotype_diagnostic_names_the_field(tmp_path):
    with pytest.raises(ValueError, match='invalid mhc_alleles'):
        load_external_ranked(manifest_args(tmp_path, mhc_alleles=['not-an-allele']))


def test_invalid_assembly_fails_before_scoring_instead_of_dropping_variants(tmp_path, monkeypatch):
    import vaxrank.external_rescoring as module
    monkeypatch.setattr(module, 'attach_per_allele_scores', forbidden)
    with pytest.raises(ValueError, match='unsupported reference_assembly'):
        load_external_ranked(manifest_args(tmp_path, reference_assembly='not-an-assembly'))


def test_configured_genome_string_is_checked(tmp_path):
    args = manifest_args(tmp_path)
    args.genome = 'GRCh37'
    with pytest.raises(ValueError, match='conflicting reference_assembly'):
        load_external_ranked(args)


def test_nonhuman_alleles_are_normalized_without_hla_string_rules(tmp_path):
    frame = pd.read_csv(INPUTS[0][1], sep='\t').iloc[[0]].copy()
    frame['allele'] = 'H2-Kb'
    path = tmp_path / 'mouse.tsv'
    frame.to_csv(path, sep='\t', index=False)
    ranked, _, _, patient, _ = load_external_ranked(manifest_args(
        tmp_path, [('lens', str(path))], reference_assembly='GRCm38', mhc_alleles=['H-2-Kb']))
    assert ranked
    assert ranked[0][0].reference_name == 'GRCm38'
    assert patient.mhc_alleles == ['H2-K*b']


@pytest.mark.parametrize('extra', ['typo: 1\n', 'patient_id: A\npatient_id: B\n'])
def test_malformed_manifest_declarations_cannot_be_silently_ignored(tmp_path, extra):
    path = tmp_path / 'bad.yaml'
    path.write_text('schema: vaxrank.input_manifest.v1\n' + extra +
                    'inputs: [{format: lens, path: report.tsv}]\n')
    with pytest.raises(ValueError):
        read_input_manifest(path)


def test_cli_reports_and_resolved_provenance_preserve_scope(tmp_path):
    args = manifest_args(tmp_path, sample_id='tumor')
    native, report = tmp_path / 'predictions.tsv', tmp_path / 'report.csv'
    main(['--input-manifest', args.input_manifest, '--output-epitopes', str(native),
          '--output-csv', str(report), '--no-processing-aware-annotation'])
    recovered = load_predictions(native)
    assert {e.input_provenance.scope.sample_id for e in recovered} == {'tumor'}
    frame = pd.read_csv(report)
    assert set(frame['input_patient_id']) == {'fixture-patient'}
    assert set(frame['input_sample_id']) == {'tumor'}
    patient = load_external_ranked(args)[3]
    args.output_dir = str(tmp_path)
    write_run_summary(args, patient, 'external')
    saved = json.loads((tmp_path / 'input_provenance.json').read_text())
    assert saved['inputs'][0]['scope']['sample_id'] == 'tumor'
    assert saved['inputs'][0]['report_declarations'] == {}
    assert saved['inputs'][0]['manifest_declarations']['patient_id'] == 'fixture-patient'
    assert 'Declared genotype:' in (tmp_path / 'run_summary.txt').read_text()
