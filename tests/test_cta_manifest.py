"""Expression admission joins new sequences and saved evidence in one run."""

from dataclasses import asdict, replace
import json

import pandas as pd
import pytest
from topiary import TopiaryPredictor
from mhctools.cli import predictors_from_args

from vaxrank import AntigenDataset
import vaxrank.cta_input as cta_input
from vaxrank.cli.entry_point import run_cli
from vaxrank.cta_expression import CTAExpressionResult
from vaxrank.cta_input import CTAExpressionSource
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.external_rescoring import prepare_reports
from vaxrank.input_scope import InputScope, read_input_manifest
from vaxrank.native_serialization import from_native_json

from . import test_antigen_input, test_cta_input
from .test_cta_expression import PRAME, MAGEA3, MAGEA4, NOT_CTA, PRAME_TRANSCRIPT, UNKNOWN_TRANSCRIPT
from .test_epitope_dataset import forbidden
from .test_external_rescoring import ALLELES

inputs = test_antigen_input.inputs
genome = test_cta_input.genome
cli_dependencies = test_cta_input.cli_dependencies


def options(**changes):
    return dict(measurement_level='gene', id_column='feature_id', value_column='patient_tpm',
                expression_unit='TPM', evidence_source='Salmon', evidence_version='1.10.3',
                assay='bulk RNA-seq', min_expression=2, **changes)


def entry(path='expression.tsv', sample='patient-001', **changes):
    return dict(format='cta_expression', path=path, sample_id=sample, cta_expression=options(**changes))


def manifest(tmp_path, entries, **scope):
    path = tmp_path / 'inputs.json'
    path.write_text(json.dumps(dict(schema='vaxrank.input_manifest.v1', patient_id='patient-001',
        reference_assembly='GRCh38', annotation='ensembl:93', mhc_alleles=ALLELES, inputs=entries, **scope)))
    return path


def expression(tmp_path):
    path = tmp_path / 'expression.tsv'
    pd.DataFrame([(PRAME, 3), (MAGEA3, 50), (MAGEA4, 4), (NOT_CTA, 10)],
                 columns=['feature_id', 'patient_tpm']).to_csv(path, sep='\t', index=False)
    return path


def cli(path, output):
    return ['--input-manifest', str(path), '--mhc-predictor', 'random', '--ensembl-release', '93',
            '--output-epitopes', str(output), '--config-text', 'epitopes.score_expr=1.0',
            '--no-processing-aware-annotation']


def test_manifest_and_alias_produce_same_cta_evidence(tmp_path, cli_dependencies):
    path = expression(tmp_path)
    alias = tmp_path / 'alias.tsv'
    run_cli(test_cta_input.argv(path, tmp_path / 'alias') + ['--output-epitopes', str(alias),
        '--vaccine-type', 'peptide', '--config-text', 'epitopes.score_expr=1.0'])
    joined = tmp_path / 'joined.tsv'
    run_cli(cli(manifest(tmp_path, [entry()]), joined))
    first, second = EpitopeDataset.load(alias), EpitopeDataset.load(joined)
    assert {e.prediction_id for e in first.epitopes} == {e.prediction_id for e in second.epitopes}
    assert first.selection['cta_expression_admission'] == second.selection['cta_expression_admission']
    assert first.antigens == second.antigens
    assert second.provenance[0].manifest_path == str(tmp_path / 'inputs.json')


@pytest.mark.parametrize('modality', ['peptide', 'mrna'])
def test_cta_viral_and_saved_mutation_share_design_and_replay(inputs, tmp_path, cli_dependencies, genome, monkeypatch, modality):
    expression(tmp_path)
    inputs.scope = replace(inputs.scope, patient_id='patient-001')
    monkeypatch.setattr('vaxrank.antigen_input.self_reference_matches', cta_input.self_reference_matches)
    viral = replace(inputs, antigens={'viral': inputs.antigens['viral']})
    viral.save(tmp_path / 'viral.json')
    mutation = replace(inputs, antigens={'mutation': inputs.antigens['mutation']})
    old = mutation.predict(TopiaryPredictor(models=[cli_dependencies]), genome=genome,
                           config=EpitopeConfig(score_expr='1.0'))
    old.save(tmp_path / 'mutation.tsv')
    path = manifest(tmp_path, [entry(), dict(format='antigens', path='viral.json'),
                              dict(format='epitopes', path='mutation.tsv')])
    factories = []
    def factory(args):
        factories.append(args.mhc_alleles)
        return [cli_dependencies]
    monkeypatch.setattr('mhctools.cli.predictors_from_args', factory)
    output = tmp_path / 'design'
    args = cli(path, output / 'native.tsv') + ['--output-dir', str(output), '--vaccine-type', modality,
        '--vaccine-peptide-length', '10', '--pdf-backend', 'weasyprint',
        '--config-text', 'vaccine_peptides.combined_score_expr=target_epitope_score']
    if modality == 'mrna':
        args += ['--mrna-antigens-per-construct', '4', '--mrna-max-constructs', '1',
                 '--mrna-no-optimize-linkers', '--mrna-codon-method', 'use_best_codon']
    run_cli(args)
    assert len(factories) == 1
    result = EpitopeDataset.load(output / 'native.tsv')
    assert {a.kind for a in result.antigens.values()} == {'CTA', 'mutation', 'viral'}
    original = {e.prediction_id: e.predictions_flat() for e in old.epitopes}
    assert {e.prediction_id: e.predictions_flat() for e in result.epitopes if e.prediction_id in original} == original
    assert len(result.selection['source_evidence']) == 1
    source = next(iter(result.selection['source_evidence'].values()))
    admission = from_native_json(source['cta_expression_admission'], CTAExpressionResult)
    assert len(admission.decisions) == 4
    assert len(admission.admitted_antigens) == 2
    assert CTAExpressionResult.load(output / 'cta_admission.json') == admission
    assert 'cta_product_audit' not in result.selection
    product = json.loads((output / 'manifest.json').read_text())
    count = len(product) if modality == 'peptide' else sum(len(p['antigens']) for p in product)
    assert count == 4
    for name in ('expression.tsv', 'viral.json', 'mutation.tsv', 'inputs.json'):
        (tmp_path / name).unlink()
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    monkeypatch.setattr('vaxrank.cta_input.admit_cta_expression', forbidden)
    monkeypatch.setattr('vaxrank.cta_input.load_hitlist_bundle', forbidden)
    replay = tmp_path / 'replay'
    run_cli(['--input-epitopes', str(output / 'native.tsv'), '--output-dir', str(replay),
             '--pdf-backend', 'weasyprint', '--no-processing-aware-annotation'])
    for fasta in output.glob('*.fasta'):
        assert fasta.read_bytes() == (replay / fasta.name).read_bytes()
    assert (output / 'manifest.json').read_bytes() == (replay / 'manifest.json').read_bytes()


def test_one_matrix_can_supply_distinct_samples_and_contracts(tmp_path, cli_dependencies):
    pd.DataFrame([(PRAME, 3, 0), (MAGEA4, 0, 7)],
                 columns=['feature_id', 'patient_tpm', 'second_tpm']).to_csv(tmp_path / 'expression.tsv', sep='\t', index=False)
    second = entry(sample='second')
    second['cta_expression'].update(value_column='second_tpm', min_expression=5)
    output = tmp_path / 'design'
    run_cli(cli(manifest(tmp_path, [entry(sample='first'), second]), output / 'native.tsv') +
            ['--output-dir', str(output), '--vaccine-type', 'peptide', '--pdf-backend', 'weasyprint'])
    result = EpitopeDataset.load(output / 'native.tsv')
    assert {p.scope.sample_id for p in result.provenance} == {'first', 'second'}
    admissions = [CTAExpressionResult.load(path) for path in (output / 'cta_admissions').glob('*.json')]
    assert len(admissions) == 2
    assert {a.input_contract.sample_id: {v.gene_id for v in a.admitted_antigens} for a in admissions} == {
        'first': {PRAME}, 'second': {MAGEA4}}
    assert len(result.selection['source_evidence']) == 2


def test_entirely_held_out_cta_source_survives_with_saved_evidence(tmp_path, cli_dependencies, genome, monkeypatch):
    expression(tmp_path)
    good = CTAExpressionSource.load(tmp_path / 'expression.tsv', options=options(),
        scope=InputScope(patient_id='patient-001', sample_id='patient-001', mhc_alleles=tuple(ALLELES)), genome=genome)
    from vaxrank.cta_input import predict_cta_source
    predict_cta_source(good, TopiaryPredictor(models=[cli_dependencies]), genome=genome,
                       epitope_config=EpitopeConfig(score_expr='1.0')).dataset.save(tmp_path / 'saved.tsv')
    pd.DataFrame([(UNKNOWN_TRANSCRIPT, 5), (PRAME_TRANSCRIPT, 0)],
                 columns=['feature_id', 'patient_tpm']).to_csv(tmp_path / 'held.tsv', sep='\t', index=False)
    held = entry('held.tsv', sample='held-sample')
    held['cta_expression']['measurement_level'] = 'transcript'
    path = manifest(tmp_path, [held, dict(format='epitopes', path='saved.tsv')])
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    output = tmp_path / 'combined.tsv'
    run_cli(cli(path, output))
    result = EpitopeDataset.load(output)
    assert len(result.selection['source_evidence']) == 2
    retained = [from_native_json(r['cta_expression_admission'], CTAExpressionResult)
                for r in result.selection['source_evidence'].values()]
    admission = next(r for r in retained if r.input_contract.sample_id == 'held-sample')
    assert not admission.admitted_antigens
    assert len(admission.decisions) == 2
    assert result.epitopes


@pytest.mark.parametrize('field,value', [('patient_id', 'other'), ('mhc_alleles', ['HLA-A*01:01']),
                                       ('reference_assembly', 'GRCh37'), ('annotation', 'ensembl:75')])
def test_conflicts_fail_before_models_and_admission_outputs(inputs, tmp_path, cli_dependencies, monkeypatch, field, value):
    expression(tmp_path)
    inputs.scope = replace(inputs.scope, patient_id='patient-001', **({field: value} if field != 'patient_id' else {}))
    if field == 'patient_id':
        inputs.scope = replace(inputs.scope, patient_id=value)
    inputs.save(tmp_path / 'new.json')
    path = manifest(tmp_path, [entry(), dict(format='antigens', path='new.json')])
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    with pytest.raises(ValueError, match='conflicting ' + field):
        run_cli(cli(path, tmp_path / 'native.tsv') + ['--output-cta-admission', str(tmp_path / 'admission.json')])
    assert not (tmp_path / 'admission.json').exists()


@pytest.mark.parametrize('change,match', [
    (lambda e: e['cta_expression'].pop('min_expression'), 'min_expression'),
    (lambda e: e['cta_expression'].update(min_expression=True), 'min_expression'),
    (lambda e: e['cta_expression'].update(min_expression=0), 'greater than zero'),
    (lambda e: e['cta_expression'].update(unknown=True), 'unknown cta_expression'),
    (lambda e: e['cta_expression'].update(exclusion_policy={'gene_patterns': 'MAGE*'}), 'lists'),
    (lambda e: e.pop('sample_id'), 'sample_id'),
    (lambda e: e['cta_expression'].update(hitlist_evidence_bundle=None), 'nonempty path'),
])
def test_invalid_contract_fails_during_manifest_loading(tmp_path, change, match):
    source = entry()
    change(source)
    with pytest.raises(ValueError, match=match):
        read_input_manifest(manifest(tmp_path, [source]))


def test_bundle_paths_resolve_relative_to_manifest(tmp_path):
    source = entry(hitlist_evidence_bundle='evidence/bundle')
    specs = read_input_manifest(manifest(tmp_path, [source]))
    assert specs[0][2]['cta_expression']['hitlist_evidence_bundle'] == str(tmp_path / 'evidence/bundle')


def test_library_source_uses_same_contract(tmp_path, cli_dependencies, genome):
    expression(tmp_path)
    specs = read_input_manifest(manifest(tmp_path, [entry()]))
    _, frame, _ = prepare_reports([(f, p) for f, p, _ in specs], EpitopeConfig(score_expr='1.0'),
        scopes=[s for _, _, s in specs], genome=genome,
        antigen_predictor_factory=lambda scope: TopiaryPredictor(models=[cli_dependencies]))
    result = frame.attrs['epitope_dataset']
    assert len(result.selection['source_evidence']) == 1
    exported = AntigenDataset.from_epitope_dataset(result)
    assert {a.gene_id for a in exported.antigens.values()} == {PRAME, MAGEA4}
    admission = from_native_json(result.selection['cta_expression_admission'], CTAExpressionResult)
    assert {d.input_identifier for d in admission.decisions} == {PRAME, MAGEA3, MAGEA4, NOT_CTA}
    assert asdict(result.provenance[0].scope)['patient_id'] == 'patient-001'


@pytest.mark.parametrize('missing', ['--ensembl-release', '--mhc-predictor'])
def test_manifest_requires_annotation_and_predictor(tmp_path, cli_dependencies, monkeypatch, missing):
    expression(tmp_path)
    args = cli(manifest(tmp_path, [entry()]), tmp_path / 'native.tsv')
    index = args.index(missing)
    del args[index:index + 2]
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    with pytest.raises(ValueError, match='CTA expression inputs require') as error:
        run_cli(args)
    assert missing in str(error.value)


@pytest.mark.parametrize('exported', [False, True])
def test_duplicate_admission_rejected_before_prediction(tmp_path, cli_dependencies, monkeypatch, exported):
    expression(tmp_path)
    path = manifest(tmp_path, [entry()])
    if exported:
        run_cli(cli(path, tmp_path / 'native.tsv'))
        extra = dict(format='epitopes', path='native.tsv')
    else:
        extra = entry()
    path = manifest(tmp_path, [entry(), extra])
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    with pytest.raises(ValueError, match='overlapping candidate identities|same external input'):
        run_cli(cli(path, tmp_path / 'duplicate.tsv'))
    assert not (tmp_path / 'duplicate.tsv').exists()


def test_multiple_sources_reject_single_admission_output(tmp_path, cli_dependencies, monkeypatch):
    expression(tmp_path)
    path = manifest(tmp_path, [entry(), entry(sample='second')])
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    admission = tmp_path / 'admission.json'
    with pytest.raises(ValueError, match='requires one CTA expression input'):
        run_cli(cli(path, tmp_path / 'native.tsv') + ['--output-cta-admission', str(admission)])
    assert not admission.exists()


def test_all_held_out_sources_write_diagnostics_without_models(tmp_path, cli_dependencies, monkeypatch):
    pd.DataFrame([(PRAME, 0)], columns=['feature_id', 'patient_tpm']).to_csv(
        tmp_path / 'expression.tsv', sep='\t', index=False)
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    output = tmp_path / 'design'
    with pytest.raises(ValueError, match='No admitted CTA proteins'):
        run_cli(cli(manifest(tmp_path, [entry()]), output / 'native.tsv') + ['--output-dir', str(output)])
    admission = CTAExpressionResult.load(output / 'cta_admission.json')
    assert not admission.admitted_antigens
    assert len(admission.decisions) == 1
    assert not (output / 'native.tsv').exists()


def test_manifest_uses_installed_mhctools_factory(tmp_path, cli_dependencies, monkeypatch):
    expression(tmp_path)
    monkeypatch.setattr('mhctools.cli.predictors_from_args', predictors_from_args)
    output = tmp_path / 'native.tsv'
    run_cli(cli(manifest(tmp_path, [entry()]), output))
    result = EpitopeDataset.load(output)
    assert result.epitopes
    assert {a.gene_id for a in result.antigens.values()} == {PRAME, MAGEA4}
