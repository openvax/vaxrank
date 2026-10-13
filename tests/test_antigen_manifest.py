"""One experimental input contract for new sequences and saved evidence."""

from dataclasses import asdict, replace
import json

import pytest
from topiary import TopiaryPredictor

from vaxrank import AntigenDataset
from vaxrank.cli.entry_point import run_cli
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.external_rescoring import prepare_reports
from vaxrank.input_scope import InputScope, read_input_manifest

from . import test_antigen_input
from .test_epitope_dataset import forbidden

inputs = test_antigen_input.inputs
dependencies = test_antigen_input.dependencies


def write_manifest(tmp_path, entries, **scope):
    path = tmp_path / 'inputs.json'
    path.write_text(json.dumps(dict(schema='vaxrank.input_manifest.v1', inputs=[
        dict(format=fmt, path=name) for fmt, name in entries], **scope)))
    return path


def design_options(output, modality='peptide'):
    result = ['--output-dir', str(output), '--vaccine-type', modality,
        '--no-processing-aware-annotation', '--pdf-backend', 'weasyprint',
        '--vaccine-peptide-length', '10', '--config-text', 'epitopes.score_expr=1.0',
        '--config-text', 'vaccine_peptides.combined_score_expr=target_epitope_score']
    if modality == 'mrna':
        result += ['--mrna-antigens-per-construct', '6', '--mrna-max-constructs', '1',
                   '--mrna-no-optimize-linkers', '--mrna-codon-method', 'use_best_codon']
    return result


@pytest.mark.parametrize('modality', ['peptide', 'mrna'])
@pytest.mark.parametrize('saved_design', [False, True])
def test_manifest_mixes_new_sequences_and_saved_evidence(inputs, dependencies, tmp_path, monkeypatch,
                                                        modality, saved_design):
    model, genome = dependencies
    old = replace(inputs, antigens={k: a for k, a in inputs.antigens.items() if k in ('mutation', 'CTA')})
    saved = old.predict(TopiaryPredictor(models=[model]), genome=genome,
                        config=EpitopeConfig(score_expr='1.0'))
    saved.save(tmp_path / 'saved.tsv')
    if saved_design:
        old.save(tmp_path / 'old.json')
        run_cli(['--input-antigens', str(tmp_path / 'old.json'), '--mhc-predictor', 'random',
                 '--ensembl-release', '93', '--output-epitopes', str(tmp_path / 'saved.tsv')]
                + design_options(tmp_path / 'old-design', modality))
        saved = EpitopeDataset.load(tmp_path / 'saved.tsv')
        # Legacy mutation designs can supply settings without window optimization.
        saved.selection['run_configuration']['effective_vaccine_peptides']['window_selection'] = None
        from vaxrank.selection_policy import configuration_digest
        saved.selection['run_configuration']['effective_sha256'] = configuration_digest(
            saved.selection['run_configuration'])
        saved.save(tmp_path / 'saved.tsv')
    new = replace(inputs, antigens={k: a for k, a in inputs.antigens.items() if k not in old.antigens},
                  scope=InputScope())
    new.save(tmp_path / 'new.json')
    manifest = write_manifest(tmp_path, [('antigens', 'new.json'), ('epitopes', 'saved.tsv')],
                              **asdict(inputs.scope))
    calls = []
    def factory(args):
        calls.append(args.mhc_alleles)
        return [model]
    monkeypatch.setattr('mhctools.cli.predictors_from_args', factory)
    model.calls.clear()
    output = tmp_path / 'design'
    run_cli(['--input-manifest', str(manifest), '--mhc-predictor', 'random',
             '--ensembl-release', '93'] + design_options(output, modality))
    assert calls == [','.join(inputs.scope.mhc_alleles)]
    combined = EpitopeDataset.load(output / 'candidate_predictions.tsv')
    assert {a.kind for a in combined.antigens.values()} == test_antigen_input.ANTIGEN_KINDS
    assert len(AntigenDataset.from_epitope_dataset(combined).antigens) == 7
    originals = {e.prediction_id: e.predictions_flat() for e in saved.epitopes}
    assert {e.prediction_id: e.predictions_flat() for e in combined.epitopes
            if e.prediction_id in originals} == originals
    fresh = next(p for p in combined.provenance if p.input_id not in {p.input_id for p in saved.provenance})
    assert fresh.manifest_path == str(manifest)
    assert fresh.manifest_declarations == json.loads(json.dumps(asdict(inputs.scope)))
    assert fresh.scope == inputs.scope
    assert 'Held-out antigen records: 1' in (output / 'run_summary.txt').read_text()
    product = json.loads((output / 'manifest.json').read_text())
    sequences = ([p['sequence'] for p in product] if modality == 'peptide' else
                 [a['aa'] for p in product for a in p['antigens']])
    assert len(sequences) == 6
    assert set(sequences) == {a.amino_acids for a in inputs.antigens.values()}
    # No upstream file, annotation or model is needed to replay the combined design.
    (tmp_path / 'new.json').unlink()
    (tmp_path / 'saved.tsv').unlink()
    manifest.unlink()
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    monkeypatch.setattr('vaxrank.antigen_input.self_reference_matches', forbidden)
    replay = tmp_path / 'replay'
    run_cli(['--input-epitopes', str(output / 'candidate_predictions.tsv'),
             '--output-dir', str(replay), '--no-processing-aware-annotation', '--pdf-backend', 'weasyprint'])
    for path in output.glob('*.fasta'):
        assert path.read_bytes() == (replay / path.name).read_bytes()
    assert (output / 'manifest.json').read_bytes() == (replay / 'manifest.json').read_bytes()


@pytest.mark.parametrize('route', ['single', 'manifest', 'external'])
def test_input_routes_preserve_same_identity(inputs, dependencies, tmp_path, route):
    path, output = tmp_path / 'new.json', tmp_path / 'native.tsv'
    inputs.save(path)
    args = {'single': ['--input-antigens', str(path)],
            'manifest': ['--input-manifest', str(write_manifest(tmp_path, [('antigens', 'new.json')]))],
            'external': ['--external-input', 'antigens=' + str(path)]}[route]
    run_cli(args + ['--mhc-predictor', 'random', '--ensembl-release', '93',
                    '--output-epitopes', str(output), '--config-text', 'epitopes.score_expr=1.0'])
    actual = EpitopeDataset.load(output)
    model, genome = dependencies
    expected = inputs.predict(TopiaryPredictor(models=[model]), genome=genome, path=path)
    assert {e.prediction_id for e in actual.epitopes} == {e.prediction_id for e in expected.epitopes}
    assert AntigenDataset.from_epitope_dataset(actual).antigens == inputs.antigens
    assert set(actual.selection['antigen_inputs']) == set(expected.selection['antigen_inputs'])
    assert actual.provenance[0].content_sha256 == expected.provenance[0].content_sha256


@pytest.mark.parametrize('field,value', [
    ('patient_id', 'another-patient'), ('reference_assembly', 'GRCh37'),
    ('annotation', 'ensembl:75'), ('mhc_alleles', ('HLA-A*01:01',)),
])
def test_conflicting_source_after_raw_input_fails_before_models(inputs, dependencies, tmp_path, monkeypatch, field, value):
    model, genome = dependencies
    saved = inputs.predict(TopiaryPredictor(models=[model]), genome=genome)
    saved.save(tmp_path / 'saved.tsv')
    replace(inputs, scope=replace(inputs.scope, **{field: value})).save(tmp_path / 'new.json')
    manifest = write_manifest(tmp_path, [('antigens', 'new.json'), ('epitopes', 'saved.tsv')])
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    with pytest.raises(ValueError, match='conflicting ' + field):
        run_cli(['--input-manifest', str(manifest), '--mhc-predictor', 'random',
                 '--ensembl-release', '93', '--output-epitopes', str(tmp_path / 'output.tsv')])
    assert not (tmp_path / 'output.tsv').exists()


@pytest.mark.parametrize('duplicate', ['copied-file', 'reformatted-file', 'native-export'])
def test_duplicate_inputs_fail_before_models(inputs, dependencies, tmp_path, monkeypatch, duplicate):
    inputs.save(tmp_path / 'new.json')
    if duplicate in ('copied-file', 'reformatted-file'):
        inputs.save(tmp_path / 'copy.json')
        if duplicate == 'reformatted-file':
            (tmp_path / 'copy.json').write_text(json.dumps(inputs.to_dict()))
        other = ('antigens', 'copy.json')
    else:
        model, genome = dependencies
        inputs.predict(TopiaryPredictor(models=[model]), genome=genome).save(tmp_path / 'saved.tsv')
        other = ('epitopes', 'saved.tsv')
    manifest = write_manifest(tmp_path, [('antigens', 'new.json'), other])
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    with pytest.raises(ValueError, match='more than once|overlapping candidate identities'):
        run_cli(['--input-manifest', str(manifest), '--mhc-predictor', 'random',
                 '--ensembl-release', '93', '--output-epitopes', str(tmp_path / 'output.tsv')])


def test_library_predicts_multiple_sources_with_one_factory(inputs, dependencies, tmp_path):
    model, genome = dependencies
    entries = []
    for kind in ('mutation', 'CTA', 'viral'):
        replace(inputs, antigens={kind: inputs.antigens[kind]}).save(tmp_path / (kind + '.json'))
        entries.append(('antigens', kind + '.json'))
    specs = read_input_manifest(write_manifest(tmp_path, entries))
    called = []
    def factory(scope):
        called.append(scope)
        return TopiaryPredictor(models=[model])
    _, frame, predictions = prepare_reports(
        [(fmt, path) for fmt, path, _ in specs], EpitopeConfig(score_expr='1.0'),
        scopes=[scope for _, _, scope in specs], genome=genome, antigen_predictor_factory=factory)
    assert len(called) == 1
    assert called[0].patient_id == 'patient'
    assert len(predictions) == 6
    assert {a.kind for a in frame.attrs['epitope_dataset'].antigens.values()} == {'mutation', 'CTA', 'viral'}


@pytest.mark.parametrize('extra,match', [
    ([], 'requires --mhc-predictor'),
    (['--mhc-predictor', 'random'], 'requires --mhc-predictor'),
    (['--mhc-predictor', 'random', '--ensembl-release', '93', '--external-predictions', 'additive'],
     'external prediction modes do not apply'),
])
def test_manifest_requires_prediction_contract(inputs, tmp_path, monkeypatch, extra, match):
    inputs.save(tmp_path / 'new.json')
    manifest = write_manifest(tmp_path, [('antigens', 'new.json')])
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    with pytest.raises(ValueError, match=match):
        run_cli(['--input-manifest', str(manifest), '--output-epitopes', str(tmp_path / 'output.tsv')] + extra)
