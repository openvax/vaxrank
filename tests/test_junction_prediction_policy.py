"""Candidate evidence and new junction queries have separate CLI policies."""

import json
from pathlib import Path

import pandas as pd
import pytest
from mhctools.pred import Prediction

from vaxrank.cli import entry_point as ep
from vaxrank.cli.arg_parser import parse_vaxrank_args
from vaxrank.epitope_io import load_predictions
from vaxrank.mrna import RNAConstruct
from vaxrank.native_serialization import from_native_json, to_native_json
from .input_scope_helpers import write_input_manifest


DATA = Path(__file__).parent / 'data/epitope_fixtures'
INPUTS = [('lens', str(DATA / 'lens_example.tsv')),
          ('pvacseq', str(DATA / 'pvacseq_example.tsv'))]
ALLELES = ['HLA-A*02:01', 'HLA-B*07:02']


def forbidden(*args, **kwargs):
    raise AssertionError('An unrequested prediction model was used')


class JunctionPredictor:
    def __init__(self, alleles):
        self.alleles = alleles
        self.calls = []

    def predict_peptides(self, peptides):
        self.calls.append(list(peptides))
        return [Prediction(
            peptide=p, allele=a, kind='pMHC_presentation', value=None, score=.1,
            percentile_rank=5.0, predictor_name='junction-test', predictor_version='2.0')
            for p in peptides for a in self.alleles]


@pytest.fixture
def junction_factory(monkeypatch):
    built = []

    def build(args):
        model = JunctionPredictor(ep.mhc_alleles_from_args(args))
        built.append((args, model))
        return model

    monkeypatch.setattr(ep, 'mhc_binding_predictor_from_args', build)
    return built


def run_external(tmp_path, name, extra=()):
    directory = tmp_path / name
    directory.mkdir()
    manifest = write_input_manifest(tmp_path / 'inputs.json', INPUTS,
                                    reference_assembly='GRCh37')
    ep.main([
        '--input-manifest', manifest, '--output-dir', str(directory),
        '--vaccine-type', 'mrna', '--ensembl-release', '75',
        '--pdf-backend', 'weasyprint', '--no-processing-aware-annotation',
        '--mrna-antigens-per-construct', '2', '--mrna-max-constructs', '1',
        '--mrna-signal-peptide', '', '--mrna-no-mitd', '--mrna-poly-a-length', '0',
        *extra])
    return directory


def test_historical_cli_scores_unchanged_and_only_junctions_are_predicted(
        tmp_path, monkeypatch, junction_factory):
    import mhctools.cli
    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', forbidden)
    written = []
    write = ep.write_mrna_outputs

    def capture(constructs, **kwargs):
        written.extend(constructs)
        return write(constructs, **kwargs)

    monkeypatch.setattr(ep, 'write_mrna_outputs', capture)
    baseline = run_external(tmp_path, 'baseline')
    assert not junction_factory
    selected = run_external(tmp_path, 'selected', ['--mrna-junction-predictor', 'junction-test'])
    assert len(junction_factory) == 1
    args, model = junction_factory[0]
    assert args.mhc_predictor == [['junction-test']]
    assert model.alleles == ALLELES
    assert model.calls
    original = load_predictions(baseline / 'candidate_predictions.tsv')
    recovered = load_predictions(selected / 'candidate_predictions.tsv')
    assert [e.to_json() for e in recovered] == [e.to_json() for e in original]
    pd.testing.assert_frame_equal(pd.read_csv(baseline / 'neoepitope_predictions.csv'),
                                  pd.read_csv(selected / 'neoepitope_predictions.csv'))
    manifest = json.loads((selected / 'manifest.json').read_text())[0]
    details = manifest['elements']['junction_swap']
    assert details['enabled'] and details['policy'] == 'explicit_junction_model'
    prediction = details['prediction']
    assert prediction['alleles'] == ALLELES
    assert prediction['peptide_lengths'] == [8, 9, 10, 11]
    assert prediction['models'] == [{'predictor_name': 'junction-test',
                                     'predictor_version': '2.0', 'kind': 'pMHC_presentation'}]
    assert [q['peptides'] for q in prediction['requests']] == model.calls
    assert prediction['scored_prediction_count'] == sum(map(len, model.calls)) * len(ALLELES)
    # Every recorded query is a chimeric junction window, not a candidate rescore.
    from vaxrank.junction_swap import junction_kmers
    from vaxrank.vaccine_library import get_linker
    antigens = manifest['antigens']
    for request in prediction['requests']:
        assert request['junction_index'] == 0
        expected = junction_kmers(antigens[0]['aa'],
                                  get_linker(request['linker_name']).amino_acids,
                                  antigens[1]['aa'], (8, 9, 10, 11))
        assert request['peptides'] == expected
    reloaded = from_native_json(to_native_json(written[-1]), RNAConstruct)
    assert reloaded.elements['junction_swap'] == details
    assert json.loads((baseline / 'manifest.json').read_text())[0]['elements'][
        'junction_swap']['policy'] == 'none'


def test_no_junction_cli_does_not_initialize_selected_model(tmp_path, monkeypatch):
    monkeypatch.setattr(ep, 'mhc_binding_predictor_from_args', forbidden)
    directory = run_external(tmp_path, 'disabled', [
        '--mrna-junction-predictor', 'unused', '--mrna-no-optimize-linkers'])
    metadata = json.loads((directory / 'manifest.json').read_text())[0]['elements']['junction_swap']
    assert metadata == {'enabled': False, 'note': 'junction prediction disabled',
                        'prediction': None, 'policy': 'none'}


def test_explicit_optimization_without_model_fails_through_cli(tmp_path, monkeypatch):
    monkeypatch.setattr(ep, 'mhc_binding_predictor_from_args', forbidden)
    with pytest.raises(ValueError, match='requires --mrna-junction-predictor'):
        run_external(tmp_path, 'missing-model', ['--mrna-optimize-linkers'])


def external_args(*extra):
    return parse_vaxrank_args(['--input-lens', INPUTS[0][1], *extra])


@pytest.mark.parametrize('input_flag', ['--input-lens', '--input-json-file'])
def test_report_coverage_does_not_supply_junction_genotype(input_flag, junction_factory):
    args = parse_vaxrank_args([input_flag, 'saved-input',
                               '--mrna-junction-predictor', 'junction-test'])
    args._inferred_mhc_alleles_from_lens = ALLELES
    with pytest.raises(ValueError, match='report allele coverage does not establish genotype'):
        ep.resolve_mhc_for_linker_optimizer(args)
    assert not junction_factory


@pytest.mark.parametrize('query', ['A0101', 'not-an-allele'])
def test_bad_junction_alleles_fail_before_model_load(query, junction_factory):
    args = external_args('--mrna-junction-predictor', 'junction-test',
                         '--mrna-junction-alleles', query)
    args._declared_mhc_alleles = ALLELES
    with pytest.raises(ValueError, match='outside declared genotype|invalid mhc_alleles'):
        ep.resolve_mhc_for_linker_optimizer(args)
    assert not junction_factory


def test_junction_subset_and_model_paths_are_independent(junction_factory):
    args = external_args('--external-predictions', 'fresh',
                         '--mhc-predictor', 'random', '--mhc-alleles', 'B0702',
                         '--mhc-predictor-path', '/candidate/tool',
                         '--mhc-predictor-models-path', '/candidate/models',
                         '--mrna-junction-predictor', 'junction-test',
                         '--mrna-junction-alleles', 'A0201',
                         '--mrna-junction-predictor-path', '/junction/tool')
    args._declared_mhc_alleles = ALLELES
    before = vars(args).copy()
    _, alleles = ep.resolve_mhc_for_linker_optimizer(args)
    assert alleles == ALLELES[:1]
    built_args = junction_factory[0][0]
    assert built_args.mhc_predictor_path == '/junction/tool'
    assert built_args.mhc_predictor_models_path is None
    assert vars(args) == before


def test_vcf_bam_auto_mode_reuses_explicit_candidate_model(junction_factory):
    args = parse_vaxrank_args(['--vcf', 'calls.vcf', '--bam', 'rna.bam',
                               '--mhc-predictor', 'random', '--mhc-alleles', 'H2-Kb'])
    _, alleles = ep.resolve_mhc_for_linker_optimizer(args)
    assert alleles == ['H2-K*b']
    assert junction_factory[0][0].mhc_predictor == args.mhc_predictor


@pytest.mark.parametrize('input_flag', ['--input-lens', '--input-json-file'])
def test_explicit_junction_model_uses_public_mhctools_factory(input_flag):
    args = parse_vaxrank_args([input_flag, 'saved-input',
                               '--mrna-junction-predictor', 'random',
                               '--mrna-junction-alleles', 'A0201'])
    predictor, alleles = ep.resolve_mhc_for_linker_optimizer(args)
    assert predictor.alleles == alleles == ALLELES[:1]
    predictions = predictor.predict_peptides(['SIINFEKLA'])
    assert predictions and all(p.allele == ALLELES[0] for p in predictions)


def test_yaml_and_cli_precedence_for_junction_policy(tmp_path, junction_factory):
    config = tmp_path / 'policy.yaml'
    config.write_text('mrna:\n  optimize_linkers: true\n'
                      '  junction_predictor: junction-test\n  junction_alleles: A0201\n')
    args = external_args('--config', str(config))
    _, alleles = ep.resolve_mhc_for_linker_optimizer(args)
    assert alleles == ALLELES[:1]
    args = external_args('--config', str(config), '--mrna-no-optimize-linkers')
    assert ep.resolve_mhc_for_linker_optimizer(args) == (None, None)
    assert len(junction_factory) == 1
    config.write_text(config.read_text().replace('optimize_linkers: true',
                                                'optimize_linkers: false'))
    assert ep.resolve_mhc_for_linker_optimizer(external_args('--config', str(config))) == (None, None)
    assert len(junction_factory) == 1
    args = external_args('--config', str(config), '--mrna-optimize-linkers')
    assert ep.resolve_mhc_for_linker_optimizer(args)[1] == ALLELES[:1]
    assert len(junction_factory) == 2
