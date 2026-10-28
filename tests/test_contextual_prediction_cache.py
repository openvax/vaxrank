"""Compose the CLI, Topiary contextual cache, and native replay contracts."""

from dataclasses import replace

import pytest
import pandas as pd
from mhctools.pred import PeptideResult
from topiary import CachedPredictor, TopiaryResult, predict_peptide_occurrences

from vaxrank.cli.entry_point import run_cli
from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.external_rescoring import rescore_candidates

from .input_scope_helpers import write_input_manifest
from .test_external_rescoring import ALLELES, DATA, ContextPredictor, candidate, config


def contextual_cache(path, *, unknown_primary=False):
    model = ContextPredictor()
    model.alleles = [*ALLELES, 'HLA-A*03:01']
    frame = predict_peptide_occurrences([
        dict(prediction_id='cached-primary', peptide='SIINFEKLA',
             n_flank='AA', c_flank='KK'),
        dict(prediction_id='cached-wt', peptide='SIINFEKLV',
             n_flank='', c_flank=''),
    ], model)
    if unknown_primary:
        mask = frame.peptide.eq('SIINFEKLA')
        frame.loc[mask, ['n_flank', 'c_flank']] = None
        frame.loc[mask, 'prediction_flanks_supplied'] = False
    TopiaryResult(frame).to_tsv(path)
    return path


def forbidden(*args, **kwargs):
    raise AssertionError('Native replay must not construct or call a predictor')


@pytest.mark.parametrize('mode', ['fresh', 'additive'])
def test_cli_contextual_cache_preserves_occurrences_and_native_replay(
        tmp_path, monkeypatch, mode):
    originals = rescore_candidates(
        [replace(candidate(name, flank), n_flank=flank, c_flank='KK')
         for name, flank in [('a', 'AA'), ('b', 'CC'), ('c', 'AA')]],
        [ContextPredictor()], ALLELES)
    source = tmp_path / 'source.tsv'
    if mode == 'fresh':
        from vaxrank.epitope_io import read_lens_report
        frame = pd.read_csv(DATA / 'lens_example.tsv', sep='\t').iloc[:3].copy()
        frame['peptide'] = 'SIINFEKLA'
        frame['pep_context'] = [epitope.source_sequence for epitope in originals]
        frame.to_csv(source, sep='\t', index=False)
        originals = read_lens_report(source, score_epitopes=False).epitopes
        options = ['--config-text', 'epitopes.score_expr=affinity[unified].score']
        fmt = 'lens'
    else:
        EpitopeDataset.from_predictions(originals, config=config()).save(source)
        options = ['--external-prediction-prefix', 'cached']
        fmt = 'epitopes'
    before = source.read_bytes()
    manifest = write_input_manifest(tmp_path / 'inputs.json', [(fmt, source)])
    cache = contextual_cache(tmp_path / 'cache.tsv')
    cache_bytes = cache.read_bytes()
    fallback = ContextPredictor()
    monkeypatch.setattr('mhctools.cli.predictors_from_args', lambda args: [fallback])
    output = tmp_path / 'predictions.tsv'
    run_cli([
        '--input-manifest', manifest, '--external-predictions', mode,
        '--mhc-predictor', 'random', '--mhc-alleles', ','.join(ALLELES),
        '--prediction-cache', str(cache), '--output-epitopes', str(output),
        '--no-processing-aware-annotation', *options,
    ])
    # Primary AA and known-terminal WT are hits; only the distinct CC context
    # is inferred. Extra cached alleles do not expand this patient's genotype.
    assert [request for batch in fallback.calls for request in batch] == [
        ('SIINFEKLA', 'CC', 'KK')]
    result = EpitopeDataset.load(output)
    assert len(result.epitopes) == 3
    assert [epitope.source_sequence for epitope in result.epitopes] == [
        epitope.source_sequence for epitope in originals]
    assert all(set(epitope.patient_alleles) == set(ALLELES) for epitope in result.epitopes)
    assert source.read_bytes() == before
    assert cache.read_bytes() == cache_bytes
    if mode == 'fresh':
        values = [epitope.predictions_flat()[0].value for epitope in result.epitopes]
        assert values[0] == values[2] != values[1]
    else:
        assert [epitope.predictions for epitope in result.epitopes] == [
            epitope.predictions for epitope in originals]
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    replay = tmp_path / 'replay.tsv'
    run_cli(['--input-epitopes', str(output), '--output-epitopes', str(replay),
             '--no-processing-aware-annotation'])
    restored = EpitopeDataset.load(replay)
    assert restored.epitopes == result.epitopes
    assert restored.provenance == result.provenance


def test_unknown_cached_flanks_do_not_cover_known_occurrence(tmp_path):
    fallback = ContextPredictor()
    cache = CachedPredictor.from_topiary_output(
        contextual_cache(tmp_path / 'cache.tsv', unknown_primary=True), fallback=fallback)
    cache.alleles = ALLELES
    fresh, = rescore_candidates([candidate('a')], [cache], ALLELES)
    assert fresh.wt.predictions_flat()
    assert [request for batch in fallback.calls for request in batch] == [
        ('SIINFEKLA', 'AA', 'KK')]


@pytest.mark.parametrize('failure', ['version', 'allele', 'kind'])
def test_contextual_cache_rejects_incompatible_or_incomplete_fallback(tmp_path, failure):
    class Invalid(ContextPredictor):
        def kind_support(self):
            support = super().kind_support()
            if failure == 'kind':
                support['pMHC_presentation'] = {'mhc_dependence': 'single_allele'}
            return support

        def predict_with_flanks(self, *args):
            results = super().predict_with_flanks(*args)
            if failure == 'version':
                return [PeptideResult(tuple(replace(p, predictor_version='2.0')
                                            for p in result.preds)) for result in results]
            if failure == 'allele':
                return [PeptideResult(result.preds[:1]) for result in results]
            return results

    cache = CachedPredictor.from_topiary_output(
        contextual_cache(tmp_path / 'cache.tsv'), fallback=Invalid())
    cache.alleles = ALLELES
    with pytest.raises(ValueError, match='version|complete .*prediction'):
        rescore_candidates([candidate('miss', 'CC')], [cache], ALLELES)
