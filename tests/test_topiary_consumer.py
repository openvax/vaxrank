"""Released Topiary APIs preserve evidence while explicit features change vaccines."""

from dataclasses import replace
import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest
from mhctools.pred import Prediction, PeptideResult
from topiary import read_exacto

from vaxrank.cli.entry_point import run_cli
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dataset import EpitopeDataset, ANTIGEN_COLUMN, dataset_ranking_result
from vaxrank.epitope_dsl import epitopes_for_ranking
from vaxrank.external_input import ExternalConstructOptions
from vaxrank.external_rescoring import (prepare_reports, read_reports, rescore_candidates,
                                       InputTablePredictor)
from vaxrank.input_scope import read_input_manifest
from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs
from vaxrank.native_serialization import to_native_json
from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
from vaxrank.vaccine_antigen import (
    VaccineAntigen, TargetableMask, AminoAcidInterval, TumorSpecificityAttestation,
)
from .input_scope_helpers import FIXTURE_SCOPE
from .test_epitope_dataset import evidence, forbidden
from .test_external_rescoring import ALLELES, INPUTS, ContextPredictor, candidate


EXACTO = Path(__file__).parent / 'data' / 'epitope_fixtures' / 'exacto'


class RerankingPredictor(ContextPredictor):
    def predict_with_flanks(self, peptides, n_flanks, c_flanks):
        self.calls.append(list(zip(peptides, n_flanks, c_flanks)))
        return [PeptideResult(tuple(Prediction(
            kind='pMHC_affinity', peptide=p, allele=a,
            value=25. if p == 'GILGFVFTL' else 2000., score=.95 if p == 'GILGFVFTL' else .1,
            predictor_name='unified', predictor_version='1.0', n_flank=n, c_flank=c)
            for a in self.alleles)) for p, n, c in zip(peptides, n_flanks, c_flanks)]


def admitted_table(path, kind):
    antigens = [VaccineAntigen(
        kind=kind, amino_acids=p, source_identifier=p,
        targetable_mask=TargetableMask((AminoAcidInterval(0, len(p)),)),
        tumor_specificity=TumorSpecificityAttestation(
            status='admitted', evidence_kind='synthetic_fixture', evidence_source='test_topiary_consumer',
            patient_specific=True, rationale_code='test_only')) for p in ('SIINFEKL', 'GILGFVFTL')]
    evidence(n_flank=['', ''], c_flank=['', ''],
             **{ANTIGEN_COLUMN: [to_native_json(a) for a in antigens]}).to_tsv(path)
    return path


def constructs(report, config):
    ranked = dataset_ranking_result(report, epitopes_for_ranking(report.epitopes, config),
                                   options=ExternalConstructOptions()).ranked
    peptides = assemble_peptide_constructs(ranked, options=PeptideConstructConfig(
        mode='minimal_epitope', min_antigen_length_aa=8))
    mrna = assemble_mrna_constructs(ranked, options=RNAConstructConfig(
        antigen_content='minimal_epitope', signal_peptide=None, include_mitd=False, optimize_linkers=False))
    return [p.sequence for p in peptides], [r.full_nt for r in mrna]


@pytest.mark.parametrize('kind', ['mutation', 'fusion', 'splice', 'CTA', 'ERV', 'viral'])
def test_additive_features_change_constructs_only_when_dsl_selects_them(kind, tmp_path, monkeypatch):
    source = admitted_table(tmp_path / 'source.tsv', kind)
    old = EpitopeConfig(score_expr='1 / affinity.value', min_epitope_score=.01)
    reports, _, _ = prepare_reports([('topiary', source)], old)
    expected = constructs(reports[0], old)
    assert expected[0] == ['SIINFEKL'] and expected[1]
    _, added_frame, added = prepare_reports([('topiary', source)], old, mode='additive',
        models=[RerankingPredictor()], prediction_prefix='fresh')
    historical = reports[0].dataset.result.df
    saved = added_frame.attrs['epitope_dataset']
    pd.testing.assert_frame_equal(saved.result.df[historical.columns], historical)
    assert [e.per_allele_scores for e in added] == [e.per_allele_scores for e in reports[0].epitopes]
    assert saved.result.extra['candidate_rescoring']['fresh']['features']
    native = tmp_path / 'native.tsv'
    saved.save(native)
    new = EpitopeConfig(score_expr='fresh__unified__pMHC_affinity__score', min_epitope_score=.5)
    changed, changed_frame, _ = prepare_reports([('epitopes', native)], new)
    actual = constructs(changed[0], new)
    assert actual[0] == ['GILGFVFTL'] and actual[1] != expected[1]
    again = tmp_path / 'again.tsv'
    changed_frame.attrs['epitope_dataset'].save(again)
    import mhctools.cli
    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', forbidden)
    replay, _, _ = prepare_reports([('epitopes', again)])
    assert constructs(replay[0], new) == actual
    with pytest.raises(ValueError, match='already exists'):
        prepare_reports([('epitopes', native)], old, mode='additive',
                        models=[RerankingPredictor()], prediction_prefix='fresh')


def test_exact_occurrences_preserve_multiple_windows_and_unknown_flanks():
    originals = [candidate('one-source', 'AA', comparator=False), candidate('one-source', 'CC', comparator=False)]
    originals[1] = replace(originals[1], offset=14, source_sequence='X' * 12 + originals[1].source_sequence)
    model = ContextPredictor()
    fresh = rescore_candidates(originals, [model], ALLELES)
    assert [e.prediction_group_key for e in fresh] == [e.prediction_group_key for e in originals]
    assert fresh[0].predictions_flat()[0].value != fresh[1].predictions_flat()[0].value
    assert len(InputTablePredictor(fresh).predict_candidates(fresh)) == 2
    unknown = replace(originals[0], source_sequence='', n_flank=None, c_flank=None)
    with pytest.raises(ValueError, match='both flanks'):
        rescore_candidates([unknown], [ContextPredictor()], ALLELES)
    peptide_only, = rescore_candidates([unknown], [ContextPredictor()], ALLELES, use_flanks=False)
    assert peptide_only.n_flank is None and peptide_only.c_flank is None


def test_cli_additive_feature_selection_and_construct_reload(tmp_path, monkeypatch):
    import mhctools.cli
    source = admitted_table(tmp_path / 'source.tsv', 'viral')
    model = RerankingPredictor()
    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', lambda args: [model])
    first = tmp_path / 'first'
    native = tmp_path / 'native.tsv'
    run_cli(['--input-topiary', str(source), '--external-predictions', 'additive',
             '--external-prediction-prefix', 'chosen', '--mhc-predictor', 'random',
             '--mhc-alleles', ','.join(ALLELES), '--output-dir', str(first),
             '--output-epitopes', str(native), '--no-processing-aware-annotation',
             '--config-value', 'epitopes.score_expr=chosen__unified__pMHC_affinity__score',
             '--min-epitope-score', '0.5'])
    assert model.calls
    saved = EpitopeDataset.load(native)
    assert [e.sequence for e in epitopes_for_ranking(saved.epitopes, saved.config)] == ['GILGFVFTL']
    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', forbidden)
    second = tmp_path / 'second'
    run_cli(['--input-epitopes', str(native), '--output-dir', str(second),
             '--no-processing-aware-annotation'])
    fasta_files = list(first.rglob('*.fasta'))
    assert len(fasta_files) >= 2 and list(first.rglob('*.pdf'))
    for path in fasta_files:
        assert path.read_bytes() == (second / path.relative_to(first)).read_bytes()


@pytest.mark.parametrize('complete', [True, False])
def test_cli_external_cache_hits_and_version_checked_fallback(tmp_path, monkeypatch, complete):
    import mhctools.cli
    from topiary import CachedPredictor, from_predictions

    source = admitted_table(tmp_path / 'source.tsv', 'viral')
    producer = RerankingPredictor()
    peptides = ['SIINFEKL', 'GILGFVFTL'] if complete else ['SIINFEKL']
    cache = CachedPredictor(from_predictions(producer.predict_dataframe(peptides)))
    cache_path = tmp_path / 'cache.tsv'
    cache.save(cache_path)
    model = RerankingPredictor()
    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', lambda args: [model])
    argv = ['--input-topiary', str(source), '--external-predictions', 'additive',
            '--external-prediction-prefix', 'cached', '--mhc-predictor', 'random',
            '--mhc-alleles', ','.join(ALLELES), '--prediction-cache', str(cache_path),
            '--output-epitopes', str(tmp_path / 'saved.tsv'), '--no-processing-aware-annotation',
            '--config-value', 'epitopes.score_expr=cached__unified__pMHC_affinity__score',
            '--min-epitope-score', '0.5']
    with pytest.raises(ValueError, match='require --external-peptide-only'):
        run_cli(argv)
    run_cli(argv + ['--external-peptide-only'])
    assert not model.calls if complete else model.calls == [[('GILGFVFTL', '', '')]]
    saved = EpitopeDataset.load(tmp_path / 'saved.tsv')
    assert saved.result.df.value.tolist() == [50., 500.]
    assert [e.sequence for e in epitopes_for_ranking(saved.epitopes, saved.config)] == ['GILGFVFTL']
    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', forbidden)
    replay, _, _ = prepare_reports([('epitopes', tmp_path / 'saved.tsv')])
    assert [e.per_allele_scores for e in replay[0].epitopes] == [e.per_allele_scores for e in saved.epitopes]


def test_additive_legacy_reports_retain_ids_historical_scores_and_features(tmp_path):
    old, _, _ = prepare_reports(INPUTS, scopes=[FIXTURE_SCOPE] * len(INPUTS))
    with pytest.raises(ValueError, match='both flanks'):
        prepare_reports(INPUTS, mode='additive', models=[ContextPredictor()],
                        prediction_prefix='comparison', scopes=[FIXTURE_SCOPE] * len(INPUTS))
    _, frame, fresh = prepare_reports(INPUTS, mode='additive', models=[ContextPredictor()],
        prediction_prefix='comparison', scopes=[FIXTURE_SCOPE] * len(INPUTS), use_flanks=False)
    assert [e.prediction_group_key for r in old for e in r.epitopes] == [e.prediction_group_key for e in fresh]
    assert [e.per_allele_scores for r in old for e in r.epitopes] == [e.per_allele_scores for e in fresh]
    dataset = frame.attrs['epitope_dataset']
    assert dataset.result.df.vaxrank_prediction_id.notna().all()
    assert dataset.result.df.comparison__unified__pMHC_affinity__value.notna().all()
    native = tmp_path / 'native.tsv'
    dataset.save(native)
    replay, _, _ = prepare_reports([('epitopes', native)])
    assert [e.per_allele_scores for e in replay[0].epitopes] == [e.per_allele_scores for e in fresh]


def exacto_manifest(path, *others, translations=False):
    item = dict(format='exacto', path=str(EXACTO / ('translations.tsv.gz' if translations else 'peptide-variants.tsv')),
                sample_id='tumor-rna', library_id='rna-library', exacto=dict(tag='exacto-fixture'))
    if not translations:
        item['exacto'].update(primary_structures=str(EXACTO / 'primary-structures.tsv'),
                            transcript_read_support=str(EXACTO / 'transcript-read-support.tsv'),
                            read_set_id='source-read-namespace')
    path.write_text(json.dumps(dict(schema='vaxrank.input_manifest.v1', **FIXTURE_SCOPE,
                                   inputs=[item, *others])))
    return path


def test_real_exacto_mixed_reports_reconcile_and_replay_without_predictors(tmp_path, monkeypatch):
    import mhctools.cli
    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', forbidden)
    manifest = exacto_manifest(tmp_path / 'inputs.json', *[
        dict(format=fmt, path=path) for fmt, path in INPUTS])
    specs = read_input_manifest(manifest)
    reports = read_reports([(fmt, path) for fmt, path, _ in specs], scopes=[s for _, _, s in specs])
    exacto = reports[0].dataset
    native_source = read_exacto(EXACTO / 'peptide-variants.tsv', sample_name='tumor-rna',
        reference_name='GRCh38', tag='exacto-fixture', primary_structures=EXACTO / 'primary-structures.tsv',
        transcript_read_support=EXACTO / 'transcript-read-support.tsv', library_id='rna-library',
        read_set_id='source-read-namespace')
    assert exacto.result.extra['combined_sources']
    assert len(exacto.epitopes) == len(native_source.df) == 228
    assert not exacto.result.df.candidate_id.notna().any()
    assert not exacto.antigens and all(not e.patient_alleles for e in exacto.epitopes)
    assert exacto.evidence_views()['rna_observations'].shape[0] > 0
    _, frame, _ = prepare_reports([], prepared_reports=reports)
    saved = frame.attrs['epitope_dataset']
    historical = saved.result.df.copy(deep=True)
    views = saved.evidence_views()
    assert len(views['candidates']) > 0
    expected_observations = (historical.source_observation_id.dropna().nunique()
                             + historical.loc[historical.source_observation_id.isna(), 'prediction_id'].nunique())
    assert views['links'].source_observation_id.nunique() == expected_observations
    assert views['rna_observations'].measurement.tolist() == exacto.evidence_views()['rna_observations'].measurement.tolist()
    pd.testing.assert_frame_equal(saved.result.df, historical)
    path = tmp_path / 'native.tsv'
    saved.save(path)
    reloaded = EpitopeDataset.load(path)
    assert reloaded.evidence_views()['links'].source_observation_id.tolist() == views['links'].source_observation_id.tolist()
    before = saved.result.df[saved.result.df.source_label.notna()]
    after = reloaded.result.df[reloaded.result.df.source_label.notna()]
    after = after[before.columns].astype(object)
    before = before.astype(object)
    pd.testing.assert_frame_equal(after.where(after.notna(), None), before.where(before.notna(), None))
    assert after.rna_observations.tolist() == exacto.result.df.rna_observations.tolist()
    out = tmp_path / 'cli.tsv'
    run_cli(['--input-manifest', str(manifest), '--output-epitopes', str(out),
             '--no-processing-aware-annotation'])
    assert len(EpitopeDataset.load(out).result.df) == len(saved.result.df)


def test_exacto_orf_only_input_is_saved_without_invented_candidates(tmp_path, monkeypatch):
    import mhctools.cli
    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', forbidden)
    manifest = exacto_manifest(tmp_path / 'orf.json', translations=True)
    native = tmp_path / 'orf.tsv'
    run_cli(['--input-manifest', str(manifest), '--output-epitopes', str(native)])
    dataset = EpitopeDataset.load(native)
    assert not dataset.epitopes and not dataset.antigens
    assert len(dataset.result.df) == 77 and dataset.result.df.protein_sequence.notna().sum() == 73
    assert dataset.evidence_views()['orfs'].shape[0] == 77


def test_historical_native_views_survive_transport_filename_changes(tmp_path):
    from .test_core_logic_config import _make_epitope

    dataset = EpitopeDataset.from_predictions([_make_epitope('SIINFEKL', 50.)])
    original = dataset.result.df.copy(deep=True)
    before = dataset.evidence_views()['links']
    path = tmp_path / 'old-name.tsv'
    dataset.save(path)
    moved = tmp_path / 'different-name.tsv'
    path.rename(moved)
    reloaded = EpitopeDataset.load(moved)
    pd.testing.assert_frame_equal(before, reloaded.evidence_views()['links'])
    pd.testing.assert_frame_equal(original, dataset.result.df)


def test_pinned_native_exacto_files_are_unmodified():
    provenance = json.loads((EXACTO / 'provenance.json').read_text())
    for name, recorded in provenance['files'].items():
        assert hashlib.sha256((EXACTO / name).read_bytes()).hexdigest() == recorded['sha256']


def test_manifest_exacto_companions_require_scope_and_reject_unknown_options(tmp_path):
    with pytest.raises(ValueError, match='sample_id'):
        prepare_reports([('exacto', EXACTO / 'peptide-variants.tsv')])
    manifest = exacto_manifest(tmp_path / 'inputs.json')
    document = json.loads(manifest.read_text())
    document['inputs'][0]['exacto']['guess_normal'] = 'true'
    manifest.write_text(json.dumps(document))
    with pytest.raises(ValueError, match='unknown Exacto options'):
        read_input_manifest(manifest)
