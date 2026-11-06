"""Installed MS data are read-only optional evidence, captured for native replay."""

import hashlib
import json
import shutil

import pandas as pd
import pytest

from vaxrank.cta_input import hitlist_annotations
from vaxrank.hitlist_evidence import capture_local_hitlist_evidence
from vaxrank.cli.entry_point import run_cli
from vaxrank.epitope_dataset import EpitopeDataset

from . import test_cta_input, test_antigen_input
from .test_cta_expression import PRAME, SEQUENCE
from .test_epitope_dataset import forbidden
from .test_external_rescoring import ALLELES

genome = test_cta_input.genome
cli_dependencies = test_cta_input.cli_dependencies
antigen_inputs = test_antigen_input.inputs
antigen_dependencies = test_antigen_input.dependencies


@pytest.fixture
def installed_data(isolate_hitlist_data, monkeypatch):
    """Real Hitlist reader, synthetic index; freshness is a producer-owned check."""
    directory = isolate_hitlist_data
    directory.mkdir()
    base = dict(peptide=SEQUENCE[:9], mhc_species='Homo sapiens', mhc_class='I',
                mhc_restriction=ALLELES[0], restriction_evidence='monoallelic',
                mhc_allele_provenance='restriction', mhc_allele_set=ALLELES[0],
                source='iedb', reference_iri='fixture-reference', pmid=123,
                assay_method='mass spectrometry', qualitative_measurement='Positive',
                source_organism='Homo sapiens', host='Homo sapiens', provenance_id='observed')
    rows = [base,
            dict(base, mhc_restriction='HLA class I', restriction_evidence='unknown',
                 mhc_allele_provenance='sample_allele_match', mhc_allele_set=';'.join(ALLELES),
                 provenance_id='donor-typed'),
            dict(base, mhc_species='Mus musculus', mhc_restriction='H2-Kb',
                 mhc_allele_set='H2-Kb', source_organism='Mus musculus', host='Mus musculus',
                 provenance_id='mouse'),
            dict(base, peptide=SEQUENCE[:10], provenance_id='longer'),
            dict(base, assay_method='fluorescence', provenance_id='not-ms'),
            dict(base, qualitative_measurement='Negative', provenance_id='negative')]
    pd.DataFrame(rows).to_parquet(directory / 'observations.parquet', index=False)
    (directory / 'observations_meta.json').write_text(json.dumps({'fixture': 'synthetic index'}))
    monkeypatch.setattr('hitlist.observations.observations_cache_is_current', lambda: True)
    return directory


def test_missing_data_never_creates_cache_or_queries(isolate_hitlist_data, monkeypatch):
    monkeypatch.setattr('hitlist.observations.load_ms_observations', forbidden)
    monkeypatch.setattr('hitlist.observations.observations_cache_is_current', forbidden)
    evidence = capture_local_hitlist_evidence([SEQUENCE[:9]], ALLELES)
    assert evidence['status'] == 'not_installed'
    assert not isolate_hitlist_data.exists()
    annotations = hitlist_annotations([SEQUENCE[:9]], None, [ALLELES[0]], local_evidence=evidence)
    assert annotations.hitlist_ms_observations.isna().all()
    assert annotations.hitlist_query_status.tolist() == ['unqueried']


@pytest.mark.parametrize('current, expected', [(False, 'stale'), (None, 'unverifiable')])
def test_unusable_index_stays_unknown(installed_data, monkeypatch, current, expected):
    monkeypatch.setattr('hitlist.observations.observations_cache_is_current', lambda: current)
    monkeypatch.setattr('hitlist.observations.load_ms_observations', forbidden)
    evidence = capture_local_hitlist_evidence([SEQUENCE[:9]], ALLELES)
    assert evidence['status'] == expected
    assert not evidence['species_observations']
    assert hitlist_annotations([SEQUENCE[:9]], None, [ALLELES[0]],
                               local_evidence=evidence).hitlist_ms_observations.isna().all()


def test_exact_species_queries_preserve_scope_and_restriction_evidence(installed_data):
    before = {p.name: p.read_bytes() for p in installed_data.iterdir()}
    evidence = capture_local_hitlist_evidence([SEQUENCE[:9], 'VVVVVVVVV'], ALLELES + ['H2-Kb'])
    assert evidence['status'] == 'captured'
    assert evidence['index_metadata_sha256'] == hashlib.sha256(before['observations_meta.json']).hexdigest()
    human = evidence['species_observations']['Homo sapiens']['rows']
    assert {row['provenance_id'] for row in human} == {'observed', 'donor-typed'}
    assert all(row['source_organism'] == 'Homo sapiens' for row in human)
    annotations = hitlist_annotations(
        [SEQUENCE[:9]] * 3 + ['VVVVVVVVV', SEQUENCE[:8], SEQUENCE[:9]], None,
        ALLELES + ['H2-Kb', ALLELES[0], ALLELES[0], 'DLA-88*501:01'], local_evidence=evidence)
    assert annotations.hitlist_ms_observations[:4].tolist() == [2, 2, 1, 0]
    assert annotations.hitlist_ms_observations[4:].isna().all()
    assert annotations.hitlist_query_status.tolist() == ['captured'] * 4 + ['unqueried'] * 2
    assert annotations.hitlist_ms_monoallelic[:3].tolist() == [1, 0, 1]
    assert annotations.hitlist_ms_donor_hla_candidate[:3].tolist() == [1, 1, 0]
    assert annotations.hitlist_ms_coarse_or_untyped[:3].tolist() == [1, 1, 0]
    for name in ('hitlist_cta_specific', 'hitlist_avoid_sequence',
                 'hitlist_tissue_review_required', 'hitlist_essential_tissue_donors'):
        assert annotations[name].isna().all()
    assert before == {p.name: p.read_bytes() for p in installed_data.iterdir()}


def test_changed_index_is_rejected(installed_data, monkeypatch):
    from hitlist.observations import load_ms_observations
    def mutate(**kwargs):
        frame = load_ms_observations(**kwargs)
        (installed_data / 'observations_meta.json').write_text('{"changed": true}')
        return frame
    monkeypatch.setattr('hitlist.observations.load_ms_observations', mutate)
    with pytest.raises(ValueError, match='data changed'):
        capture_local_hitlist_evidence([SEQUENCE[:9]], ALLELES)


@pytest.mark.parametrize('column', ['hitlist_ms_observations', 'hitlist_ms_monoallelic'])
def test_local_evidence_cli_and_portable_replay(
        installed_data, tmp_path, monkeypatch, cli_dependencies, column):
    expression = tmp_path / 'expression.tsv'
    pd.DataFrame([(PRAME, 3)], columns=['feature_id', 'patient_tpm']).to_csv(
        expression, sep='\t', index=False)
    first, native = tmp_path / 'first', tmp_path / 'native.tsv'
    run_cli(test_cta_input.argv(expression, first) + [
        '--output-epitopes', str(native), '--vaccine-type', 'peptide',
        '--config-text', 'epitopes.filter_expr=' + column + ' > 0'])
    saved = EpitopeDataset.load(native)
    assert saved.selection['hitlist_local_evidence']['status'] == 'captured'
    assert saved.result.df.hitlist_ms_observations.max() == 2
    assert saved.epitopes
    assert len({e.prediction_id for e in saved.epitopes}) == len(saved.epitopes)
    assert all(e.source_class == 'self' for e in saved.epitopes)
    observed = [e for e in saved.epitopes if e.sequence == SEQUENCE[:9]]
    assert observed and all({p.allele for p in e.predictions_flat()} == set(ALLELES) for e in observed)
    retained = {ALLELES[0]} if column == 'hitlist_ms_monoallelic' else set(ALLELES)
    assert all({allele for allele, score in e.per_allele_scores.items() if score > 0} == retained
               for e in observed)
    expression.unlink()
    shutil.rmtree(installed_data)
    monkeypatch.setattr('vaxrank.cta_input.capture_local_hitlist_evidence', forbidden)
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    replay, restored_path = tmp_path / 'replay', tmp_path / 'restored.tsv'
    run_cli(['--input-epitopes', str(native), '--output-epitopes', str(restored_path),
             '--output-dir', str(replay), '--vaccine-type', 'peptide', '--no-processing-aware-annotation'])
    restored = EpitopeDataset.load(restored_path)
    assert restored.selection['hitlist_local_evidence'] == saved.selection['hitlist_local_evidence']
    assert restored.epitopes == saved.epitopes
    assert (first / 'cta_design.json').read_bytes() == (replay / 'cta_design.json').read_bytes()
    assert (first / 'vaccine.fasta').read_bytes() == (replay / 'vaccine.fasta').read_bytes()


def test_global_topiary_cli_preserves_occurrences_and_offline_replay(
        installed_data, tmp_path, monkeypatch):
    from topiary import TopiaryResult
    source, native, replay = (tmp_path / name for name in ('source.tsv', 'native.tsv', 'replay.tsv'))
    peptide = SEQUENCE[:9]
    rows = [dict(peptide=sequence, allele=allele, prediction_id=identity,
                 kind='pMHC_affinity', value=20., prediction_method_name='model')
            for identity, sequence in [('one', peptide), ('two', peptide), ('absent', 'VVVVVVVVV')]
            for allele in ALLELES]
    TopiaryResult(pd.DataFrame(rows)).to_tsv(source)
    run_cli(['--input-topiary', str(source), '--output-epitopes', str(native),
             '--hitlist-local-evidence', '--config-text', 'epitopes.score_expr=1.0',
             '--config-text', 'epitopes.filter_expr=hitlist_ms_monoallelic > 0'])
    saved = EpitopeDataset.load(native)
    assert len(saved.epitopes) == 3
    assert len({e.prediction_id for e in saved.epitopes}) == 3
    for candidate in saved.epitopes:
        expected = {ALLELES[0]: 1.} if candidate.sequence == peptide else {}
        assert {a: s for a, s in candidate.per_allele_scores.items() if s > 0} == expected
    captures = saved.selection['source_evidence']
    assert len(captures) == 1
    query = next(iter(captures.values()))['hitlist_local_queries'][0]
    assert query['status'] == 'captured'
    assert query['query']['peptides'] == sorted([peptide, 'VVVVVVVVV'])
    assert len(query['species_observations']['Homo sapiens']['rows']) == 2
    shutil.rmtree(installed_data)
    source.unlink()
    monkeypatch.setattr('vaxrank.hitlist_evidence.capture_local_hitlist_evidence', forbidden)
    run_cli(['--input-epitopes', str(native), '--output-epitopes', str(replay),
             '--hitlist-local-evidence'])
    again = EpitopeDataset.load(replay)
    assert again.selection['source_evidence'] == captures
    assert [e.per_allele_scores for e in again.epitopes] == [e.per_allele_scores for e in saved.epitopes]
    assert again.result.df.hitlist_ms_observations.tolist() == saved.result.df.hitlist_ms_observations.tolist()


@pytest.mark.parametrize('enabled', [False, True])
def test_global_opt_in_and_unavailable_replay(tmp_path, monkeypatch, enabled, isolate_hitlist_data):
    from .test_epitope_dataset import evidence
    from vaxrank.external_rescoring import prepare_reports
    source, native = tmp_path / 'input.tsv', tmp_path / 'native.tsv'
    evidence().to_tsv(source)
    if not enabled:
        monkeypatch.setattr('vaxrank.hitlist_evidence.capture_local_hitlist_evidence', forbidden)
    _, frame, _ = prepare_reports([('topiary', source)], hitlist_local_evidence=enabled)
    dataset = frame.attrs['epitope_dataset']
    assert ('hitlist_ms_observations' in dataset.result.df) == enabled
    if enabled:
        assert dataset.result.df.hitlist_ms_observations.isna().all()
        captured = next(iter(dataset.selection['source_evidence'].values()))['hitlist_local_queries']
        assert captured[0]['status'] == 'not_installed'
        dataset.save(native)
        monkeypatch.setattr('vaxrank.hitlist_evidence.capture_local_hitlist_evidence', forbidden)
        _, replay, _ = prepare_reports([('epitopes', native)], hitlist_local_evidence=True)
        assert replay.attrs['epitope_dataset'].selection['source_evidence'] == dataset.selection['source_evidence']
    assert not isolate_hitlist_data.exists()


@pytest.mark.parametrize('mode', ['input', 'fresh', 'additive'])
def test_global_legacy_table_modes_retain_capture(installed_data, mode, tmp_path):
    from .test_external_rescoring import INPUTS, ContextPredictor, scoped_reports
    from vaxrank.epitope_config import EpitopeConfig
    options = {} if mode == 'input' else dict(models=[ContextPredictor()], alleles=ALLELES)
    if mode == 'additive':
        options['prediction_prefix'] = 'fresh'
        options['use_flanks'] = False
    reports, frame, candidates = scoped_reports(
        INPUTS, EpitopeConfig(score_expr='1.0', min_epitope_score=0), mode=mode,
        hitlist_local_evidence=True, **options)
    dataset = frame.attrs['epitope_dataset']
    assert candidates and len(reports) == 2
    assert set(dataset.result.df.hitlist_query_status) == {'captured'}
    assert (dataset.result.df.hitlist_ms_observations == 0).all()
    assert len(dataset.selection['source_evidence']) == 2
    assert len(dataset.construction_reports) == 2
    dataset.save(tmp_path / 'native.tsv')
    assert EpitopeDataset.load(tmp_path / 'native.tsv').selection['source_evidence'] == dataset.selection['source_evidence']


def test_query_collector_keeps_scope_and_nondefault_index(installed_data):
    from vaxrank.hitlist_evidence import LocalHitlistLookup
    lookup = LocalHitlistLookup()
    human = pd.DataFrame({'peptide': [SEQUENCE[:9]], 'allele': [ALLELES[0]]}, index=[7])
    mouse = pd.DataFrame({'peptide': ['VVVVVVVVV'], 'allele': ['H2-Kb']}, index=[11])
    lookup(human)
    lookup(mouse)
    assert human.hitlist_ms_observations.tolist() == [2]
    assert mouse.hitlist_ms_observations.tolist() == [0]
    exported = pd.DataFrame({'peptide': [SEQUENCE[:9]] * 2 + ['VVVVVVVVV'],
                             'allele': [ALLELES[0], 'H2-Kb', 'H2-Kb']}, index=[4, 5, 6])
    lookup.annotate_captured(exported)
    assert exported.hitlist_ms_observations.loc[4] == 2
    assert pd.isna(exported.hitlist_ms_observations.loc[5])  # Never queried this pair.
    assert exported.hitlist_ms_observations.loc[6] == 0
    assert len(lookup.queries) == 2


def test_direct_prediction_uses_hitlist_before_filtering(installed_data, monkeypatch):
    from unittest.mock import MagicMock
    from pyensembl import EnsemblRelease
    from topiary import TopiaryPredictor
    from vaxrank.epitope_logic import predict_epitopes
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.hitlist_evidence import LocalHitlistLookup
    from vaxrank.vaccine_antigen import (
        VaccineAntigen, SelfReferenceMatch, TargetableMask, TumorSpecificityAttestation,
        AminoAcidInterval)
    from .test_external_rescoring import ContextPredictor
    monkeypatch.setattr('vaxrank.epitope_logic.ReferenceProteome', MagicMock())
    monkeypatch.setattr('vaxrank.epitope_logic.self_reference_matches',
        lambda peptides, antigen, genome: {p: SelfReferenceMatch(p, False, antigen.kind) for p in peptides})
    lookup = LocalHitlistLookup()
    candidates = predict_epitopes(
        TopiaryPredictor(models=[ContextPredictor()]),
        antigen=VaccineAntigen(kind='mutation', amino_acids=SEQUENCE, gene_name='synthetic',
            targetable_mask=TargetableMask((AminoAcidInterval(0, len(SEQUENCE)),)),
            tumor_specificity=TumorSpecificityAttestation(
                status='admitted', evidence_kind='synthetic_fixture', evidence_source=__name__)),
        genome=EnsemblRelease(93), annotate=lookup,
        epitope_config=EpitopeConfig(filter_expr='hitlist_ms_monoallelic > 0',
                                     score_expr='1.0', min_epitope_score=0))
    assert len(candidates) == 2  # The exact sequence occurs at two source positions.
    assert all(e.sequence == SEQUENCE[:9] for e in candidates)
    assert all(e.per_allele_scores == {ALLELES[0]: 1.} for e in candidates)
    assert len(lookup.queries) == 1


def test_global_lookup_rejects_legacy_ranked_json(tmp_path):
    with pytest.raises(ValueError, match='legacy ranked JSON cannot be rescored'):
        run_cli(['--input-json-file', str(tmp_path / 'old.json'), '--hitlist-local-evidence',
                 '--output-csv', str(tmp_path / 'out.csv')])


def test_global_pooled_native_adds_only_unqueried_sources(installed_data, tmp_path, monkeypatch):
    from .test_epitope_dataset import evidence
    from .input_scope_helpers import FIXTURE_SCOPE
    from vaxrank.external_rescoring import prepare_reports
    a, b, native, mixed = (tmp_path / name for name in ('a.tsv', 'b.tsv', 'a-native.tsv', 'mixed.tsv'))
    evidence(peptide=[SEQUENCE[:9], 'VVVVVVVVV']).to_tsv(a)
    evidence(peptide=['AAAAAAAAA', 'LLLLLLLLL']).to_tsv(b)
    _, frame, _ = prepare_reports([('topiary', a)], scopes=[FIXTURE_SCOPE], hitlist_local_evidence=True)
    original = frame.attrs['epitope_dataset']
    original.save(native)
    _, frame, _ = prepare_reports([('epitopes', native), ('topiary', b)], scopes=[FIXTURE_SCOPE] * 2)
    frame.attrs['epitope_dataset'].save(mixed)
    queries = []
    def capture(peptides, alleles):
        queries.append(set(peptides))
        return capture_local_hitlist_evidence(peptides, alleles)
    monkeypatch.setattr('vaxrank.hitlist_evidence.capture_local_hitlist_evidence', capture)
    _, frame, _ = prepare_reports([('epitopes', mixed)], hitlist_local_evidence=True)
    result = frame.attrs['epitope_dataset']
    assert queries == [{'AAAAAAAAA', 'LLLLLLLLL'}]
    assert len(result.selection['source_evidence']) == 2
    source_ids = {e.input_provenance.input_id for e in result.epitopes}
    assert set(result.selection['source_evidence']) == source_ids
    for key, value in original.selection['source_evidence'].items():
        assert result.selection['source_evidence'][key] == value
    assert set(result.result.df.hitlist_query_status) == {'captured'}




def test_global_antigen_cli_and_source_evidence_replay(
        antigen_inputs, antigen_dependencies, installed_data, tmp_path, monkeypatch):
    from dataclasses import replace
    from vaxrank.antigen_input import AntigenDataset
    records = {key: replace(value, amino_acids=SEQUENCE[:10])
               for key, value in antigen_inputs.antigens.items() if key in ('viral', 'mutation')}
    source, native, replay = (tmp_path / name for name in ('antigens.json', 'native.tsv', 'replay.tsv'))
    AntigenDataset(records, antigen_inputs.scope).save(source)
    run_cli(['--input-antigens', str(source), '--mhc-predictor', 'random',
             '--ensembl-release', '93', '--output-epitopes', str(native),
             '--hitlist-local-evidence', '--config-text', 'epitopes.score_expr=1.0',
             '--config-text', 'epitopes.filter_expr=hitlist_ms_monoallelic > 0'])
    saved = EpitopeDataset.load(native)
    assert {a.kind for a in saved.antigens.values()} == {'viral', 'mutation'}
    assert len(saved.epitopes) == 4  # Two occurrences per source, two alleles each.
    observed = [e for e in saved.epitopes if e.sequence == SEQUENCE[:9]]
    assert len(observed) == 2
    assert all(e.per_allele_scores == {ALLELES[0]: 1.} for e in observed)
    source.unlink()
    shutil.rmtree(installed_data)
    monkeypatch.setattr('vaxrank.hitlist_evidence.capture_local_hitlist_evidence', forbidden)
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    run_cli(['--input-epitopes', str(native), '--output-epitopes', str(replay),
             '--hitlist-local-evidence'])
    again = EpitopeDataset.load(replay)
    assert again.antigens == saved.antigens
    assert again.selection['source_evidence'] == saved.selection['source_evidence']
    assert [e.per_allele_scores for e in again.epitopes] == [e.per_allele_scores for e in saved.epitopes]


def test_collector_can_query_after_allele_free_rows(installed_data):
    from vaxrank.hitlist_evidence import LocalHitlistLookup
    lookup = LocalHitlistLookup()
    unknown = pd.DataFrame({'peptide': [SEQUENCE[:9]], 'allele': ['']})
    measured = pd.DataFrame({'peptide': [SEQUENCE[:9]], 'allele': [ALLELES[0]]})
    lookup(unknown)
    lookup(measured)
    assert [q['status'] for q in lookup.queries] == ['no_candidates', 'captured']
    assert unknown.hitlist_ms_observations.isna().all()
    assert measured.hitlist_ms_observations.tolist() == [2]


def test_collector_rejects_index_change_between_queries(installed_data):
    from vaxrank.hitlist_evidence import LocalHitlistLookup
    lookup = LocalHitlistLookup()
    lookup(pd.DataFrame({'peptide': [SEQUENCE[:9]], 'allele': [ALLELES[0]]}))
    (installed_data / 'observations_meta.json').write_text('{"fixture": "changed index"}')
    with pytest.raises(ValueError, match='changed between queries'):
        lookup(pd.DataFrame({'peptide': ['VVVVVVVVV'], 'allele': [ALLELES[0]]}))
