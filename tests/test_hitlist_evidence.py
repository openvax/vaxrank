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

from . import test_cta_input
from .test_cta_expression import PRAME, SEQUENCE
from .test_epitope_dataset import forbidden
from .test_external_rescoring import ALLELES

genome = test_cta_input.genome
cli_dependencies = test_cta_input.cli_dependencies


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
