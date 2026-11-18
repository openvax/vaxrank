"""Common RNA quantities survive scoring, slicing, native files and reports."""
from dataclasses import replace
import math
from types import SimpleNamespace

import pytest
from isovar import RNAFeature, RNASamplingGeometry, RNASamplingProfile, measure_rna_feature
from isovar.rna_path_linkage import OriginalPathMatch

from vaxrank.antigen_input import AntigenDataset
from vaxrank.antigen_report import antigen_evidence_sections
from vaxrank.combined_score_dsl import evaluate_combined_score
from vaxrank.native_serialization import to_native_json, from_native_json
from vaxrank.rna_support import RNAFeatureEvidence, rna_support_bindings
from vaxrank.vaccine_antigen import ANTIGEN_KINDS, VaccineAntigen
from vaxrank.vaccine_peptide import VaccinePeptide
from .test_vaccine_antigen import cta_antigen
from .test_combined_score_dsl import _make_vp


def support(*, count=4, sample='S', profile=True, scope='window', interval=(0, 10)):
    measurement = measure_rna_feature(
        RNAFeature('feature', 'path', ((100, 130),), 'coding_interval'),
        [OriginalPathMatch(('RG', str(i), 0), 'placement', 70, 170) for i in range(count)],
        evidence_scope=(sample, 'library'),
        profile=RNASamplingProfile('original-library', (RNASamplingGeometry(((0, 100),), 100),),
                                  1000000, 1000000, 'whole_library') if profile else None)
    return RNAFeatureEvidence(measurement, scope, interval if scope == 'window' else None)


@pytest.mark.parametrize('kind', sorted(ANTIGEN_KINDS))
def test_all_categories_use_identical_bindings_and_formula(kind):
    antigen = replace(cta_antigen(), kind=kind, rna_support=(support(),))
    vp = VaccinePeptide(antigen=antigen, epitopes=[], combined_score_expr='sqrt(rna_fragments) * target_epitope_score')
    vp.target_epitope_score = 3
    assert vp.combined_score == 6
    assert evaluate_combined_score('sqrt(rna_reads_per_million) * target_epitope_score', vp) == 6
    assert evaluate_combined_score('rna_fragments_per_effective_kb_per_million', vp) == pytest.approx(4000 / 71)
    assert evaluate_combined_score('rna_effective_length_reads', vp) == 71


def test_exact_window_support_is_not_reused_for_subwindows():
    parent = replace(cta_antigen('ACDEFGHIKLMNPQRSTVWY'), rna_support=(support(interval=(3, 13)),))
    assert rna_support_bindings(parent)['rna_fragments'] is None
    window = parent.sliced(3, 13)
    assert rna_support_bindings(window)['rna_fragments'] == 4
    assert not window.sliced(1, 9).rna_support
    assert parent.sliced(1, 16).sliced(2, 12) == window


def test_event_evidence_retains_its_declared_scope_when_sliced():
    parent = replace(cta_antigen(), rna_support=(support(scope='event'),))
    assert parent.sliced(2, 8).rna_support == parent.rna_support
    assert rna_support_bindings(parent.sliced(2, 8))['rna_fragments'] == 4


def test_missing_zero_and_multiple_libraries_are_distinct():
    for records in ((), (support(profile=False),), (support(), support(sample='T2'))):
        antigen = replace(cta_antigen(), rna_support=records)
        vp = SimpleNamespace(antigen=antigen, mutant_protein_fragment=None, target_epitope_score=3)
        with pytest.raises(ValueError, match='Unavailable RNA score bindings'):
            evaluate_combined_score('rna_fragments_per_million', vp)
        assert evaluate_combined_score('target_epitope_score', vp) == 3
    vp.antigen = replace(cta_antigen(), rna_support=(support(count=0),))
    assert evaluate_combined_score('rna_fragments_per_million', vp) == 0


def test_json_native_hash_and_report_roundtrip(tmp_path):
    antigen = replace(cta_antigen(), rna_support=(support(),))
    data = AntigenDataset({'cta': antigen})
    path = tmp_path / 'antigens.json'
    data.save(path)
    restored = AntigenDataset.load(path)
    assert restored == data
    assert {restored.antigens['cta'], antigen} == {antigen}
    assert from_native_json(to_native_json(antigen), VaccineAntigen) == antigen
    sections = antigen_evidence_sections(antigen)
    text = repr(sections)
    assert 'Common RNA feature support' in text
    assert 'fragments per effective kb per million' in text
    assert 'original-library' in text
    assert 'Libraries remain separate' in text


def test_legacy_mutation_score_stays_frozen_with_new_measurements():
    vp = _make_vp(n_alt_fragments=9)
    before = vp.combined_score
    fragment = replace(vp.mutant_protein_fragment, rna_support=(support(scope='event'),))
    after = VaccinePeptide(mutant_protein_fragment=fragment, epitopes=[])
    after.target_epitope_score = vp.target_epitope_score
    assert after.combined_score == before == math.sqrt(9) * vp.target_epitope_score
    assert evaluate_combined_score('rna_fragments', after) == 4
    _, sliced = next(fragment.generate_subsequences(8))
    assert sliced.rna_support == fragment.rna_support


def test_invalid_and_duplicate_attachments_fail():
    with pytest.raises(ValueError, match='beyond'):
        replace(cta_antigen(), rna_support=(support(interval=(0, 11)),))
    with pytest.raises(ValueError, match='Duplicate'):
        replace(cta_antigen(), rna_support=(support(), support()))
    with pytest.raises(ValueError, match='scope'):
        RNAFeatureEvidence(support().measurement, 'pooled')
    fragment = _make_vp().mutant_protein_fragment
    with pytest.raises(ValueError, match='beyond'):
        replace(fragment, rna_support=(support(interval=(0,len(fragment.amino_acids)+1)),))


def test_common_rna_cli_prediction_report_and_offline_replay(tmp_path, monkeypatch):
    from .test_antigen_input import inputs, dependencies
    from topiary import TopiaryPredictor
    from vaxrank.cli.entry_point import run_cli
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dataset import EpitopeDataset
    # Use the suite's deterministic predictor and reference adapters.
    data = inputs.__wrapped__()
    model, genome = dependencies.__wrapped__(monkeypatch)
    data.antigens = {name: replace(a, rna_support=(support(),)) for name, a in data.antigens.items()}
    predicted = data.predict(TopiaryPredictor(models=[model]), genome=genome,
                             config=EpitopeConfig(score_expr='1.0'))
    assert set(predicted.scoring_frame().rna_fragments.dropna()) == {4}
    from vaxrank.epitope_dsl import score_predictions
    scored = score_predictions(predicted.epitopes, EpitopeConfig(score_expr='column(rna_fragments).sqrt()'),
                               topiary_df=predicted.scoring_frame())
    assert set(scored) == {2}
    source = tmp_path / 'input.json'
    data.save(source)
    output = tmp_path / 'design'
    expression = 'sqrt(rna_fragments_per_million) * target_epitope_score'
    run_cli(['--input-antigens', str(source), '--mhc-predictor', 'random',
        '--ensembl-release', '93', '--output-dir', str(output), '--vaccine-type', 'peptide',
        '--no-processing-aware-annotation', '--pdf-backend', 'weasyprint', '--report-candidate-details',
        '--vaccine-peptide-length', '10', '--config-text', 'epitopes.score_expr=1.0',
        '--config-text', 'vaccine_peptides.combined_score_expr=' + expression])
    native = output / 'candidate_predictions.tsv'
    saved = EpitopeDataset.load(native)
    assert all(a.rna_support for a in saved.antigens.values())
    report = (output / 'vaccine_report.txt').read_text()
    assert 'Common RNA feature support' in report
    assert 'original-library' in report
    source.unlink()
    def no_prediction(*args, **kwargs):
        raise AssertionError('Replay must use captured evidence')
    monkeypatch.setattr('mhctools.cli.predictors_from_args', no_prediction)
    replay = tmp_path / 'replay'
    run_cli(['--input-epitopes', str(native), '--output-dir', str(replay),
             '--pdf-backend', 'weasyprint', '--no-processing-aware-annotation'])
    assert (output / 'manifest.json').read_bytes() == (replay / 'manifest.json').read_bytes()


def _large_profile(n_geometries=300):
    geometries = tuple(
        RNASamplingGeometry(((0, 80), (120 + i, 200 + i)), 200 + i, weight=10 + (i % 7))
        for i in range(n_geometries)
    )
    return RNASamplingProfile(
        'shared-whole-library', geometries, 5_000_000, 2_500_000, 'whole_library')


def test_native_export_deduplicates_profiles_measurements_and_shared_records(tmp_path, monkeypatch):
    import json
    import vaxrank.rna_support as rna_mod
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dataset import ANTIGEN_COLUMN, DATASET_METADATA, EpitopeDataset, read_table
    from .test_core_logic_config import _make_epitope

    rna_mod._validate_canonical_measurement.cache_clear()
    profile = _large_profile(300)
    measurement_a = measure_rna_feature(
        RNAFeature('feat-a', 'path-a', ((100, 130),), 'coding_interval'),
        [OriginalPathMatch(('RG', f'read-{i}', 0), 'p', 50, 250) for i in range(6)],
        evidence_scope=('sample-1', 'lib-1'),
        profile=profile,
    )
    measurement_b = measure_rna_feature(
        RNAFeature('feat-b', 'path-b', ((200, 230),), 'coding_interval'),
        [OriginalPathMatch(('RG', f'other-{i}', 0), 'p', 150, 350) for i in range(3)],
        evidence_scope=('sample-1', 'lib-1'),
        profile=profile,
    )
    seq_a = 'M' + 'ACDEFGHIKLMNPQRSTVWY' * 6
    seq_b = 'M' + 'WVTQSRPNMLKIHGFEDCAY' * 6
    antigen_a = replace(cta_antigen(seq_a), source_identifier='cta-a',
                        rna_support=(RNAFeatureEvidence(measurement_a, 'event'),))
    antigen_b = replace(cta_antigen(seq_b), source_identifier='cta-b',
                        rna_support=(RNAFeatureEvidence(measurement_b, 'event'),))

    def build_dataset(n_epitopes, prof):
        meas_a = measure_rna_feature(
            RNAFeature('feat-a', 'path-a', ((100, 130),), 'coding_interval'),
            [OriginalPathMatch(('RG', f'read-{i}', 0), 'p', 50, 250) for i in range(6)],
            evidence_scope=('sample-1', 'lib-1'),
            profile=prof,
        )
        meas_b = measure_rna_feature(
            RNAFeature('feat-b', 'path-b', ((200, 230),), 'coding_interval'),
            [OriginalPathMatch(('RG', f'other-{i}', 0), 'p', 150, 350) for i in range(3)],
            evidence_scope=('sample-1', 'lib-1'),
            profile=prof,
        )
        ant_a = replace(cta_antigen(seq_a), source_identifier='cta-a',
                        rna_support=(RNAFeatureEvidence(meas_a, 'event'),))
        ant_b = replace(cta_antigen(seq_b), source_identifier='cta-b',
                        rna_support=(RNAFeatureEvidence(meas_b, 'event'),))
        frag_a = replace(_make_vp().mutant_protein_fragment, amino_acids=seq_a,
                         rna_support=(RNAFeatureEvidence(meas_a, 'event'),))
        epitopes, antigens, fragments = [], {}, {}
        for idx in range(n_epitopes):
            offset = idx % (len(seq_a) - 9)
            pep = seq_a[offset:offset + 9]
            pid = f'occ-a-{idx}'
            epitopes.append(replace(_make_epitope(pep, ic50=50.0, source_sequence=seq_a, offset=offset),
                                    prediction_id=pid))
            antigens[pid] = ant_a
            fragments[pid] = frag_a
        pid_b = 'occ-b-0'
        epitopes.append(replace(_make_epitope(seq_b[1:10], ic50=80.0, source_sequence=seq_b, offset=1),
                                prediction_id=pid_b))
        antigens[pid_b] = ant_b
        ds = EpitopeDataset.from_predictions(epitopes, config=EpitopeConfig(score_expr='1.0'))
        ds.antigens = antigens
        ds.mutation_fragments = fragments
        return ds, meas_a, meas_b

    tiny_profile = _large_profile(2)
    two_tiny = tmp_path / 'two_tiny.tsv'
    large_tiny = tmp_path / 'large_tiny.tsv'
    two_big = tmp_path / 'two_big.tsv'
    large_path = tmp_path / 'large_big.tsv'
    build_dataset(2, tiny_profile)[0].save(two_tiny)
    build_dataset(120, tiny_profile)[0].save(large_tiny)
    build_dataset(2, profile)[0].save(two_big)
    large_ds, measurement_a, measurement_b = build_dataset(120, profile)
    large_ds.save(large_path)

    raw_payload = read_table(large_path).extra[DATASET_METADATA]
    assert len(raw_payload['rna_sampling_profiles']) == 1
    assert len(raw_payload['rna_measurements']) == 2
    assert len(raw_payload['shared_antigens']) == 1
    assert len(raw_payload['shared_fragments']) == 1
    assert ANTIGEN_COLUMN not in read_table(large_path).df.columns

    profile_delta = len(json.dumps(measurement_a['sampling_profile'])) - len(
        json.dumps(measure_rna_feature(
            RNAFeature('feat-a', 'path-a', ((100, 130),), 'coding_interval'),
            [], evidence_scope=('sample-1', 'lib-1'), profile=tiny_profile)['sampling_profile'])
    )
    assert profile_delta > 20_000
    # Marginal cost of 118 additional epitopes is identical regardless of profile size.
    assert (large_path.stat().st_size - two_big.stat().st_size) == (
        large_tiny.stat().st_size - two_tiny.stat().st_size
    )
    # Increasing the library profile across 241 occurrences adds the profile only once.
    assert (large_path.stat().st_size - large_tiny.stat().st_size) < 2 * profile_delta

    # Validation runs once per unique measurement during a cold load.
    rna_mod._validate_canonical_measurement.cache_clear()
    calls = []
    real_validate = rna_mod.validate_rna_feature_support

    def counting_validate(record):
        calls.append(record['measurement_id'])
        return real_validate(record)

    monkeypatch.setattr(rna_mod, 'validate_rna_feature_support', counting_validate)
    loaded = EpitopeDataset.load(large_path)
    assert len(calls) == 2
    assert set(calls) == {measurement_a['measurement_id'], measurement_b['measurement_id']}
    assert loaded.antigens == large_ds.antigens
    assert {k: f.rna_support for k, f in loaded.mutation_fragments.items()} == {
        k: f.rna_support for k, f in large_ds.mutation_fragments.items()
    }

    # AntigenDataset also deduplicates shared whole-library profiles in-file.
    antigen_file = tmp_path / 'antigens.json'
    antigen_ds = AntigenDataset({'a': antigen_a, 'b': antigen_b})
    antigen_ds.save(antigen_file)
    antigen_doc = json.loads(antigen_file.read_text())
    assert len(antigen_doc['rna_sampling_profiles']) == 1
    assert len(antigen_doc['rna_measurements']) == 2
    assert AntigenDataset.load(antigen_file) == antigen_ds


def test_post_load_dict_mutation_cannot_corrupt_siblings_or_bypass_validation(tmp_path):
    from vaxrank.epitope_dataset import EpitopeDataset
    from .test_core_logic_config import _make_epitope

    ev = support()
    seq = 'ACDEFGHIKLMNPQRSTVWY'
    antigen = replace(cta_antigen(seq), rna_support=(ev,))
    other_antigen = replace(cta_antigen(seq), source_identifier='cta-2', rna_support=(ev,))
    e1 = replace(_make_epitope(seq[:9], ic50=50.0, source_sequence=seq, offset=0), prediction_id='p1')
    e2 = replace(_make_epitope(seq[1:10], ic50=60.0, source_sequence=seq, offset=1), prediction_id='p2')
    ds = EpitopeDataset.from_predictions([e1, e2])
    ds.antigens = {'p1': antigen, 'p2': other_antigen}
    path = tmp_path / 'native.tsv'
    ds.save(path)
    loaded = EpitopeDataset.load(path)

    # Mutating one decoded measurement dict does not alter another antigen's measurement.
    loaded.antigens['p1'].rna_support[0].measurement['values']['fragments'] = 999
    assert loaded.antigens['p2'].rna_support[0].measurement['values']['fragments'] == 4
    with pytest.raises(ValueError, match='disagree|mismatch'):
        replace(loaded.antigens['p1'], source_identifier='revalidated')
    with pytest.raises(ValueError, match='disagree|mismatch'):
        loaded.save(tmp_path / 'tampered-save.tsv')


@pytest.mark.parametrize('tamper,match', [
    (lambda p: next(iter(p['rna_sampling_profiles'].values())).update(library_reads=999), 'sampling profile checksum'),
    (lambda p: next(iter(p['rna_measurements'].values()))['values'].update(fragments=999), 'measurement checksum'),
    (lambda p: p.update(rna_measurements={'0' * 64: next(iter(p['rna_measurements'].values()))}), 'measurement checksum'),
    (lambda p: p['shared_antigens'].update({next(iter(p['shared_antigens'])): '{}'}), 'shared_antigens checksum'),
    (lambda p: p['antigens'].update({'p1': '{"__vaxrank_shared__":"' + 'f' * 64 + '"}'}), 'Unknown native shared record'),
    (lambda p: p['antigens'].update({'p1': '{"__vaxrank_shared__":"not-hex"}'}), 'Malformed native shared record'),
])
def test_tampered_deduplicated_native_payloads_rejected(tmp_path, tamper, match):
    import hashlib
    import json
    from vaxrank.epitope_dataset import DATASET_METADATA, EpitopeDataset, read_table
    from .test_core_logic_config import _make_epitope

    seq = 'ACDEFGHIKLMNPQRSTVWY'
    antigen = replace(cta_antigen(seq), rna_support=(support(),))
    e1 = replace(_make_epitope(seq[:9], ic50=50.0, source_sequence=seq, offset=0), prediction_id='p1')
    e2 = replace(_make_epitope(seq[1:10], ic50=60.0, source_sequence=seq, offset=1), prediction_id='p2')
    ds = EpitopeDataset.from_predictions([e1, e2])
    ds.antigens = {'p1': antigen, 'p2': antigen}
    path = tmp_path / 'native.tsv'
    ds.save(path)

    table = read_table(path)
    tamper(table.extra[DATASET_METADATA])
    table.to_tsv(path)
    with pytest.raises(ValueError, match=match):
        EpitopeDataset.load(path)

    # Even if an attacker recomputes the outer rna_measurements SHA-256 key, Isovar validation rejects it.
    ds.save(path)
    table = read_table(path)
    payload = table.extra[DATASET_METADATA]
    old_digest, meas = next(iter(payload['rna_measurements'].items()))
    profile = next(iter(payload['rna_sampling_profiles'].values()))
    meas['values']['fragments'] = 999
    full = {**meas, 'sampling_profile': profile}
    rehashed = hashlib.sha256(
        json.dumps(full, sort_keys=True, separators=(',', ':')).encode()
    ).hexdigest()
    payload['rna_measurements'] = {rehashed: meas}
    old_shared_digest, shared_json = next(iter(payload['shared_antigens'].items()))
    updated_shared = shared_json.replace(old_digest, rehashed)
    new_shared_digest = hashlib.sha256(updated_shared.encode()).hexdigest()
    payload['shared_antigens'] = {new_shared_digest: updated_shared}
    for k in payload['antigens']:
        payload['antigens'][k] = json.dumps({'__vaxrank_shared__': new_shared_digest}, separators=(',', ':'))
    table.to_tsv(path)
    with pytest.raises(ValueError, match='mismatch|disagree'):
        EpitopeDataset.load(path)


def test_legacy_inline_native_and_antigen_files_still_load(tmp_path):
    from dataclasses import asdict
    import json
    from vaxrank.epitope_dataset import DATASET_METADATA, DATASET_SCHEMA, EpitopeDataset
    from topiary import TopiaryResult
    from .test_core_logic_config import _make_epitope

    seq = 'ACDEFGHIKLMNPQRSTVWY'
    antigen = replace(cta_antigen(seq), rna_support=(support(),))
    epitope = replace(_make_epitope(seq[:9], ic50=50.0, source_sequence=seq, offset=0), prediction_id='legacy-1')
    base = EpitopeDataset.from_predictions([epitope])

    # Write legacy inline vaxrank.epitope_dataset.v1 without deduplication tables.
    legacy_tsv = tmp_path / 'legacy.tsv'
    result = TopiaryResult(base.result.df.copy(), metadata=base.result.metadata)
    result.extra = {
        DATASET_METADATA: {
            'schema': DATASET_SCHEMA,
            'epitopes': [to_native_json(epitope)],
            'provenance': [],
            'config': None,
            'antigens': {'legacy-1': to_native_json(antigen)},
            'selection': {},
            'mutation_fragments': {},
            'construction_reports': [],
            'direct_sources': [],
            'policy_evidence': [],
        }
    }
    result.to_tsv(legacy_tsv)
    loaded = EpitopeDataset.load(legacy_tsv)
    assert loaded.antigens['legacy-1'] == antigen

    # Write legacy inline vaxrank.antigen_dataset.v1 without rna_measurements / rna_sampling_profiles.
    legacy_json = tmp_path / 'legacy_antigens.json'
    legacy_json.write_text(json.dumps({
        'schema': 'vaxrank.antigen_dataset.v1',
        'scope': {},
        'antigens': {'cta': asdict(antigen)},
        'provenance': [],
    }))
    assert AntigenDataset.load(legacy_json).antigens['cta'] == antigen


def test_shared_custom_reference_fragments_relocate_with_rna_support(tmp_path, monkeypatch):
    import shutil
    import pyensembl.download_cache
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dataset import EpitopeDataset
    from vaxrank.mutant_protein_fragment import MutantProteinFragment
    from .test_core_logic_config import _make_epitope
    from .test_epitope_dataset import forbidden
    from .test_isovar_support_contract import MUTANT, reconstruct, reference as _support_reference
    from .test_native_references import delete_original_reference

    ref = _support_reference.__wrapped__(tmp_path, SimpleNamespace(param='+'))
    result = reconstruct(tmp_path, ref, [(MUTANT, [(len(MUTANT), 'M')], 0)] * 2)
    ev = support(scope='event')
    fragment = replace(MutantProteinFragment.from_isovar_result(result), rna_support=(ev,))
    antigen = replace(VaccineAntigen.from_mutant_protein_fragment(fragment), rna_support=(ev,))
    e1 = replace(_make_epitope(fragment.amino_acids[0:9], ic50=100., source_sequence=fragment.amino_acids, offset=0),
                 prediction_id='occ-1')
    e2 = replace(_make_epitope(fragment.amino_acids[1:10], ic50=120., source_sequence=fragment.amino_acids, offset=1),
                 prediction_id='occ-2')
    ds = EpitopeDataset.from_predictions([e1, e2], config=EpitopeConfig(score_expr='1 / affinity.value'))
    ds.antigens = {'occ-1': antigen, 'occ-2': antigen}
    ds.mutation_fragments = {'occ-1': fragment, 'occ-2': fragment}
    archive_path = tmp_path / 'archive' / 'shared_native.tsv'
    ds.save(archive_path)

    moved = tmp_path / 'relocated'
    shutil.move(str(archive_path.parent), str(moved))
    delete_original_reference(tmp_path)
    monkeypatch.setattr(pyensembl.download_cache.DownloadCache, '_fetch', forbidden)

    loaded = EpitopeDataset.load(moved / 'shared_native.tsv')
    assert loaded.antigens['occ-1'] == loaded.antigens['occ-2'] == antigen
    assert loaded.mutation_fragments['occ-1'] is not loaded.mutation_fragments['occ-2']
    assert loaded.mutation_fragments['occ-1'].rna_support == (ev,)
    loaded.index_native_references()
    assert loaded.mutation_fragments['occ-2'].supporting_reference_transcripts[0].gene.name == 'G'
