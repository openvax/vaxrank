"""Published-policy reference parity, evidence requirements and real CLI replay."""

import json
from pathlib import Path

import msgspec
import numpy as np
import pandas as pd
import pytest
from topiary import read_lens, read_pvacseq

from vaxrank.cli.arg_parser import make_vaxrank_arg_parser
from vaxrank.cli.entry_point import run_cli
from vaxrank.config.policies import list_policies, show_policy
from vaxrank.config.loader import extract_epitope_config_kwargs
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.epitope_dsl import score_predictions, attach_per_allele_scores
from vaxrank.selection_policy import decode_evaluations, encode_evaluations
from .test_epitope_dsl import _predictions_df as _frame

DATA = Path(__file__).parent / 'data' / 'published_policies'


def _predictions_df(rows):
    return _frame([dict(allele='HLA-A*02:01', **row) for row in rows])


def config(name):
    return msgspec.convert(extract_epitope_config_kwargs(show_policy(name)), EpitopeConfig)


def test_policy_catalog_and_inspection_without_inputs(capsys):
    parser = make_vaxrank_arg_parser()
    with pytest.raises(SystemExit) as result:
        parser.parse_args(['--list-policies'])
    assert result.value.code == 0
    assert 'lens-v1\t' in capsys.readouterr().out
    for policy in list_policies():
        definition = show_policy(policy['name'])
        assert definition['name'] == policy['name']
        assert definition['epitopes']['selection_policy']['expanded']
    with pytest.raises(SystemExit) as result:
        parser.parse_args(['--show-policy', 'lens-v1'])
    assert result.value.code == 0
    displayed = msgspec.yaml.decode(capsys.readouterr().out)
    assert displayed == show_policy('lens-v1')
    assert displayed['policy_metadata']['source_version'].startswith('LENS 1.9')
    with pytest.raises(SystemExit) as result:
        parser.parse_args(['--show-policy', '../unknown'])
    assert result.value.code == 2
    assert '--list-policies' in capsys.readouterr().err


def test_lens_real_report_uses_actual_coding_support_and_optional_ccf():
    path = Path(__file__).parent / 'data/epitope_fixtures/real_lens_subsets/lens_v1.9_real_subset.tsv'
    report = read_lens(path)
    dataset = EpitopeDataset.from_topiary(report, sample_name='test-patient')
    frame = dataset.scoring_frame()
    cfg = config('lens-v1')
    scores = score_predictions(dataset.epitopes, cfg, topiary_df=frame)
    raw = pd.read_csv(path, sep='\t')
    support = np.log2(raw.rna_reads_covering_genomic_origin_with_peptide_cds)
    expected = np.maximum(0, 1 - raw['mhcflurry_2.1.1.aff'] / 1000) * support / support.max() * raw.ccf.fillna(1)
    decisions = scores.attrs['policy_evaluation'].occurrences
    from mhcgnomes import parse
    raw['allele'] = raw.allele.map(lambda allele: parse(allele).to_string())
    actual = decisions.merge(raw[['peptide', 'allele']].assign(expected=expected),
                             on=['peptide', 'allele'])
    np.testing.assert_allclose(actual.score, actual.expected)
    assert show_policy('lens-v1')['vaccine_peptides']['combined_score_expr'] == 'target_epitope_score'
    replay, = decode_evaluations(encode_evaluations([scores.attrs['policy_evaluation']]))
    pd.testing.assert_frame_equal(replay.occurrences, decisions)


@pytest.mark.parametrize('counts', [[0., 0.], [1., 1.], [None, None]])
def test_lens_undefined_normalization_stays_unscorable(counts):
    frame = _predictions_df([
        dict(peptide='SIINFEKL', value=100., peptide_offset=0),
        dict(peptide='GILGFVFTL', value=200., peptide_offset=1)])
    frame['lens_rna_reads_covering_genomic_origin_with_peptide_cds'] = counts
    scores = score_predictions([], config('lens-v1'), topiary_df=frame)
    assert scores.empty
    assert scores.attrs['policy_evaluation'].occurrences.raw_score.isna().all()


def test_lens_optional_missing_ccf_and_support_scale_have_observable_effect():
    frame = _predictions_df([
        dict(peptide='SIINFEKL', value=100.),
        dict(peptide='GILGFVFTL', value=400.)])
    frame['lens_rna_reads_covering_genomic_origin_with_peptide_cds'] = [4., 16.]
    actual = score_predictions([], config('lens-v1'), topiary_df=frame)
    assert actual.tolist() == pytest.approx([.45, .6])
    frame['ccf'] = [.5, None]
    weighted = score_predictions([], config('lens-v1'), topiary_df=frame)
    assert weighted.tolist() == pytest.approx([.225, .6])
    frame['lens_rna_reads_covering_genomic_origin_with_peptide_cds'] = [64., 16.]
    changed = score_predictions([], config('lens-v1'), topiary_df=frame)
    assert changed.tolist() == pytest.approx([.45, .4])


def test_pvacseq_matches_versioned_upstream_sort_and_native_replay(tmp_path):
    path = DATA / 'pvacseq-aggregate.tsv'
    expected = json.loads((DATA / 'pvacseq-expected.json').read_text())['expected_ids']
    result = read_pvacseq(path)
    dataset = EpitopeDataset.from_topiary(result, sample_name='test-patient')
    partial = dataset.scoring_frame().query('value.isna() and percentile_rank.notna()')
    assert len(partial) > 0
    candidates = {candidate.prediction_group_source: candidate for candidate in dataset.epitopes}
    for row in partial.itertuples():
        assert row.allele in candidates[row.prediction_id].patient_alleles
        assert not candidates[row.prediction_id].predictions
    dataset.config = config('pvacseq-aggregate-v1')
    dataset.epitopes = tuple(attach_per_allele_scores(
        dataset.epitopes, dataset.config, topiary_df=dataset.scoring_frame(),
        policy_evaluations=dataset.policy_evaluations))
    evaluation, = dataset.policy_evaluations
    decisions = evaluation.occurrences
    ids = dataset.scoring_frame().drop_duplicates('prediction_id').set_index('prediction_id').variant
    assert decisions.sort_values('score', ascending=False, kind='stable').prediction_id.map(ids).tolist() == expected
    native = tmp_path / 'native.tsv'
    dataset.save(native)
    loaded = EpitopeDataset.load(native)
    assert loaded.config == dataset.config
    pd.testing.assert_frame_equal(loaded.policy_evaluations[0].occurrences, decisions)
    scores = score_predictions(loaded.epitopes, loaded.config, topiary_df=loaded.scoring_frame())
    assert scores.tolist() == decisions.score.tolist()


@pytest.mark.parametrize('name', ['tesla-presentation-v1', 'tesla-recognition-v1'])
@pytest.mark.parametrize('affinity,abundance,stability,eligible', [
    (33., 34., 1.5, True), (34., 34., 1.5, False),
    (33., 33., 1.5, False), (33., 34., 1.4, False),
    (33., None, 1.5, False), (33., 34., None, False),
])
def test_tesla_presentation_strict_thresholds_and_missing_evidence(
        name, affinity, abundance, stability, eligible):
    frame = _predictions_df([dict(peptide='SIINFEKL', value=affinity, wt_value=1000.)])
    frame['tumor_abundance_tpm'] = abundance
    frame['foreignness'] = None
    stable = frame.assign(kind='pMHC_stability', value=stability)
    scores = score_predictions([], config(name), topiary_df=pd.concat([frame, stable], ignore_index=True))
    assert (not scores.empty) == eligible
    if eligible:
        assert scores.iloc[0] == pytest.approx(1 / (1 + affinity))


@pytest.mark.parametrize('wt,foreignness,eligible', [
    (1000., None, True), (None, 1.01e-16, True),
    (1000., 'absent', True), (None, 'absent', False),
    (300., 1e-16, False), (None, None, False), (0., 0., False),
])
def test_tesla_recognition_or_branch_and_unknowns(wt, foreignness, eligible):
    frame = _predictions_df([dict(peptide='SIINFEKL', value=30., wt_value=wt)])
    frame['tumor_abundance_tpm'] = 34.
    if foreignness != 'absent':
        frame['foreignness'] = foreignness
    frame = pd.concat([frame, frame.assign(kind='pMHC_stability', value=2.)], ignore_index=True)
    assert (not score_predictions([], config('tesla-recognition-v1'), topiary_df=frame).empty) == eligible


def test_published_table_cli_replays_without_predictors(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('Historical policy evaluation must not run predictors')
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    normalized = tmp_path / 'source.tsv'
    result = read_pvacseq(DATA / 'pvacseq-aggregate.tsv')
    result.df['sample_name'] = 'test-patient'
    result.to_tsv(normalized)
    native, first, second = (tmp_path / name for name in ('native.tsv', 'first.csv', 'second.csv'))
    replayed = tmp_path / 'replayed.tsv'
    run_cli(['--input-topiary', str(normalized), '--config', 'builtin:pvacseq-aggregate-v1',
             '--output-epitopes', str(native), '--output-csv', str(first)])
    run_cli(['--input-epitopes', str(native), '--output-csv', str(second),
             '--output-epitopes', str(replayed)])
    a, b = pd.read_csv(first), pd.read_csv(second)
    pd.testing.assert_frame_equal(a[['rank', 'variant', 'vaxrank_score']],
                                  b[['rank', 'variant', 'vaxrank_score']])
    assert (tmp_path / 'native.tsv.policy/selection_policy.json').exists()
    saved = EpitopeDataset.load(replayed).selection['run_configuration']
    assert saved['configuration']['policy_metadata'] == show_policy('pvacseq-aggregate-v1')['policy_metadata']
    assert saved['effective_vaccine_peptides']['combined_score_expr'] == 'target_epitope_score'
    run_cli(['--input-epitopes', str(native), '--output-epitopes', str(replayed),
             '--config-text', 'vaccine_peptides.combined_score_expr=2 * target_epitope_score'])
    saved = EpitopeDataset.load(replayed).selection['run_configuration']
    assert saved['effective_vaccine_peptides']['combined_score_expr'] == '2 * target_epitope_score'


def test_neoepitope_report_ranks_before_rounding(tmp_path):
    from vaxrank.epitope_io import write_neoepitope_report
    evidence = _predictions_df([
        dict(peptide='SIINFEKL', value=2001., prediction_id='lower'),
        dict(peptide='GILGFVFTL', value=2000., prediction_id='higher')])
    report = evidence[['prediction_id', 'peptide', 'peptide_offset', 'allele']].rename(columns={
        'prediction_id': 'Prediction identity', 'peptide': 'Mutant peptide sequence',
        'peptide_offset': 'Peptide offset', 'allele': 'Allele'})
    cfg = EpitopeConfig(selection_policy=dict(name='precision', score_by='1 / affinity.value'))
    path = tmp_path / 'ranks.csv'
    write_neoepitope_report(report, [], topiary_df=evidence,
                            epitope_config=cfg, csv_report_path=path)
    actual = pd.read_csv(path)
    assert actual['Prediction identity'].tolist() == ['higher', 'lower']
    assert actual.vaxrank_score.tolist() == [.0005, .0005]


def test_published_construct_scores_replay_and_explicit_override(tmp_path, monkeypatch):
    from .test_epitope_dataset import evidence
    from vaxrank.epitope_dataset import ANTIGEN_COLUMN
    from vaxrank.native_serialization import to_native_json
    from vaxrank.vaccine_antigen import (
        AminoAcidInterval, TargetableMask, TumorSpecificityAttestation, VaccineAntigen)
    import vaxrank.external_rescoring as external
    antigens = [VaccineAntigen(
        kind='viral', amino_acids=peptide,
        targetable_mask=TargetableMask((AminoAcidInterval(0, len(peptide)),)),
        tumor_specificity=TumorSpecificityAttestation(
            status='admitted', evidence_kind='synthetic_fixture',
            evidence_source='test_published_policies', patient_specific=True,
            rationale_code='test_only'), source_identifier=peptide)
        for peptide in ('SIINFEKL', 'GILGFVFTL')]
    source = tmp_path / 'source.tsv'
    evidence(**{ANTIGEN_COLUMN: [to_native_json(a) for a in antigens],
                'lens_rna_reads_covering_genomic_origin_with_peptide_cds': [4., 16.]}).to_tsv(source)
    outcomes = []
    original = external.load_unified_external
    def capture(*args, **kwargs):
        result = original(*args, **kwargs)
        outcomes.append([peptides[0].combined_score for _, peptides in result[0]])
        return result
    monkeypatch.setattr(external, 'load_unified_external', capture)
    native = tmp_path / 'native.tsv'
    run_cli(['--input-topiary', str(source), '--config', 'builtin:lens-v1',
             '--output-epitopes', str(native)])
    run_cli(['--input-epitopes', str(native), '--output-csv', str(tmp_path / 'replay.csv')])
    run_cli(['--input-epitopes', str(native), '--output-csv', str(tmp_path / 'derived.csv'), '--config-text',
             'vaccine_peptides.combined_score_expr=2 * target_epitope_score'])
    assert outcomes[0] == pytest.approx([.5, .475])
    assert outcomes[1] == pytest.approx(outcomes[0])
    assert outcomes[2] == pytest.approx([1., .95])
