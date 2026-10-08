"""Consumer contracts for portable Topiary policies and the frozen baseline."""

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path

import msgspec
import pandas as pd
import pytest
from topiary import replay_selection_policy
from topiary import read_tsv
from topiary.ranking import Column

from vaxrank.config.loader import (
    load_vaxrank_config, extract_epitope_config_kwargs, resolve_config_path,
)
from vaxrank.cli.epitope_config_args import epitope_config_from_args
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dsl import score_predictions
from vaxrank.selection_policy import (
    encode_evaluations, decode_evaluations, save_policy_evaluations,
)
from vaxrank.native_serialization import from_native_json, to_native_json
from vaxrank.peptide import PeptideConstruct
from .test_epitope_dsl import _predictions_df


def frozen_config():
    return msgspec.convert(extract_epitope_config_kwargs(
        load_vaxrank_config(config_path='builtin:openvax-v1')), EpitopeConfig)


def test_frozen_policy_matches_legacy_scores_and_replays(tmp_path):
    frame = _predictions_df([
        dict(peptide='SIINFEKL', allele='HLA-A*02:01', peptide_offset=i, value=value)
        for i, value in enumerate([0., 50., 350., 500., 1000., 4999., 5000., None])])
    legacy = score_predictions([], EpitopeConfig(), topiary_df=frame)
    actual = score_predictions([], frozen_config(), topiary_df=frame)
    pd.testing.assert_series_equal(actual, legacy[legacy >= 1e-5], check_names=False)
    evaluation = actual.attrs['policy_evaluation']
    assert len(evaluation.occurrences) == 8
    assert len(evaluation.evidence.df) == len(frame)
    decoded, = decode_evaluations(encode_evaluations([evaluation]))
    pd.testing.assert_frame_equal(evaluation.occurrences, decoded.occurrences)
    record, = save_policy_evaluations([evaluation], tmp_path)
    restored = replay_selection_policy(read_tsv(tmp_path / record['evidence']))
    pd.testing.assert_frame_equal(evaluation.occurrences, restored.occurrences)


@pytest.mark.parametrize('renderer', ['repr', 'to_expr_string'])
@pytest.mark.parametrize('expression, measurements, threshold, expected_scores, eligible', [
    pytest.param(
        Column('review_score') - (Column('penalty') - Column('credit')),
        dict(review_score=[10., 8.], penalty=[3., 3.], credit=[1., 1.]),
        7., [8., 6.], [True, False], id='grouped-subtraction'),
    pytest.param(
        Column('signed_score').clip(-1, 1), dict(signed_score=[-2., 0., 2.]),
        -1., [-1., 0., 1.], [True, True, True], id='negative-clip-bound'),
    pytest.param(
        Column('review score') + 10 * Column('review label').eq("approved 'quoted'\\line\nnext"),
        {'review score': [1., 9.], 'review label': ["approved 'quoted'\\line\nnext", 'rejected']},
        10., [11., 9.], [True, False], id='categorical-quoted-equality'),
    pytest.param(
        10 * ~Column('review label').isin(['rejected', 'pending']),
        {'review label': ['approved', 'rejected', 'pending']},
        1., [10., 0., 0.], [True, False, False], id='negated-membership'),
])
def test_rendered_policy_yaml_keeps_scores_selection_and_evidence_replay(
        tmp_path, renderer, expression, measurements, threshold, expected_scores, eligible):
    """Saved policies preserve arithmetic, transform parameters and categories."""
    def render(node):
        return repr(node) if renderer == 'repr' else node.to_expr_string()

    policy = dict(
        name='saved-expression', score_by=render(expression), min_score=None,
        criteria=[dict(name='accepted', expression=render(expression >= threshold),
                       role='eligibility')], filter_by='criterion("accepted")')
    path = tmp_path / 'policy.yaml'
    path.write_bytes(msgspec.yaml.encode({'epitopes': {'selection_policy': policy}}))
    cfg = msgspec.convert(extract_epitope_config_kwargs(
        load_vaxrank_config(config_path=str(path))), EpitopeConfig)
    frame = _predictions_df([
        dict(peptide='SIINFEKL', allele='HLA-A*02:01', peptide_offset=i, value=50.,
             **{name: values[i] for name, values in measurements.items()})
        for i in range(len(expected_scores))])

    scores = score_predictions([], cfg, topiary_df=frame)
    assert scores.tolist() == [score for score, keep in zip(expected_scores, eligible) if keep]
    evaluation = scores.attrs['policy_evaluation']
    assert evaluation.occurrences.eligible.tolist() == eligible
    # Topiary evaluates scores after eligibility; rejected occurrences stay
    # in the audit with a missing score rather than an evaluated value.
    assert evaluation.occurrences.loc[eligible, 'score'].tolist() == scores.tolist()
    assert evaluation.occurrences.loc[[not keep for keep in eligible], 'score'].isna().all()
    assert len(evaluation.evidence.df) == len(frame)
    decoded, = decode_evaluations(encode_evaluations([evaluation]))
    pd.testing.assert_frame_equal(evaluation.occurrences, decoded.occurrences)
    record, = save_policy_evaluations([evaluation], tmp_path / 'evidence')
    restored = replay_selection_policy(read_tsv(tmp_path / 'evidence' / record['evidence']))
    assert restored.policy.sha256 == evaluation.policy.sha256
    pd.testing.assert_frame_equal(evaluation.occurrences, restored.occurrences)


def test_policy_overlay_preserves_explicit_null_and_cli_minimum(tmp_path):
    overlay = tmp_path / 'custom.yaml'
    overlay.write_text('name: custom\nepitopes:\n  selection_policy:\n    name: custom\n    min_score: null\n')
    merged = load_vaxrank_config(config_path=['builtin:openvax-v1', str(overlay)])
    cfg = msgspec.convert(extract_epitope_config_kwargs(merged), EpitopeConfig)
    assert cfg.selection_policy['min_score'] is None
    args = argparse.Namespace(min_epitope_score=.4)
    cfg = epitope_config_from_args(args, merged)
    assert cfg.min_epitope_score == cfg.selection_policy['min_score'] == .4
    with pytest.raises(ValueError, match='selection_policy'):
        EpitopeConfig(selection_policy=cfg.selection_policy, score_expr='1.0')


def test_named_criteria_rejections_and_explicit_model_choice():
    frame = _predictions_df([
        dict(peptide='SIINFEKL', allele='HLA-A*02:01', value=50., prediction_method_name='mhcflurry'),
        dict(peptide='SIINFEKL', allele='HLA-A*02:01', value=800., prediction_method_name='netmhcpan')])
    cfg = EpitopeConfig(selection_policy=dict(
        name='test', score_by='affinity.value',
        default_methods={'pMHC_affinity': 'netmhcpan'},
        criteria=[dict(name='binder', expression='affinity.value < 500', role='eligibility')],
        filter_by='criterion("binder")'))
    scored = score_predictions([], cfg, topiary_df=frame)
    assert scored.empty
    assert not scored.attrs['policy_evaluation'].audit.empty


def test_saved_policy_can_be_derived_with_yaml_or_cli_without_stale_expansion(tmp_path):
    saved = tmp_path / 'saved.yaml'
    saved.write_bytes(msgspec.yaml.encode({'epitopes': {'selection_policy': frozen_config().selection_policy}}))
    overlay = tmp_path / 'override.yaml'
    overlay.write_text('epitopes:\n  selection_policy:\n    score_by: "0.5"\n')
    for kwargs in (dict(config_path=[str(saved), str(overlay)]),
                   dict(config_path=str(saved), expr_overrides=['epitopes.selection_policy.score_by=0.5'])):
        merged = load_vaxrank_config(**kwargs)
        cfg = msgspec.convert(extract_epitope_config_kwargs(merged), EpitopeConfig)
        assert cfg.selection_policy['score_by'] == '0.5'
        assert cfg.selection_policy['expanded']['score_by']['expression'] == '0.5'
    corrupted = frozen_config().selection_policy
    corrupted['expanded']['score_by'] = '0'
    saved.write_bytes(msgspec.yaml.encode({'epitopes': {'selection_policy': corrupted}}))
    with pytest.raises(ValueError, match='expanded'):
        load_vaxrank_config(config_path=[str(saved), str(overlay)])


def test_peptide_construct_native_codec():
    construct = PeptideConstruct('one', 'SIINFEKL', ['gene'],
                                 components={'n_terminal_acetylation': True})
    assert from_native_json(to_native_json(construct), PeptideConstruct) == construct
    assert from_native_json(to_native_json(replace(construct, name='two')), PeptideConstruct).name == 'two'


def test_named_ranking_criteria_select_and_replay_representatives(tmp_path):
    from topiary import TopiaryResult, combine_sources
    from .test_epitope_dataset import evidence
    from vaxrank.external_rescoring import prepare_reports
    from vaxrank.epitope_dataset import EpitopeDataset
    first = evidence().long_df.iloc[:1].copy()
    second = first.assign(new_feature=100.)
    result = combine_sources({'first': TopiaryResult(first), 'second': TopiaryResult(second)},
                              sample_name='patient')
    source = tmp_path / 'input.tsv'
    result.to_tsv(source)
    cfg = EpitopeConfig(selection_policy=dict(name='rank-test', score_by='1.0', duplicates='best',
        criteria=[dict(name='support', expression='new_feature', role='ranking')],
        ranking_by=[dict(expression='criterion("support")', ascending=False)]))
    reports, frame, _ = prepare_reports([('topiary', source)], cfg)
    dataset = frame.attrs['epitope_dataset']
    selected = dataset.select_representatives()
    assert len(selected) == 1
    expected = result.df.loc[result.df.new_feature == 100., 'source_observation_id'].iloc[0]
    assert next(iter(selected))[0] == expected
    assert reports[0].dataset.select_representatives() == selected
    saved = tmp_path / 'native.tsv'
    dataset.save(saved)
    restored = EpitopeDataset.load(saved)
    assert restored.select_representatives() == selected
    assert restored.config == cfg
    assert len(restored.result.df) == 2


@pytest.mark.parametrize('name', ['openvax-v1', 'presentation-v1', 'presentation-ba-v1',
                                  'presentation-ba-cterm-v1', 'self-trim-v1', 'peptide-serum-v1'])
def test_bundled_configs_compose(name):
    from vaxrank.config.loader import extract_vaccine_config_kwargs, extract_construct_kwargs
    from vaxrank.vaccine_config import VaccineConfig
    from vaxrank.peptide import PeptideConstructConfig
    merged = load_vaxrank_config(config_path=['builtin:openvax-v1', 'builtin:' + name])
    msgspec.convert(extract_epitope_config_kwargs(merged), EpitopeConfig)
    msgspec.convert(extract_vaccine_config_kwargs(merged), VaccineConfig)
    kwargs = extract_construct_kwargs(merged, 'peptide')
    # YAML terminal-chemistry names are translated by the CLI constructor.
    kwargs['n_terminal_acetylation'] = kwargs.pop('n_terminal_acetyl')
    kwargs['c_terminal_amidation'] = kwargs.pop('c_terminal_amide')
    PeptideConstructConfig(**kwargs)


def test_frozen_definition_digest():
    # This golden hash catches accidental edits to the released definition.
    digest = hashlib.sha256(resolve_config_path('builtin:openvax-v1').read_bytes()).hexdigest()
    assert digest == '22f13665dcc23c71eb3973e6dff41fb16361a2dbb2780f497a76f99aae4cf9da'


def test_real_sid_model_evidence_replays_all_named_scores():
    directory = Path(__file__).parent / 'data' / 'selection_policy'
    frame = read_tsv(directory / 'sid-evidence.tsv').df
    expected = json.loads((directory / 'sid-expected.json').read_text())
    for name, records in expected.items():
        merged = load_vaxrank_config(config_path=['builtin:openvax-v1', 'builtin:' + name])
        cfg = msgspec.convert(extract_epitope_config_kwargs(merged), EpitopeConfig)
        scored = score_predictions([], cfg, topiary_df=frame).rename('score').reset_index()
        pd.testing.assert_frame_equal(scored, pd.DataFrame(records), check_dtype=False)


def test_named_report_combines_sources_without_copying_evidence_or_reapplying_minimum(tmp_path):
    from vaxrank.epitope_io import write_neoepitope_report

    class OpaqueEvidence:
        def __deepcopy__(self, memo):
            raise AssertionError('Report rows must not copy their evidence graph')

    first = _predictions_df([dict(peptide='SIINFEKL', allele='HLA-A*02:01', value=100.)]).assign(prediction_id='first')
    second = first.assign(prediction_id='second', source_sequence_name='second')
    evidence = pd.concat([first, second], ignore_index=True)
    report = evidence[['prediction_id', 'peptide', 'peptide_offset', 'allele']].rename(columns={
        'prediction_id': 'Prediction identity', 'peptide': 'Mutant peptide sequence',
        'peptide_offset': 'Peptide offset', 'allele': 'Allele'})
    report.attrs = {'input_score_frames': [first, second], 'evidence': OpaqueEvidence()}
    cfg = EpitopeConfig(selection_policy=dict(name='zero-is-eligible', score_by='0.0', min_score=None))
    output = tmp_path / 'report.csv'
    write_neoepitope_report(report, [], topiary_df=evidence, epitope_config=cfg,
                            csv_report_path=output)
    actual = pd.read_csv(output)
    assert len(actual) == 2
    assert actual.vaxrank_rank_eligible.all()
    assert (actual.vaxrank_score == 0.).all()
