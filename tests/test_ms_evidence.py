"""Evidence identity and the published MS power-transform contract."""

import math

import pandas as pd
import pytest
from topiary import resolve_selection_policy, TopiaryResult, evaluate_selection_policy

from vaxrank.ms_evidence import MSParameters, ms_criteria, support_features

A = 'HLA-A*02:01'
B = 'HLA-B*07:02'


def observation(**values):
    return dict(dict(observation_id='o1', peptide='SIINFEKLA', eligible=True,
                     restriction=A, restriction_basis='unknown', references=['pmid:1']), **values)


def adjusted(base, features, parameters=None):
    parameters = parameters or MSParameters()
    frame = pd.DataFrame([dict(
        prediction_id='one', peptide='SIINFEKLA', allele=A, peptide_offset=0,
        base=base, **{'ms_' + k: features[k] for k in ('resolved', 'presence', 'reference')})])
    policy = resolve_selection_policy(dict(name='ms-evidence-v1',
        criteria=ms_criteria(parameters, 'base'), score_by='criterion("ms_adjusted")'))
    result = evaluate_selection_policy(TopiaryResult(frame, form='long'), policy,
        group_keys=['prediction_id', 'peptide', 'peptide_offset', 'allele'])
    return result.occurrences.score.iloc[0]


def test_recurrence_is_resolved_lineage_not_rows_or_labels():
    one = observation(typed_alleles=[A], unit_id='donor:1', unit_kind='native_donor')
    features = support_features([one], A, prediction_qualifies=True)
    assert features['resolved'] == features['presence'] == 3
    assert adjusted(.2, features) == pytest.approx(.2 ** (1 / 3))
    assert support_features([one] * 20, A, prediction_qualifies=True) == features
    donors = [dict(one, observation_id=str(i), unit_id=f'donor:{i}') for i in range(20)]
    twenty = support_features(donors, A, prediction_qualifies=True)
    assert twenty['resolved'] == 60
    assert adjusted(.2, twenty) == pytest.approx(.793, abs=.0005)
    unresolved = [dict(one, unit_id=None, sample_label=str(i)) for i in range(20)]
    assert support_features(unresolved, A, prediction_qualifies=True)['resolved'] == 0
    assert support_features(unresolved, A, prediction_qualifies=True)['presence'] == 3


def test_filtered_bridge_and_shared_donor_prevent_paper_bonuses():
    rows = [observation(references=[f'pmid:{i}']) for i in range(20)]
    assert support_features(rows, A)['reference'] == 5
    connected = [dict(o, unit_id='donor:1', unit_kind='native_donor') for o in rows]
    assert support_features(connected, A)['reference'] == .25
    bridge = observation(references=[f'pmid:{i}' for i in range(20)], quality=0)
    assert support_features(rows + [bridge], A)['reference'] == .25


def test_exact_name_and_inferred_singleton_do_not_establish_restriction():
    row = observation()
    f = support_features([row], A)
    assert f['presence'] == 1 and not f['direct']
    assert support_features([dict(row, restriction='HLA-A*02')], A)['presence'] == .25
    mono = dict(row, monoallelic_verified=True)
    assert support_features([mono], A)['presence'] == 8
    assert support_features([mono], A)['direct']
    assert support_features([mono], B)['presence'] == .25
    assert not support_features([mono], B)['direct']
    assert not support_features([dict(mono, quality=0)], A)['direct']


def test_qualification_and_typing_stay_allele_specific():
    row = observation(restriction=None, typed_alleles=[A])
    assert support_features([row], A)['presence'] == 1
    assert support_features([row], A, prediction_qualifies=True)['presence'] == 3
    assert support_features([row], B, prediction_qualifies=True)['presence'] == .5


def test_uncapped_transform_and_explicit_cap_keep_zero_one_and_missing():
    low = dict(resolved=255, presence=3, reference=.75)
    high = dict(resolved=3000, presence=3, reference=.75)
    assert adjusted(.2, high) > adjusted(.2, low)
    assert adjusted(.2, high) == pytest.approx(.880, abs=.0005)
    assert adjusted(0, high) == 0
    assert adjusted(1, high) == 1
    assert math.isnan(adjusted(float('nan'), high))
    assert adjusted(.2, high, MSParameters(gain=0)) == .2
    assert adjusted(.2, high, MSParameters(max_relaxation=1)) == .2
    assert adjusted(.2, high, MSParameters(max_relaxation=3)) == pytest.approx(.2 ** (1/3))
    assert adjusted(.05, high) > .1 > adjusted(.05, dict(resolved=0, presence=0, reference=0))


@pytest.mark.parametrize('options', [dict(gain=-1), dict(gain=float('nan')),
    dict(max_relaxation=.5), dict(max_relaxation=float('inf')), dict(reference_discount=2)])
def test_invalid_parameters_are_rejected(options):
    with pytest.raises(ValueError):
        MSParameters(**options)


@pytest.mark.parametrize('allele', ['HLA-A*02', 'HLA-A*02:01N', 'HLA-A*02:01/A*02:02'])
def test_patient_allele_boundary_matches_topiary(allele):
    from vaxrank.allele_validation import parse_prediction_allele
    with pytest.raises(ValueError):
        parse_prediction_allele(allele)


def test_frame_policy_crosses_soft_gate_but_preserves_hard_filters():
    import json
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.ms_evidence import score_frame
    keys = ['prediction_id', 'peptide', 'peptide_offset', 'allele']
    row = dict(prediction_id='candidate', peptide='SIINFEKLA', peptide_offset=0,
               allele=A, kind='pMHC_affinity', prediction_method_name='fixture',
               predictor_version='1', value=900., score=None, percentile_rank=None,
               base=.05, allowed=1)
    obs = observation(mhc_species='Homo sapiens', mhc_class='I', monoallelic_verified=True)
    frame = pd.DataFrame([dict(row, hitlist_ms_records_json=json.dumps([obs]))])
    cfg = EpitopeConfig(ms_evidence=dict(base_score='base', soft_threshold=.1))
    values = score_frame(frame, cfg, group_keys=keys)
    assert len(values) == 1 and values.iloc[0] > .1
    trace = values.attrs['ms_scoring'][('candidate', 'SIINFEKLA', 0)][A]
    assert trace['baseline'] == 0 and trace['base'] == .05
    assert score_frame(frame, EpitopeConfig(ms_evidence=cfg.ms_evidence,
        filter_expr='allowed == 0'), group_keys=keys).empty
    no_gain = dict(cfg.ms_evidence, parameters={'gain': 0})
    assert score_frame(frame, EpitopeConfig(ms_evidence=no_gain), group_keys=keys).empty
    frame['hitlist_ms_records_json'] = '[]'
    assert score_frame(frame, cfg, group_keys=keys).empty


def test_default_profile_without_ms_matches_existing_scores():
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import score_predictions
    frame = pd.DataFrame([dict(prediction_id=str(i), peptide='SIINFEKLA', peptide_offset=0,
        allele=A, kind='pMHC_affinity', prediction_method_name='fixture',
        predictor_version='1', value=value, affinity=value, score=None, percentile_rank=None)
        for i, value in enumerate([1., 50., 350., 500., 1500., 5000., 50000.])])
    cfg = EpitopeConfig()
    old = score_predictions([], cfg, topiary_df=frame)
    old = old.loc[old >= cfg.min_epitope_score]
    new = score_predictions([], EpitopeConfig(ms_evidence={}), topiary_df=frame)
    pd.testing.assert_series_equal(old, new, check_names=False, check_exact=False)


def test_ms_only_score_requires_direct_evidence_and_explicit_policy():
    import json
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.ms_evidence import score_frame
    keys = ['prediction_id', 'peptide', 'peptide_offset', 'allele']
    obs = observation(mhc_species='Homo sapiens', mhc_class='I', monoallelic_verified=True)
    frame = pd.DataFrame([dict(prediction_id='candidate', peptide='SIINFEKLA', peptide_offset=0,
        allele=A, kind='pMHC_affinity', prediction_method_name='fixture', predictor_version='1',
        value=None, score=None, percentile_rank=None, base=float('nan'),
        hitlist_ms_records_json=json.dumps([obs]))])
    options = dict(base_score='base', soft_threshold=.1)
    assert score_frame(frame, EpitopeConfig(ms_evidence=options), group_keys=keys).empty
    options['ms_only_score'] = .7
    result = score_frame(frame, EpitopeConfig(ms_evidence=options), group_keys=keys)
    assert result.iloc[0] == .7
    frame['base'] = 0.
    assert score_frame(frame, EpitopeConfig(ms_evidence=options), group_keys=keys).iloc[0] == 0.


def test_window_credits_ms_once_across_alleles_and_repeated_positions(tmp_path):
    from dataclasses import replace
    from .test_vaccine_peptide_config import _make_fragment, _make_mutant_epitope
    from vaxrank.vaccine_peptide import VaccinePeptide
    from vaxrank.epitope_io import save_predictions, load_predictions
    e = replace(_make_mutant_epitope(), per_allele_scores={A: .2, B: .1},
        ms_scoring={A: dict(base=.2, baseline=.2, adjusted=.7),
                    B: dict(base=.1, baseline=.1, adjusted=.8)})
    assert e.epitope_score == pytest.approx(1.)
    repeated = replace(e, offset=10)
    window = VaccinePeptide(mutant_protein_fragment=_make_fragment(), epitopes=[e, repeated])
    assert window.target_epitope_score == pytest.approx(.6 + .7)
    assert replace(e, per_allele_scores={A: .2}).epitope_score == pytest.approx(.7)
    assert e.sliced(0, 20).ms_scoring == e.ms_scoring
    path = tmp_path / 'epitopes.tsv'
    save_predictions([e], path)
    saved = load_predictions(path)[0]
    assert saved.ms_scoring == e.ms_scoring
    assert saved.epitope_score == e.epitope_score


def test_hitlist_normalization_preserves_uncertainty_and_chemical_form():
    from vaxrank.ms_lookup import normalize_observation
    row = dict(peptide='SIINFEKLA', has_ptm=False, pmid='123', assay_iri='iedb:1',
        mhc_species='Homo sapiens', mhc_class='I', mhc_restriction=A,
        is_monoallelic=True, restriction_evidence='unknown', mhc_allele_set=[A],
        mhc_allele_provenance='restriction', sample='same descriptive label')
    n = normalize_observation(row)
    assert n['eligible'] and not n['monoallelic_verified']
    assert not n['typed_alleles'] and 'unit_id' not in n
    assert not normalize_observation(dict(row, has_ptm=None))['eligible']
    assert not normalize_observation(dict(row, has_ptm=True))['eligible']
    assert normalize_observation(dict(row, restriction_evidence='monoallelic'))['monoallelic_verified']
    unknown = normalize_observation(dict(row, mhc_class=None, restriction_evidence='monoallelic'))
    assert not unknown['monoallelic_verified'] and not unknown['restriction']


def test_ms_config_cli_and_unknown_fields():
    from vaxrank.cli.arg_parser import make_vaxrank_arg_parser
    from vaxrank.cli.epitope_config_args import epitope_config_from_args
    parser = make_vaxrank_arg_parser()
    cfg = epitope_config_from_args(parser.parse_args(['--bam', 'fixture.bam', '--mhc-predictor', 'random', '--ms-evidence-package', '/tmp/test',
        '--config-value', 'evidence.ms.parameters.gain=2']))
    assert cfg.ms_evidence == {'package': '/tmp/test', 'parameters': {'gain': 2}}


def test_enumeration_retains_missing_candidates_without_invented_predictions():
    from types import SimpleNamespace
    from vaxrank.ms_lookup import complete_candidate_frame
    from vaxrank.epitope_dataset import EpitopeDataset
    model = SimpleNamespace(alleles=[A, B], default_peptide_lengths=[9])
    frame = complete_candidate_frame(pd.DataFrame(), {'protein': 'SIINFEKLAA'}, model)
    assert len(frame) == 4
    assert frame[['value', 'score', 'kind', 'prediction_method_name']].isna().all().all()
    frame['sample_name'] = 'patient'
    frame['source_sequence'] = 'SIINFEKLAA'
    data = EpitopeDataset.from_topiary(TopiaryResult(frame, form='long'), label='sample')
    assert len(data.epitopes) == 2
    assert all(not e.predictions_flat() for e in data.epitopes)
    assert all(e.patient_alleles == (A, B) for e in data.epitopes)
    assert len(complete_candidate_frame(frame, {'protein': 'SIINFEKLAA'}, model)) == 4
