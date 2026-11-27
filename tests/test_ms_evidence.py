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


@pytest.mark.parametrize('minimum', [0., 1e-5])
def test_default_profile_without_ms_matches_existing_scores(minimum):
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import score_predictions
    frame = pd.DataFrame([dict(prediction_id=str(i), peptide='SIINFEKLA', peptide_offset=0,
        allele=A, kind='pMHC_affinity', prediction_method_name='fixture',
        predictor_version='1', value=value, affinity=value, score=None, percentile_rank=None)
        for i, value in enumerate([1., 50., 350., 500., 1500., 5000., 50000.])])
    cfg = EpitopeConfig(min_epitope_score=minimum)
    old = score_predictions([], cfg, topiary_df=frame)
    old = old.loc[old >= cfg.min_epitope_score]
    new = score_predictions([], EpitopeConfig(ms_evidence={}, min_epitope_score=minimum), topiary_df=frame)
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


def test_normalize_observation_gates_hitlist_perturbed_and_excluded_elutions():
    from vaxrank.ms_lookup import normalize_observation

    base = dict(
        peptide='SIINFEKL',
        has_ptm=False,
        is_non_peptide_ligand=False,
        is_potential_contaminant=False,
        mhc_species='Homo sapiens',
        mhc_class='I',
        mhc_restriction=A,
        restriction_evidence='monoallelic',
        is_monoallelic=True,
        pmid=33633747,
        provenance_id='obs-1',
        source='iedb',
    )
    # Unperturbed baseline row remains eligible.
    unperturbed = normalize_observation(dict(base, condition_category='unperturbed', perturbation='', apm_perturbed=False))
    assert unperturbed['eligible'] is True
    assert unperturbed['exclusion_reasons'] == []
    assert unperturbed['condition_category'] == 'unperturbed'

    # Drug-treated (e.g. flucloxacillin-only) arm is excluded from baseline presentation boost.
    drug = normalize_observation(dict(base, condition_category='drug_treatment', perturbation='flucloxacillin'))
    assert drug['eligible'] is False
    assert 'perturbed_condition' in drug['exclusion_reasons']
    assert drug['perturbation'] == 'flucloxacillin'

    # APC machinery perturbation (e.g. HLA-DM / ERAP1 knockout) is excluded.
    apm = normalize_observation(dict(base, condition_category='HLA_DM_perturbation', apm_perturbed=True))
    assert apm['eligible'] is False
    assert 'apm_perturbed' in apm['exclusion_reasons']
    assert 'perturbed_condition' in apm['exclusion_reasons']

    # Explicit control arm from a perturbation study remains eligible.
    control = normalize_observation(dict(base, condition_category='drug_treatment', is_control_arm=True))
    assert control['eligible'] is True

    # Excluded elution condition (e.g. blank / filler assay) is gated out.
    excluded_elution = normalize_observation(dict(base, sample_attribution='elution_conditions_excluded'))
    assert excluded_elution['eligible'] is False
    assert 'excluded_elution_condition' in excluded_elution['exclusion_reasons']

    # Perturbed rows still do not award positive MS boost in support_features.
    feat = support_features([drug, apm, excluded_elution], A)
    assert feat['resolved'] == 0 and feat['presence'] == 0 and feat['reference'] == 0 and feat['direct'] is False

    # Hitlist 2.3+ flat intervention columns gate drug/APM/cytokine arms even
    # when coarse condition_category is 'unperturbed' or is_control_arm is True.
    flat_drug = normalize_observation(dict(
        base,
        condition_category='unperturbed',
        condition_drugs='flucloxacillin',
        elution_condition_id='elc-fluclox-1',
    ))
    assert flat_drug['eligible'] is False
    assert 'perturbed_condition' in flat_drug['exclusion_reasons']
    assert flat_drug['condition_drugs'] == 'flucloxacillin'
    assert flat_drug['source_context']['elution_condition_id'] == 'elc-fluclox-1'

    flat_apm = normalize_observation(dict(
        base,
        condition_category='unperturbed',
        is_control_arm=True,
        condition_overexpression_genes='HLA-DM',
    ))
    assert flat_apm['eligible'] is False
    assert 'apm_perturbed' in flat_apm['exclusion_reasons']
    assert 'perturbed_condition' in flat_apm['exclusion_reasons']
    assert flat_apm['condition_overexpression_genes'] == 'HLA-DM'

    flat_cytokine = normalize_observation(dict(
        base,
        perturbation_type='cytokine',
        perturbation_agent='IFN-gamma',
    ))
    assert flat_cytokine['eligible'] is False
    assert 'apm_perturbed' in flat_cytokine['exclusion_reasons']
    assert 'perturbed_condition' in flat_cytokine['exclusion_reasons']

    # Positive absence tokens ('none') in flat intervention columns stay eligible.
    flat_none = normalize_observation(dict(
        base,
        condition_category='unperturbed',
        condition_drugs='none',
        condition_cytokines='none',
    ))
    assert flat_none['eligible'] is True
    assert flat_none['exclusion_reasons'] == []


class _FakeEvidenceStore:
    """Minimal verified EvidenceStore test double exercising the v2 typed contract."""

    def __init__(self, manifest, tables):
        self.manifest = manifest
        self._tables = tables
        self.bindings = {name: object() for name in tables}

    def query(self, evidence_type, *, where=None, limit=None):
        df = self._tables[evidence_type].copy()
        if where:
            for key, value in where.items():
                if isinstance(value, (list, tuple, set)):
                    df = df.loc[df[key].isin(value)]
                else:
                    df = df.loc[df[key] == value]
        return df


def _write_fake_store_dir(tmp_path, manifest):
    import json
    pkg = tmp_path / 'evidence_store_pkg'
    pkg.mkdir(parents=True, exist_ok=True)
    (pkg / 'manifest.json').write_text(json.dumps(manifest, sort_keys=True), encoding='utf-8')
    return pkg


def test_evidence_store_reanalysis_counts_once_and_independent_origins_count_separately(tmp_path, monkeypatch):
    """#820: one reviewed origin across original + reanalysis counts once; two independent origins count separately."""
    import vaxrank.ms_lookup as ms_lookup

    manifest = {
        'schema_version': 6,
        'package_family': 'supporting',
        'files': {'typed/mhc_observations.parquet': 'deadbeef'},
    }
    pkg = _write_fake_store_dir(tmp_path, manifest)

    tables = {
        'mhc_positive_ms': pd.DataFrame([
            # Abelin 2017 original monoallelic observation
            dict(
                observation_id='obs-abelin-1',
                peptide_sequence='SIINFEKLA',
                has_ptm=False,
                is_non_peptide_ligand=False,
                is_potential_contaminant=False,
                source_dataset='iedb',
                source_assay_id='iedb:101',
                source_pmid='28228285',
                source_reference_iri='https://doi.org/10.1016/abelin2017',
                peptide_mhc_restriction=A,
                restriction_evidence='monoallelic',
            ),
            # Sarkizova 2020 reanalysis of the exact same origin (different PMID, different acquisition)
            dict(
                observation_id='obs-sarkizova-reanalysis',
                peptide_sequence='SIINFEKLA',
                has_ptm=False,
                is_non_peptide_ligand=False,
                is_potential_contaminant=False,
                source_dataset='iedb',
                source_assay_id='iedb:202',
                source_pmid='31844290',
                source_reference_iri='https://doi.org/10.1038/sarkizova2020',
                peptide_mhc_restriction=A,
                restriction_evidence='monoallelic',
            ),
            # Independent second monoallelic origin
            dict(
                observation_id='obs-independent-origin',
                peptide_sequence='SIINFEKLA',
                has_ptm=False,
                is_non_peptide_ligand=False,
                is_potential_contaminant=False,
                source_dataset='iedb',
                source_assay_id='iedb:303',
                source_pmid='33333333',
                source_reference_iri='https://doi.org/10.1038/other2021',
                peptide_mhc_restriction=A,
                restriction_evidence='monoallelic',
            ),
        ]),
        'assay_context': pd.DataFrame([
            dict(observation_id='obs-abelin-1', arm_id='arm-abelin', condition_id='721_221_a0201',
                 arm_link_status='linked', cell_line_name='721.221', sample_attribution='reviewed'),
            dict(observation_id='obs-sarkizova-reanalysis', arm_id='arm-sarkizova', condition_id='721_221_a0201',
                 arm_link_status='linked', cell_line_name='721.221', sample_attribution='reviewed'),
            dict(observation_id='obs-independent-origin', arm_id='arm-other', condition_id='k562_a0201',
                 arm_link_status='linked', cell_line_name='K562', sample_attribution='reviewed'),
        ]),
        'reviewed_lineage': pd.DataFrame([
            dict(claim_id='c1', arm_id='arm-abelin', lineage_axis='experimental_origin',
                 lineage_status='resolved', identity_ids=['origin:MSV000080527:721_221_hla_a_02_01']),
            dict(claim_id='c2', arm_id='arm-sarkizova', lineage_axis='experimental_origin',
                 lineage_status='resolved', identity_ids=['origin:MSV000080527:721_221_hla_a_02_01']),
            dict(claim_id='c3', arm_id='arm-other', lineage_axis='experimental_origin',
                 lineage_status='resolved', identity_ids=['origin:MSV000099999:k562_hla_a_02_01']),
        ]),
        'arm_mhc_claim': pd.DataFrame([
            dict(claim_id='m1', arm_id='arm-abelin', mhc_claim_role='monoallelic_host',
                 evidence_basis='monoallelic', reported_mhc=A, cellular_mhc_genotype=[A],
                 genotype_complete_loci=['HLA-A']),
            dict(claim_id='m2', arm_id='arm-sarkizova', mhc_claim_role='monoallelic_host',
                 evidence_basis='monoallelic', reported_mhc=A, cellular_mhc_genotype=[A],
                 genotype_complete_loci=['HLA-A']),
            dict(claim_id='m3', arm_id='arm-other', mhc_claim_role='monoallelic_host',
                 evidence_basis='monoallelic', reported_mhc=A, cellular_mhc_genotype=[A],
                 genotype_complete_loci=['HLA-A', 'HLA-B', 'HLA-C']),
        ]),
    }

    fake_store = _FakeEvidenceStore(manifest, tables)
    monkeypatch.setattr(ms_lookup, 'verified_evidence_store', lambda d, h: fake_store)

    evidence = ms_lookup.capture_package(pkg, ['SIINFEKLA'], [A])
    assert evidence['source_contract'] == 'verified_hitlist_evidence_store'
    assert evidence['package_schema_version'] == 6
    assert evidence['index_metadata_sha256']

    rows = evidence['species_observations']['Homo sapiens']['rows']
    assert len(rows) == 3
    normalized = [ms_lookup.normalize_observation(r) for r in rows]

    # Abelin + Sarkizova alone: repeated reanalysis of one origin counts ONCE in R and ONCE in J.
    reanalysis_pair = [n for n in normalized if n['observation_id'] in ('obs-abelin-1', 'obs-sarkizova-reanalysis')]
    feat_pair = support_features(reanalysis_pair, A)
    assert feat_pair['resolved'] == 8.0
    assert feat_pair['presence'] == 8.0
    assert feat_pair['reference'] == pytest.approx(0.25 * 8.0)
    assert len(feat_pair['units']) == 1
    assert feat_pair['direct'] is True

    # Adding the independent second origin counts separately in R (8 + 8 = 16).
    feat_all = support_features(normalized, A)
    assert feat_all['resolved'] == 16.0
    assert len(feat_all['units']) == 2
    assert feat_all['reference'] == pytest.approx(0.25 * 16.0)

    # Coverage diagnostics preserve resolved unit counts and unavailable tables.
    diag = evidence['coverage_diagnostics']
    assert diag['total_observations'] == 3
    assert diag['eligible_observations'] == 3
    assert diag['resolved_engineered_origin_observations'] == 3
    assert diag['unique_resolved_units'] == {'native_donor': 0, 'engineered_origin': 2}
    assert 'arm_factor' in diag['unavailable_facets']
    assert 'reviewed_lineage' in diag['available_tables']


def test_evidence_store_keeps_unknown_and_ambiguous_lineage_unresolved_and_gates_perturbations(tmp_path, monkeypatch):
    """#820: unresolved/pooled lineage never inflates R; perturbed arm factors are gated while control arms stay eligible."""
    import json
    import vaxrank.ms_lookup as ms_lookup

    manifest = {'schema_version': 7, 'package_family': 'reference'}
    pkg = _write_fake_store_dir(tmp_path, manifest)

    tables = {
        'mhc_positive_ms': pd.DataFrame([
            # Unresolved native donor with two pooled donor IDs -> must stay unresolved (R=0)
            dict(
                observation_id='obs-pooled',
                peptide_sequence='SIINFEKLA',
                has_ptm=False,
                is_non_peptide_ligand=False,
                is_potential_contaminant=False,
                source_dataset='hla_ligand_atlas',
                source_pmid='33858848',
                peptide_mhc_restriction=A,
                restriction_evidence='predicted',
            ),
            # Resolved native donor in perturbed drug arm -> excluded from baseline boost
            dict(
                observation_id='obs-drug-arm',
                peptide_sequence='SIINFEKLA',
                has_ptm=False,
                is_non_peptide_ligand=False,
                is_potential_contaminant=False,
                source_dataset='iedb',
                source_pmid='44444444',
                peptide_mhc_restriction=A,
                restriction_evidence='experimental',
            ),
            # Resolved native donor in explicit control arm -> remains eligible and increments R
            dict(
                observation_id='obs-control-arm',
                peptide_sequence='SIINFEKLA',
                has_ptm=False,
                is_non_peptide_ligand=False,
                is_potential_contaminant=False,
                source_dataset='iedb',
                source_pmid='44444444',
                peptide_mhc_restriction=A,
                restriction_evidence='experimental',
            ),
        ]),
        'assay_context': pd.DataFrame([
            dict(observation_id='obs-pooled', arm_id='arm-pooled', condition_id='donor_pool',
                 arm_link_status='linked', cell_type='PBMC'),
            dict(observation_id='obs-drug-arm', arm_id='arm-drug', condition_id='fluclox_treated',
                 arm_link_status='linked', cell_type='B-LCL'),
            dict(observation_id='obs-control-arm', arm_id='arm-ctrl', condition_id='vehicle_ctrl',
                 arm_link_status='linked', cell_type='B-LCL'),
        ]),
        'reviewed_lineage': pd.DataFrame([
            dict(claim_id='l1', arm_id='arm-pooled', lineage_axis='donor',
                 lineage_status='resolved', identity_ids=['donor:D1', 'donor:D2']),
            dict(claim_id='l2', arm_id='arm-drug', lineage_axis='donor',
                 lineage_status='resolved', identity_ids=['donor:D3']),
            dict(claim_id='l3', arm_id='arm-ctrl', lineage_axis='donor',
                 lineage_status='resolved', identity_ids=['donor:D4']),
        ]),
        'arm_mhc_claim': pd.DataFrame([
            dict(claim_id='mc1', arm_id='arm-pooled', mhc_claim_role='genotype',
                 evidence_basis='predicted', cellular_mhc_genotype=[A, B],
                 genotype_complete_loci=['HLA-A', 'HLA-B']),
            dict(claim_id='mc2', arm_id='arm-drug', mhc_claim_role='restriction',
                 evidence_basis='experimental', cellular_mhc_genotype=[A],
                 genotype_complete_loci=['HLA-A', 'HLA-B', 'HLA-C']),
            dict(claim_id='mc3', arm_id='arm-ctrl', mhc_claim_role='restriction',
                 evidence_basis='experimental', cellular_mhc_genotype=[A],
                 genotype_complete_loci=['HLA-A', 'HLA-B', 'HLA-C']),
        ]),
        'arm_factor': pd.DataFrame([
            dict(claim_id='f1', arm_id='arm-drug', factor_type='drug_treatment',
                 factor_value='flucloxacillin', claim_status='present'),
            dict(claim_id='f2', arm_id='arm-ctrl', factor_type='drug_treatment',
                 factor_value='flucloxacillin', claim_status='present'),
        ]),
        'arm_control_link': pd.DataFrame([
            dict(claim_id='cl1', control_arm_id='arm-ctrl', comparison_arm_id='arm-drug'),
        ]),
    }

    fake_store = _FakeEvidenceStore(manifest, tables)
    monkeypatch.setattr(ms_lookup, 'verified_evidence_store', lambda d, h: fake_store)

    lookup = ms_lookup.MSLookup({'package': str(pkg), 'lookup': 'required'})
    frame = pd.DataFrame([dict(peptide='SIINFEKLA', allele=A)])
    lookup(frame)
    records = json.loads(frame[ms_lookup.RECORDS_COLUMN].iloc[0])
    by_id = {r['observation_id']: r for r in records}

    # Pooled donors -> unresolved (no unit_id), genotype uncertainty='partial', restriction_basis='predicted'.
    assert 'unit_id' not in by_id['obs-pooled']
    assert by_id['obs-pooled']['restriction_basis'] == 'predicted'
    assert by_id['obs-pooled']['genotype_uncertainty'] == 'partial'
    assert by_id['obs-pooled']['eligible'] is True
    assert support_features([by_id['obs-pooled']], A, prediction_qualifies=True)['resolved'] == 0.0
    assert support_features([by_id['obs-pooled']], A, prediction_qualifies=True)['presence'] == 3.0

    # Drug-treated arm is excluded; explicit control arm remains eligible with native_donor unit_id.
    assert by_id['obs-drug-arm']['eligible'] is False
    assert 'perturbed_condition' in by_id['obs-drug-arm']['exclusion_reasons']
    assert by_id['obs-control-arm']['eligible'] is True
    assert by_id['obs-control-arm']['unit_id'] == 'donor:D4'
    assert by_id['obs-control-arm']['unit_kind'] == 'native_donor'

    # Combined support: only the eligible control arm contributes to resolved R (experimental=6.0).
    feat = support_features(records, A, prediction_qualifies=True)
    assert feat['resolved'] == 6.0
    assert feat['direct'] is True

    # Replay query preserves package identity & coverage diagnostics.
    assert len(lookup.queries) == 1
    q = lookup.queries[0]
    assert q['source_contract'] == 'verified_hitlist_evidence_store'
    assert q['coverage_diagnostics']['excluded_observations'] == 1
    assert q['coverage_diagnostics']['unresolved_lineage_observations'] == 1
    assert q['coverage_diagnostics']['unique_resolved_units'] == {'native_donor': 2, 'engineered_origin': 0}


def test_evidence_store_cross_species_and_missing_manifest_rejection(tmp_path, monkeypatch):
    """#820: EvidenceStore lineage resolution works cross-species via mhcgnomes and rejects invalid packages."""
    import vaxrank.ms_lookup as ms_lookup

    with pytest.raises(ValueError, match='missing manifest.json'):
        ms_lookup.capture_package(tmp_path / 'empty_dir', ['SIINFEKL'], ['H2-Kb'])

    manifest = {'schema_version': 6, 'package_family': 'supporting'}
    pkg = _write_fake_store_dir(tmp_path, manifest)
    tables = {
        'mhc_positive_ms': pd.DataFrame([
            dict(
                observation_id='obs-mouse-1',
                peptide_sequence='SIINFEKL',
                has_ptm=False,
                is_non_peptide_ligand=False,
                is_potential_contaminant=False,
                source_dataset='iedb',
                source_pmid='11111111',
                peptide_mhc_restriction='H2-Kb',
                restriction_evidence='experimental',
            ),
        ]),
        'assay_context': pd.DataFrame([
            dict(observation_id='obs-mouse-1', arm_id='arm-m1', condition_id='c57bl6_spleen',
                 arm_link_status='linked', cell_type='splenocyte'),
        ]),
        'reviewed_lineage': pd.DataFrame([
            dict(claim_id='lm1', arm_id='arm-m1', lineage_axis='donor',
                 lineage_status='resolved', identity_ids=['donor:mouse-c57bl6-1']),
        ]),
    }
    monkeypatch.setattr(ms_lookup, 'verified_evidence_store', lambda d, h: _FakeEvidenceStore(manifest, tables))

    ev = ms_lookup.capture_package(pkg, ['SIINFEKL'], ['H2-Kb'])
    rows = ev['species_observations']['Mus musculus']['rows']
    assert len(rows) == 1
    norm = ms_lookup.normalize_observation(rows[0])
    assert norm['mhc_species'] == 'Mus musculus'
    assert norm['mhc_class'] == 'I'
    assert norm['unit_id'] == 'donor:mouse-c57bl6-1'
    assert norm['unit_kind'] == 'native_donor'
    feat = support_features([norm], 'H2-Kb')
    assert feat['resolved'] == 6.0
    assert feat['direct'] is True
