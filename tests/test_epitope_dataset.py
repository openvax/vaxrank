"""File and CLI regression coverage for evidence-preserving inputs."""
from dataclasses import replace

import pandas as pd
import pytest
from topiary import TopiaryResult, combine_sources

from vaxrank.cli.entry_point import run_cli
from vaxrank.candidate_epitope import CandidateEpitope
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dataset import EpitopeDataset, ANTIGEN_COLUMN
from vaxrank.epitope_dsl import attach_per_allele_scores
from vaxrank.epitope_io import load_predictions, save_predictions
from vaxrank.external_input import ExternalConstructOptions
from vaxrank.external_rescoring import prepare_reports
from vaxrank.native_serialization import to_native_json
from vaxrank.vaccine_antigen import (
    AminoAcidInterval, TargetableMask, TumorSpecificityAttestation, VaccineAntigen,
)

ALLELE = 'HLA-A*02:01'


def evidence(**columns):
    return TopiaryResult(pd.DataFrame({
        'peptide': ['SIINFEKL', 'GILGFVFTL'], 'allele': [ALLELE] * 2,
        'kind': ['pMHC_affinity'] * 2, 'value': [50., 500.],
        'prediction_method_name': ['model'] * 2, 'predictor_version': ['01'] * 2,
        'score': [None, 0.], 'n_flank': [None, ''], 'c_flank': ['', None],
        'prediction_id': ['original-one', 'original-two'],
        'details': [{'count': 0, 'known': False}, {'names': ['001', 'NA']}],
        'new_feature': [1., 100.], **columns,
    }), extra={'producer': {'run': '001', 'parameters': [False, None, 0]}})


def forbidden(*args, **kwargs):
    raise AssertionError('Original-only input must not instantiate a predictor')


@pytest.mark.parametrize('suffix', ['csv', 'tsv'])
def test_typed_evidence_and_config_roundtrip(suffix, tmp_path):
    dataset = EpitopeDataset.from_topiary(evidence(), sample_name='patient')
    dataset.config = EpitopeConfig(score_expr='new_feature / affinity.value', min_epitope_score=0)
    dataset.epitopes = tuple(attach_per_allele_scores(
        dataset.epitopes, dataset.config, topiary_df=dataset.scoring_frame()))
    path = tmp_path / ('native.' + suffix)
    dataset.save(path)
    loaded = EpitopeDataset.load(path)
    pd.testing.assert_frame_equal(
        loaded.result.df[dataset.result.df.columns], dataset.result.df, check_dtype=False)
    assert loaded.result.extra == dataset.result.extra
    assert loaded.config == dataset.config
    assert [to_native_json(e) for e in loaded.epitopes] == [to_native_json(e) for e in dataset.epitopes]
    assert loaded.result.df.prediction_id.tolist() == ['original-one', 'original-two']
    assert loaded.scoring_frame().prediction_id.tolist() == loaded.result.df.source_observation_id.tolist()
    rescored = attach_per_allele_scores(loaded.epitopes, loaded.config, topiary_df=loaded.scoring_frame())
    assert [e.per_allele_scores for e in rescored] == [e.per_allele_scores for e in dataset.epitopes]
    assert load_predictions(path) == list(loaded.epitopes)


@pytest.mark.parametrize('missing', [float('nan'), pd.NA])
def test_missing_flanks_become_native_nulls_without_losing_empty_strings(missing, tmp_path):
    source = evidence(n_flank=[missing, ''], c_flank=['', missing],
                      wt_peptide=['SIINFEKA', 'GILGFVFTA'], wt_value=[1000., 2000.],
                      wt_n_flank=[missing, ''], wt_c_flank=['', missing])
    dataset = EpitopeDataset.from_topiary(source, sample_name='patient')
    for candidate, expected in zip(dataset.epitopes, [(None, ''), ('', None)]):
        for peptide in (candidate, candidate.wt):
            assert (peptide.n_flank, peptide.c_flank) == expected
            assert [(p.n_flank, p.c_flank) for p in peptide.predictions_flat()] == [expected]
    path = tmp_path / 'native.tsv'
    dataset.save(path)
    assert load_predictions(path) == list(dataset.epitopes)


def test_simple_table_and_native_cli_reload(tmp_path, monkeypatch):
    import mhctools.cli
    import vaxrank.cli.entry_point
    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', forbidden)
    monkeypatch.setattr(vaxrank.cli.entry_point, 'annotate_predictions_with_processing', forbidden)
    source = tmp_path / 'table.tsv'
    evidence().to_tsv(source)
    native = tmp_path / 'native.tsv'
    first_report = tmp_path / 'first.csv'
    run_cli(['--input-topiary', str(source), '--output-epitopes', str(native),
             '--output-csv', str(first_report),
             '--config-value', 'epitopes.score_expr=new_feature / affinity.value'])
    second_report = tmp_path / 'second.csv'
    run_cli(['--input-epitopes', str(native), '--output-csv', str(second_report)])
    first, second = pd.read_csv(first_report), pd.read_csv(second_report)
    assert first.vaxrank_score.tolist() == second.vaxrank_score.tolist() == [0.2, 0.02]
    assert first['Construction limitation'].str.contains('missing antigen').all()
    assert len(load_predictions(native)) == 2


def test_legacy_native_payload_uses_shared_cli(tmp_path):
    candidates = EpitopeDataset.from_topiary(evidence(), sample_name='patient').epitopes
    candidates = tuple(replace(e, prediction_id='') for e in candidates)
    source = tmp_path / 'legacy.tsv'
    save_predictions(candidates, source)
    output = tmp_path / 'report.csv'
    run_cli(['--input-epitopes', str(source), '--output-csv', str(output),
             '--no-processing-aware-annotation'])
    assert len(pd.read_csv(output)) == 2


@pytest.mark.parametrize('with_predictions', [False, True])
@pytest.mark.parametrize('enriched', [False, True])
def test_native_candidates_without_predictions_survive_reload(
        with_predictions, enriched, tmp_path, monkeypatch):
    import mhctools.cli
    import vaxrank.cli.entry_point
    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', forbidden)
    monkeypatch.setattr(vaxrank.cli.entry_point, 'annotate_predictions_with_processing', forbidden)
    unpredicted = CandidateEpitope(sequence='AAAAAAAA', prediction_id='payload-only')
    candidates = [unpredicted]
    if with_predictions:
        candidates.extend(EpitopeDataset.from_topiary(
            evidence(), sample_name='patient').epitopes)
    source = tmp_path / 'input.tsv'
    if enriched:
        EpitopeDataset.from_predictions(candidates).save(source)
    else:
        save_predictions(candidates, source)
    native = tmp_path / 'output.tsv'
    run_cli(['--input-epitopes', str(source), '--output-epitopes', str(native)])
    loaded = EpitopeDataset.load(native)
    assert len(loaded.epitopes) == len(candidates)
    assert replace(loaded.epitopes[0], input_provenance=None) == unpredicted
    assert len(loaded.result.df) == (2 if with_predictions else 0)
    again = tmp_path / 'again.tsv'
    run_cli(['--input-epitopes', str(native), '--output-epitopes', str(again)])
    assert load_predictions(again) == list(loaded.epitopes)


@pytest.mark.parametrize('kind', ['mutation', 'fusion', 'splice', 'CTA', 'ERV', 'viral'])
def test_explicit_antigens_survive_file_reload_and_build(kind, tmp_path):
    from vaxrank.epitope_dataset import dataset_ranking_result
    antigens = [VaccineAntigen(
        kind=kind, amino_acids=peptide,
        targetable_mask=TargetableMask((AminoAcidInterval(0, len(peptide)),)),
        tumor_specificity=TumorSpecificityAttestation(
            status='admitted', evidence_kind='synthetic_fixture',
            evidence_source='test_epitope_dataset', patient_specific=True,
            rationale_code='test_only'), source_identifier=peptide)
        for peptide in ('SIINFEKL', 'GILGFVFTL')]
    source = tmp_path / 'source.tsv'
    evidence(**{ANTIGEN_COLUMN: [to_native_json(a) for a in antigens]}).to_tsv(source)
    reports, frame, _ = prepare_reports([('topiary', source)],
        EpitopeConfig(score_expr='new_feature / affinity.value'))
    result = dataset_ranking_result(reports[0], reports[0].epitopes,
                                    options=ExternalConstructOptions())
    assert len(result.ranked) == 2
    assert result.ranked[0][1][0].target_epitope_score == pytest.approx(.2)
    native = tmp_path / 'native.tsv'
    frame.attrs['epitope_dataset'].save(native)
    reloaded, _, _ = prepare_reports([('epitopes', native)])
    again = dataset_ranking_result(reloaded[0], reloaded[0].epitopes,
                                   options=ExternalConstructOptions())
    assert [(source, p[0].target_epitope_score) for source, p in again.ranked] == [
        (source, p[0].target_epitope_score) for source, p in result.ranked]


def test_combined_sources_retain_all_observations_and_reject_missing_scope(tmp_path):
    first = tmp_path / 'a.tsv'
    second = tmp_path / 'b.tsv'
    evidence().to_tsv(first)
    evidence(new_feature=[5., 10.]).to_tsv(second)
    with pytest.raises(ValueError, match='requires declared'):
        prepare_reports([('topiary', first), ('topiary', second)])
    combined = combine_sources({'one': evidence(), 'two': evidence()}, sample_name='patient')
    dataset = EpitopeDataset.from_topiary(combined)
    assert len(dataset.epitopes) == 4
    assert len(dataset.result.df) == 4


def test_legacy_native_keeps_saved_custom_scores(tmp_path):
    candidates = EpitopeDataset.from_topiary(evidence(), sample_name='patient').epitopes
    candidates = tuple(replace(e, per_allele_scores={ALLELE: 123.}) for e in candidates)
    source = tmp_path / 'old.tsv'
    save_predictions(candidates, source)
    _, frame, loaded = prepare_reports([('epitopes', source)])
    assert [e.per_allele_scores for e in loaded] == [{ALLELE: 123.}] * 2
    native = tmp_path / 'new.tsv'
    frame.attrs['epitope_dataset'].save(native)
    _, _, again = prepare_reports([('epitopes', native)])
    assert [e.per_allele_scores for e in again] == [e.per_allele_scores for e in loaded]


def test_source_local_model_choices_survive_combined_reload(tmp_path):
    scope = dict(patient_id='patient', reference_assembly='GRCh38', mhc_alleles=[ALLELE])
    inputs = []
    for method, values in [('mhcflurry', [20., 40.]), ('netmhcpan', [200., 400.])]:
        path = tmp_path / (method + '.tsv')
        evidence(prediction_method_name=[method] * 2, value=values).to_tsv(path)
        inputs.append(('topiary', path))
    _, frame, before = prepare_reports(inputs, scopes=[scope, scope])
    dataset = frame.attrs['epitope_dataset']
    assert dataset.result.df.prediction_id.tolist() == ['original-one', 'original-two'] * 2
    path = tmp_path / 'combined.tsv'
    dataset.save(path)
    _, _, after = prepare_reports([('epitopes', path)])
    assert [e.per_allele_scores for e in after] == [e.per_allele_scores for e in before]
    assert len(after) == 4


def test_native_patient_conflict_is_rejected_before_scoring(tmp_path, monkeypatch):
    source = tmp_path / 'input.tsv'
    evidence().to_tsv(source)
    scope = dict(patient_id='patient', reference_assembly='GRCh38', mhc_alleles=[ALLELE])
    _, frame, _ = prepare_reports([('topiary', source)], scopes=[scope])
    native = tmp_path / 'native.tsv'
    frame.attrs['epitope_dataset'].save(native)
    import vaxrank.external_rescoring as module
    monkeypatch.setattr(module, 'attach_per_allele_scores', forbidden)
    with pytest.raises(ValueError, match='conflicting patient_id'):
        prepare_reports([('epitopes', native)], scopes=[{**scope, 'patient_id': 'another'}])


def test_normalized_wildtype_predictions_keep_their_own_provenance(tmp_path):
    dataset = EpitopeDataset.from_topiary(evidence(
        wt_peptide=['SIINFEKA', 'GILGFVFTA'], wt_value=[1000., 2000.],
        wt_score=[None, 0.], wt_prediction_method_name=['wt-model'] * 2,
        wt_predictor_version=['02'] * 2), sample_name='patient')
    path = tmp_path / 'native.tsv'
    dataset.save(path)
    wt = load_predictions(path)[0].wt
    assert wt.sequence == 'SIINFEKA'
    assert wt.predictions_flat()[0].predictor_name == 'wt-model'
    assert wt.predictions_flat()[0].predictor_version == '02'
    assert wt.predictions_flat()[0].score is None


def test_repeated_discoveries_select_once_without_summing_support(tmp_path):
    from types import SimpleNamespace
    from vaxrank.external_input import load_external_ranked
    from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
    from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs
    from .input_scope_helpers import write_input_manifest
    antigen = VaccineAntigen(
        kind='fusion', amino_acids='SIINFEKL',
        targetable_mask=TargetableMask((AminoAcidInterval(0, 8),)),
        tumor_specificity=TumorSpecificityAttestation(
            status='admitted', evidence_kind='synthetic_fixture',
            evidence_source='test_epitope_dataset', patient_specific=True,
            rationale_code='test_only'), source_identifier='fusion')
    paths = []
    for i, value in enumerate([50., 100.]):
        result = evidence(**{ANTIGEN_COLUMN: [to_native_json(antigen), None]})
        result.df = result.df.iloc[:1].copy()
        result.df['value'] = value
        result.df['n_rna_alt'] = 3
        path = tmp_path / f'input-{i}.tsv'
        result.to_tsv(path)
        paths.append(('topiary', str(path)))
    args = SimpleNamespace(input_manifest=write_input_manifest(tmp_path / 'inputs.json', paths),
                           duplicate_candidates='best')
    ranked, frame, candidates, _, _ = load_external_ranked(args)
    assert len(candidates) == 2
    assert len(ranked) == 1
    assert ranked[0][1][0].target_epitope_score == max(e.epitope_score for e in candidates)
    assert frame.attrs['epitope_dataset'].result.df.n_rna_alt.tolist() == [3, 3]
    assert len(frame.attrs['epitope_dataset'].selection['representatives']) == 1
    assert assemble_peptide_constructs(ranked, options=PeptideConstructConfig(
        mode='minimal_epitope', min_antigen_length_aa=8))
    assert assemble_mrna_constructs(ranked, options=RNAConstructConfig(
        antigen_content='minimal_epitope', signal_peptide=None, include_mitd=False,
        optimize_linkers=False))
    native = tmp_path / 'native.tsv'
    frame.attrs['epitope_dataset'].save(native)
    again = load_external_ranked(SimpleNamespace(input_epitopes=str(native)))
    assert len(again[0]) == 1
    assert again[0][0][1][0].target_epitope_score == ranked[0][1][0].target_epitope_score


def test_original_allele_spelling_is_preserved_with_canonical_consumer_identity():
    dataset = EpitopeDataset.from_topiary(evidence(allele=['A0201', ALLELE]), sample_name='patient')
    assert dataset.result.df.allele.tolist() == ['A0201', ALLELE]
    assert dataset.scoring_frame().allele.tolist() == [ALLELE, ALLELE]
    assert all(e.patient_alleles == (ALLELE,) for e in dataset.epitopes)


def test_rna_only_rows_and_string_filters_survive_reload(tmp_path):
    combined = combine_sources({
        'predictions': evidence(status=['keep', 'drop']),
        'rna': pd.DataFrame({'protein_sequence': ['MAAASIINFEKL'], 'n_rna_alt': [0]})
    }, sample_name='patient')
    source = tmp_path / 'combined.tsv'
    combined.to_tsv(source)
    cfg = EpitopeConfig(filter_expr="status == 'keep'", score_expr='1 / affinity.value')
    _, frame, candidates = prepare_reports([('topiary', source)], cfg)
    assert [bool(e.per_allele_scores) for e in candidates] == [True, False]
    native = tmp_path / 'native.tsv'
    frame.attrs['epitope_dataset'].save(native)
    _, reloaded_frame, reloaded = prepare_reports([('epitopes', native)])
    assert [e.per_allele_scores for e in reloaded] == [e.per_allele_scores for e in candidates]
    assert len(reloaded_frame.attrs['epitope_dataset'].result.df) == 3


def test_topiary_additive_rescoring_retains_original_values_and_changes_only_chosen_policy(tmp_path):
    from topiary import rescore_candidates

    class Model:
        alleles = [ALLELE]
        uses_flanking_sequences = False
        calls = 0

        def kind_support(self):
            return {'pMHC_affinity': {'mhc_dependence': 'single_allele'}}

        def predict_dataframe(self, peptides, **kwargs):
            self.calls += 1
            return pd.DataFrame([
                dict(peptide=p, allele=ALLELE, kind='pMHC_affinity',
                     value=800. if p == 'SIINFEKL' else 10., score=.5,
                     prediction_method_name='testmodel', predictor_version='2') for p in peptides])

    model = Model()
    original = combine_sources({'original': evidence()}, sample_name='patient')
    result = rescore_candidates(original, model, prefix='new', use_flanks=False)
    calls = model.calls
    path = tmp_path / 'rescored.tsv'
    result.to_tsv(path)
    _, frame, candidates = prepare_reports([('topiary', path)], EpitopeConfig(score_expr='1 / affinity.value'))
    native = tmp_path / 'native.tsv'
    frame.attrs['epitope_dataset'].save(native)
    _, _, fresh = prepare_reports([('epitopes', native)], EpitopeConfig(
        score_expr='1 / new__testmodel__pMHC_affinity__value'))
    assert max(candidates, key=lambda e: e.epitope_score).sequence == 'SIINFEKL'
    assert max(fresh, key=lambda e: e.epitope_score).sequence == 'GILGFVFTL'
    assert frame.attrs['epitope_dataset'].result.df.value.tolist() == [50., 500.]
    assert model.calls == calls
