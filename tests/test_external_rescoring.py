"""Multi-source prediction evidence must preserve source and model identity."""
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from mhctools.base_predictor import BasePredictor
from mhctools.pred import Prediction, PeptideResult

from vaxrank.candidate_epitope import CandidateEpitope, Peptide
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_io import load_predictions, read_lens_report, write_neoepitope_report
from vaxrank.external_input import ExternalConstructOptions, load_external_ranked
from vaxrank.external_rescoring import (
    InputTablePredictor, external_inputs, prepare_reports, rescore_candidates,
)

DATA = Path(__file__).parent / 'data/epitope_fixtures'
INPUTS = [('lens', str(DATA / 'lens_example.tsv')),
          ('pvacseq', str(DATA / 'pvacseq_example.tsv'))]
ALLELES = ['HLA-A*02:01', 'HLA-B*07:02']


class ContextPredictor(BasePredictor):
    uses_flanking_sequences = True
    flank_length = 2
    predictor_version = '1.0'

    def __init__(self):
        super().__init__(alleles=ALLELES, default_peptide_lengths=[9])
        self.calls = []

    def predict_with_flanks(self, peptides, n_flanks, c_flanks):
        self.calls.append(list(zip(peptides, n_flanks, c_flanks)))
        return [PeptideResult(tuple(
            Prediction(kind='pMHC_affinity', peptide=p, allele=a,
                       value=20 + sum(map(ord, n + c)) / 100,
                       score=0.8, percentile_rank=0.5,
                       predictor_name='unified', predictor_version='1.0',
                       n_flank=n, c_flank=c)
            for a in self.alleles))
            for p, n, c in zip(peptides, n_flanks, c_flanks)]


def candidate(name, context='AA', comparator=True):
    sequence = 'SIINFEKLA'
    return CandidateEpitope(
        sequence=sequence, source_sequence=context + sequence + 'KK', offset=2,
        prediction_id=name, patient_alleles=tuple(ALLELES),
        overlaps_mutation=True,
        comparators={'wt': Peptide(sequence='SIINFEKLV')} if comparator else {})


def config():
    return EpitopeConfig(filter_expr='affinity[unified].value < 500',
                         score_expr='affinity[unified].score')


def test_fresh_contexts_are_batched_without_merging_occurrences():
    model = ContextPredictor()
    originals = [candidate('a', 'AA'), candidate('b', 'CC'), candidate('c', 'AA')]
    fresh = rescore_candidates(originals, [model], ALLELES)
    assert len(model.calls) == 2
    requests = [q for batch in model.calls for q in batch]
    assert requests.count(('SIINFEKLA', 'AA', 'KK')) == 1
    assert requests.count(('SIINFEKLA', 'CC', 'KK')) == 1
    assert requests.count(('SIINFEKLV', '', '')) == 1
    assert fresh[0].predictions_flat()[0].value == fresh[2].predictions_flat()[0].value
    assert fresh[0].predictions_flat()[0].value != fresh[1].predictions_flat()[0].value
    assert [e.prediction_id for e in fresh] == ['a', 'b', 'c']
    assert all(e.wt.predictions_flat()[0].predictor_name == 'unified' for e in fresh)
    assert all(not e.predictions_flat() for e in originals)


def test_anonymous_wildtype_cannot_keep_old_predictions_after_rescoring():
    old = Prediction(kind='pMHC_affinity', peptide='', allele=ALLELES[0],
                     predictor_name='old', score=0.9, value=10)
    e = replace(candidate('a'), comparators={'wt': Peptide('', predictions=[old])})
    fresh, = rescore_candidates([e], [ContextPredictor()], ALLELES)
    assert fresh.wt.sequence == ''
    assert not fresh.wt.predictions_flat()
    assert e.wt.predictions_flat() == (old,)


def test_input_table_predictor_replays_unknown_versions_without_models():
    reports, _, _ = prepare_reports([INPUTS[1]])
    originals = reports[0].epitopes
    cache = InputTablePredictor(originals)
    replay = cache.predict_candidates(originals)
    assert [e.predictions for e in replay] == [e.predictions for e in originals]
    assert {p.predictor_version for e in replay for p in e.predictions_flat()} == {''}
    with pytest.raises(KeyError):
        cache.predict_candidates([replace(originals[0], prediction_id='absent')])
    with pytest.raises(ValueError, match='context'):
        cache.predict_candidates([replace(originals[0], n_flank='wrong')])


def test_fresh_reports_preserve_original_values_and_allow_new_alleles(tmp_path):
    saved = tmp_path / 'input.tsv'
    reports, frame, epitopes = prepare_reports(
        INPUTS, config(), mode='fresh', models=[ContextPredictor()],
        alleles=ALLELES, input_predictions_path=saved)
    assert {r.source_format for r in reports} == {'lens', 'pvacseq'}
    assert all(e.per_allele_scores == dict.fromkeys(ALLELES, 0.8) for e in epitopes)
    assert set(frame['Prediction evidence']) == {'fresh'}
    assert set(frame['Allele']) == set(ALLELES)
    assert 'Input Predicted mutant pMHC affinity' in frame
    assert 'unified[1.0] pMHC_affinity value' in frame
    # Original table alleles and unknown predictor versions remain historical.
    recovered = load_predictions(saved)
    assert len(recovered) == len(epitopes)
    assert any(p.predictor_name == 'pvacseq' and not p.predictor_version
               for e in recovered for p in e.predictions_flat())
    assert all(p.predictor_name == 'unified' for e in epitopes for p in e.predictions_flat())
    assert frame['Predicted mutant pMHC affinity'].eq('').all()
    assert all(e.prediction_id.startswith('{') for e in recovered)
    assert not frame['Prediction identity'].isna().any()
    assert frame.attrs['topiary_df']['prediction_id'].nunique() == len(epitopes)
    assert frame['Input source'].nunique() == 2


def test_combined_report_emits_global_and_source_ranks(tmp_path):
    _, frame, epitopes = prepare_reports(INPUTS, config(), mode='fresh',
                                         models=[ContextPredictor()], alleles=ALLELES)
    path = tmp_path / 'report.csv'
    write_neoepitope_report(frame, epitopes, csv_report_path=path, epitope_config=config())
    observed = pd.read_csv(path)
    assert list(observed['rank']) == list(range(1, len(observed) + 1))
    for _, source in observed.groupby('Input source', sort=False):
        assert list(source['source_rank']) == list(range(1, len(source) + 1))
    assert observed['vaxrank_score'].eq(0.8).all()
    assert list(observed['Prediction identity']) == list(frame['Prediction identity'])


def test_identity_namespaces_isolate_two_reports_of_the_same_candidate(tmp_path):
    second = tmp_path / 'second.tsv'
    df = pd.read_csv(INPUTS[0][1], sep='\t')
    for metric in ('netmhcpan_4.1b.aff_nm', 'mhcflurry_2.1.1.aff'):
        df[metric] = 9999
    df.to_csv(second, sep='\t', index=False)
    reports, _, _ = prepare_reports([INPUTS[0], ('lens', str(second))])
    assert {e.prediction_id for e in reports[0].epitopes}.isdisjoint(
        e.prediction_id for e in reports[1].epitopes)
    assert any(e.epitope_score > 0 for e in reports[0].epitopes)
    assert all(e.epitope_score == 0 for e in reports[1].epitopes)


def test_duplicate_input_is_rejected_even_at_a_different_path(tmp_path):
    copy = tmp_path / 'copy.tsv'
    copy.write_bytes(Path(INPUTS[0][1]).read_bytes())
    with pytest.raises(ValueError, match='more than once'):
        prepare_reports([INPUTS[0], ('lens', str(copy))])


def test_fresh_scoring_does_not_validate_live_model_against_old_table():
    with pytest.raises(ValueError):
        read_lens_report(INPUTS[0][1], epitope_config=config())
    reports, _, eps = prepare_reports([INPUTS[0]], config(), mode='fresh',
                                      models=[ContextPredictor()], alleles=ALLELES)
    assert eps and reports[0].epitopes


def test_predictions_missing_a_requested_peptide_fail():
    class Incomplete(ContextPredictor):
        def predict_with_flanks(self, peptides, n_flanks, c_flanks):
            return []
    with pytest.raises(ValueError, match='omitted requested peptides'):
        rescore_candidates([candidate('a')], [Incomplete()], ALLELES)


def test_cli_parses_repeated_sources_and_conditional_model_configuration():
    from vaxrank.cli.arg_parser import parse_vaxrank_args
    args = parse_vaxrank_args(['--external-input', 'lens=one.tsv', '--external-input',
                              'pvacseq=two.tsv', '--external-predictions', 'fresh',
                              '--mhc-predictor', 'random', '--mhc-alleles', ALLELES[0]])
    assert external_inputs(args) == [('lens', 'one.tsv'), ('pvacseq', 'two.tsv')]
    assert args.external_predictions == 'fresh'
    assert parse_vaxrank_args(['--external-input', 'lens=one.tsv']).mhc_predictor is None
    with pytest.raises(ValueError, match='requires lens=PATH'):
        external_inputs(SimpleNamespace(external_input=['exacto=unsupported.json']))


@pytest.mark.parametrize('fmt,path', INPUTS)
def test_single_input_alias_has_the_same_scores_provenance_and_constructs(fmt, path):
    from vaxrank.cli.arg_parser import parse_vaxrank_args

    alias = load_external_ranked(parse_vaxrank_args(['--input-' + fmt, path]))
    repeated = load_external_ranked(parse_vaxrank_args(['--external-input', fmt + '=' + path]))
    assert alias[0] and repeated[0]
    assert [(str(source), [(v.amino_acids, v.combined_score) for v in peptides])
            for source, peptides in alias[0]] == [
        (str(source), [(v.amino_acids, v.combined_score) for v in peptides])
        for source, peptides in repeated[0]]
    pd.testing.assert_frame_equal(alias[1], repeated[1])
    assert [e.to_dict() for e in alias[2]] == [e.to_dict() for e in repeated[2]]
    assert alias[3].inputs == repeated[3].inputs
    assert alias[3].mhc_alleles == repeated[3].mhc_alleles


@pytest.mark.parametrize('source', [
    ['--input-lens', 'missing.tsv'],
    ['--input-pvacseq', 'missing.tsv'],
    ['--external-input', 'lens=missing.tsv'],
])
@pytest.mark.parametrize('extra', [
    ['--mhc-predictor', 'random'],
    ['--mhc-alleles', 'HLA-A*01:01'],
    ['--mhc-alleles-file', 'missing-alleles.txt'],
])
@pytest.mark.parametrize("export_originals", [False, True])
def test_input_mode_rejects_ignored_options_before_reading_or_exporting(
        source, extra, export_originals, tmp_path):
    from vaxrank.cli.arg_parser import parse_vaxrank_args

    output = tmp_path / 'originals.tsv'
    export_args = ['--output-input-predictions', str(output)] if export_originals else []
    args = parse_vaxrank_args(source + extra + export_args)
    with pytest.raises(ValueError, match='--external-predictions fresh'):
        load_external_ranked(args)
    assert not output.exists()


def test_unified_candidates_reach_both_vaccine_designs(monkeypatch):
    import mhctools.cli
    from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
    from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs
    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', lambda args: [ContextPredictor()])
    monkeypatch.setattr(mhctools.cli, 'mhc_alleles_from_args', lambda args: ALLELES)
    args = SimpleNamespace(
        external_input=[fmt + '=' + path for fmt, path in INPUTS],
        external_predictions='fresh', mhc_predictor=['test'], genome=None)
    ranked, _, _, patient, _ = load_external_ranked(args, epitope_config=config())
    assert ranked
    assert patient.mhc_alleles == ALLELES
    assert len(patient.inputs) == 2
    assert assemble_peptide_constructs(ranked, options=PeptideConstructConfig(mode='slp'))
    assert assemble_mrna_constructs(ranked, options=RNAConstructConfig(
        signal_peptide=None, include_mitd=False, poly_a_length=10,
        antigens_per_construct=2, max_constructs=1, optimize_linkers=False))


def test_mixed_cli_reranks_with_fresh_values_and_exports_originals(monkeypatch, tmp_path):
    import mhctools.cli
    from vaxrank.cli.entry_point import main

    inputs = [arg for fmt, path in INPUTS for arg in ('--external-input', fmt + '=' + path)]
    original_csv = tmp_path / 'original.csv'
    original_native = tmp_path / 'original.tsv'
    main(inputs + ['--output-csv', str(original_csv),
                   '--output-input-predictions', str(original_native)])
    original_frame = pd.read_csv(original_csv).sort_values('rank')
    promoted = original_frame.iloc[-1]['Mutant peptide sequence']

    class RerankingPredictor(ContextPredictor):
        def predict_with_flanks(self, *args):
            return [PeptideResult(tuple(
                replace(p, score=0.9 if p.peptide == promoted else 0.1)
                for p in result.preds))
                for result in super().predict_with_flanks(*args)]

    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', lambda args: [RerankingPredictor()])
    fresh_csv = tmp_path / 'fresh.csv'
    fresh_native = tmp_path / 'fresh.tsv'
    original_copy = tmp_path / 'original-copy.tsv'
    main(inputs + [
        '--external-predictions', 'fresh', '--mhc-predictor', 'random',
        '--mhc-alleles', ','.join(ALLELES),
        '--config-text', 'epitopes.filter_expr=affinity[unified].value < 500',
        '--config-text', 'epitopes.score_expr=affinity[unified].score',
        '--output-csv', str(fresh_csv), '--output-epitopes', str(fresh_native),
        '--output-input-predictions', str(original_copy)])
    fresh_frame = pd.read_csv(fresh_csv).sort_values('rank')
    assert original_frame.iloc[0]['Mutant peptide sequence'] != promoted
    assert fresh_frame.iloc[0]['Mutant peptide sequence'] == promoted
    assert fresh_frame.iloc[0]['vaxrank_score'] == 0.9
    assert original_copy.read_bytes() == original_native.read_bytes()
    fresh = load_predictions(fresh_native)
    assert {p.predictor_name for e in fresh for p in e.predictions_flat()} == {'unified'}
    assert {e.prediction_id for e in fresh} == {
        e.prediction_id for e in load_predictions(original_native)}
    assert set(fresh_frame['Input format']) == {'lens', 'pvacseq'}


def test_same_variant_is_not_double_counted_for_construct_selection(tmp_path):
    second = tmp_path / 'second.tsv'
    df = pd.read_csv(INPUTS[0][1], sep='\t')
    df['additional_annotation'] = 'independent report; same candidate identities'
    df.to_csv(second, sep='\t', index=False)
    args = SimpleNamespace(external_input=['lens=' + INPUTS[0][1], 'lens=' + str(second)])
    ranked, frame, _, patient, _ = load_external_ranked(args)
    assert frame['Input source'].nunique() == 2
    original = read_lens_report(INPUTS[0][1])
    from vaxrank.external_input import lens_ranking_result
    baseline = lens_ranking_result(original, original.epitopes, options=ExternalConstructOptions())
    assert len(ranked) == len(baseline.ranked)
    assert patient.num_somatic_variants == baseline.input_summary.num_somatic_variants


def test_historical_scores_remain_source_scoped_through_report_export(tmp_path):
    reports, frame, eps = prepare_reports(INPUTS)
    for report, (fmt, path) in zip(reports, INPUTS):
        from vaxrank.external_rescoring import READERS
        baseline = READERS[fmt](path)
        assert [e.per_allele_scores for e in report.epitopes] == [
            e.per_allele_scores for e in baseline.epitopes]
        assert any(e.epitope_score > 0 for e in report.epitopes)
    path = tmp_path / 'historical.csv'
    write_neoepitope_report(frame, eps, csv_report_path=path)
    observed = pd.read_csv(path)
    assert observed.groupby('Input format')['vaxrank_score'].max().gt(0).all()


def test_fresh_scoring_namespaces_historical_comparator_metrics():
    _, frame, _ = prepare_reports(INPUTS, config(), mode='fresh',
                                  models=[ContextPredictor()], alleles=ALLELES)
    scoring = frame.attrs['topiary_df']
    assert 'input_wt_value' in scoring
    assert 'wt_value' not in scoring
    assert 'input_mhcflurry_affinity_value' in scoring
    assert 'mhcflurry_affinity_value' not in scoring


def test_wildtype_table_sequence_does_not_claim_mutant_source_context():
    from vaxrank.epitope_io import read_pvacseq_report
    for epitope in read_pvacseq_report(INPUTS[1][1]).epitopes:
        if epitope.wt:
            assert epitope.wt.source_sequence == ''
            assert epitope.wt.n_flank == epitope.wt.c_flank == ''


def test_missing_allele_coverage_is_not_fabricated():
    class Incomplete(ContextPredictor):
        def predict_with_flanks(self, *args):
            return [PeptideResult(r.preds[:1]) for r in super().predict_with_flanks(*args)]
    with pytest.raises(ValueError, match='allele coverage'):
        rescore_candidates([candidate('a')], [Incomplete()], ALLELES)


def test_haplotype_model_requires_explicit_transport_support():
    class Haplotype(ContextPredictor):
        def kind_support(self):
            return {'pMHC_presentation': {'mhc_dependence': 'haplotype'}}
    with pytest.raises(ValueError, match='haplotype'):
        rescore_candidates([candidate('a')], [Haplotype()], ALLELES)


def test_input_mode_does_not_construct_a_live_predictor(monkeypatch):
    import mhctools.cli
    def forbidden(*args, **kwargs):
        raise AssertionError('Historical table replay must not initialize a model')
    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', forbidden)
    ranked, _, _, _, _ = load_external_ranked(SimpleNamespace(
        external_input=[fmt + '=' + path for fmt, path in INPUTS]))
    assert ranked


def test_sequence_matches_link_observations_without_merging_context_or_abundance(tmp_path):
    second = tmp_path / 'second.tsv'
    alternative = pd.read_csv(INPUTS[0][1], sep='\t')
    alternative.loc[0, 'tpm'] = 99
    alternative.loc[1, 'pep_context'] = 'AA' + alternative.loc[1, 'peptide'] + 'CC'
    alternative.to_csv(second, sep='\t', index=False)
    reports, frame, _ = prepare_reports([INPUTS[0], ('lens', str(second)), INPUTS[1]])
    a, b = (r.report_df for r in reports[:2])
    assert a.iloc[0]['Context sequence SHA256'] == b.iloc[0]['Context sequence SHA256']
    assert a.iloc[0]['Prediction identity'] != b.iloc[0]['Prediction identity']
    assert a.iloc[1]['Peptide sequence SHA256'] == b.iloc[1]['Peptide sequence SHA256']
    assert a.iloc[1]['Context sequence SHA256'] != b.iloc[1]['Context sequence SHA256']
    assert {r.rows[0]['gene_tpm'] for r in reports[:2]} == {42.5, 99}
    assert set(frame.loc[frame['Input format'] == 'pvacseq', 'Context extent']) == {'epitope only'}
    assert set(a['Context extent']) == {'reported window'}


def test_run_summary_lists_combined_sources_and_configured_genotype(tmp_path):
    from vaxrank.cli.entry_point import write_run_summary
    args = SimpleNamespace(
        external_input=[fmt + '=' + path for fmt, path in INPUTS],
        external_predictions='fresh', output_dir=str(tmp_path),
        mhc_alleles=','.join(ALLELES), vaccine_type=['peptide', 'mrna'],
        _inferred_mhc_alleles_from_lens=ALLELES)
    write_run_summary(args, None, source='external')
    summary = (tmp_path / 'run_summary.txt').read_text()
    assert all(path in summary for _, path in INPUTS)
    assert 'Prediction evidence: fresh' in summary
    assert 'inferred from report' not in summary
    assert 'full pipeline' not in summary


def test_allele_free_prediction_keeps_its_scope():
    class Processing(ContextPredictor):
        def kind_support(self):
            return {'antigen_processing': {'mhc_dependence': 'none'}}

        def predict_with_flanks(self, peptides, n_flanks, c_flanks):
            return [PeptideResult((Prediction(
                kind='antigen_processing', peptide=p, score=0.7,
                predictor_name='processing', predictor_version='1'),))
                for p in peptides]

    fresh, = rescore_candidates([candidate('a')], [Processing()], ALLELES)
    assert len(fresh.predictions_flat()) == 1
    assert not fresh.predictions_flat()[0].allele
    assert fresh.patient_alleles == tuple(ALLELES)


def test_empty_candidate_input_is_explicit(tmp_path):
    path = tmp_path / 'empty.tsv'
    pd.read_csv(INPUTS[0][1], sep='\t').iloc[:0].to_csv(path, sep='\t', index=False)
    with pytest.raises(ValueError, match='no candidate peptides'):
        prepare_reports([('lens', str(path))])
