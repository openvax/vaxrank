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
