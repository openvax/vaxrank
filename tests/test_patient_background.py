"""Patient comparisons use exact nucleotide contexts and retain uncertainty."""

from copy import deepcopy
import gzip
import json
from pathlib import Path

import pytest

from vaxrank.patient_background import PatientBackgroundPairs
from . import test_antigen_input
from . import test_isovar_support_contract

dependencies = test_antigen_input.dependencies
reference = test_isovar_support_contract.reference


@pytest.fixture
def cases():
    path = Path(__file__).parent / 'data/isovar/patient_background.json.gz'
    return json.loads(gzip.decompress(path.read_bytes()))


def matched(case):
    source = PatientBackgroundPairs(case['pairs'])
    event = case['hypotheses']['events'][0]
    protein = event['protein_hypotheses'][0]
    return source.match_hypothesis(case['hypotheses'], event['event_id'], protein['hypothesis_id'])


def test_focal_baseline_keeps_patient_background_instead_of_annotation(cases):
    case = cases['focal']
    result = matched(case).compare(1, 4)
    assert result['state'] == 'available'
    assert result['sequence'] == 'AKD'  # Annotation would give AKG.
    assert result['comparison'] == 'amino_acid_divergent'
    pair, = case['pairs']['events'][0]['pairs']
    assert {e['origin'] for e in pair['retained_edits']} == {
        'germline_variant', 'co_somatic_variant', 'unexplained'}
    assert result['alternatives'][0]['baseline_nucleotide_interval'] == [5, 14]


def test_explicit_haplotype_policy_preserves_actual_reversions(cases):
    case = cases['haplotype']
    result = matched(case).compare(0, 4)
    assert result['baseline_policy'] == 'supplied_somatic_haplotype'
    pair, = case['pairs']['events'][0]['pairs']
    assert {e['origin'] for e in pair['removed_edits']} == {'nominated_variant', 'co_somatic_variant'}
    assert all(e['origin'] != 'germline_variant' for e in pair['removed_edits'])


def test_germline_indel_and_identical_peptide_remain_evidence(cases):
    case = cases['germline_deletion']
    result = matched(case).compare(0, 4)
    assert result['state'] == 'available'
    assert result['comparison'] == 'identical'
    assert result['sequence'] == 'MARN'
    retained = case['pairs']['events'][0]['pairs'][0]['retained_edits']
    assert any(e['origin'] == 'germline_variant' and e['alt_bases'] == '' for e in retained)


def test_inframe_reversion_has_explicit_unequal_peptide_lengths(cases):
    result = matched(cases['inframe_deletion']).compare(0, 3)
    assert result['state'] == 'available'
    assert result['sequence'] == 'MAKD'
    assert result['same_length'] is False
    assert result['alternatives'][0]['baseline_protein_interval'] == [0, 4]


def test_frameshift_does_not_take_same_offset_reference_slice(cases):
    result = matched(cases['frameshift']).compare(0, 4)
    assert result['state'] == 'unavailable'
    assert result['sequence'] is None
    assert result['comparison'] == 'unavailable'
    assert result['alternatives'][0]['reason'] == 'peptide_boundaries_do_not_share_baseline_frame'


def test_original_sid_ntf3_exact_transcripts_reconcile_without_losing_germline(cases):
    case = cases['sid_ntf3']
    result = matched(case).compare(0, 20)
    assert result['sequence'] == 'NKLSKQMVDVNENYQSTLPK'
    assert {a['transcript_id'] for a in result['alternatives']} == {'ENST00000331010', 'ENST00000423158'}
    assert len({a['pair_id'] for a in result['alternatives']}) == 2
    pairs = case['pairs']['events'][0]['pairs']
    assert all(any(e['origin'] == 'germline_variant' for e in p['retained_edits']) for p in pairs)


def test_tiled_phase_is_retained_without_claiming_observation(cases):
    source = PatientBackgroundPairs(cases['tiled']['pairs'])
    saved = source.source_evidence[source.evidence_id]['result']
    pair, = saved['events'][0]['pairs']
    assert pair['phase']['state'] == 'assembly_only'
    assert pair['baseline']['observation'] == 'computed_only'
    assert pair['baseline']['normal_tissue_observed'] is False
    assert pair['baseline']['translation_observed'] is False
    assert saved == cases['tiled']['pairs']
    assert matched(cases['tiled']).compare(0, 4)['state'] == 'available'


def test_missing_compatible_pair_is_not_dropped_from_reconciliation(cases):
    case = cases['two_transcripts']
    case['pairs']['events'][0]['pairs'].pop()
    result = matched(case).compare(0, 4)
    assert result['state'] == 'ambiguous'
    assert result['sequence'] is None
    assert len(result['alternatives']) == 2
    assert {a['state'] for a in result['alternatives']} == {'available', 'unavailable'}
    assert 'exact_pair_not_recorded' in {a['reason'] for a in result['alternatives']}


def test_missing_event_is_explicitly_unavailable(cases):
    case = cases['focal']
    case['pairs']['events'] = []
    result = matched(case).compare(0, 4)
    assert result['state'] == 'unavailable'
    assert result['alternatives'][0]['reason'] == 'exact_pair_not_recorded'


def test_same_amino_acids_do_not_override_nucleotide_identity(cases):
    case = cases['focal']
    case['pairs']['events'][0]['pairs'][0]['translation_id'] = 'a-different-nucleotide-context'
    with pytest.raises(ValueError, match='translation/transcript provenance'):
        matched(case)


def test_same_pair_identity_cannot_override_conflicting_event(cases):
    case = cases['focal']
    case['pairs']['events'][0]['variant']['start'] += 1
    with pytest.raises(ValueError, match='event provenance'):
        matched(case)


@pytest.mark.parametrize('start,end', [(-1, 2), (0, 5), (2, 2), (True, 2), (0, 2.0)])
def test_invalid_peptide_bounds_are_rejected(cases, start, end):
    with pytest.raises(ValueError, match='interval'):
        matched(cases['focal']).compare(start, end)


def test_join_does_not_mutate_original_evidence(cases):
    original = deepcopy(cases['two_transcripts'])
    match = matched(cases['two_transcripts'])
    first = match.compare(0, 4)
    second = match.compare(1, 3)
    assert len(first['alternatives']) == len(second['alternatives']) == 2
    assert cases['two_transcripts'] == original


@pytest.mark.parametrize('field', ['sample_id', 'source'])
def test_join_rejects_other_sample_or_read_set(cases, field):
    case = cases['focal']
    case['pairs'][field] = 'different'
    case['pairs']['evidence_scope'] = [case['pairs']['sample_id'], case['pairs']['source']]
    for support in case['pairs']['evidence_sets'].values():
        support['evidence_scope'] = list(case['pairs']['evidence_scope'])
    with pytest.raises(ValueError, match='sample/read-set provenance'):
        matched(case)


def test_identical_repeated_pairs_do_not_create_extra_contexts(cases):
    case = cases['sid_ntf3']
    case['pairs']['events'][0]['pairs'] *= 3
    source = PatientBackgroundPairs(case['pairs'])
    assert source.source_evidence[source.evidence_id]['result'] == case['pairs']
    assert len(matched(case).compare(0, 20)['alternatives']) == 2


def test_unrequested_comparator_length_remains_available_but_unmeasured(cases):
    from topiary import TopiaryPredictor
    from vaxrank.candidate_epitope import CandidateEpitope
    from vaxrank.patient_background import attach_patient_background
    from .test_external_rescoring import ContextPredictor
    context = matched(cases['inframe_deletion'])
    candidate = CandidateEpitope(sequence=context.sequence[:3], source_sequence=context.sequence,
                                offset=0, overlaps_mutation=True)
    model = ContextPredictor()
    result, = attach_patient_background([candidate], TopiaryPredictor(models=[model]), context)
    assert result.reference_comparator.sequence == 'MAKD'
    assert not result.reference_comparator.predictions_flat()
    assert result.comparison_evidence['prediction_state'] == 'outside_requested_peptide_lengths'
    assert result.comparison_evidence['same_length'] is False
    assert not model.calls
    from vaxrank.report_display import epitope_display
    from .test_report_clarity import ranked_pair, prediction
    antigen = ranked_pair()[0][0]
    assert 'length outside the requested prediction lengths' in epitope_display(
        result, antigen, {}, prediction())['comparison_note']


def sid_candidate(cases):
    from vaxrank.candidate_epitope import CandidateEpitope
    context = matched(cases['sid_ntf3'])
    return context, CandidateEpitope(sequence=context.sequence[4:13], source_sequence=context.sequence,
        offset=4, overlaps_mutation=True, overlaps_targetable=True)


def test_comparator_prediction_uses_topiary_and_preserves_separate_reference_role(cases):
    from dataclasses import replace
    from topiary import TopiaryPredictor
    from vaxrank.candidate_epitope import Peptide
    from vaxrank.patient_background import attach_patient_background
    from .test_external_rescoring import ContextPredictor

    context, candidate = sid_candidate(cases)
    candidate = replace(candidate, comparators={'wt': Peptide('AAAAAAAAA')})
    model = ContextPredictor()
    result, = attach_patient_background([candidate], TopiaryPredictor(models=[model]), context)
    assert result.wt.sequence == 'AAAAAAAAA'
    assert result.reference_comparator.sequence == 'KQMVDVNEN'
    assert result.reference_comparator.predictions_flat()
    assert {p.peptide for p in result.reference_comparator.predictions_flat()} == {'KQMVDVNEN'}
    assert result.comparison_evidence['prediction_state'] == 'measurements_recorded'
    assert len(result.comparison_evidence['alternatives']) == 2
    assert candidate.comparison_evidence is None


def test_missing_patient_context_never_uses_annotation_fallback(cases):
    from dataclasses import replace
    from topiary import TopiaryPredictor
    from vaxrank.candidate_epitope import Peptide
    from vaxrank.patient_background import attach_patient_background
    from .test_external_rescoring import ContextPredictor

    _, candidate = sid_candidate(cases)
    model = ContextPredictor()
    result, = attach_patient_background([replace(candidate, comparators={'wt': Peptide('AAAAAAAAA')})],
                                        TopiaryPredictor(models=[model]), None)
    assert result.wt is not None and result.reference_comparator is None
    assert result.comparison_evidence['state'] == 'unavailable'
    assert not model.calls


def test_comparison_records_survive_native_save_load_and_slicing(cases, tmp_path):
    from topiary import TopiaryPredictor
    from vaxrank.epitope_dataset import EpitopeDataset
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.patient_background import attach_patient_background
    from .test_external_rescoring import ContextPredictor

    context, candidate = sid_candidate(cases)
    result, = attach_patient_background([candidate], TopiaryPredictor(models=[ContextPredictor()]), context)
    source = PatientBackgroundPairs(cases['sid_ntf3']['pairs'])
    dataset = EpitopeDataset.from_predictions([result],
        config=EpitopeConfig(comparison_reference='patient_background'),
        selection={'source_evidence': source.source_evidence})
    path = tmp_path / 'native.tsv'
    dataset.save(path)
    restored = EpitopeDataset.load(path)
    assert restored.epitopes == (result,)
    assert restored.selection['source_evidence'] == source.source_evidence
    sliced = restored.epitopes[0].sliced(2, 18)
    assert sliced.offset == 2
    assert sliced.comparison_evidence['mutant_protein_interval'] == [4, 13]
    assert sliced.comparators == result.comparators


def test_report_comparisons_remain_specific_to_each_source_occurrence(cases):
    from dataclasses import replace
    from topiary import TopiaryPredictor
    from vaxrank import EpitopeConfig
    from vaxrank.epitope_dataset import EpitopeDataset
    from vaxrank.patient_background import attach_patient_background
    from .test_external_rescoring import ContextPredictor, ALLELES
    context, first = sid_candidate(cases)
    second = replace(first, offset=5, sequence=context.sequence[5:14])
    compared = attach_patient_background([first, second], TopiaryPredictor(models=[ContextPredictor()]), context)
    candidates = [replace(e, prediction_id='shared-source', patient_alleles=tuple(ALLELES),
        predictions=[replace(p, peptide=e.sequence) for p in e.reference_comparator.predictions_flat()])
        for e in compared]
    dataset = EpitopeDataset.from_predictions(candidates,
        config=EpitopeConfig(comparison_reference='patient_background'))
    frame = dataset.report_frame()
    for row in frame.to_dict('records'):
        assert row['Comparison sequence'] == {4: 'KQMVDVNEN', 5: 'QMVDVNENY'}[row['Peptide offset']]
        assert row['Predicted comparison pMHC affinity']


def test_fresh_predictions_update_comparator_state_without_inventing_flanks(cases):
    from dataclasses import replace
    from topiary import TopiaryPredictor
    from vaxrank.external_rescoring import rescore_candidates
    from vaxrank.patient_background import attach_patient_background
    from .test_external_rescoring import ContextPredictor, ALLELES
    context, candidate = sid_candidate(cases)
    initial = ContextPredictor()
    initial.default_peptide_lengths = [10]  # Comparator is known but was not requested.
    compared, = attach_patient_background([candidate], TopiaryPredictor(models=[initial]), context)
    assert compared.comparison_evidence['prediction_state'] == 'outside_requested_peptide_lengths'
    model = ContextPredictor()
    fresh, = rescore_candidates([replace(compared, prediction_id='source')], [model], ALLELES)
    assert fresh.comparison_evidence['prediction_state'] == 'measurements_recorded'
    assert compared.comparison_evidence['prediction_state'] == 'outside_requested_peptide_lengths'
    assert all((not p.n_flank and not p.c_flank) for p in fresh.reference_comparator.predictions_flat())
    assert fresh.comparison_evidence['prediction_context'] == 'peptide_only'


@pytest.mark.parametrize('fault', ['source', 'interval', 'comparator', 'policy'])
def test_native_replay_rejects_contradictory_comparison_evidence(cases, tmp_path, fault):
    from dataclasses import replace
    from topiary import TopiaryPredictor
    from vaxrank import EpitopeConfig
    from vaxrank.epitope_dataset import EpitopeDataset
    from vaxrank.patient_background import attach_patient_background
    from .test_external_rescoring import ContextPredictor
    context, candidate = sid_candidate(cases)
    compared, = attach_patient_background([candidate], TopiaryPredictor(models=[ContextPredictor()]), context)
    source = PatientBackgroundPairs(cases['sid_ntf3']['pairs'])
    dataset = EpitopeDataset.from_predictions([compared],
        config=EpitopeConfig(comparison_reference='patient_background'),
        selection={'source_evidence': source.source_evidence})
    if fault == 'source':
        dataset.selection['source_evidence'][source.evidence_id]['result']['isovar_version'] = 'changed'
    elif fault == 'interval':
        compared.comparison_evidence['mutant_protein_interval'] = [5, 14]
    elif fault == 'comparator':
        dataset.epitopes = (replace(compared, comparators={}),)
    else:
        dataset.config = EpitopeConfig()
    path = tmp_path / 'contradictory.tsv'
    dataset.save(path)
    with pytest.raises(ValueError, match='patient-background'):
        EpitopeDataset.load(path)


def test_direct_pair_file_supplies_native_scope_and_conflicts_are_rejected(cases, tmp_path):
    from vaxrank import EpitopeConfig
    from vaxrank.cli.arg_parser import parse_vaxrank_args
    from vaxrank.construct_replay import pipeline_native_dataset
    from vaxrank.patient_background import direct_patient_background
    pair_file = tmp_path / 'pairs.json'
    pair_file.write_text(json.dumps(cases['sid_ntf3']['pairs']))
    vcf, bam = tmp_path / 'input.vcf', tmp_path / 'input.bam'
    vcf.write_text('synthetic input for scope test')
    bam.write_bytes(b'synthetic input for scope test')
    argv = ['--vcf', str(vcf), '--bam', str(bam), '--mhc-predictor', 'random',
            '--mhc-alleles', 'HLA-A*02:01', '--comparison-reference', 'patient_background',
            '--patient-background-pairs', str(pair_file)]
    args = parse_vaxrank_args(argv)
    config = EpitopeConfig(comparison_reference='patient_background')
    background = direct_patient_background(args, config)
    dataset = pipeline_native_dataset(args, [], config, patient_background=background)
    assert dataset.provenance[0].scope.sample_id == 'patient'
    assert dataset.provenance[0].scope.library_id == 'RNA'
    with pytest.raises(ValueError, match='matching patient-background'):
        direct_patient_background(parse_vaxrank_args(argv + ['--rna-source', 'different']), config)


@pytest.mark.parametrize('extension', ['html', 'txt'])
def test_report_names_patient_baseline_and_retains_both_contexts(cases, extension):
    from dataclasses import replace
    from topiary import TopiaryPredictor
    from vaxrank.patient_background import attach_patient_background
    from vaxrank.vaccine_antigen import VaccineAntigen, TargetableMask, AminoAcidInterval
    from vaxrank.vaccine_peptide import VaccinePeptide
    from .test_external_rescoring import ContextPredictor
    from .test_isovar_inputs import attestation
    from .test_product_reports import template, render

    context, candidate = sid_candidate(cases)
    compared, = attach_patient_background([candidate], TopiaryPredictor(models=[ContextPredictor()]), context)
    # Use measured comparator leaves as deterministic test-only target predictions.
    leaves = [replace(p, peptide=compared.sequence) for p in compared.reference_comparator.predictions_flat()]
    compared = replace(compared, predictions=leaves, per_allele_scores={'HLA-A*02:01': 1.})
    antigen = VaccineAntigen(kind='mutation', amino_acids=context.sequence,
        targetable_mask=TargetableMask((AminoAcidInterval(4, 13),)),
        tumor_specificity=attestation(), gene_name='NTF3', source_identifier=context.hypothesis_id)
    peptide = VaccinePeptide(antigen=antigen, epitopes=[compared], combined_score_expr='target_epitope_score')
    source = PatientBackgroundPairs(cases['sid_ntf3']['pairs'])
    report = render(template([(antigen, [peptide])], source_evidence=source.source_evidence,
                             candidate_details=True), extension)
    assert 'Patient background (focal)' in report
    assert 'Patient-background comparison 1' in report
    assert 'ENST00000331010' in report and 'ENST00000423158' in report
    assert 'germline variant' in report
    assert 'protein expression' in report


def test_cli_and_yaml_resolve_the_same_comparator_configuration():
    from vaxrank.cli.arg_parser import make_vaxrank_arg_parser
    from vaxrank.cli.epitope_config_args import epitope_config_from_args
    from vaxrank.config.loader import extract_epitope_config_kwargs
    from vaxrank.epitope_config import EpitopeConfig

    parser = make_vaxrank_arg_parser()
    args = parser.parse_args(['--vcf', 'input.vcf', '--bam', 'input.bam',
        '--mhc-predictor', 'random',
        '--comparison-reference', 'patient_background',
        '--patient-background-baseline', 'supplied_somatic_haplotype'])
    expected = EpitopeConfig(**extract_epitope_config_kwargs({'epitopes': {
        'comparison_reference': 'patient_background',
        'patient_background_baseline': 'supplied_somatic_haplotype'}}))
    assert epitope_config_from_args(args, {}) == expected


def test_generic_isovar_cli_retains_pairs_and_replays_without_inputs(cases, tmp_path, dependencies, monkeypatch):
    from vaxrank import antigens_from_isovar
    from vaxrank.input_scope import InputScope
    from vaxrank.cli.entry_point import run_cli
    from vaxrank.epitope_dataset import EpitopeDataset
    from .test_isovar_inputs import attestation
    from .test_epitope_dataset import forbidden
    from .test_external_rescoring import ALLELES

    case = cases['sid_ntf3']
    source = antigens_from_isovar(case['hypotheses'], species='homo_sapiens',
        scope=InputScope(patient_id='test', sample_id='patient', reference_assembly='GRCh38',
                         annotation='ensembl:93', mhc_alleles=tuple(ALLELES)),
        tumor_specificity=attestation(), options={'admission': {'allow_partial': True}})
    path, pairs = tmp_path / 'antigens.json', tmp_path / 'pairs.json'
    source.save(path)
    pairs.write_text(json.dumps(case['pairs']))
    native, report = tmp_path / 'native.tsv', tmp_path / 'report.txt'
    run_cli(['--input-antigens', str(path), '--mhc-predictor', 'random', '--ensembl-release', '93',
        '--comparison-reference', 'patient_background', '--patient-background-pairs', str(pairs),
        '--output-epitopes', str(native), '--output-ascii-report', str(report),
        '--report-candidate-details', '--no-processing-aware-annotation',
        '--vaccine-peptide-length', '20', '--config-text', 'epitopes.score_expr=1.0',
        '--config-text', 'vaccine_peptides.combined_score_expr=target_epitope_score'])
    saved = EpitopeDataset.load(native)
    compared = [e for e in saved.epitopes if e.comparison_evidence is not None]
    assert compared and any(e.reference_comparator is not None for e in compared)
    assert 'Patient background (focal)' in report.read_text()
    record = saved.selection['source_evidence'][compared[0].comparison_evidence['evidence_id']]
    assert record['result'] == case['pairs']
    path.unlink()
    pairs.unlink()
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    replay = tmp_path / 'replay.html'
    run_cli(['--input-epitopes', str(native), '--output-html-report', str(replay),
             '--report-candidate-details', '--no-processing-aware-annotation'])
    assert 'Patient background (focal)' in replay.read_text()
    assert 'ENST00000331010' in replay.read_text()
    assert EpitopeDataset.load(native).epitopes == saved.epitopes


def test_live_rna_export_binds_selected_isovar_object_and_core_predictions(tmp_path, reference, monkeypatch):
    from types import SimpleNamespace
    from vaxrank import EpitopeConfig, VaccineConfig, run_vaxrank
    from vaxrank.patient_background import direct_patient_background
    from .test_external_rescoring import ContextPredictor

    monkeypatch.setenv('VAXRANK_REF_PEPTIDES_DIR', str(tmp_path / 'vaxrank-cache'))
    sequence = test_isovar_support_contract.MUTANT
    reads = [(sequence, [(len(sequence), 'M')], 0)] * 2
    result = test_isovar_support_contract.reconstruct(tmp_path, reference, reads)
    config = EpitopeConfig(comparison_reference='patient_background', score_expr='1.0')
    args = SimpleNamespace(rna_sample_id='sample', rna_source='synthetic-bam')
    assert direct_patient_background(args, config) is None  # Valid preflight, no RNA work.
    pairs = direct_patient_background(args, config, [result])
    assert pairs.export['sample_id'] == 'sample'
    context = pairs.contexts_for_results([result])[id(result)]
    assert context.sequence == result.top_protein_sequence.amino_acids
    design = run_vaxrank([result], ContextPredictor(), patient_background=pairs,
        epitope_config=config, vaccine_config=VaccineConfig(
            preferred_peptide_length=10, min_peptide_length=10, max_peptide_length=10,
            padding_around_mutation=0, combined_score_expr='target_epitope_score'))
    assert design.patient_background_pairs == pairs.export
    assert design.ranked_vaccine_peptides
    epitopes = [e for _, peptides in design.ranked_vaccine_peptides for p in peptides for e in p.epitopes]
    assert epitopes and all(e.reference_comparator is not None for e in epitopes)
    assert {e.reference_comparator.sequence for e in epitopes} == {'MAAAAKGGG', 'AAAAKGGGG'}
    assert all(e.wt is None for e in epitopes)
