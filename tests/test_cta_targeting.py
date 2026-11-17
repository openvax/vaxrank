"""Current owner sequence roles affect targets and survive native replay."""

from dataclasses import replace
from types import SimpleNamespace

import pytest
from pyensembl import EnsemblRelease
from orfdata import ReferenceContext

from vaxrank.cta_targeting import admit_expression, apply_cta_targeting, peptide_containment, targeting_evidence
from vaxrank.cta_admission import CTAAdmissionPolicy
from vaxrank.cta_expression import CTAExpressionResult
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.input_scope import InputScope, InputProvenance
from vaxrank.vaccine_antigen import VaccineAntigen, TargetableMask, AminoAcidInterval
from vaxrank.vaccine_peptide import VaccinePeptide

from .test_cta_expression import contract, PRAME, PRAME_TRANSCRIPT
from .test_core_logic_config import _make_epitope
from .test_isovar_inputs import attestation


@pytest.fixture
def bundle():
    context = ReferenceContext.from_genome(EnsemblRelease(115))
    reference = dict(taxon_id=9606, assembly_accession=context.assembly_accession,
        annotation_source='Ensembl', annotation_release='115',
        identifier_namespaces={'gene': 'Ensembl gene', 'transcript': 'Ensembl transcript',
                               'protein': 'Ensembl protein'}, reference_id='fixture-reference')
    rows = []
    def add(name, sequence, role, qualified=False, segments=None):
        rows.append(dict(record_id=name, sequence=sequence, shared_sequence=role,
            direct_target_candidate=qualified, valid_segments=segments or [dict(start=0, end=len(sequence))],
            occurrences=[dict(gene_id=PRAME if qualified else 'ENSG00000111111',
                gene_name='PRAME' if qualified else name, transcript_id=PRAME_TRANSCRIPT + '.3',
                direct_target_candidate=qualified)]))
    add('representative', 'SIINFEKLGILGFVFTL', 'allow', True)
    add('unsafe', 'AASIINFEKLGKK', 'exclude')
    add('uncertain', 'ACDEFGHIK', 'review')
    add('interrupted', 'GILGX FVF'.replace(' ', ''), 'exclude',
        segments=[dict(start=0, end=4), dict(start=5, end=8)])
    def segments(**kwargs):
        for row in rows:
            for interval in row['valid_segments']:
                yield dict(**interval, sequence=row['sequence'][interval['start']:interval['end']],
                    record_id=row['record_id'], shared_sequence=row['shared_sequence'],
                    occurrences=row['occurrences'])
    return SimpleNamespace(identity='oncoref:fixture', manifest={'reference': reference},
                           coverage={'fixture': True}, records=lambda: iter(rows), segments=segments)


def test_containment_respects_segment_boundaries_and_keeps_every_source(bundle):
    matches = peptide_containment(['SIINFEKL', 'GILGFVFTL', 'GILGFVF'], bundle)
    assert [m['role'] for m in matches['SIINFEKL']] == ['allow', 'exclude']
    assert [m['record_id'] for m in matches['GILGFVF']] == ['representative']
    assert matches['SIINFEKL'][1]['start'] == 2


def dataset():
    sequence = 'SIINFEKLGILGFVFTLACDEFGHIKLMNPQRSTV'
    antigen = VaccineAntigen(kind='CTA', amino_acids=sequence, source_identifier='cta', gene_name='PRAME',
        targetable_mask=TargetableMask((AminoAcidInterval(0, len(sequence)),)),
        tumor_specificity=attestation())
    epitopes = [replace(_make_epitope(p, ic50=20., offset=offset, source_sequence=sequence),
                       prediction_id=str(index), overlaps_targetable=True)
                for index, (p, offset) in enumerate([('SIINFEKL', 0), ('GILGFVFTL', 8),
                                                     ('ACDEFGHIK', 17), ('LMNPQRSTV', 26)])]
    result = EpitopeDataset.from_predictions(epitopes, config=EpitopeConfig(score_expr='1.0'),
        antigens={e.prediction_id: antigen for e in epitopes})
    result.provenance = (InputProvenance('input', 'antigens', '', '', InputScope(
        patient_id='patient', reference_assembly='GRCh38', annotation='ensembl:115',
        host_taxon_id=9606, mhc_alleles=('HLA-A*02:01',))),)
    return result


def test_target_decisions_drive_shared_selection_and_offline_replay(bundle, tmp_path):
    data = dataset()
    apply_cta_targeting(data, bundle)
    assert [targeting_evidence(e)['status'] for e in data.epitopes] == [
        'excluded', 'eligible', 'review', 'unassessed']
    assert [e.is_targetable for e in data.epitopes] == [False, True, False, False]
    peptide = VaccinePeptide(antigen=data.antigens['0'], epitopes=data.epitopes)
    assert [e.sequence for e in peptide.target_epitopes] == ['GILGFVFTL']
    assert peptide.self_epitopes[0].self_reference_match.sources[0].protein_id == 'unsafe'
    path = tmp_path / 'native.tsv'
    data.save(path)
    again = EpitopeDataset.load(path)
    assert again.epitopes == data.epitopes
    assert again.selection['source_evidence'] == data.selection['source_evidence']
    from vaxrank.antigen_biology import antigen_biology
    view = antigen_biology(peptide.antigen, peptides=[peptide])
    assert any('Unsafe' in row['status'] and 'unsafe' in row['sources'] for row in view['targeting'])
    assert 'mutant' not in view['sequence_html']


def test_target_reference_scope_mismatch_fails(bundle):
    data = dataset()
    data.provenance = (replace(data.provenance[0], scope=replace(data.provenance[0].scope,
                                                               annotation='ensembl:93')),)
    with pytest.raises(ValueError, match='annotation'):
        apply_cta_targeting(data, bundle)


def test_replay_preserves_cta_reference_under_score_overrides():
    from vaxrank.external_rescoring import report_config
    from msgspec.structs import replace as replace_config
    data = dataset()
    data.config = replace_config(data.config, cta_sequence_bundle='current')
    report = SimpleNamespace(dataset=data, source_format='epitopes')
    resolved = report_config([report], EpitopeConfig(score_expr='2.0'))
    assert resolved.score_expr == '2.0' and resolved.cta_sequence_bundle == 'current'
    with pytest.raises(ValueError, match='Saved CTA targeting'):
        report_config([report], EpitopeConfig(cta_sequence_bundle='/different/reference'))


def test_current_expression_admission_uses_qualified_sequence_and_preserves_missing(bundle, tmp_path):
    path = tmp_path / 'expression.tsv'
    path.write_text('feature_id\tpatient_tpm\n' + PRAME + '\t12\nENSG00000111111\tNaN\n')
    admission = admit_expression(path, contract(), CTAAdmissionPolicy(2, 'TPM'), bundle,
                                 dataset().provenance[0].scope)
    assert [d.status for d in admission.decisions] == ['admitted', 'unsafe_reference']
    assert admission.decisions[1].value is None
    antigen, = admission.admitted_antigens
    assert antigen.amino_acids == 'SIINFEKLGILGFVFTL' and antigen.protein_ids == ('representative',)
    saved = tmp_path / 'admission.json'
    admission.save(saved)
    assert CTAExpressionResult.load(saved) == admission


def test_identical_qualified_cta_sequences_keep_all_accessions_in_one_candidate(bundle, tmp_path):
    from copy import deepcopy
    records = list(bundle.records())
    duplicate = deepcopy(records[0])
    duplicate['record_id'] = 'equivalent-protein'
    duplicate['occurrences'][0]['transcript_id'] = 'ENST00000398742.2'
    records.append(duplicate)
    bundle.records = lambda: iter(records)
    path = tmp_path / 'expression.tsv'
    path.write_text('feature_id\tpatient_tpm\n' + PRAME + '\t12\n')
    admission = admit_expression(path, contract(), CTAAdmissionPolicy(2, 'TPM'), bundle,
                                 dataset().provenance[0].scope)
    antigen, = admission.admitted_antigens
    assert antigen.protein_ids == ('equivalent-protein', 'representative')
    assert len(antigen.transcript_ids) == len(admission.decisions[0].assessment.equivalent_occurrences) == 2


def test_current_cta_cli_reports_unsafe_regions_and_replays_without_reference(bundle, tmp_path, monkeypatch):
    from .test_cta_input import argv
    from .test_external_rescoring import ContextPredictor
    from .test_epitope_dataset import forbidden
    from vaxrank.cli.entry_point import run_cli
    monkeypatch.setattr('vaxrank.cta_targeting.load_bundle', lambda path: bundle)
    monkeypatch.setattr('mhctools.cli.predictors_from_args', lambda args: [ContextPredictor()])
    monkeypatch.setattr('vaxrank.cli.entry_point.resolve_ensembl_release',
                        lambda args: setattr(args, 'genome', EnsemblRelease(115)))
    monkeypatch.setattr('vaxrank.cta_input.capture_local_hitlist_evidence', lambda *args: {})
    path = tmp_path / 'expression.tsv'
    path.write_text('feature_id\tpatient_tpm\n' + PRAME + '\t12\n')
    output, saved = tmp_path / 'output', tmp_path / 'native.tsv'
    args = argv(path, output)
    args[args.index('--ensembl-release') + 1] = '115'
    # Explicit report paths exercise the normal route without extra product audits.
    index = args.index('--output-dir')
    args[index:index + 2] = ['--output-html-report', str(output / 'report.html')]
    run_cli(args + ['--cta-sequence-bundle', 'current', '--output-epitopes', str(saved),
        '--vaccine-peptide-length', '16', '--config-text', 'epitopes.score_expr=1.0',
        '--config-text', 'vaccine_peptides.combined_score_expr=target_epitope_score'])
    native = EpitopeDataset.load(saved)
    assert any(targeting_evidence(e)['status'] == 'excluded' for e in native.epitopes)
    assert any(e.is_targetable for e in native.epitopes)
    assert 'Do not target: Unsafe containment' in (output / 'report.html').read_text()
    monkeypatch.setattr('vaxrank.cta_targeting.load_bundle', forbidden)
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    run_cli(['--input-epitopes', str(saved), '--output-html-report', str(output / 'replay.html'),
             '--no-processing-aware-annotation'])
    assert 'Do not target: Unsafe containment' in (output / 'replay.html').read_text()
