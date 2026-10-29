"""Admission, sample support and catalog evidence remain separate in reports."""

from dataclasses import replace
import json

import pytest

from vaxrank import antigens_from_isovar
from vaxrank.antigen_report import antigen_evidence_sections
from vaxrank.vaccine_antigen import TumorSpecificityEvidence

from .test_isovar_inputs import attestation, result
from .test_isovar_inputs import dependencies  # noqa: F401 -- shared deterministic predictor fixture
from .test_product_reports import render, template
from .test_report_clarity import ranked_pair


def rows(sections):
    return {label: value for section in sections for label, value in section['rows']}


@pytest.mark.parametrize('extension', ['html', 'txt'])
def test_research_admission_and_review_flag_do_not_claim_validation(extension):
    antigen, peptides = ranked_pair()[0]
    assumption = replace(antigen.tumor_specificity,
        evidence_kind='explicit_research_inclusion_assumption',
        rationale_code='research_only_not_clinically_adjudicated', requires_review=False)
    antigen = replace(antigen, tumor_specificity=assumption)
    peptides[0].antigen = antigen
    text = render(template([(antigen, peptides)], candidate_details=True), extension)
    assert 'Construct admission' in text and 'explicit_research_inclusion_assumption' in text
    assert 'does not by itself establish tumor specificity' in text
    assert 'completed scientific review is not recorded by this flag' in text
    assert 'Tumor specificity status' not in text and 'Requires review' not in text


def test_colo_partialness_support_and_catalog_have_distinct_subjects():
    data = antigens_from_isovar(result('colo829_dynll2'), tumor_specificity=attestation(),
        options={'admission': {'allow_partial': True, 'allow_ambiguous': True}})
    antigen = next(a for a in data.antigens.values() if a.tumor_specificity.admits_construct)
    candidate = json.loads(dict(antigen.source_metadata)['isovar_candidate'])
    sections = antigen_evidence_sections(antigen, data.source_evidence)
    displayed = rows(sections)
    assert displayed['Partial coding candidate'] == 'yes'
    assert displayed['Ambiguous coding interpretation'] == 'yes'
    assert displayed['Allow partial coding candidates'] == 'yes'
    assert displayed['Allow ambiguous coding candidates'] == 'yes'
    assert displayed['RNA reconstruction support — fragments'] == candidate['rna_support']['fragments']
    assert displayed['Minimum coding-interval fragment depth'] == candidate['coding_support']['minimum_fragment_depth']
    reference = next(s for s in sections if s['title'] == 'Reference ORF evidence')
    assert 'not evidence observed in this sample' in reference['note']
    assert dict(reference['rows'])['Reference catalog'] == 'gencode_phase1'
    assert all('peptide-occurrence support' not in label for label in displayed)


def test_historical_candidate_does_not_invent_ambiguity_or_zero_support():
    data = antigens_from_isovar(result('ncorf'), tumor_specificity=attestation())
    antigen, = data.antigens.values()
    values = dict(antigen.source_metadata)
    candidate = json.loads(values['isovar_candidate'])
    for name in ('partial_coding_candidate', 'ambiguous_coding_candidate', 'rna_support'):
        candidate.pop(name)
    values['isovar_candidate'] = json.dumps(candidate)
    antigen = replace(antigen, source_metadata=tuple(values.items()))
    displayed = rows(antigen_evidence_sections(antigen))
    assert displayed['Partial coding candidate'] == 'no'  # explicit historical sequence_extent
    assert displayed['Ambiguous coding interpretation'] == 'Not recorded'
    assert displayed['Candidate RNA support'] == 'Not recorded'
    assert displayed['Allow partial coding candidates'] == 'Not recorded'
    assert displayed['Coding extent'] == 'complete_coding_candidate'
    candidate.pop('sequence_extent')
    values['isovar_candidate'] = json.dumps(candidate)
    antigen = replace(antigen, source_metadata=tuple(values.items()))
    assert rows(antigen_evidence_sections(antigen))['Partial coding candidate'] == 'Not recorded'


def test_sv_occurrence_does_not_borrow_union_counts():
    data = antigens_from_isovar(result('sv_occurrences'), tumor_specificity=attestation(),
        options={'admission': {'allow_ambiguous': True, 'allow_partial': True}})
    for antigen in data.antigens.values():
        candidate = json.loads(dict(antigen.source_metadata)['isovar_candidate'])
        displayed = rows(antigen_evidence_sections(antigen, data.source_evidence))
        assert displayed['ORF candidate union support — fragments'] == candidate['rna_support']['fragments']
        support = candidate['occurrence_rna_support'] or {}
        assert displayed['ORF occurrence support — fragments'] == support.get('fragments', 'Not recorded')


def test_each_attestation_record_keeps_scope_units_and_provenance():
    antigen = ranked_pair()[0][0]
    records = (TumorSpecificityEvidence('catalog', 'reference-study', subject_id='ORF-A',
                                      evidence_version='v1', patient_specific=False),
               TumorSpecificityEvidence('expression', 'tumor-assay', subject_id='gene-A',
                                        sample_id='sample-A', patient_specific=True,
                                        numeric_value=12., unit='TPM', threshold=3., passed=True))
    antigen = replace(antigen, tumor_specificity=replace(antigen.tumor_specificity,
        patient_specific=True, evidence_records=records))
    sections = [s for s in antigen_evidence_sections(antigen)
                if s['title'].startswith('Admission evidence')]
    assert len(sections) == 2
    first, second = (dict(s['rows']) for s in sections)
    assert first['Source'] == 'reference-study' and first['Declared patient-specific'] == 'no'
    assert second['Source'] == 'tumor-assay' and second['Sample'] == 'sample-A'
    assert second['Value'] == 12 and second['Units'] == 'TPM' and second['Threshold'] == 3


def test_multiple_origins_are_retained_without_guessing_a_reference_pair():
    data = antigens_from_isovar(result('ncorf'), tumor_specificity=attestation())
    antigen, = data.antigens.values()
    evidence, = data.source_evidence.values()
    antigen = replace(antigen, source_metadata=antigen.source_metadata + (('isovar_evidence_id', 'other-origin'),))
    catalog = dict(data.source_evidence, **{'other-origin': dict(evidence,
        options={'admission': {'allow_partial': True, 'allow_ambiguous': True, 'min_fragments': 9}})})
    sections = antigen_evidence_sections(antigen, catalog)
    sources = [dict(s['rows']) for s in sections if s['title'].startswith('Recorded Isovar source')]
    assert len(sources) == 2
    assert {s['Minimum admission fragments'] for s in sources} == {2, 9}
    assert not any(s['title'] == 'Reference ORF evidence' for s in sections)
    assert sum('pairing is not uniquely recorded' in s['note'] for s in sections) == 2


@pytest.mark.parametrize('value', ['[1, 2]', 'invalid JSON'])
def test_malformed_reserved_candidate_metadata_has_an_actionable_error(value):
    antigen = replace(ranked_pair()[0][0], source_metadata=(('isovar_candidate', value),))
    with pytest.raises(ValueError, match='Invalid isovar_candidate report metadata'):
        antigen_evidence_sections(antigen)


def test_isoform_zero_support_does_not_claim_observed_coding_sequence():
    data = antigens_from_isovar(result('isoforms'))
    antigen = next(iter(data.antigens.values()))
    displayed = rows(antigen_evidence_sections(antigen, data.source_evidence))
    assert displayed['RNA reconstruction support — fragments'] == 0
    assert displayed['RNA sequence status'] == 'not_observed'
    assert displayed['All coding bases observed'] == 'no'
    assert displayed['RNA sequence origin'] == 'reference_model'


def test_viral_reference_attribution_is_separate_from_sample_support():
    data = antigens_from_isovar(result('viral'))
    antigen = next(iter(data.antigens.values()))
    sections = antigen_evidence_sections(antigen, data.source_evidence)
    reference = next(s for s in sections if s['title'].startswith('Coding reference'))
    assert dict(reference['rows'])['Source species'] == 'human_papillomavirus_16'
    assert dict(reference['rows'])['Reference catalog'] == 'NCBI RefSeq NC_001526.4'
    assert dict(reference['rows'])['Protein matches reference'] == 'Not recorded'
    assert 'separate from RNA evidence' in reference['note']
    assert rows(sections)['RNA reconstruction support — fragments'] == 2


@pytest.mark.parametrize('family', ['fusion', 'sv'])
def test_path_junction_support_and_partner_context_are_visible(family):
    if family == 'fusion':
        from isovar import reconstruct_fusion
        from .test_isovar_fusion import example
        producer = reconstruct_fusion(*example())
    else:
        producer = result('sv')
    data = antigens_from_isovar(producer)
    antigen = next(iter(data.antigens.values()))
    displayed = rows(antigen_evidence_sections(antigen, data.source_evidence))
    expected = producer['paths'][0]['junctions'][0]['direct_support']['fragments']
    assert displayed['Junction 1 direct support — fragments'] == expected
    if family == 'fusion':
        assert displayed['Donor compatible transcripts'] == 'donor'
        assert displayed['Acceptor compatible transcripts'] == 'acceptor'


@pytest.mark.parametrize('family', ['ncorf', 'sv', 'fusion', 'isoforms', 'viral'])
@pytest.mark.usefixtures('dependencies')
def test_source_evidence_survives_fresh_and_replayed_full_report(family, tmp_path, monkeypatch):
    from vaxrank.cli.entry_point import run_cli
    from vaxrank.epitope_dataset import EpitopeDataset
    if family == 'fusion':
        from isovar import reconstruct_fusion
        from .test_isovar_fusion import example
        producer = reconstruct_fusion(*example())
    else:
        producer = result(family)
    dataset = antigens_from_isovar(producer, tumor_specificity=attestation(),
        options={'admission': {'allow_partial': True, 'allow_ambiguous': True}})
    # Bind these test records to the deterministic host/predictor fixture.
    dataset.scope = replace(dataset.scope, patient_id='patient', reference_assembly='GRCh38',
        annotation='ensembl:93', mhc_alleles=('HLA-A*02:01', 'HLA-B*07:02'))
    source, native = tmp_path / 'source.json', tmp_path / 'predictions.tsv'
    dataset.save(source)
    run_cli(['--input-antigens', str(source), '--mhc-predictor', 'random', '--ensembl-release', '93',
        '--output-epitopes', str(native), '--output-html-report', str(tmp_path / 'fresh.html'),
        '--report-candidate-details', '--no-processing-aware-annotation', '--vaccine-peptide-length', '10',
        '--config-text', 'epitopes.score_expr=1.0',
        '--config-text', 'vaccine_peptides.combined_score_expr=target_epitope_score'])
    assert EpitopeDataset.load(native).selection['source_evidence'] == dataset.source_evidence
    def forbidden(*args, **kwargs):
        raise AssertionError('Native report replay must not reconstruct or predict')
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    monkeypatch.setattr('vaxrank.isovar_input.reconstruct_isovar_antigens', forbidden)
    run_cli(['--input-epitopes', str(native), '--output-html-report', str(tmp_path / 'replay.html'),
        '--report-candidate-details', '--no-processing-aware-annotation'])
    fresh, replay = ((tmp_path / (name + '.html')).read_text() for name in ('fresh', 'replay'))
    for text in (fresh, replay):
        assert 'Candidate RNA and coding evidence' in text
        assert 'Recorded Isovar source' in text
        assert 'Construct admission' in text
        assert 'completed scientific review is not recorded' in text
    # Compare the complete candidate sections, including all evidence rows.
    def candidates(text):
        return text.split('<h2 id="candidate-details">', 1)[1].split('<h2 id="ranking">', 1)[0]
    assert candidates(fresh) == candidates(replay)
