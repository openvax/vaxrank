"""Reports distinguish source windows from the products actually emitted."""

from dataclasses import replace
from io import StringIO

import pytest

from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs
from vaxrank.patient_info import PatientInfo
from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
from vaxrank.report import TemplateDataCreator, render_report

from .test_assembly_antigens import candidate


def template(ranked=(), products=None, **kwargs):
    return TemplateDataCreator(
        ranked, PatientInfo('product-report'), final_review='', reviewers='',
        args_for_report={'manufacturability': False, 'wt_epitopes': True},
        input_json_file=None, products=products, **kwargs).compute_template_data()


def render(data, extension):
    output = StringIO()
    render_report(data, output, 'templates/template.' + extension)
    return output.getvalue()


@pytest.mark.parametrize('extension', ['txt', 'html'])
def test_reports_use_actual_refined_packed_products_and_preserve_candidate_ranking(extension):
    ranked = [candidate(0), candidate(1, prefix='GGG'), candidate(2, 'GILGFVFTL')]
    options = dict(min_antigen_length_aa=9, max_antigen_length_aa=9, max_constructs=1)
    peptide = assemble_peptide_constructs(ranked, PeptideConstructConfig(**options))
    mrna = assemble_mrna_constructs(ranked, RNAConstructConfig(
        **options, antigens_per_construct=2, optimize_linkers=False, codon_method='use_best_codon'),
        host_taxon_id=9606)
    data = template(ranked, {'peptide': peptide, 'mrna': mrna})
    peptide_group, mrna_group = data['product_groups']
    assert peptide_group['products'][0]['sequence'] == 'SIINFEKLA'
    assert peptide_group['products'][0]['length'] == 9
    assert peptide_group['counts']['antigen_slots'] == 1
    assert mrna_group['counts']['antigen_slots'] == 2
    assert peptide_group['products'][0]['antigens'][0]['sources'] == 2
    assert mrna_group['products'][0]['sequence'] == mrna[0].sequence
    rows = data['ranked_candidates']
    assert len(rows) == 3 and data['variants'] == []
    assert rows[0]['window'] == 'AASIINFEKLAKK'
    assert peptide[0].name in rows[0]['products'] and peptide[0].name in rows[1]['products']
    assert peptide[0].name not in rows[2]['products']
    assert mrna[0].name in rows[2]['products']
    text = render(data, extension)
    assert 'Emitted products' in text and 'Candidate ranking' in text
    assert '9 aa' in text and 'AASIINFEKLAKK' in text
    assert 'Candidate window details' not in text
    assert 'Full nucleotide sequence, including polyA' in text
    if extension == 'html':
        assert '<link' not in text and 'yui.yahooapis.com' not in text
        assert '<style' in text and 'border-collapse' in text


@pytest.mark.parametrize('extension', ['txt', 'html'])
def test_no_assembly_and_successful_empty_assembly_are_distinct(extension):
    assert 'No construct files were emitted' in render(template(), extension)
    text = render(template(products={'peptide': [], 'mrna': []}), extension)
    assert 'No peptide products were emitted' in text
    assert 'No mrna products were emitted' in text
    assert 'No candidates passed ranking' in text
    assert 'request product assembly' not in text


def test_historical_product_does_not_claim_zero_source_evidence():
    ranked = [candidate(0)]
    product, = assemble_mrna_constructs(ranked, RNAConstructConfig(
        optimize_linkers=False, codon_method='use_best_codon'), host_taxon_id=9606)
    old_product = replace(product, antigens=[{k: v for k, v in a.items() if k != 'sources'}
                                           for a in product.antigens])
    data = template(ranked, {'mrna': [old_product]})
    mapping = data['product_groups'][0]['products'][0]['antigens'][0]
    assert mapping['sources'] is None and mapping['transcripts'] is None
    assert data['ranked_candidates'][0]['products'] == 'Selection mappings unavailable'
    assert 'unavailable candidates' in render(data, 'txt')


@pytest.mark.parametrize('extension', ['txt', 'html'])
def test_detailed_appendix_uses_target_terms_and_compact_provenance(extension):
    antigen, peptides = candidate(0)
    antigen = replace(antigen, source_metadata=antigen.source_metadata + (
        ('isovar_family', 'ncorf'), ('isovar_source_metadata_json', '{"raw_internal_key": "private_raw_payload"}')))
    peptides[0].antigen = antigen
    data = template([(antigen, peptides)], candidate_details=True)
    text = render(data, extension)
    assert 'Candidate window details' in text and 'Predicted target epitopes' in text
    assert 'Mutant epitopes' not in text and 'Wildtype epitopes' not in text
    assert 'private_raw_payload' not in text and 'raw_internal_key' not in text
    assert 'shared-reconstruction' in text


def test_mixed_historical_and_current_products_do_not_claim_unselected_sources():
    ranked = [candidate(0), candidate(1, 'GILGFVFTL')]
    old, = assemble_peptide_constructs(ranked[:1], PeptideConstructConfig())
    old = replace(old, antigens=[])
    current = assemble_mrna_constructs(ranked[1:], RNAConstructConfig(
        optimize_linkers=False, codon_method='use_best_codon'), host_taxon_id=9606)
    data = template(ranked, {'peptide': [old], 'mrna': current})
    assert data['ranked_candidates'][0]['products'] == 'Selection mappings incomplete'
    assert data['ranked_candidates'][1]['products'] == current[0].name


@pytest.mark.parametrize('extension', ['txt', 'html'])
def test_shared_product_reports_attribute_score_and_separate_candidate_contexts(extension):
    first, second = candidate(0), candidate(1, prefix='GGG')
    second[1][0].combined_score_expr = '2 * window_epitope_score'
    ranked = [first, second]
    [product] = assemble_peptide_constructs(ranked, PeptideConstructConfig(
        min_antigen_length_aa=9, max_antigen_length_aa=9))
    data = template(ranked, {'peptide': [product]})
    row = data['product_groups'][0]['products'][0]['antigens'][0]
    assert row['candidate_contexts'] == 2
    assert row['representative']['source_ids'] == ['S1']
    text = render(data, extension)
    assert 'Allocation source: GENE0' in text
    assert '2 distinct candidate contexts' in text
    assert 'Complete protein identity is not established by this window' in text
    assert 'Scores and RNA measurements remain source-specific' in text


@pytest.mark.parametrize('extension', ['txt', 'html'])
@pytest.mark.parametrize('modality', ['peptide', 'mrna'])
def test_shared_product_reports_distinguish_identical_vs_distinct_complete_proteins_and_separate_rna(
        extension, modality):
    from isovar import RNAFeature, RNASamplingGeometry, RNASamplingProfile, measure_rna_feature
    from isovar.rna_path_linkage import OriginalPathMatch
    from vaxrank.rna_support import RNAFeatureEvidence
    from vaxrank.vaccine_antigen import TumorSpecificityEvidence
    from .test_assembly_antigens import _full_protein_candidate

    window = 'SIINFEKLAKKVV'
    same_protein = 'MAAAA' + window + 'GGGGGGGGGG'
    diff_protein_a = 'MVVVV' + window + 'HHHHHHHH'
    diff_protein_b = 'MLLLL' + window + 'YYYYYYYYYYYY'

    def _rna_ev(feature_id, sample_id, n_reads, eff_template_len):
        obs = [OriginalPathMatch(('rg1', f'{feature_id}_q{i}', 0), 'placement', 0, 80)
               for i in range(n_reads)]
        p = RNASamplingProfile(
            'original-library',
            (RNASamplingGeometry(((0, 80),), eff_template_len),),
            1000000, 1000000, 'whole_library',
        )
        m = measure_rna_feature(
            RNAFeature(feature_id, 'path-' + feature_id, ((10, 49),), 'coding_interval'),
            obs,
            evidence_scope=(sample_id, 'lib-' + sample_id),
            profile=p,
        )
        return RNAFeatureEvidence(measurement=m, scope='window', protein_interval=(0, len(window)))

    # 1) Identical complete proteins across two paralogs with different scores, RNA, and temporal facts
    src_a, vps_a = _full_protein_candidate(
        'CTAG1A', 'ENSG00000184033', 'ENST_CTAG1A', 'ENSP_CTAG1A', same_protein, 5, 5 + len(window), score=1.6)
    src_b, vps_b = _full_protein_candidate(
        'CTAG1B', 'ENSG00000184034', 'ENST_CTAG1B', 'ENSP_CTAG1B', same_protein, 5, 5 + len(window), score=0.9)
    ag_a = replace(
        vps_a[0].antigen,
        rna_support=(_rna_ev('ctag1a', 'T0', 7, 100),),
        tumor_specificity=replace(
            vps_a[0].antigen.tumor_specificity,
            evidence_records=(TumorSpecificityEvidence(
                evidence_kind='patient_tumor_expression', evidence_source='rnaseq',
                sample_id='T0', patient_specific=True, passed=True,
                numeric_value=42.5, unit='TPM', threshold=5.0),),
        ),
    )
    ag_b = replace(
        vps_b[0].antigen,
        rna_support=(_rna_ev('ctag1b', 'T1', 3, 140),),
        tumor_specificity=replace(
            vps_b[0].antigen.tumor_specificity,
            status='overridden', override_reason='paralog review',
            evidence_records=(TumorSpecificityEvidence(
                evidence_kind='patient_tumor_expression', evidence_source='rnaseq',
                sample_id='T1', patient_specific=True, passed=False,
                numeric_value=2.0, unit='TPM', threshold=5.0),),
        ),
    )
    ranked_same = [(src_a, [replace(vps_a[0], antigen=ag_a)]), (src_b, [replace(vps_b[0], antigen=ag_b)])]
    opts = dict(min_antigen_length_aa=len(window), max_antigen_length_aa=len(window), antigens_per_construct=1)
    products_same = (
        assemble_peptide_constructs(ranked_same, PeptideConstructConfig(**opts))
        if modality == 'peptide'
        else assemble_mrna_constructs(ranked_same, RNAConstructConfig(**opts, optimize_linkers=False), host_taxon_id=9606)
    )
    data_same = template(ranked_same, {modality: products_same})
    text_same = render(data_same, extension)
    assert f'Identical complete supplied source protein ({len(same_protein)} aa) across recorded origins' in text_same
    assert 'Source-protein groups (1): CTAG1A, CTAG1B (ENSP_CTAG1A, ENSP_CTAG1B)' in text_same
    assert 'Displayed representative: CTAG1A (CTA)' in text_same
    assert 'Alternative origin: CTAG1B (CTA)' in text_same
    assert '7 reads, 7 fragments' in text_same and '3 reads, 3 fragments' in text_same
    assert '10 reads' not in text_same and '10 fragments' not in text_same
    assert 'T0: 42.5 TPM' in text_same and 'T1: 2 TPM' in text_same
    assert 'admission overridden' in text_same
    # Source-annotation ranking stays separate per origin
    assert len(data_same['ranked_candidates']) == 2

    # 2) Distinct complete source proteins sharing one emitted window
    ranked_diff = [
        _full_protein_candidate('PARALOG1', 'ENSG00000100001', 'ENST_P1', 'ENSP_P1',
                                diff_protein_a, 5, 5 + len(window), score=1.5),
        _full_protein_candidate('PARALOG2', 'ENSG00000100002', 'ENST_P2', 'ENSP_P2',
                                diff_protein_b, 5, 5 + len(window), score=1.1),
    ]
    products_diff = (
        assemble_peptide_constructs(ranked_diff, PeptideConstructConfig(**opts))
        if modality == 'peptide'
        else assemble_mrna_constructs(ranked_diff, RNAConstructConfig(**opts, optimize_linkers=False), host_taxon_id=9606)
    )
    data_diff = template(ranked_diff, {modality: products_diff})
    text_diff = render(data_diff, extension)
    assert '2 distinct complete source proteins share this window' in text_diff
    assert 'Source-protein groups (2):' in text_diff
    assert f'PARALOG1 (ENSP_P1) [complete, {len(diff_protein_a)} aa]' in text_diff
    assert f'PARALOG2 (ENSP_P2) [complete, {len(diff_protein_b)} aa]' in text_diff
