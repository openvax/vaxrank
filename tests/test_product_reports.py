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
        **options, antigens_per_construct=2, optimize_linkers=False, codon_method='use_best_codon'))
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
        optimize_linkers=False, codon_method='use_best_codon'))
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
        optimize_linkers=False, codon_method='use_best_codon'))
    data = template(ranked, {'peptide': [old], 'mrna': current})
    assert data['ranked_candidates'][0]['products'] == 'Selection mappings incomplete'
    assert data['ranked_candidates'][1]['products'] == current[0].name
