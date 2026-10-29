"""Event kind and coding reference survive conversion, products and reports."""

from dataclasses import replace
import json

import pandas as pd
import pytest

from vaxrank import AntigenDataset, antigens_from_isovar
from vaxrank.coding_context import antigen_source_label, coding_contexts
from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
from vaxrank.report import make_csv_report, make_pdf_report

from .test_assembly_antigens import candidate
from .test_isovar_inputs import attestation, result
from .test_product_reports import render, template


def with_context(index, kind, orf_type):
    antigen, peptides = candidate(index)
    context = dict(kind=kind, orf_type=orf_type, orf_id='orf-%d' % index,
                   source='fixture', reference_key='reference-1', source_transcript_id='tx-%d' % index)
    antigen = replace(antigen, kind='mutation', source_metadata=antigen.source_metadata + (
        ('isovar_family', 'variant_orfs'),
        ('isovar_candidate', json.dumps(dict(coding_reference=context)))))
    peptides[0].antigen = antigen
    return antigen, peptides


@pytest.mark.parametrize('kind,orf_type,label', [
    ('extra_orf', 'dORF', 'noncanonical dORF'),
    ('annotated', None, 'annotated CDS'),
    (None, None, 'coding context not recorded'),
    ('future_kind', 'future_type', 'future_kind / future_type'),
])
def test_producer_context_survives_native_antigen_round_trip(tmp_path, kind, orf_type, label):
    producer = result('ncorf')
    producer['models'][0]['reference'].update(kind=kind, orf_type=orf_type)
    dataset = antigens_from_isovar(producer, tumor_specificity=attestation())
    path = tmp_path / 'antigens.json'
    dataset.save(path)
    antigen, = AntigenDataset.load(path).antigens.values()
    assert antigen.kind == 'mutation'
    assert antigen_source_label(antigen) == 'mutation / ' + label
    context, = coding_contexts(antigen)
    assert context['kind'] == kind and context['orf_type'] == orf_type
    assert context['orf_id'] == 'noncoding_orf' and context['source'] == 'synthetic'
    assert context['reference_key'] == producer['models'][0]['reference']['reference_key']


def test_historical_catalog_recovery_requires_unambiguous_model_pair():
    dataset = antigens_from_isovar(result('ncorf'), tumor_specificity=attestation())
    antigen, = dataset.antigens.values()
    metadata = dict(antigen.source_metadata)
    value = json.loads(metadata['isovar_candidate'])
    value.pop('coding_reference')
    metadata['isovar_candidate'] = json.dumps(value)
    historical = replace(antigen, source_metadata=tuple(metadata.items()))
    assert coding_contexts(historical) == []
    assert 'not recorded' in antigen_source_label(historical)
    assert coding_contexts(historical, dataset.source_evidence) == coding_contexts(antigen)
    merged = replace(historical, source_metadata=historical.source_metadata + (
        ('isovar_evidence_id', 'other-origin'),))
    assert coding_contexts(merged, dataset.source_evidence) == []


def test_merged_explicit_contexts_retain_both_origins():
    first, _ = with_context(0, 'extra_orf', 'dORF')
    second, _ = with_context(1, 'annotated', None)
    merged = replace(first, source_metadata=first.source_metadata + second.source_metadata)
    assert {json.dumps(c, sort_keys=True) for c in coding_contexts(merged)} == {
        json.dumps(c, sort_keys=True) for c in coding_contexts(first) + coding_contexts(second)}
    label = antigen_source_label(merged)
    assert label.startswith('mutation / ') and 'noncanonical dORF' in label and 'annotated CDS' in label


def test_shared_product_retains_mixed_contexts_in_reports_and_manifest(tmp_path):
    ranked = [with_context(0, 'extra_orf', 'dORF'), with_context(1, 'annotated', None)]
    products = assemble_peptide_constructs(ranked, PeptideConstructConfig())
    product, = products
    sources = product.antigens[0]['sources']
    assert len(sources) == 2
    assert {s['coding_contexts'][0]['kind'] for s in sources} == {'extra_orf', 'annotated'}
    data = template(ranked, {'peptide': products}, candidate_details=True)
    for extension in ('html', 'txt'):
        text = render(data, extension)
        assert 'mutation / noncanonical dORF' in text and 'mutation / annotated CDS' in text
        assert 'orf-0' in text and 'orf-1' in text and 'reference-1' in text
    pdf = tmp_path / 'report.pdf'
    make_pdf_report(data, pdf, backend='weasyprint')
    assert pdf.read_bytes().startswith(b'%PDF-')

    # Old saved product mappings can display recovered context without edits.
    old = replace(product, antigens=[dict(a, sources=[{k: v for k, v in s.items()
        if k != 'coding_contexts'} for s in a['sources']]) for a in product.antigens])
    before = json.dumps(old.antigens, sort_keys=True)
    old_data = template(ranked, {'peptide': [old]})
    assert old_data['product_groups'][0]['products'][0]['antigens'][0]['context_labels'] == [
        'mutation / noncanonical dORF', 'mutation / annotated CDS']
    assert json.dumps(old.antigens, sort_keys=True) == before
    assert old.name == product.name and old.sequence == product.sequence


def test_tabular_reports_preserve_all_context_fields(tmp_path):
    ranked = [with_context(0, 'extra_orf', 'dORF'), with_context(1, 'annotated', None)]
    csv, xlsx = tmp_path / 'candidates.csv', tmp_path / 'candidates.xlsx'
    make_csv_report(ranked, csv_report_path=csv, excel_report_path=xlsx)
    for frame in (pd.read_csv(csv), pd.read_excel(xlsx, sheet_name='All')):
        assert list(frame.source_kind) == ['mutation', 'mutation']
        assert list(frame.coding_context) == ['noncanonical dORF', 'annotated CDS']
        contexts = [json.loads(value)[0] for value in frame.coding_references]
        assert [value['orf_id'] for value in contexts] == ['orf-0', 'orf-1']
        assert all(value['reference_key'] == 'reference-1' for value in contexts)
