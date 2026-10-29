"""Shared report tables preserve source identities and differing evidence."""

from dataclasses import replace
from html.parser import HTMLParser
import json

import pandas as pd
import pytest

from vaxrank.candidate_epitope import Peptide
from vaxrank.input_coverage import report_source_key
from vaxrank.native_serialization import to_native_json
from vaxrank.patient_info import PatientInfo
from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
from vaxrank.report import TemplateDataCreator, make_csv_report
from vaxrank.vaccine_antigen import AminoAcidInterval, TargetableMask
from vaxrank.vaccine_peptide import VaccinePeptide

from .test_product_reports import render, template
from .test_report_clarity import ranked_pair


def duplicate_pair():
    first, = ranked_pair()
    antigen, peptides = first
    other = replace(antigen, source_identifier='second-source')
    epitope = replace(peptides[0].epitopes[0], prediction_id='second-prediction',
        source_name='second-epitope-origin', predictions=[replace(p, source_sequence_name='second-predictor-origin')
                                                       for p in peptides[0].epitopes[0].predictions_flat()])
    return [first, (other, [VaccinePeptide(antigen=other, epitopes=[epitope])])]


def tables(data, kind):
    return [table for table in data['evidence_tables'] if table['kind'] == kind]


class Links(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.identities, self.targets = [], set()
        self.feed(html)

    def handle_starttag(self, tag, attributes):
        attributes = dict(attributes)
        if 'id' in attributes:
            self.identities.append(attributes['id'])
        if attributes.get('href', '').startswith('#'):
            self.targets.add(attributes['href'][1:])


@pytest.mark.parametrize('details', [False, True])
def test_identical_evidence_shares_tables_and_preserves_all_sources_and_products(tmp_path, details):
    ranked = duplicate_pair()
    products = assemble_peptide_constructs(ranked, PeptideConstructConfig())
    before = to_native_json(dict(antigens=[a for a, _ in ranked],
        epitopes=[p.epitopes for _, ps in ranked for p in ps], products=products))
    data = template(ranked, {'peptide': products}, candidate_details=details)
    assert len(data['ranked_candidates']) == len(data['source_directory']) == 2
    assert data['product_groups'][0]['products'][0]['antigens'][0]['source_ids'] == ['S1', 'S2']
    if details:
        group, = data['candidate_groups']
        assert [s['identifier'] for s in group['sources']] == ['S1', 'S2']
        assert len(tables(data, 'predictions')) == len(tables(data, 'window')) == 1
        assert all(t['sources'] == ['S1', 'S2'] for t in data['evidence_tables'])
    else:
        assert data['candidate_groups'] == data['evidence_tables'] == []
    for extension in ('html', 'txt'):
        report = render(data, extension)
        assert 'second-source' in report and ranked[0][0].source_identifier in report
        if details:
            assert 'second-epitope-origin' in report and 'second-predictor-origin' in report
    links = Links(render(data, 'html'))
    assert len(links.identities) == len(set(links.identities))
    assert links.targets <= set(links.identities)
    csv = tmp_path / 'candidates.csv'
    make_csv_report(ranked, csv_report_path=csv)
    assert len(pd.read_csv(csv)) == 2
    assert before == to_native_json(dict(antigens=[a for a, _ in ranked],
        epitopes=[p.epitopes for _, ps in ranked for p in ps], products=products))


@pytest.mark.parametrize('difference', [
    'flanks', 'parent_sequence', 'target_mask', 'raw_prediction', 'comparator',
    'support', 'admission', 'gene', 'input_evidence',
])
def test_different_biological_evidence_is_not_merged(difference):
    ranked = duplicate_pair()
    antigen, peptides = ranked[1]
    epitope = peptides[0].epitopes[0]
    if difference == 'flanks':
        epitope = replace(epitope, n_flank='GG', c_flank='WW')
    elif difference == 'parent_sequence':
        epitope = replace(epitope, source_sequence='GG' + epitope.source_sequence, offset=epitope.offset + 2)
    elif difference == 'target_mask':
        antigen = replace(antigen, targetable_mask=TargetableMask((AminoAcidInterval(7, 8),)))
    elif difference == 'raw_prediction':
        # Both display as 78.82 nM. Formatting must not decide equivalence.
        epitope = replace(epitope, predictions=[replace(p, value=78.821) if p.value is not None else p
                                              for p in epitope.predictions_flat()])
    elif difference == 'comparator':
        epitope = replace(epitope, comparators={'wt': Peptide(sequence='KRFHATASF')})
    elif difference == 'support':
        antigen = replace(antigen, source_metadata=antigen.source_metadata + (
            ('isovar_candidate', json.dumps(dict(rna_support=dict(fragments=7, reads=9)))),))
    elif difference == 'admission':
        antigen = replace(antigen, tumor_specificity=replace(antigen.tumor_specificity,
            requires_review=True, rationale_code='manual_review_pending'))
    elif difference == 'gene':
        antigen = replace(antigen, gene_name='OTHER_GENE', gene_id='other-gene')
    else:
        epitope = replace(epitope, input_evidence=(('recorded_support', 'separate_observation'),))
    ranked[1] = (antigen, [VaccinePeptide(antigen=antigen, epitopes=[epitope])])
    data = template(ranked, candidate_details=True)
    assert len(data['candidate_groups']) == 2
    if difference in {'raw_prediction', 'comparator', 'flanks'}:
        assert len(tables(data, 'predictions')) == 2
    if difference == 'parent_sequence':
        assert len(tables(data, 'window')) == 2
        contexts = [t['value']['contexts'][0] for t in tables(data, 'window')]
        assert contexts[0]['sequence'] != contexts[1]['sequence']
        assert contexts[0]['occurrences'][0]['offset'] != contexts[1]['occurrences'][0]['offset']
    report = render(data, 'html')
    if difference == 'support':
        assert 'RNA reconstruction support — fragments</th><td>7' in report
    if difference == 'admission':
        assert 'manual_review_pending' in report
    if difference == 'gene':
        assert 'OTHER_GENE' in report


def test_held_out_sources_remain_in_directory_and_share_only_exact_outcomes():
    ranked = duplicate_pair()
    sources = [dict(identity=a.display_identifier, source_key=report_source_key(a),
        gene=a.gene_name, family='variant_orfs', context='noncanonical dORF', inputs=[],
        rank=i, admitted=True, stage='Ranked', reasons=[], retained_epitopes=1)
        for i, (a, _) in enumerate(ranked, 1)]
    sources += [dict(sources[0], identity='held-out', source_key='held-out-key', rank=None,
                     admitted=False, stage='Admission', reasons=['Explicitly held out'], retained_epitopes=0)]
    data = TemplateDataCreator(ranked,
        PatientInfo('held-out-report', input_coverage={'sources': sources}),
        final_review='', reviewers='', args_for_report={'manufacturability': False, 'wt_epitopes': True},
        input_json_file=None,
        candidate_details=True).compute_template_data()
    assert len(data['source_directory']) == 3
    outcomes = data['input_coverage']['source_groups']
    assert [len(row['members']) for row in outcomes] == [2, 1]
    assert outcomes[-1]['stage'] == 'Admission' and outcomes[-1]['members'][0]['rank'] is None
    for extension in ('html', 'txt'):
        report = render(data, extension)
        assert 'held-out' in report and 'Explicitly held out' in report


def test_empty_report_has_no_dangling_source_links():
    data = template(products={'peptide': []}, candidate_details=True)
    assert data['source_directory'] == data['candidate_groups'] == data['evidence_tables'] == []
    links = Links(render(data, 'html'))
    assert links.targets <= set(links.identities)


def test_one_source_can_link_multiple_distinct_evidence_groups():
    antigen, peptides = ranked_pair()[0]
    other = VaccinePeptide(antigen=antigen, epitopes=[replace(peptides[0].epitopes[0], n_flank='GG')])
    data = template([(antigen, peptides), (antigen, [other])], candidate_details=True)
    source, = data['source_directory']
    assert source['groups'] == ['G1', 'G2']
    assert len(tables(data, 'predictions')) == 2
    links = Links(render(data, 'html'))
    assert links.targets <= set(links.identities)
