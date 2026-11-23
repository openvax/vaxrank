"""Scientific comparisons remain legible and attributable in the report."""

from dataclasses import replace

import pytest
from mhctools.pred import Prediction

from vaxrank.candidate_epitope import CandidateEpitope, Peptide
from vaxrank.processing_prediction import ProcessingPrediction
from vaxrank.patient_info import PatientInfo
from vaxrank.report import TemplateDataCreator
from vaxrank.report_display import epitope_display, marked_sequence, model_inventory
from vaxrank.vaccine_antigen import AminoAcidInterval, TargetableMask
from vaxrank.vaccine_peptide import VaccinePeptide

from .test_assembly_antigens import candidate
from .test_product_reports import render, template


def prediction(model='netMHCpan', version='4.2c', kind='pMHC_affinity', value=78.82, score=0.99):
    return Prediction(kind=kind, predictor_name=model, predictor_version=version,
                      peptide='KRFHATISF', allele='HLA-B*27:05', value=value,
                      score=score, percentile_rank=0.004)


def ranked_pair(wt='KRFHATVSF', include_flurry=False):
    antigen, _ = candidate(0, sequence='KRFHATISF', prefix='HG')
    antigen = replace(antigen, targetable_mask=TargetableMask((AminoAcidInterval(8, 9),)))
    predictions = [prediction(), prediction(kind='pMHC_presentation', value=None)]
    if include_flurry:
        predictions += [prediction('mhcflurry', '2.3', value=50.),
                        prediction('mhcflurry', '2.3', 'pMHC_presentation', None, .8),
                        replace(prediction('mhcflurry', '2.3', 'antigen_processing', None, .4), allele='')]
    epitope = CandidateEpitope(sequence='KRFHATISF', source_sequence=antigen.amino_acids,
                              offset=2, predictions=predictions, occurs_in_reference=False,
                              overlaps_targetable=True, per_allele_scores={'HLA-B*27:05': .94},
                              comparators={'wt': Peptide(sequence=wt, predictions=[prediction(value=89.19)])})
    return [(antigen, [VaccinePeptide(antigen=antigen, epitopes=[epitope])])]


def test_full_report_leads_with_coverage_and_products_and_attributes_scores_and_comparator():
    data = template(ranked_pair(), candidate_details=True)
    html = render(data, 'html')
    assert (html.index('id="input-coverage"') < html.index('id="products"')
            < html.index('id="candidate-details"') < html.index('id="run-details"'))
    assert 'KRFHAT<span class="mutant">I</span>SF' in html
    assert 'KRFHAT<span class="reference-change">V</span>SF' in html
    assert 'NetMHCpan Score_EL: <strong>0.990</strong>' in html
    assert 'Presentation rank: 0.004%' in html
    assert 'Version 4.2c' in html
    assert 'Target IC50: <strong>78.82 nM</strong>' in html
    assert 'Comparator IC50: 89.19 nM' in html
    assert '1 candidate window(s)' not in html
    assert 'No results in the reported candidates.' in html
    assert 'Annotation and computed patient backgrounds do not establish observed patient-normal tissue or protein expression.' in html
    text = render(data, 'txt')
    assert text.index('Emitted products') < text.index('=== Candidate window details')
    assert 'NetMHCpan: Score_EL' in text


def test_flurry_processing_and_pepsickle_are_separate_recorded_evidence():
    ranked = ranked_pair(include_flurry=True)
    epitope = ranked[0][1][0].epitopes[0]
    pp = ProcessingPrediction(epitope.sequence, epitope.source_sequence, 'pepsickle', 2,
                              predictor_version='1.0', c_term_cleavage_prob=.8,
                              max_internal_cut_prob=.2, processing_score=.8)
    html = render(template(ranked, candidate_details=True,
                           processing_predictions_by_key={pp.key(): pp}), 'html')
    assert 'mhcflurry presentation: <strong>0.800</strong>' in html
    assert 'mhcflurry processing: 0.400' in html
    assert 'Pepsickle C-term: 0.80' in html
    assert 'Pepsickle max internal: 0.20' in html
    assert 'Cleavage composite: 0.80' in html
    assert 'Recorded cleavage annotations' in html


@pytest.mark.parametrize('wt,note', [('', 'Sequence unavailable'), ('KRFHATVS', 'Different lengths; no residue alignment inferred.')])
def test_missing_and_unequal_comparators_are_explicit(wt, note):
    html = render(template(ranked_pair(wt), candidate_details=True), 'html')
    assert note in html
    assert 'class="reference-change">' not in html.split('<table class="pure-table target-epitopes">', 1)[1].split('</table>', 1)[0]


def test_repeated_sequence_does_not_guess_a_target_location():
    ranked = ranked_pair()
    antigen, peptides = ranked[0]
    epitope = replace(peptides[0].epitopes[0], source_sequence='other source', offset=0)
    display = epitope_display(epitope, antigen, {}, prediction())
    assert not display['location_known']
    assert str(display['target']) == epitope.sequence


def test_disjoint_target_intervals_and_junction_are_preserved_and_escaped():
    assert str(marked_sequence('ACDEFG', [(1, 2), (4, 5)])) == (
        'A<span class="mutant">C</span>DE<span class="mutant">F</span>G')
    text = str(marked_sequence('ACDE', [(2, 2)]))
    assert 'AC<span class="junction" title="Target junction">|</span>DE' == text
    assert '<script>' not in str(marked_sequence('<script>'))


def test_disabled_and_unavailable_cleavage_are_distinct():
    assert model_inventory([], {}, False)[-1]['status'].startswith('Disabled for this run')
    assert model_inventory([], {}, None)[-1]['status'].startswith('No cleavage results recorded')


def test_wt_affinity_matches_predictor_version_without_fallback():
    creator = TemplateDataCreator.__new__(TemplateDataCreator)
    creator.processing_predictions_by_key = {}
    epitope = ranked_pair()[0][1][0].epitopes[0]
    epitope = replace(epitope, comparators={'wt': Peptide(sequence='KRFHATVSF', predictions=[
        prediction(version='4.1', value=200.), prediction(version='4.2c', value=89.19)])})
    assert creator.epitope_data(epitope, prediction())['WT IC50'] == '89.19 nM'
    assert creator.epitope_data(epitope, prediction(version='4.0'))['WT IC50'] == 'No prediction'


@pytest.mark.parametrize('values,expected', [([89.19, 89.19], '89.19 nM'),
                                             ([89.19, 200.0], 'Conflicting predictions')])
@pytest.mark.parametrize('role', ['wt', 'patient_background'])
def test_comparator_affinity_does_not_choose_between_conflicting_measurements(values, expected, role):
    from vaxrank.epitope_io import neoepitope_comparison_row

    creator = TemplateDataCreator.__new__(TemplateDataCreator)
    creator.processing_predictions_by_key = {}
    original = ranked_pair()[0][1][0].epitopes[0]
    evidence = ({'kind': 'patient_background_comparison', 'state': 'available',
                 'baseline_policy': 'focal', 'comparison': 'amino_acid_divergent'},) if role != 'wt' else ()
    for ordered in (values, list(reversed(values))):
        epitope = replace(original, input_evidence=evidence, comparators={role: Peptide(
            sequence='KRFHATVSF', predictions=[prediction(value=value) for value in ordered])})
        row = creator.epitope_data(epitope, prediction())
        assert row['WT IC50' if role == 'wt' else 'Comparison IC50'] == expected
        if role == 'patient_background':
            exported = neoepitope_comparison_row(epitope, 'HLA-B*27:05', 'netMHCpan', '4.2c')
            assert bool(exported['Predicted comparison pMHC affinity']) == (expected != 'Conflicting predictions')
            assert exported['Comparison affinity state'] == ('conflicting' if len(set(values)) > 1 else 'recorded')
            assert neoepitope_comparison_row(epitope, 'HLA-B*27:05')['Comparison affinity state'] == 'model_not_selected'


def test_product_highlighting_uses_consistent_mapped_masks_without_changing_products():
    from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
    from vaxrank.report_display import product_target_intervals

    ranked = ranked_pair()
    products = assemble_peptide_constructs(ranked, PeptideConstructConfig())
    before = [(p.sequence, p.antigens) for p in products]
    data = template(ranked, {'peptide': products})
    view = data['product_groups'][0]['products'][0]
    assert 'class="mutant">I</span>' in view['sequence_html']
    assert view['display_name'].startswith('GENE0')
    assert view['antigens'][0]['genes'] == ''
    assert [(p.sequence, p.antigens) for p in products] == before
    # Same source sequence with differing target masks cannot acquire a union highlight.
    antigen, peptides = ranked[0]
    alternative = replace(antigen, targetable_mask=TargetableMask((AminoAcidInterval(2, 3),)))
    ranked.append((alternative, [VaccinePeptide(antigen=alternative, epitopes=peptides[0].epitopes)]))
    masks = product_target_intervals(products[0], ranked)
    assert list(masks.values()) == [None]


@pytest.mark.parametrize('interval', [(1, 2), (11, 12)])
def test_adjacent_target_residue_does_not_become_a_deletion_junction(interval):
    antigen, peptides = ranked_pair()[0]
    antigen = replace(antigen, targetable_mask=TargetableMask((AminoAcidInterval(*interval),)))
    display = epitope_display(peptides[0].epitopes[0], antigen, {}, prediction())
    assert str(display['target']) == 'KRFHATISF'


def test_product_with_an_unavailable_source_mask_is_not_partially_highlighted():
    from copy import deepcopy
    from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
    from vaxrank.report_display import product_target_intervals

    ranked = ranked_pair()
    product, = assemble_peptide_constructs(ranked, PeptideConstructConfig())
    antigens = deepcopy(product.antigens)
    missing = dict(antigens[0]['sources'][0], source_identifier='unavailable-source')
    antigens[0]['sources'].append(missing)
    product = replace(product, antigens=antigens)
    assert list(product_target_intervals(product, ranked).values()) == [None]


def test_input_coverage_lists_families_without_evidence_once():
    data = template(ranked_pair())
    families = data["input_coverage"]["families"]
    empty = [f["name"] for f in families if f["status"] != "Supplied"]
    assert len(empty) >= 5
    listing = "No supplied evidence recorded: " + ", ".join(empty) + "."
    for kind in ("txt", "html"):
        text = render(data, kind)
        assert listing in text
        assert text.count("No supplied evidence recorded") <= 3
        assert "Recorded sources: 0" not in text
    # It starts its own line instead of running on from the sentence before.
    assert "\n" + listing in render(data, "txt")


def test_family_with_evidence_keeps_its_own_block_in_both_formats():
    data = template(ranked_pair())
    data["input_coverage"]["families"][0].update(status="Supplied", records=3, ranked=2)
    name = data["input_coverage"]["families"][0]["name"]
    assert "%s: Supplied\n  Recorded sources: 3" % name in render(data, "txt")
    html = render(data, "html")
    assert "<td>%s</td>\n<td>3</td>" % name in html or "<td>%s</td><td>3</td>" % name in html
    assert "<th>Input evidence</th>" not in html
    listing = [line for line in render(data, "txt").splitlines()
               if line.startswith("No supplied evidence recorded:")]
    assert listing and name not in listing[0]


def test_report_leads_with_survey_then_enriched_ranking_before_products_and_appendices():
    from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs

    first = ranked_pair()[0]
    second_antigen, _ = candidate(1, sequence='GILGFVFTL', prefix='AA')
    second_ep = CandidateEpitope(
        sequence='GILGFVFTL', source_sequence=second_antigen.amino_acids,
        offset=2, predictions=[prediction(value=120.0)], occurs_in_reference=False,
        overlaps_targetable=True, per_allele_scores={'HLA-B*27:05': 0.4},
    )
    second = (second_antigen, [VaccinePeptide(antigen=second_antigen, epitopes=[second_ep])])
    ranked = [first, second]
    products = {'peptide': assemble_peptide_constructs(
        ranked, PeptideConstructConfig(max_constructs=1))}
    data = template(ranked, products=products, candidate_details=True)
    html = render(data, 'html')
    # Survey -> Candidate ranking -> Emitted products -> Antigen biology -> Selection rules -> Details -> Methods
    assert (
        html.index('id="input-coverage"')
        < html.index('id="ranking"')
        < html.index('id="products"')
        < html.index('id="antigen-biology"')
        < html.index('id="selection-rules"')
        < html.index('id="candidate-details"')
        < html.index('id="run-details"')
    )
    assert '<th>Input evidence</th>' not in html
    assert 'Selected candidates (1)' in html
    assert 'Additional ranked candidates not emitted in products (1)' in html
    assert 'KRFHATISF · HLA-B*27:05 (78.8 nM)' in html
    assert html.count('Full antigen labels and source identities are in') == 1
    txt = render(data, 'txt')
    assert txt.index('=== Input coverage ===') < txt.index('=== Candidate ranking ===') < txt.index('=== Emitted products ===')


def test_blank_patient_id_is_omitted_not_printed_empty():
    def data_for(patient_id):
        return TemplateDataCreator(
            ranked_pair(), PatientInfo(patient_id), final_review='', reviewers='',
            args_for_report={'manufacturability': False, 'wt_epitopes': True},
            input_json_file=None, products=None).compute_template_data()

    for kind in ("txt", "html"):
        assert "Patient ID" not in render(data_for(""), kind)
        assert "Patient ID" in render(data_for("SID-001"), kind)


def test_variant_report_includes_overlapping_genes_only_when_multiple_or_different():
    from dataclasses import dataclass
    from types import SimpleNamespace

    @dataclass(frozen=True)
    class _VariantStub:
        contig: str
        start: int
        gene_names: tuple[str, ...]

    creator = TemplateDataCreator.__new__(TemplateDataCreator)
    creator.dna_vaf_by_variant = {}
    creator.source_vaf_by_variant = {}
    fragment = SimpleNamespace(
        gene_name='GENE_A',
        rna_vaf=0.5,
        dna_vaf=0.4,
        n_alt_reads=10,
        n_ref_reads=10,
        n_dna_alt=None,
        n_dna_ref=None,
        n_dna_overlapping=None,
        dna_evidence_method=None,
        n_other_reads=None,
        n_alt_fragments=None,
        n_ref_fragments=None,
        n_other_fragments=None,
        rna_evidence_method=None,
        rna_evidence_subject=None,
        sequence_source=None,
    )
    top_vp = SimpleNamespace(mutant_protein_fragment=fragment, combined_score=1.25)

    single_gene_variant = _VariantStub(contig='1', start=100, gene_names=('GENE_A',))
    single_data = creator._variant_data(single_gene_variant, top_vp)
    assert 'Overlapping genes' not in single_data

    multi_gene_variant = _VariantStub(contig='1', start=100, gene_names=('GENE_A', 'GENE_B'))
    multi_data = creator._variant_data(multi_gene_variant, top_vp)
    assert multi_data['Gene name'] == 'GENE_A'
    assert multi_data['Overlapping genes'] == 'GENE_A, GENE_B'
    keys = list(multi_data.keys())
    assert keys.index('Overlapping genes') == keys.index('Gene name') + 1


def test_cta_and_mutation_highlighting_strictly_exclude_unsafe_self_overlaps_and_callout_somatic():
    from vaxrank.antigen_biology import antigen_biology
    from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
    from vaxrank.report_display import effective_target_intervals
    from vaxrank.vaccine_antigen import (
        SelfReferenceMatch,
        SelfReferenceSource,
        TumorSpecificityAttestation,
        VaccineAntigen,
    )

    cta_seq = 'MKYSEKISYVYMKRNYKAMTKLGFK'
    cta = VaccineAntigen(
        kind='CTA',
        amino_acids=cta_seq,
        targetable_mask=TargetableMask((AminoAcidInterval(0, len(cta_seq)),)),
        tumor_specificity=TumorSpecificityAttestation(
            status='admitted', evidence_kind='cta_expression', evidence_source='test',
        ),
        gene_name='SSX1',
        gene_id='ENSG00000126752',
        source_identifier='SSX1:test',
        source_protein_offset=39,
    )
    eligible_ep = CandidateEpitope(
        sequence='KYSEKISYV', source_sequence=cta_seq, offset=1,
        predictions=[prediction(value=40.0)], occurs_in_reference=True,
        occurs_in_non_CTA_reference=False, overlaps_targetable=True,
        per_allele_scores={'HLA-B*27:05': 0.9},
        self_reference_match=SelfReferenceMatch('KYSEKISYV', False, 'CTA', source_provenance_complete=True),
    )
    unsafe_ep = CandidateEpitope(
        sequence='EKISYVYMK', source_sequence=cta_seq, offset=4,
        predictions=[prediction(value=25.0)], occurs_in_reference=True,
        occurs_in_non_CTA_reference=True, overlaps_targetable=False,
        per_allele_scores={'HLA-B*27:05': 0.95},
        self_reference_match=SelfReferenceMatch(
            'EKISYVYMK', True, 'CTA',
            sources=(SelfReferenceSource(gene_id='ENSG00000188167', gene_name='SSX7'),),
            source_provenance_complete=True,
        ),
    )
    eff = effective_target_intervals(cta, [eligible_ep, unsafe_ep])
    # Eligible [1, 10) minus Unsafe [4, 13) -> [1, 4) strictly non-overlapping with excluded [4, 13)
    assert eff['target'] == ((1, 4),)
    assert eff['excluded'] == ((4, 13),)
    assert eff['excluded_genes'] == ('SSX7',)

    cta_vp = VaccinePeptide(antigen=cta, epitopes=[eligible_ep, unsafe_ep])
    mut_source, mut_peptides = ranked_pair()[0]
    ranked = [(cta, [cta_vp]), (mut_source, mut_peptides)]
    products = {'peptide': assemble_peptide_constructs(ranked, PeptideConstructConfig())}
    data = template(ranked, products=products)
    html = render(data, 'html')

    # Mixed CTA + somatic run highlights top somatic targets separately AND marks windows
    assert 'Top somatic &amp; non-CTA targets (1 of 1 ranked)' in html
    assert 'class="targetable">KYS</span><span class="self-excluded"' in html
    assert 'Unsafe/self overlap: 9 aa (SSX7)' in html
    assert 'KRFHAT<span class="mutant">I</span>SF' in html

    # Full protein context renders boxed vaccine window and Unsafe/self exclusion strictly outside targetable
    full_seq = ('A' * 39) + cta_seq + ('G' * 20)
    bio = antigen_biology(
        cta,
        peptides=[cta_vp],
        protein_context=dict(
            full_sequence=full_seq,
            full_excluded_intervals=((43, 52),),
            full_excluded_genes=('SSX7',),
            window_excluded_intervals=((4, 13),),
            window_excluded_genes=('SSX7',),
        ),
    )
    assert bio['full_protein'] is not None
    assert bio['full_protein']['length'] == len(full_seq)
    assert bio['full_protein']['full_excluded_aa'] == 9
    assert 'class="vaccine-window"' in bio['full_protein']['sequence_html']
    assert 'class="self-excluded"' in bio['full_protein']['sequence_html']
    assert 'class="targetable"' in bio['full_protein']['sequence_html']


def test_direct_coverage_and_variant_counts_handle_none_coding_effect():
    from dataclasses import dataclass
    from types import SimpleNamespace
    from vaxrank.input_coverage import record_direct_coverage
    from vaxrank.vaxrank_results import VaxrankResults

    @dataclass(frozen=True)
    class _Var:
        contig: str
        start: int
        ref: str
        alt: str

    v1 = _Var(contig='1', start=100, ref='A', alt='T')
    v2 = _Var(contig='1', start=200, ref='G', alt='C')
    iso_coding = SimpleNamespace(
        variant=v1, top_gene_name='TP53',
        predicted_effect_modifies_protein_sequence=True,
        has_mutant_protein_sequence_from_rna=True,
        num_alt_fragments=10, num_ref_fragments=5,
        filter_values={'min_num_alt_reads': True},
    )
    iso_unannotatable = SimpleNamespace(
        variant=v2, top_gene_name='',
        predicted_effect_modifies_protein_sequence=None,
        has_mutant_protein_sequence_from_rna=False,
        num_alt_fragments=0, num_ref_fragments=0,
        filter_values={'min_num_alt_reads': False},
    )
    results = VaxrankResults(
        isovar_results=[iso_coding, iso_unannotatable],
        variant_to_vaccine_peptides_dict={v1: ['vp']},
        ranked_vaccine_peptides=[(v1, ['vp'])],
    )
    counts = results.variant_counts()
    assert counts['num_total_variants'] == 2
    assert counts['num_coding_effect_variants'] == 1
    assert counts['num_variants_with_rna_support'] == 1
    assert counts['num_variants_with_vaccine_peptides'] == 1

    dataset = SimpleNamespace(
        provenance=(SimpleNamespace(input_id='vcf_bam:test'),),
        selection={},
    )
    obs = record_direct_coverage(dataset, results, [v1, v2], window_selection_performed=True)
    assert obs['counts']['loaded'] == 2
    assert obs['counts']['coding_effect'] == 1
    assert obs['counts']['coding_with_rna_protein'] == 1
    assert obs['counts']['with_candidate_windows'] == 1
    assert obs['variant_outcomes'][1]['is_coding_nonsynonymous'] is None


def test_compact_variant_labels_and_single_window_biology_alignment():
    from vaxrank.antigen_biology import antigen_biology
    from vaxrank.report_display import compact_allele, compact_variant_label

    long_del = 'G' + 'A' * 736
    assert compact_allele('A') == 'A'
    assert compact_allele('TAGTGC') == 'TAGTGC'
    assert compact_allele(long_del) == 'GAAAAA…(737bp)'
    assert compact_variant_label('chr14 g.101980529G>A') == 'chr14 g.101980529G>A'
    assert compact_variant_label('chr2 g.189010890_189011626del' + long_del) == (
        'chr2 g.189010890_189011626del (737 bp)'
    )

    # Multiple overlapping epitopes on one SNV window must merge into ONE
    # contiguous alignment without triggering annotation queries.
    antigen, peptides = ranked_pair()[0]
    antigen = replace(antigen, amino_acids='HGKRFHATISFAK')
    ep1 = replace(peptides[0].epitopes[0], source_sequence=antigen.amino_acids, offset=2)
    ep2 = replace(
        ep1, sequence='RFHATISFA', offset=3,
        per_allele_scores={'HLA-B*27:05': 0.5},
        comparators={'wt': Peptide(sequence='RFHATVSFA', predictions=[prediction(value=120.0)])},
    )
    ep3 = replace(
        ep1, sequence='FHATISFAK', offset=4,
        per_allele_scores={'HLA-B*27:05': 0.3},
        comparators={'wt': Peptide(sequence='FHATVSFAK', predictions=[prediction(value=150.0)])},
    )
    vp_multi = VaccinePeptide(antigen=antigen, epitopes=[ep1, ep2, ep3])
    bio_merged = antigen_biology(antigen, peptides=[vp_multi])
    assert len(bio_merged['comparisons']) == 1
    assert bio_merged['comparisons'][0]['tumor'] == 'KRFHATISFAK'
    assert bio_merged['comparisons'][0]['reference'] == 'KRFHATVSFAK'
    assert bio_merged['comparisons'][0]['tumor_start'] == 3
    assert bio_merged['comparisons'][0]['tumor_end'] == 13


def test_compact_simple_peptide_table_and_collapsible_large_ledgers():
    from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
    from vaxrank.report_display import compact_manufacturability

    assert compact_manufacturability({
        'cterm_7mer_gravy_score': '-0.55714',
        'max_7mer_gravy_score': '0.9',
        'cysteine_count': '0',
        'n_terminal_methionine': '0',
        'aspartate_proline_bond_count': '0',
    }) == (('C-term GRAVY', '-0.56'), ('Max GRAVY', '0.9'), ('Flags', 'all 0'))
    assert compact_manufacturability({
        'cterm_7mer_gravy_score': 0.1,
        'max_7mer_gravy_score': 1.5,
        'cysteine_count': '2',
        'n_terminal_methionine': '0',
    }) == (('C-term GRAVY', '0.1'), ('Max GRAVY', '1.5'), ('Cys', '2'))

    # Distinct sequences so each peptide has 1 unique source mapping.
    seqs = ['SIINFEKLA', 'GILGFVFTL', 'NLVPMVATV']
    ranked = [candidate(i, sequence=seqs[i], prefix='AA') for i in range(3)]
    products = assemble_peptide_constructs(ranked, PeptideConstructConfig())
    data = TemplateDataCreator(
        ranked, PatientInfo('compact-test'), final_review='', reviewers='',
        args_for_report={'manufacturability': True, 'wt_epitopes': True},
        input_json_file=None, products={'peptide': products},
    ).compute_template_data()
    html = render(data, 'html')
    assert 'class="pure-table pure-table-bordered compact-peptides"' in html
    assert 'id="product-' in html
    assert 'C-term GRAVY:' in html
    # Manufacturability metrics have explicit separators so pdfkit never jams them together.
    assert '</span> &middot; <span class="mfg-item">' in html

    # Large direct variant_outcomes (>20) collapse into <details class="audit-details">.
    data['input_coverage']['direct'] = [{
        'input_id': 'vcf_bam:test',
        'producer_version': '3.68.0',
        'unit': 'parsed variant records',
        'stage': 'original direct pipeline',
        'counts': {'loaded': 25, 'annotatable': 25, 'coding_effect': 25,
                   'coding_with_rna_protein': 0, 'with_candidate_windows': None},
        'isovar_filter_thresholds': {},
        'variant_outcomes': [
            {
                'gene_name': f'G{i}', 'contig': '1', 'start': 100 + i,
                'ref': 'A' * 40 if i == 0 else 'A', 'alt': 'T',
                'is_coding_nonsynonymous': True, 'rna_support': False,
                'has_vaccine_peptide': None, 'isovar_filter_results': {'rna_support': False},
            }
            for i in range(25)
        ],
        'annotation_unavailable': [],
    }]
    html_big = render(data, 'html')
    assert '<summary>Original variant outcomes: vcf_bam:test (25 variants)</summary>' in html_big
    assert 'AAAAAA…(40bp)&gt;T' in html_big
