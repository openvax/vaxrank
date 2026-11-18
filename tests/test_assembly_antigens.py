"""Product budgets count emitted sequences, not alternative source identifiers."""

import json

import pytest

from vaxrank.assembly_antigens import antigens_for_assembly, summarize_products
from vaxrank.candidate_epitope import CandidateEpitope
from vaxrank.core_logic import SOURCE_AGNOSTIC_RANKING_RULES
from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs, write_mrna_outputs
from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs, write_peptide_outputs
from vaxrank.vaccine_antigen import (
    AminoAcidInterval, TargetableMask, TumorSpecificityAttestation, VaccineAntigen,
)
from vaxrank.vaccine_peptide import VaccinePeptide


def candidate(index, sequence='SIINFEKLA', prefix='AA'):
    antigen = VaccineAntigen(
        kind='mutation', amino_acids=prefix + sequence + 'KK',
        targetable_mask=TargetableMask((AminoAcidInterval(len(prefix), len(prefix) + len(sequence)),)),
        tumor_specificity=TumorSpecificityAttestation(
            status='admitted', evidence_kind='synthetic_fixture', evidence_source='test', patient_specific=True),
        gene_name='GENE%d' % index, gene_id='gene-%d' % index, transcript_ids=('transcript-%d' % index,),
        protein_ids=('protein-%d' % index,), species='homo_sapiens', source_identifier='candidate-%d' % index,
        source_metadata=(('isovar_evidence_id', 'shared-reconstruction'),))
    epitope = CandidateEpitope(sequence=sequence, source_sequence=antigen.amino_acids, offset=len(prefix),
        prediction_id='prediction-%d' % index, occurs_in_reference=False, overlaps_targetable=True,
        per_allele_scores={'HLA-A*02:01': 1.})
    return antigen, [VaccinePeptide(antigen=antigen, epitopes=[epitope],
        combined_score_expr='window_epitope_score', ranking_rules=SOURCE_AGNOSTIC_RANKING_RULES)]


@pytest.mark.parametrize('modality', ['peptide', 'mrna'])
@pytest.mark.parametrize('content', ['mutation_spanning', 'minimal_epitope'])
@pytest.mark.parametrize('per_construct', [1, 2])
@pytest.mark.parametrize('alleles', [[], ['HLA-A*02:01']])
def test_distinct_products_preserve_every_source_before_caps(tmp_path, modality, content, per_construct, alleles):
    # Different candidate windows can yield the same final sequence. A later
    # duplicate must still map to the first product even after the cap is full.
    ranked = [candidate(0), candidate(1, prefix='GGG'), candidate(2, 'GILGFVFTL'), candidate(3, prefix='VV')]
    event = dict(reference_name='GRCh38', contig='1', start=123, ref='G', alt='A')
    evidence = {'shared-reconstruction': {'isovar_result': {'variant': event}}}
    kwargs = dict(antigen_content=content, min_antigen_length_aa=9, max_antigen_length_aa=9,
                  antigens_per_construct=per_construct, max_constructs=2 // per_construct)
    if modality == 'peptide':
        options = PeptideConstructConfig(**kwargs)
        products = assemble_peptide_constructs(ranked, options, target_alleles=alleles, source_evidence=evidence)
        write_peptide_outputs(products, tmp_path / 'products.fasta', tmp_path / 'manifest.json', options=options)
    else:
        options = RNAConstructConfig(**kwargs, optimize_linkers=False, codon_method='use_best_codon')
        products = assemble_mrna_constructs(ranked, options, target_alleles=alleles,
                                           source_evidence=evidence, host_taxon_id=9606)
        write_mrna_outputs(products, tmp_path, tmp_path / 'manifest.json')
        assert [a['aa'] for p in products for a in p.antigens] == ['SIINFEKLA', 'GILGFVFTL']
    assert len(products) == 2 // per_construct
    assert sum(len(p.antigen_names) for p in products) == 2
    mappings = [sources for p in products for sources in p.antigen_sources.values()]
    assert [len(sources) for sources in mappings] == [3, 1]
    assert {s['source_identifier'] for sources in mappings for s in sources} == {
        'candidate-%d' % i for i in range(4)}
    assert {s['gene_ids'][0] for s in mappings[0]} == {'gene-0', 'gene-1', 'gene-3'}
    assert {s['transcript_ids'][0] for s in mappings[0]} == {'transcript-0', 'transcript-1', 'transcript-3'}
    assert {s['protein_ids'][0] for s in mappings[0]} == {'protein-0', 'protein-1', 'protein-3'}
    assert all(s['source_event'] == event and s['evidence_group'] == 'shared-reconstruction'
               and s['species'] == 'homo_sapiens' for sources in mappings for s in sources)
    manifest = json.loads((tmp_path / 'manifest.json').read_text())
    assert [p['antigens'] for p in manifest] == [p.antigens for p in products]
    for product in products:
        sequence = product.sequence if modality == 'peptide' else product.cds_aa
        for antigen in product.antigens:
            assert sequence[antigen['start_aa']:antigen['end_aa']] == antigen['aa']
    assert all(p.target_epitope_score == 1. for _, peptides in ranked for p in peptides)
    summary = summarize_products(products)
    assert summary['unique_antigen_sequences'] == 2
    assert summary['source_candidates'] == 4
    assert summary['source_events'] == 1
    assert summary['source_candidates_without_event'] == 0


def test_requested_alternative_windows_keep_distinct_sequences():
    first, duplicate, different = candidate(0), candidate(1), candidate(2, 'GILGFVFTL')
    ranked = [(first[0], first[1] + duplicate[1] + different[1])]
    products = assemble_peptide_constructs(ranked, PeptideConstructConfig(
        candidates_per_slot=3, max_constructs=2, min_antigen_length_aa=9, max_antigen_length_aa=9))
    assert [p.sequence for p in products] == ['SIINFEKLA', 'GILGFVFTL']
    assert [len(next(iter(p.antigen_sources.values()))) for p in products] == [2, 1]


def test_colliding_source_names_do_not_merge_distinct_sequences():
    first, second = candidate(0), candidate(0, 'GILGFVFTL')
    # Same named source with two candidate windows and no unique producer IDs.
    groups = antigens_for_assembly([first, second], PeptideConstructConfig(), 'peptide')
    assert len(groups) == 2 and len({g.name for g in groups}) == 2
    assert groups[0].amino_acids != groups[1].amino_acids


@pytest.mark.parametrize('modality', ['peptide', 'mrna'])
def test_shared_sequence_labels_and_fasta_name_every_origin(tmp_path, modality):
    ranked = [candidate(0), candidate(1, prefix='GGG')]
    settings = dict(min_antigen_length_aa=9, max_antigen_length_aa=9, antigens_per_construct=1)
    if modality == 'peptide':
        products = assemble_peptide_constructs(ranked, PeptideConstructConfig(**settings))
        path = tmp_path / 'vaccine.fasta'
        write_peptide_outputs(products, path)
    else:
        products = assemble_mrna_constructs(ranked, RNAConstructConfig(**settings, optimize_linkers=False),
                                            host_taxon_id=9606)
        write_mrna_outputs(products, tmp_path)
        path = tmp_path / 'full.fasta'
    [product] = products
    assert all(gene in product.antigen_names[0] for gene in ('GENE0', 'GENE1'))
    assert product.antigens[0]['gene_names'] == ['GENE0', 'GENE1']
    header = path.read_text().splitlines()[0]
    assert 'GENE1' in header.split('genes=', 1)[1]


def test_legacy_product_fasta_description_remains_reproducible(tmp_path):
    from dataclasses import replace
    first, second = candidate(0), candidate(1)
    [product] = assemble_peptide_constructs([first, second], PeptideConstructConfig())
    old_name = 'GENE0 (mutation)'
    legacy = replace(product, antigen_names=[old_name], antigens=[
        dict({k: v for k, v in product.antigens[0].items() if k != 'gene_names'}, name=old_name)])
    path = tmp_path / 'legacy.fasta'
    write_peptide_outputs([legacy], path)
    assert path.read_text().splitlines()[0].endswith('genes=GENE0')


def _full_protein_candidate(
        gene_name, gene_id, transcript_id, protein_id, full_sequence,
        window_start, window_end, *, score=1.0, completeness='complete',
        source_identifier=None):
    full_antigen = VaccineAntigen(
        kind='CTA',
        amino_acids=full_sequence,
        targetable_mask=TargetableMask((AminoAcidInterval(0, len(full_sequence)),)),
        tumor_specificity=TumorSpecificityAttestation(
            status='admitted', evidence_kind='synthetic_fixture',
            evidence_source='test', patient_specific=True),
        gene_name=gene_name,
        gene_id=gene_id,
        transcript_ids=(transcript_id,) if transcript_id else (),
        protein_ids=(protein_id,) if protein_id else (),
        species='homo_sapiens',
        source_identifier=gene_id if source_identifier is None else source_identifier,
        source_protein_completeness=completeness,
    )
    window_antigen = full_antigen.sliced(window_start, window_end)
    epitope_seq = window_antigen.amino_acids[:9]
    epitope = CandidateEpitope(
        sequence=epitope_seq,
        source_sequence=window_antigen.amino_acids,
        offset=0,
        prediction_id='pred-%s-%d' % (gene_id or transcript_id or protein_id, window_start),
        occurs_in_reference=False,
        overlaps_targetable=True,
        per_allele_scores={'HLA-A*02:01': score},
    )
    return full_antigen, [VaccinePeptide(
        antigen=window_antigen,
        epitopes=[epitope],
        combined_score_expr='window_epitope_score',
        ranking_rules=SOURCE_AGNOSTIC_RANKING_RULES,
    )]


@pytest.mark.parametrize('modality', ['peptide', 'mrna'])
def test_identical_complete_proteins_group_together_while_distinct_proteins_stay_separate(
        tmp_path, modality):
    shared_window = 'SIINFEKLAKKVV'
    same_full = 'MAAAA' + shared_window + 'GGGGGGGGGG'
    diff_full_1 = 'MVVVV' + shared_window + 'HHHHHHHH'
    diff_full_2 = 'MLLLL' + shared_window + 'YYYYYYYYYYYY'
    diff_full_3 = 'MNNNN' + shared_window + 'DDDDDDDDDD'

    identical_ranked = [
        _full_protein_candidate('CTAG1A', 'ENSG00000184033', 'ENST_CTAG1A', 'ENSP_CTAG1A',
                                same_full, 5, 5 + len(shared_window), score=1.5),
        _full_protein_candidate('CTAG1B', 'ENSG00000184034', 'ENST_CTAG1B', 'ENSP_CTAG1B',
                                same_full, 5, 5 + len(shared_window), score=1.1),
    ]
    distinct_ranked = [
        _full_protein_candidate('PARALOG1', 'ENSG00000100001', 'ENST_P1', 'ENSP_P1',
                                diff_full_1, 5, 5 + len(shared_window), score=1.8),
        _full_protein_candidate('PARALOG2', 'ENSG00000100002', 'ENST_P2', 'ENSP_P2',
                                diff_full_2, 5, 5 + len(shared_window), score=1.4),
        _full_protein_candidate('PARALOG3', 'ENSG00000100003', 'ENST_P3', 'ENSP_P3',
                                diff_full_3, 5, 5 + len(shared_window), score=1.0),
    ]

    opts = dict(min_antigen_length_aa=len(shared_window), max_antigen_length_aa=len(shared_window),
                antigens_per_construct=1)
    for label, ranked, expected_groups, expected_lengths in (
            ('same', identical_ranked, 1, [len(same_full)]),
            ('distinct', distinct_ranked, 3, [len(diff_full_1), len(diff_full_2), len(diff_full_3)])):
        sub = tmp_path / label
        sub.mkdir()
        if modality == 'peptide':
            products = assemble_peptide_constructs(ranked, PeptideConstructConfig(**opts))
            write_peptide_outputs(products, sub / 'vaccine.fasta', sub / 'manifest.json')
            fasta_text = (sub / 'vaccine.fasta').read_text()
        else:
            products = assemble_mrna_constructs(
                ranked, RNAConstructConfig(**opts, optimize_linkers=False), host_taxon_id=9606)
            write_mrna_outputs(products, sub, sub / 'manifest.json')
            fasta_text = (sub / 'full.fasta').read_text()
        assert len(products) == 1
        manifest = json.loads((sub / 'manifest.json').read_text())
        seg = manifest[0]['antigens'][0]
        assert seg['aa'] == shared_window
        assert [s['is_representative'] for s in seg['sources']] == [True] + [False] * (len(ranked) - 1)
        groups = seg['source_protein_groups']
        assert len(groups) == expected_groups
        assert [g['length_aa'] for g in groups] == expected_lengths
        assert all(g['completeness'] == 'complete' and g['protein_sha256'] for g in groups)
        for source_pair in ranked:
            assert source_pair[0].gene_name in fasta_text


@pytest.mark.parametrize('modality', ['peptide', 'mrna'])
def test_duplicate_origins_overlapping_windows_and_missing_identifiers(tmp_path, modality):
    # Missing gene_name, gene_id and source_identifier; distinguished by transcript_ids.
    unnamed_a = VaccineAntigen(
        kind='mutation', amino_acids='SIINFEKLA',
        targetable_mask=TargetableMask((AminoAcidInterval(0, 9),)),
        tumor_specificity=TumorSpecificityAttestation(
            status='admitted', evidence_kind='synthetic_fixture', evidence_source='test', patient_specific=True),
        transcript_ids=('ENST_UNNAMED_1',), protein_ids=('ENSP_UNNAMED_1',),
    )
    unnamed_b = VaccineAntigen(
        kind='mutation', amino_acids='SIINFEKLA',
        targetable_mask=TargetableMask((AminoAcidInterval(0, 9),)),
        tumor_specificity=TumorSpecificityAttestation(
            status='admitted', evidence_kind='synthetic_fixture', evidence_source='test', patient_specific=True),
        transcript_ids=('ENST_UNNAMED_2',), protein_ids=('ENSP_UNNAMED_2',),
        source_protein_completeness='partial',
    )
    overlapping_window = VaccineAntigen(
        kind='mutation', amino_acids='SIINFEKLAKK',
        targetable_mask=TargetableMask((AminoAcidInterval(0, 11),)),
        tumor_specificity=TumorSpecificityAttestation(
            status='admitted', evidence_kind='synthetic_fixture', evidence_source='test', patient_specific=True),
        transcript_ids=('ENST_UNNAMED_1',), protein_ids=('ENSP_UNNAMED_1',),
    )

    def _vp(ag, offset, pred_id):
        ep = CandidateEpitope(sequence='SIINFEKLA', source_sequence=ag.amino_acids, offset=offset,
                              prediction_id=pred_id, occurs_in_reference=False, overlaps_targetable=True,
                              per_allele_scores={'HLA-A*02:01': 1.0})
        return VaccinePeptide(antigen=ag, epitopes=[ep], combined_score_expr='window_epitope_score',
                              ranking_rules=SOURCE_AGNOSTIC_RANKING_RULES)

    # Duplicate occurrence of unnamed_a + distinct unnamed_b (same 9-aa window) + distinct overlapping 11-aa window
    ranked = [
        (unnamed_a, [_vp(unnamed_a, 0, 'pred-u1')]),
        (unnamed_a, [_vp(unnamed_a, 0, 'pred-u1')]),
        (unnamed_b, [_vp(unnamed_b, 0, 'pred-u2')]),
        (overlapping_window, [_vp(overlapping_window, 0, 'pred-overlap')]),
    ]
    opts = dict(min_antigen_length_aa=9, max_antigen_length_aa=11, antigens_per_construct=1, max_constructs=5)
    if modality == 'peptide':
        products = assemble_peptide_constructs(ranked, PeptideConstructConfig(**opts))
    else:
        products = assemble_mrna_constructs(
            ranked, RNAConstructConfig(**opts, optimize_linkers=False), host_taxon_id=9606)
    # Distinct overlapping windows ('SIINFEKLA' vs 'SIINFEKLAKK') remain separate constructs!
    emitted_aas = [a['aa'] for p in products for a in p.antigens]
    assert emitted_aas == ['SIINFEKLA', 'SIINFEKLAKK']
    shared_seg = products[0].antigens[0]
    # Duplicate unnamed_a collapsed cleanly -> 2 distinct sources (unnamed_a, unnamed_b)
    assert len(shared_seg['sources']) == 2
    # Partial and unassessed sources never merge into one protein group
    assert [g['completeness'] for g in shared_seg['source_protein_groups']] == ['unassessed', 'partial']


def test_source_protein_completeness_validation():
    attestation = TumorSpecificityAttestation(
        status='admitted', evidence_kind='synthetic_fixture', evidence_source='test', patient_specific=True)
    mask = TargetableMask((AminoAcidInterval(0, 9),))
    with pytest.raises(ValueError, match='Unknown source_protein_completeness'):
        VaccineAntigen(kind='CTA', amino_acids='SIINFEKLA', targetable_mask=mask,
                       tumor_specificity=attestation, source_protein_completeness='guessed')
    with pytest.raises(ValueError, match='Only complete source proteins carry'):
        VaccineAntigen(kind='CTA', amino_acids='SIINFEKLA', targetable_mask=mask,
                       tumor_specificity=attestation, source_protein_completeness='partial',
                       source_protein_length=100)
    with pytest.raises(ValueError, match='nonzero offset requires explicit'):
        VaccineAntigen(kind='CTA', amino_acids='SIINFEKLA', targetable_mask=mask,
                       tumor_specificity=attestation, source_protein_offset=5,
                       source_protein_completeness='complete')
