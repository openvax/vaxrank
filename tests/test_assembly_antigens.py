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
        products = assemble_mrna_constructs(ranked, options, target_alleles=alleles, source_evidence=evidence)
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
