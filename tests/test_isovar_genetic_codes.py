"""The native consumer validates the producer's effective code, not source nulls."""

import json

import pytest
from isovar import AlleleRead, ProteinSequenceCreator, ReadEvidence, reconstruct_variant_orfs
from isovar.dna import reverse_complement_dna
from orfdata import TranscriptORF
from pyensembl import EnsemblRelease
from varcode import Variant

from vaxrank import AntigenDataset, antigens_from_isovar
from .test_isovar_inputs import variant_inputs as variant_inputs


@pytest.mark.parametrize('name,table', [('MT-ND5-201', 2), ('MT-ND6-201', 2), ('TP53-001', 1)])
def test_native_annotated_effective_code_and_replay(name, table, tmp_path):
    genome = EnsemblRelease(87)
    transcript = genome.transcripts_by_name(name)[0]
    reference = TranscriptORF.from_annotation(transcript)
    positions = [p for start, end in sorted(transcript.exon_intervals)
                 for p in range(start, end + 1)]
    reverse = transcript.strand == '-'
    if reverse:
        positions.reverse()
    index = reference.start + 201
    sequence = reference.transcript_sequence
    ref, alt = sequence[index], 'C' if sequence[index] != 'C' else 'A'
    prefix, suffix = sequence[index - 72:index], sequence[index + 1:index + 94]
    if reverse:
        ref, alt = reverse_complement_dna(ref), reverse_complement_dna(alt)
        prefix, suffix = reverse_complement_dna(suffix), reverse_complement_dna(prefix)
    variant = Variant(transcript.contig, positions[index], ref, alt, genome)
    evidence = ReadEvidence.from_variant_and_allele_reads(
        variant, [AlleleRead(prefix, alt, suffix, name='synthetic-mito-control')])
    result = reconstruct_variant_orfs(
        variant, evidence, [reference], sample_id='synthetic', source='synthetic-RNA',
        protein_sequence_creator=ProteinSequenceCreator(
            protein_sequence_length=49, protein_sequence_preference='support', min_variant_sequence_coverage=1))
    assert result['schema'] == 'isovar.variant_rna_orfs.v2'
    result = json.loads(json.dumps(result))
    dataset = antigens_from_isovar(result)
    assert dataset.antigens
    assert {a.amino_acids for a in dataset.antigens.values()} == {
        c['protein_sequence'] for c in result['protein_candidates']}
    for antigen in dataset.antigens.values():
        candidate = json.loads(dict(antigen.source_metadata)['isovar_candidate'])
        assert candidate['translation_table'] == table
        assert candidate['translation_table_source'] == 'reference_context'
        assert candidate['rna_support']['fragments'] == 1
    stored, = dataset.source_evidence.values()
    assert stored['isovar_result']['models'][0]['reference']['translation_table'] is None
    assert stored['isovar_result'] == result
    path = tmp_path / 'native-genetic-code.antigens.json'
    dataset.save(path)
    assert AntigenDataset.load(path) == dataset


@pytest.mark.parametrize('field,value', [
    ('translation_table', None), ('translation_table', True), ('translation_table', '1'),
    ('translation_table', 999), ('translation_table', 2),
    ('translation_table_source', None), ('translation_table_source', 'reference_context'),
])
def test_native_v2_rejects_missing_or_conflicting_code(variant_inputs, field, value):
    result = reconstruct_variant_orfs(**variant_inputs)
    result['protein_candidates'][0][field] = value
    with pytest.raises(ValueError, match='translation.table|provenance'):
        antigens_from_isovar(result)


def test_native_v2_does_not_fall_back_when_code_is_missing(variant_inputs):
    result = reconstruct_variant_orfs(**variant_inputs)
    del result['protein_candidates'][0]['translation_table']
    with pytest.raises(ValueError, match='requires an integer'):
        antigens_from_isovar(result)
