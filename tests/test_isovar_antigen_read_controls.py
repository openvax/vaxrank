"""Live antigen sources use Isovar's read policy and retain it on replay."""

import json

import pysam
import pytest
from isovar import FusionReference, IsoformTarget

from vaxrank import AntigenDataset, reconstruct_isovar_antigens


@pytest.fixture
def antigen_inputs(tmp_path):
    sequence = 'ATG' + 'GCTGAC' * 20 + 'TAA'
    fasta = tmp_path / 'reference.fa'
    fasta.write_text('>test\n' + sequence + '\n')
    bam_path = tmp_path / 'reads.bam'
    header = {'SQ': [{'SN': 'test', 'LN': len(sequence)}],
              'RG': [{'ID': 'library', 'SM': 'sample'}]}
    with pysam.AlignmentFile(bam_path, 'wb', header=header) as bam:
        for name, flag, mapq, qualities in [
            ('primary', 0, 60, True), ('ambiguous', 0, 0, True),
            ('secondary', 256, 60, False), ('duplicate', 1024, 60, True),
        ]:
            read = pysam.AlignedSegment(bam.header)
            read.query_name, read.query_sequence = name, sequence
            read.flag, read.mapping_quality = flag, mapq
            read.reference_id, read.reference_start = 0, 0
            read.cigarstring = f'{len(sequence)}M'
            if qualities:
                read.query_qualities = pysam.qualitystring_to_array('I' * len(sequence))
            read.set_tag('RG', 'library')
            bam.write(read)
    reference = FusionReference('T', 'test', 'test:1', 'test', '+',
                                ((0, len(sequence)),), sequence, 0, len(sequence))
    return bam_path, fasta, IsoformTarget('G', 'TEST', reference)


@pytest.mark.parametrize('family', ['viral', 'isoforms'])
@pytest.mark.parametrize('reads,expected', [
    ({}, 3),
    ({'min_mapping_quality': 10, 'use_secondary_alignments': False}, 1),
    ({'use_reads_without_base_qualities': False}, 2),
])
def test_live_antigen_read_policy_changes_evidence_and_survives_replay(
        antigen_inputs, family, reads, expected, tmp_path):
    bam_path, fasta, target = antigen_inputs
    options = {'reads': reads}
    with pysam.AlignmentFile(bam_path) as bam:
        inputs = dict(bam=bam, sample_id='sample', source='library', scratch_dir=tmp_path)
        if family == 'viral':
            inputs.update(genome_fasta=fasta, contigs=['test'], reference_cache=tmp_path)
            options['reconstruction'] = {'reference_profile': None}
        else:
            inputs['targets'] = [target]
        dataset = reconstruct_isovar_antigens(family, inputs=inputs, options=options)
    evidence, = dataset.source_evidence.values()
    assert evidence['isovar_result']['support']['target_fragments'] == expected
    policy = evidence['reconstruction']['controls']['reads']
    assert policy['min_mapping_quality'] == reads.get('min_mapping_quality', 0)
    assert policy['use_secondary_alignments'] == reads.get('use_secondary_alignments', True)
    assert policy['use_reads_without_base_qualities'] == reads.get('use_reads_without_base_qualities', True)
    assert not policy['use_duplicate_reads']
    assert not policy['use_soft_clipped_bases']
    saved = tmp_path / 'antigens.json'
    dataset.save(saved)
    assert AntigenDataset.load(saved) == dataset
    # The policy is portable JSON, not a serialized collector or callback.
    assert json.loads(json.dumps(policy)) == policy


@pytest.mark.parametrize('family', ['sv', 'viral', 'isoforms'])
def test_input_collector_is_rejected_instead_of_silently_replaced(family):
    from isovar import ReadCollector

    with pytest.raises(ValueError, match='options.reads'):
        reconstruct_isovar_antigens(family, inputs={'read_collector': ReadCollector()})


@pytest.mark.parametrize('family', ['viral', 'isoforms'])
def test_invalid_read_control_fails_before_acquisition(family):
    with pytest.raises(ValueError, match='unknown or inapplicable'):
        reconstruct_isovar_antigens(family, inputs={}, options={'reads': {'unsupported': 1}})
