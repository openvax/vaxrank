# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import absolute_import, print_function, division

from varcode import load_vcf_fast as load_vcf
from pysam import AlignmentFile
from isovar.protein_sequence import variants_to_protein_sequences
from isovar.default_parameters import (
    MIN_TRANSCRIPT_PREFIX_LENGTH,
    MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
    PROTEIN_SEQUENCE_LENGTH,
    MAX_PROTEIN_SEQUENCES_PER_VARIANT,
    MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
    MIN_READ_MAPPING_QUALITY,
)

from .vaccine_peptide import VaccinePeptide

def vaccine_peptides(
        vcf_path,
        rna_bam_path,
        vaccine_peptide_length=25,
        num_windows_around_mutation=None,
        min_reads_supporting_rna_sequence=MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
        max_transcript_mismatches=MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
        max_protein_sequences_per_variant=1,
        min_mapping_quality=MIN_READ_MAPPING_QUALITY):
    variants = load_vcf(vcf_path)
    rna_reads = AlignmentFile(rna_bam_path)
    if num_windows_around_mutation is None:
        # allow mutation at any position other than first or last amino
        # acid in the sequence
        num_windows_around_mutation = vaccine_peptide_length - 2
    protein_sequence_length = vaccine_peptide_length + num_windows_around_mutation - 1
    protein_sequence_generator = variants_to_protein_sequences(
        variants,
        rna_reads,
        transcript_id_whitelist=None,
        protein_sequence_length=protein_sequence_length,
        min_reads_supporting_rna_sequence=min_reads_supporting_rna_sequence,
        min_transcript_prefix_length=min_transcript_prefix_length,
        max_transcript_mismatches=max_transcript_mismatches,
        max_protein_sequences_per_variant=max_protein_sequences_per_variant,
        min_mapping_quality=min_mapping_quality)

    for variant, protein_sequences in protein_sequence_generator:
        for protein_sequence_record in protein_sequences:
            amino_acid_sequence = protein_sequence_record.amino_acid_sequence
            n_amino_acids = len(amino_acid_sequence)
            for i in range(n_amino_acids - vaccine_peptide_length + 1):
                subsequence = amino_acid_sequence[i:i + vaccine_peptide_length]
                yield VaccinePeptide(
                    amino_acid_sequence=subsequence,
                    mutation_start_offset=mutation_start_offset,
                    mutation_end_offset=mutation_end_offset,
                    n_mutant_residues=mutation_end_offset - mutation_start_offset,
                    deletion=is_deletion,
                    mutation_distance_from_edge=mutation_distance_from_edge,
                    epitope_score=epitope_score,
                    start_offset_in_protein=start_offset_in_protein + i,
                    variant_chr=variant.original_contig,
                    variant_pos=variant.original_start,
                    variant_ref=variant.original_ref,
                    variant_alt=variant.original_alt)

def vaccine_peptides_dataframe(vcf_path, rna_bam_path, **kwargs):
    columns = OrderedDict()
    for field_name in VaccinePeptide._fields:
        columns[field_name] = []
    """
        ("chr", []),
        ("pos", []),
        ("ref", []),
        ("alt", []),
        ("amino_acid_sequence", []),
        # offsets of codons containing mutant nucleotides, even if
        # synonymous
        ("mutant_codon_start_offset", []),
        ("mutant_codon_end_offset", []),
        ("num_mutant_codons", [])
        # first non-synonymous codons resulting from the mutation
        ("mutant_amino_acid_start_offset", []),
        ("mutant_amino_acid_end_offset", []),

    """
    vaccine_peptide_generator = vaccine_peptides(
        vcf_path=vcf_path,
        rna_bam_path=rna_bam_path, **kwargs)

    for vaccine_peptide in vaccine_peptide_generator:
        for field_name in VaccinePeptide._fields:
            columns[field_name].append(getattr(vaccine_peptide, field_name))

def vaccine_peptides_dataframe_from_args(args):
    return vaccine_peptides_dataframe(
        vcf_path=args.vcf,
        rna_bam_path=args.bam,
        vaccine_peptide_length=args.vaccine_peptide_length,
        min_reads_supporting_rna_sequence=args.min_reads,
        max_transcript_mismatches=args.max_reference_transcript_mismatches,
        max_protein_sequences_per_variant=args.max_protein_sequences_per_variant,
        min_mapping_quality=args.min_mapping_quality)
