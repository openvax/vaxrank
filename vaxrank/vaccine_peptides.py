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

from collections import OrderedDict

from varcode import load_vcf_fast as load_vcf
from pysam import AlignmentFile
from isovar.protein_sequence import variants_to_protein_sequences
from isovar.default_parameters import (
    MIN_TRANSCRIPT_PREFIX_LENGTH,
    MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
    MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
    MIN_READ_MAPPING_QUALITY,
)
from topiary.commandline_args import mhc_binding_predictor_from_args

from .vaccine_peptide import VaccinePeptide

def vaccine_peptides(
        vcf_path,
        rna_bam_path,
        mhc_binding_predictor,
        vaccine_peptide_length=25,
        max_offset_from_center=0,
        min_reads_supporting_rna_sequence=MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
        max_transcript_mismatches=MAX_REFERENCE_TRANSCRIPT_MISMATCHES,

        min_mapping_quality=MIN_READ_MAPPING_QUALITY):
    variants = load_vcf(vcf_path)
    rna_reads = AlignmentFile(rna_bam_path)

    protein_fragment_length = vaccine_peptide_length + 2 * max_offset_from_center

    protein_fragment_generator = variants_to_protein_sequences(
        variants,
        rna_reads,
        transcript_id_whitelist=None,
        protein_sequence_length=protein_fragment_length,
        min_reads_supporting_rna_sequence=min_reads_supporting_rna_sequence,
        min_transcript_prefix_length=min_transcript_prefix_length,
        max_transcript_mismatches=max_transcript_mismatches,
        max_protein_sequences_per_variant=1,
        min_mapping_quality=min_mapping_quality)

    full_sequence_dictionary = {}
    for variant, protein_sequence_records in protein_fragment_generator:
        for protein_sequence_record in protein_sequence_records:
            seq = protein_sequence_record.amino_acids
            full_sequence_dictionary[(variant, protein_sequence_record)] = seq

    epitopes_dict = mhc_binding_predictor.predict(full_sequence_dictionary)
    for (variant, protein_sequence_record), epitopes in epitopes_dict.items():
        print(variant)
        print(protein_sequence_record)
        print(epitopes)
        print("---")
        """
        epitope_predictions = mhc_binding_predictor.predict([amino_acids])

            n_amino_acids = len(amino_acids)
            midpoint = n_amino_acids // 2
            mutant_aa_start = protein_sequence_record.variant_aa_interval_start=24
            mutant_aa_end = protein_sequence_record.variant_aa_interval_end = 33

            for i in range(n_amino_acids - vaccine_peptide_length + 1):
                subsequence = amino_acids[i:i + vaccine_peptide_length]
                yield VaccinePeptide(
                    variant=variant,
                    isovar_protein_sequence=protein_sequence_record,
                    frameshift=protein_sequence_record.frameshift,
                    gene=protein_sequence_record.gene,
                    amino_acids=subsequence,)
        """
        """
        variant_chr=variant.original_contig,
        variant_pos=variant.original_start,
        variant_ref=variant.original_ref,
        variant_alt=variant.original_alt,
        amino_acid_sequence=subsequence,
        mutation_start_offset=mutation_start_offset,
        mutation_end_offset=mutation_end_offset,
        n_mutant_residues=mutation_end_offset - mutation_start_offset,
        deletion=is_deletion,
        mutation_distance_from_edge=mutation_distance_from_edge,
        epitope_score=epitope_score,
        start_offset_in_protein=start_offset_in_protein + i)
        """

def vaccine_peptides_dataframe(
        vcf_path, rna_bam_path, mhc_binding_predictor, **kwargs):
    columns = OrderedDict()
    for field_name in VaccinePeptide._fields:
        columns[field_name] = []

    vaccine_peptide_generator = vaccine_peptides(
        vcf_path=vcf_path,
        rna_bam_path=rna_bam_path,
        mhc_binding_predictor=mhc_binding_predictor,
        **kwargs)

    for vaccine_peptide in vaccine_peptide_generator:
        for field_name in VaccinePeptide._fields:
            columns[field_name].append(getattr(vaccine_peptide, field_name))

def vaccine_peptides_dataframe_from_args(args):
    print(args)
    mhc_binding_predictor = mhc_binding_predictor_from_args(args)
    return vaccine_peptides_dataframe(
        vcf_path=args.vcf,
        rna_bam_path=args.bam,
        mhc_binding_predictor=mhc_binding_predictor,
        vaccine_peptide_length=args.vaccine_peptide_length,
        min_reads_supporting_rna_sequence=args.min_reads,
        max_transcript_mismatches=args.max_reference_transcript_mismatches,
        max_protein_sequences_per_variant=1,
        min_mapping_quality=args.min_mapping_quality)
