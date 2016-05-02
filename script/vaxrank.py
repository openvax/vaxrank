#!/usr/bin/env python
#
# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function, division, absolute_import
import logging
import argparse

import varcode
from pysam import AlignmentFile

from isovar.protein_sequence import variants_to_protein_sequences_dataframe
from isovar.default_parameters import (
    MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
    MIN_TRANSCRIPT_PREFIX_LENGTH,
    MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
    PROTEIN_SEQUENCE_LEGNTH,
    MAX_PROTEIN_SEQUENCES_PER_VARIANT,
)

from isovar import args

args.extend_parser_with_necessary_args(parser)

parser = argparse.ArgumentParser()

parser.add_argument(
    "--vcf",
    required=True,
    help="VCF file containing somatic variants")

parser.add_argument(
    "--bam",
    required=True,
    help="BAM file containing RNAseq reads")

parser.add_argument(
    "--genome",
    default=None)

parser.add_argument(
    "--min-reads",
    type=int,
    default=MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE)

parser.add_argument(
    "--protein-sequence-length",
    default=PROTEIN_SEQUENCE_LEGNTH,
    type=int)

parser.add_argument(
    "--max-sequences-per-variant",
    type=int,
    default=MAX_PROTEIN_SEQUENCES_PER_VARIANT)

parser.add_argument(
    "--max-reference-transcript-mismatches",
    type=int,
    default=MAX_REFERENCE_TRANSCRIPT_MISMATCHES)

parser.add_argument(
    "--min-transcript-prefix-length",
    type=int,
    default=MIN_TRANSCRIPT_PREFIX_LENGTH,
    help=(
        "Number of nucleotides before the variant we try to match against "
        "a reference transcript. Values greater than zero exclude variants "
        "near the start codon of transcripts without 5' UTRs."))

parser.add_argument(
    "--min-mapping-quality",
    type=int,
    default=0,
    help="Minimum MAPQ value to allow for a read")

parser.add_argument(
    "--output",
    default="isovar-translate-variants-results.csv",
    help="Name of CSV file which contains predicted sequences")

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    args = parser.parse_args()

    print(args)

    variants = varcode.load_vcf(
        args.vcf,
        genome=args.genome)

    samfile = AlignmentFile(args.bam)
    df = variants_to_protein_sequences_dataframe(
        variants=variants,
        samfile=samfile,
        protein_sequence_length=args.protein_sequence_length,
        min_reads_supporting_rna_sequence=args.min_reads,
        min_transcript_prefix_length=args.min_transcript_prefix_length,
        max_transcript_mismatches=args.max_reference_transcript_mismatches,
        max_protein_sequences_per_variant=args.max_sequences_per_variant,
        min_mapping_quality=args.min_mapping_quality)
    print(df)
    df.to_csv(args.output)
