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
import sys
import logging

from isovar.args.variant_sequences import make_variant_sequences_arg_parser
from isovar.args.rna_reads import allele_reads_generator_from_args
from isovar.args.variants import variants_from_args

from topiary.commandline_args.mhc import (
    add_mhc_args,
    mhc_alleles_from_args,
    mhc_binding_predictor_from_args,
)

from .core_logic import ranked_vaccine_peptides, dataframe_from_ranked_list
from .report import ascii_report_from_ranked_vaccine_peptides

# inherit all commandline options from Isovar
arg_parser = make_variant_sequences_arg_parser(
    prog="vaxrank",
    description=(
        "Select personalized vaccine peptides from cancer variants, "
        "expression data, and patient HLA type."),
)

add_mhc_args(arg_parser)

arg_parser.add_argument(
    "--output-csv",
    default="vaccine-peptides.csv",
    help="Name of CSV file which contains predicted sequences")

arg_parser.add_argument(
    "--output-ascii-report",
    default="vaccine-peptides-report.txt",
    help="Path to ASCII vaccine peptide report")

vaccine_peptide_group = arg_parser.add_argument_group("Vaccine peptide options")
vaccine_peptide_group.add_argument(
    "--vaccine-peptide-length",
    default=25,
    type=int,
    help="Number of amino acids in the vaccine peptides (default %(default)s)")

vaccine_peptide_group.add_argument(
    "--padding-around-mutation",
    default=0,
    type=int,
    help=(
        "Number of off-center windows around the mutation to consider "
        "as vaccine peptides (default %(default)s)"
    ))

vaccine_peptide_group.add_argument(
    "--max-vaccine-peptides-per-mutation",
    default=1,
    type=int,
    help="Number of vaccine peptides to generate for each mutation")

vaccine_peptide_group.add_argument(
    "--max-mutations-in-report",
    default=10,
    type=int,
    help="Number of mutations to report")

def main(args_list=None):
    """
    Script to generate vaccine peptide predictions from somatic cancer variants,
    patient HLA type, and tumor RNA-seq data.

    Example usage:
        vaxrank
            --vcf somatic.vcf \
            --bam rnaseq.bam \
            --vaccine-peptide-length 25 \
            --output-csv vaccine-peptides.csv
    """
    if args_list is None:
        args_list = sys.argv[1:]

    logging.basicConfig(level=logging.DEBUG)
    args = arg_parser.parse_args(args_list)
    print(args)

    variants = variants_from_args(args)
    print(variants)

    mhc_alleles = mhc_alleles_from_args(args)
    print("MHC alleles: %s" % (mhc_alleles,))
    mhc_predictor = mhc_binding_predictor_from_args(args)

    # generator that for each variant gathers all RNA reads, both those
    # supporting the variant and reference alleles
    reads_generator = allele_reads_generator_from_args(args)

    ranked_list = ranked_vaccine_peptides(
        reads_generator=reads_generator,
        mhc_predictor=mhc_predictor,
        vaccine_peptide_length=args.vaccine_peptide_length,
        padding_around_mutation=args.padding_around_mutation,
        max_vaccine_peptides_per_variant=args.max_vaccine_peptides_per_mutation,
        min_reads_supporting_cdna_sequence=args.min_reads_supporting_variant_sequence)

    df = dataframe_from_ranked_list(ranked_list)
    print(df)
    df.to_csv(args.output_csv, index=False)

    ascii_report = ascii_report_from_ranked_vaccine_peptides(
        ranked_variants_with_vaccine_peptides=ranked_list,
        mhc_alleles=mhc_alleles,
        variants=variants,
        bam_path=args.bam)

    print(ascii_report)

    with open(args.output_ascii_report, "w") as f:
        f.write(ascii_report)
