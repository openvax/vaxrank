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

from isovar.args.protein_sequences import (
    make_protein_sequences_arg_parser,
    protein_sequences_generator_from_args,
)

# inherit all commandline options from Isovar
arg_parser = make_protein_sequences_arg_parser(
    prog="vaxrank",
    description=(
        "Select personalized vaccine peptides from cancer variants, "
        "expression data, and patient HLA type."),
)

arg_parser.add_argument(
    "--output-csv",
    default="vaccine_peptides.csv",
    help="Name of CSV file which contains predicted sequences")

arg_parser.add_argument(
    "--vaccine-peptide-length",
    default=25,
    type=int)


def main(args_list=None):
    """
    Script to generate vaccine peptide predictions from somatic cancer variants,
    patient HLA type, and tumor RNA-seq data.

    Example usage:
        vaxrank
            --vcf somatic.vcf \
            --bam rnaseq.bam \
            --min-vaccine-peptide-length 23 \
            --max-vaccine-peptide-length 27 \
            --output-csv vaccine-peptides.csv
    """
    if args_list is None:
        args_list = sys.argv[1:]

    logging.basicConfig(level=logging.DEBUG)
    args = arg_parser.parse_args(args_list)

    for variant, protein_sequences in protein_sequences_generator_from_args(args):
        print(variant)
        print(protein_sequences)
    """
    df_vaccine_peptides = vaccine_peptides_dataframe_from_args(args)
    print(df_vaccine_peptides)
    if args.output_csv:
        df_vaccine_peptides.to_csv(args.output_csv)
    """
