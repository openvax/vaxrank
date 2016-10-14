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
import logging.config
import pkg_resources

from varcode.cli import variant_collection_from_args
from mhctools.cli import (
    add_mhc_args,
    mhc_alleles_from_args,
    mhc_binding_predictor_from_args,
)
from isovar.cli.variant_sequences import make_variant_sequences_arg_parser
from isovar.cli.rna_reads import allele_reads_generator_from_args

from .core_logic import ranked_vaccine_peptides, dataframe_from_ranked_list
from .report import (
    compute_template_data,
    make_ascii_report,
    make_html_report,
    make_pdf_report,
)


logging.config.fileConfig(pkg_resources.resource_filename(__name__, 'logging.conf'))
logger = logging.getLogger(__name__)


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
    default="",
    help="Name of CSV file which contains predicted sequences")

arg_parser.add_argument(
    "--output-ascii-report",
    default="",
    help="Path to ASCII vaccine peptide report")

arg_parser.add_argument(
    "--output-html-report",
    default="",
    help="Path to HTML vaccine peptide report")

arg_parser.add_argument(
    "--output-pdf-report",
    default="",
    help="Path to PDF vaccine peptide report")

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


    args = arg_parser.parse_args(args_list)
    logger.info(args)

    if (len(args.output_csv) == 0 and
            len(args.output_ascii_report) == 0 and
            len(args.output_html_report) == 0 and
            len(args.output_pdf_report) == 0):
        raise ValueError(
            "Must specify at least one of: --output-csv, "
            "--output-ascii-report, "
            "--output-html-report, "
            "--output-pdf-report")

    variants = variant_collection_from_args(args)
    logger.info(variants)

    mhc_alleles = mhc_alleles_from_args(args)
    logger.info("MHC alleles: %s", mhc_alleles)
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
    logger.debug(df)

    if args.output_csv:
        df.to_csv(args.output_csv, index=False)

    template_data = compute_template_data(
        ranked_variants_with_vaccine_peptides=ranked_list,
        mhc_alleles=mhc_alleles,
        variants=variants,
        bam_path=args.bam)

    if args.output_ascii_report:
        make_ascii_report(
            template_data=template_data,
            ascii_report_path=args.output_ascii_report)

    if args.output_html_report:
        make_html_report(
            template_data=template_data,
            html_report_path=args.output_html_report)

    if args.output_pdf_report:
        make_pdf_report(
            template_data=template_data,
            pdf_report_path=args.output_pdf_report)
