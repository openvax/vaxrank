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

from argparse import ArgumentParser
from isovar.cli.rna_reads import allele_reads_generator_from_args
from isovar.cli.translation import add_translation_args
from isovar.cli.variant_sequences import make_variant_sequences_arg_parser
from mhctools.cli import (
    add_mhc_args,
    mhc_alleles_from_args,
    mhc_binding_predictor_from_args,
)

import serializable
from varcode.cli import variant_collection_from_args

from .core_logic import ranked_vaccine_peptides
from .report import (
    make_ascii_report,
    make_html_report,
    make_pdf_report,
    make_csv_report,
    make_minimal_neoepitope_report,
    TemplateDataCreator,
    PatientInfo,
)

logger = logging.getLogger(__name__)


def new_run_arg_parser():
    # inherit commandline options from Isovar
    arg_parser = make_variant_sequences_arg_parser(
        prog="vaxrank",
        description=(
            "Select personalized vaccine peptides from cancer variants, "
            "expression data, and patient HLA type."),
    )
    add_translation_args(arg_parser)
    add_mhc_args(arg_parser)
    add_vaccine_peptide_args(arg_parser)
    add_output_args(arg_parser)
    add_optional_output_args(arg_parser)
    add_supplemental_report_args(arg_parser)
    return arg_parser


def cached_run_arg_parser():
    arg_parser = ArgumentParser(
        prog="vaxrank",
        description=(
            "Select personalized vaccine peptides from cancer variants, "
            "expression data, and patient HLA type."),
    )
    arg_parser.add_argument(
        "--input-json-file",
        default="",
        help="Path to JSON file containing results of vaccine peptide report")
    add_output_args(arg_parser)
    add_optional_output_args(arg_parser)
    add_supplemental_report_args(arg_parser)
    return arg_parser


# Lets the user specify whether they want to see particular sections in the report.
def add_optional_output_args(arg_parser):
    manufacturability_args = arg_parser.add_mutually_exclusive_group(required=False)
    manufacturability_args.add_argument(
        "--include-manufacturability-in-report",
        dest="manufacturability",
        action="store_true")
    manufacturability_args.add_argument(
        "--no-manufacturability-in-report",
        dest="manufacturability",
        action="store_false")
    arg_parser.set_defaults(manufacturability=True)

    wt_epitope_args = arg_parser.add_mutually_exclusive_group(required=False)
    wt_epitope_args.add_argument(
        "--include-non-overlapping-epitopes-in-report",
        dest="wt_epitopes",
        action="store_true",
        help="Set to true to include a report section for each vaccine peptide containing "
             "strong binders that do not overlap the mutation")

    wt_epitope_args.add_argument(
        "--no-non-overlapping-epitopes-in-report",
        dest="wt_epitopes",
        action="store_false",
        help="Set to false to exclude report information for each vaccine peptide about "
             "strong binders that do not overlap the mutation")
    arg_parser.set_defaults(wt_epitopes=True)


def add_output_args(arg_parser):
    output_args_group = arg_parser.add_argument_group("Output options")

    output_args_group.add_argument(
        "--output-patient-id",
        default="",
        help="Patient ID to use in report")

    output_args_group.add_argument(
        "--output-csv",
        default="",
        help="Name of CSV file which contains predicted sequences")

    output_args_group.add_argument(
        "--output-ascii-report",
        default="",
        help="Path to ASCII vaccine peptide report")

    output_args_group.add_argument(
        "--output-html-report",
        default="",
        help="Path to HTML vaccine peptide report")

    output_args_group.add_argument(
        "--output-pdf-report",
        default="",
        help="Path to PDF vaccine peptide report")

    output_args_group.add_argument(
        "--output-json-file",
        default="",
        help="Path to JSON file containing results of vaccine peptide report")

    output_args_group.add_argument(
        "--output-xlsx-report",
        default="",
        help="Path to XLSX vaccine peptide report worksheet, one sheet per variant. This is meant "
             "for use by the vaccine manufacturer.")

    output_args_group.add_argument(
        "--output-neoepitope-report",
        default="",
        help="Path to XLSX neoepitope report, containing information focusing on short peptide "
             "sequences.")

    output_args_group.add_argument(
        "--num-epitopes-per-peptide",
        type=int,
        help="Number of top-ranking epitopes for each vaccine peptide to include in the "
             "neoepitope report.")

    output_args_group.add_argument(
        "--output-reviewed-by",
        default="",
        help="Comma-separated list of reviewer names")

    output_args_group.add_argument(
        "--output-final-review",
        default="",
        help="Name of final reviewer of report")

    output_args_group.add_argument(
        "--log-path",
        default="python.log",
        help="File path to write the vaxrank Python log to")

    output_args_group.add_argument(
        "--max-mutations-in-report",
        type=int,
        help="Number of mutations to report")


def add_vaccine_peptide_args(arg_parser):
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
        "--min-epitope-score",
        default=1e-10,
        type=float,
        help="Ignore epitopes whose normalized score falls below this threshold")


def add_supplemental_report_args(arg_parser):
    report_args_group = arg_parser.add_argument_group("Supplemental report options")
    report_args_group.add_argument(
        "--cosmic_vcf_filename",
        default="",
        help="Local path to COSMIC vcf")


def check_args(args):
    if not (args.output_csv or
            args.output_ascii_report or
            args.output_html_report or
            args.output_pdf_report or
            args.output_json_file or
            args.output_xlsx_report or
            args.output_neoepitope_report):
        raise ValueError(
            "Must specify at least one of: --output-csv, "
            "--output-xlsx-report, "
            "--output-ascii-report, "
            "--output-html-report, "
            "--output-pdf-report, "
            "--output-neoepitope-report, "
            "--output-json-file")

def ranked_variant_list_with_metadata(args):
    """
    Computes all the data needed for report generation.

    Parameters
    ----------
    args : Namespace
      Parsed user args from this run

    Returns a dictionary containing 3 items:
    - ranked variant/vaccine peptide list
    - a dictionary of command-line arguments used to generate it
    - patient info object
    """
    if hasattr(args, 'input_json_file'):
        with open(args.input_json_file) as f:
            data = serializable.from_json(f.read())
            # the JSON data from the previous run will have the older args saved, which may need to
            # be overridden with args from this run (which all be output related)
            data['args'].update(vars(args))

            # if we need to truncate the variant list based on max_mutations_in_report, do that here
            if len(data['variants']) > args.max_mutations_in_report:
                data['variants'] = data['variants'][:args.max_mutations_in_report]
            return data

    # get various things from user args
    mhc_alleles = mhc_alleles_from_args(args)
    logger.info("MHC alleles: %s", mhc_alleles)
    variants = variant_collection_from_args(args)
    logger.info("Variants: %s", variants)
    # generator that for each variant gathers all RNA reads, both those
    # supporting the variant and reference alleles
    reads_generator = allele_reads_generator_from_args(args)
    mhc_predictor = mhc_binding_predictor_from_args(args)

    ranked_list, variants_count_dict = ranked_vaccine_peptides(
        reads_generator=reads_generator,
        mhc_predictor=mhc_predictor,
        vaccine_peptide_length=args.vaccine_peptide_length,
        padding_around_mutation=args.padding_around_mutation,
        max_vaccine_peptides_per_variant=args.max_vaccine_peptides_per_mutation,
        min_alt_rna_reads=args.min_alt_rna_reads,
        min_variant_sequence_coverage=args.min_variant_sequence_coverage,
        min_epitope_score=args.min_epitope_score,
        num_mutant_epitopes_to_keep=args.num_epitopes_per_peptide,
        variant_sequence_assembly=args.variant_sequence_assembly)

    ranked_list_for_report = ranked_list[:args.max_mutations_in_report]

    patient_info = PatientInfo(
        patient_id=args.output_patient_id,
        vcf_paths=variants.sources,
        bam_path=args.bam,
        mhc_alleles=mhc_alleles,
        num_somatic_variants=len(variants),
        num_coding_effect_variants=variants_count_dict['num_coding_effect_variants'],
        num_variants_with_rna_support=variants_count_dict['num_variants_with_rna_support'],
        num_variants_with_vaccine_peptides=variants_count_dict['num_variants_with_vaccine_peptides']
    )

    # return variants, patient info, and command-line args
    data = {
        'variants': ranked_list_for_report,
        'patient_info': patient_info,
        'args': vars(args),
    }
    logger.info('About to save args: %s', data['args'])

    # save JSON data if necessary. as of time of writing, vaxrank takes ~25 min to run,
    # most of which is core logic. the formatting is super fast, and it can
    # be useful to save the data to be able to iterate just on the formatting
    if args.output_json_file:
        with open(args.output_json_file, 'w') as f:
            f.write(serializable.to_json(data))
            logger.info('Wrote JSON report data to %s', args.output_json_file)

    return data


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

    if "--input-json-file" in args_list:
        arg_parser = cached_run_arg_parser()
    else:
        arg_parser = new_run_arg_parser()

    args = arg_parser.parse_args(args_list)
    logging.config.fileConfig(
        pkg_resources.resource_filename(__name__, 'logging.conf'),
            defaults={'logfilename': args.log_path})

    logger.info(args)
    check_args(args)

    data = ranked_variant_list_with_metadata(args)
    ranked_variant_list = data['variants']
    patient_info = data['patient_info']
    args_for_report = data['args']

    ###################
    # CSV-based reports
    ###################
    if args.output_csv or args.output_xlsx_report:
        make_csv_report(
            ranked_variant_list,
            excel_report_path=args.output_xlsx_report,
            csv_report_path=args.output_csv)

    if args.output_neoepitope_report:
        make_minimal_neoepitope_report(
            ranked_variant_list,
            num_epitopes_per_peptide=args.num_epitopes_per_peptide,
            excel_report_path=args.output_neoepitope_report)

    ########################
    # Template-based reports
    ########################

    if not (args.output_ascii_report or args.output_html_report or args.output_pdf_report):
        return

    input_json_file = args.input_json_file if hasattr(args, 'input_json_file') else None
    template_data_creator = TemplateDataCreator(
        ranked_variants_with_vaccine_peptides=ranked_variant_list,
        patient_info=patient_info,
        final_review=args.output_final_review,
        reviewers=args.output_reviewed_by,
        args_for_report=args_for_report,
        input_json_file=input_json_file,
        cosmic_vcf_filename=args.cosmic_vcf_filename)
    template_data = template_data_creator.compute_template_data()

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
