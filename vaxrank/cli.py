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

from isovar import isovar_results_to_dataframe
from isovar.cli import (make_isovar_arg_parser, run_isovar_from_parsed_args,)
from mhctools.cli import (
    add_mhc_args,
    mhc_alleles_from_args,
    mhc_binding_predictor_from_args,
)

import pandas as pd
import serializable
from varcode.cli import variant_collection_from_args

from . import __version__
from .core_logic import run_vaxrank
from .gene_pathway_check import GenePathwayCheck
from .report import (
    make_ascii_report,
    make_html_report,
    make_pdf_report,
    make_csv_report,
    make_minimal_neoepitope_report,
    TemplateDataCreator,
)
from .patient_info import PatientInfo

logger = logging.getLogger(__name__)


def make_vaxrank_arg_parser():
    # create common parser with the --version flag
    parent_parser = ArgumentParser('parent', add_help=False)
    parent_parser.add_argument('--version', action='version', version='Vaxrank %s' % (__version__,))

    # inherit commandline options from Isovar
    arg_parser = make_isovar_arg_parser(
        prog="vaxrank",
        description=(
            "Select personalized vaccine peptides from cancer variants, "
            "expression data, and patient HLA type."),
        parents=[parent_parser],
    )
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
        "--output-reviewed-by",
        default="",
        help="Comma-separated list of reviewer names")

    output_args_group.add_argument(
        "--output-final-review",
        default="",
        help="Name of final reviewer of report")

    output_args_group.add_argument(
        "--output-passing-variants-csv",
        default="",
        help="Path to CSV file containing some metadata about every variant that has passed all "
             "variant caller filters")

    output_args_group.add_argument(
        "--output-isovar-csv",
        default="",
        help="Path to CSV file containing raw RNA counts and filtering metadata "
             "for all variants (generated by Isovar)")

    output_args_group.add_argument(
        "--log-path",
        default="python.log",
        help="File path to write the vaxrank Python log to")

    output_args_group.add_argument(
        "--max-mutations-in-report",
        default=None,
        type=int,
        help="Number of mutations to report")


def add_vaccine_peptide_args(arg_parser):
    vaccine_peptide_group = arg_parser.add_argument_group("Vaccine peptide options")
    vaccine_peptide_group.add_argument(
        "--vaccine-peptide-length",
        default=25,
        type=int,
        help="Number of amino acids in the vaccine peptides. (default: %(default)s)")

    vaccine_peptide_group.add_argument(
        "--padding-around-mutation",
        default=5,
        type=int,
        help=(
            "Number of off-center windows around the mutation to consider "
            "as vaccine peptides. (default: %(default)s)"
        ))

    vaccine_peptide_group.add_argument(
        "--max-vaccine-peptides-per-mutation",
        default=1,
        type=int,
        help=(
            "Number of vaccine peptides to generate for each mutation. "
            "(default: %(default)s)"
        ))

    vaccine_peptide_group.add_argument(
        "--min-epitope-score",
        default=1e-10,
        type=float,
        help=(
            "Ignore predicted MHC ligands whose normalized binding score "
            "falls below this threshold. (default: %(default)s)"))

    vaccine_peptide_group.add_argument(
        "--num-epitopes-per-vaccine-peptide",
        type=int,
        help=(
            "Maximum number of mutant epitopes to consider when scoring "
            "each vaccine peptide. (default: %(default)s)"))


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
            args.output_neoepitope_report or
            args.output_passing_variants_csv or
            args.output_isovar_csv):
        raise ValueError(
            "Must specify at least one of: --output-csv, "
            "--output-xlsx-report, "
            "--output-ascii-report, "
            "--output-html-report, "
            "--output-pdf-report, "
            "--output-neoepitope-report, "
            "--output-json-file, "
            "--output-passing-variants-csv, "
            "--output-isovar-csv")

def run_vaxrank_from_parsed_args(args):
    mhc_predictor = mhc_binding_predictor_from_args(args)

    args.protein_sequence_length = (
            args.vaccine_peptide_length + 2 * args.padding_around_mutation
    )

    # Vaxrank is going to evaluate multiple vaccine peptides containing
    # the same mutation so need a longer sequence from Isovar
    isovar_results = run_isovar_from_parsed_args(args)

    if args.output_isovar_csv:
        df = isovar_results_to_dataframe(isovar_results)
        df.to_csv(args.output_isovar_csv, index=False)

    return run_vaxrank(
        isovar_results=isovar_results,
        mhc_predictor=mhc_predictor,
        vaccine_peptide_length=args.vaccine_peptide_length,
        max_vaccine_peptides_per_variant=args.max_vaccine_peptides_per_mutation,
        min_epitope_score=args.min_epitope_score,
        num_mutant_epitopes_to_keep=args.num_epitopes_per_vaccine_peptide)

def ranked_vaccine_peptides_with_metadata_from_parsed_args(args):
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

    vaxrank_results = run_vaxrank_from_parsed_args(args)

    variants_count_dict = vaxrank_results.variant_counts()
    assert len(variants) == variants_count_dict['num_total_variants'], \
        "Len(variants) is %d but variants_count_dict came back with %d" % (
            len(variants), variants_count_dict['num_total_variants'])

    if args.output_passing_variants_csv:
        variant_metadata_dicts = vaxrank_results.variant_properties(
            gene_pathway_check=GenePathwayCheck())
        df = pd.DataFrame(variant_metadata_dicts)
        df.to_csv(args.output_passing_variants_csv, index=False)

    ranked_variants_with_vaccine_peptides = vaxrank_results.ranked_vaccine_peptides
    ranked_variants_with_vaccine_peptides_for_report = \
        ranked_variants_with_vaccine_peptides[:args.max_mutations_in_report]
    patient_info = PatientInfo(
        patient_id=args.output_patient_id,
        vcf_paths=variants.sources,
        bam_path=args.bam,
        mhc_alleles=mhc_alleles,
        num_somatic_variants=variants_count_dict['num_total_variants'],
        num_coding_effect_variants=variants_count_dict['num_coding_effect_variants'],
        num_variants_with_rna_support=variants_count_dict['num_variants_with_rna_support'],
        num_variants_with_vaccine_peptides=variants_count_dict['num_variants_with_vaccine_peptides']
    )
    # return variants, patient info, and command-line args
    data = {
        # TODO:
        #  change this field to 'ranked_variants_with_vaccine_peptides'
        #  but figure out how to do it in a backwards compatible way
        'variants': ranked_variants_with_vaccine_peptides_for_report,
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

def configure_logging(args):
    logging.config.fileConfig(
        pkg_resources.resource_filename(
            __name__,
            'logging.conf'),
        defaults={'logfilename': args.log_path})

def choose_arg_parser(args_list):
    # TODO: replace this with a command sub-parser
    if "--input-json-file" in args_list:
        return cached_run_arg_parser()
    else:
        return make_vaxrank_arg_parser()

def parse_vaxrank_args(args_list):
    arg_parser = choose_arg_parser(args_list)
    return arg_parser.parse_args(args_list)

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

    args = parse_vaxrank_args(args_list)
    configure_logging(args)
    logger.info(args)
    check_args(args)

    data = ranked_vaccine_peptides_with_metadata_from_parsed_args(args)

    ranked_variants_with_vaccine_peptides = data['variants']
    patient_info = data['patient_info']
    args_for_report = data['args']

    ###################
    # CSV-based reports
    ###################
    if args.output_csv or args.output_xlsx_report:
        make_csv_report(
            ranked_variants_with_vaccine_peptides,
            excel_report_path=args.output_xlsx_report,
            csv_report_path=args.output_csv)

    if args.output_neoepitope_report:
        make_minimal_neoepitope_report(
            ranked_variants_with_vaccine_peptides,
            num_epitopes_per_peptide=args.num_epitopes_per_vaccine_peptide,
            excel_report_path=args.output_neoepitope_report)

    ########################
    # Template-based reports
    ########################

    if not (args.output_ascii_report or args.output_html_report or args.output_pdf_report):
        return

    input_json_file = args.input_json_file if hasattr(args, 'input_json_file') else None
    template_data_creator = TemplateDataCreator(
        ranked_variants_with_vaccine_peptides=ranked_variants_with_vaccine_peptides,
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
