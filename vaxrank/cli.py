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
import sys
import logging
import logging.config
import pkg_resources

from argparse import ArgumentParser
from isovar.cli.rna_reads import allele_reads_generator_from_args
from isovar.cli.variant_sequences import make_variant_sequences_arg_parser
from mhctools.cli import (
    add_mhc_args,
    mhc_alleles_from_args,
    mhc_binding_predictor_from_args,
)
from six.moves import cPickle
from varcode.cli import variant_collection_from_args

from .core_logic import ranked_vaccine_peptides
from .report import (
    make_ascii_report,
    make_html_report,
    make_pdf_report,
    make_csv_report,
    TemplateDataCreator,
)


logging.config.fileConfig(pkg_resources.resource_filename(__name__, 'logging.conf'))
logger = logging.getLogger(__name__)


def new_run_arg_parser():
    # inherit all commandline options from Isovar
    arg_parser = make_variant_sequences_arg_parser(
        prog="vaxrank",
        description=(
            "Select personalized vaccine peptides from cancer variants, "
            "expression data, and patient HLA type."),
    )
    add_mhc_args(arg_parser)
    add_vaccine_peptide_args(arg_parser)
    add_output_args(arg_parser)
    return arg_parser


def cached_run_arg_parser():
    arg_parser = ArgumentParser(
        prog="vaxrank",
        description=(
            "Select personalized vaccine peptides from cancer variants, "
            "expression data, and patient HLA type."),
    )
    arg_parser.add_argument(
        "--input-pickle-file",
        default="",
        help="Path to pickle file containing results of vaccine peptide report")
    add_output_args(arg_parser)
    return arg_parser


def add_output_args(arg_parser):
    output_args_group = arg_parser.add_argument_group("Output options")

    output_args_group.add_argument(
        "--output-patient-id",
        default="UNKNOWN",
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
        "--output-pickle-file",
        default="",
        help="Path to pickle file containing results of vaccine peptide report")

    output_args_group.add_argument(
        "--output-csv-report-dir",
        default="",
        help="Path to CSV vaccine peptide report dir, one file per variant")

    output_args_group.add_argument(
        "--output-reviewed-by",
        default="",
        help="Comma-separated list of reviewer names")

    output_args_group.add_argument(
        "--output-final-review",
        default="",
        help="Name of final reviewer of report")


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
        "--max-mutations-in-report",
        default=10,
        type=int,
        help="Number of mutations to report")

    vaccine_peptide_group.add_argument(
        "--min-epitope-score",
        default=0.001,
        type=float,
        help="Ignore epitopes whose normalized score falls below this threshold")


def check_args(args):
    if not (args.output_csv or
            args.output_ascii_report or
            args.output_html_report or
            args.output_pdf_report or
            args.output_pickle_file or
            args.output_csv_report_dir):
        raise ValueError(
            "Must specify at least one of: --output-csv, "
            "--output-csv-report-dir, "
            "--output-ascii-report, "
            "--output-html-report, "
            "--output-pdf-report, "
            "--output-pickle-file")
    if args.output_patient_id == "UNKNOWN":
        logger.warn("Please specify --output-patient-id if possible; defaulting to unknown")


def ranked_variant_list_with_saved_args(args):
    """
    Returns a dictionary containing 2 items: a ranked variant/vaccine peptide list,
    and command-line arguments used to generate it.
    """
    if hasattr(args, 'input_pickle_file'):
        with open(args.input_pickle_file, 'rb') as f:
            data = cPickle.load(f)
            logger.info('Loaded pickle data from %s', args.input_pickle_file)
            return data

    # generator that for each variant gathers all RNA reads, both those
    # supporting the variant and reference alleles
    reads_generator = allele_reads_generator_from_args(args)
    mhc_predictor = mhc_binding_predictor_from_args(args)

    ranked_list = ranked_vaccine_peptides(
        reads_generator=reads_generator,
        mhc_predictor=mhc_predictor,
        vaccine_peptide_length=args.vaccine_peptide_length,
        padding_around_mutation=args.padding_around_mutation,
        max_vaccine_peptides_per_variant=args.max_vaccine_peptides_per_mutation,
        min_reads_supporting_cdna_sequence=args.min_reads_supporting_variant_sequence,
        min_epitope_score=args.min_epitope_score)

    ranked_list_for_report = ranked_list[:args.max_mutations_in_report]
    data = {
        'variants': ranked_list_for_report,
        'saved_args': args,
    }
    logger.info('About to save args: %s', data['saved_args'])

    # save pickle data if necessary. as of time of writing, vaxrank takes ~25 min to run,
    # most of which is core logic. the formatting is super fast, and it can
    # be useful to save the data to be able to iterate just on the formatting
    if args.output_pickle_file:
        with open(args.output_pickle_file, 'wb') as f:
            cPickle.dump(data, f, protocol=2)
            logger.info('Wrote pickle file to %s', args.output_pickle_file)

    return data


def output_data_frames(ranked_list_for_report, args):
    """
    Outputs the ranked list to dataframes.

    Parameters
    ----------
    ranked_list_for_report : list(Variant, [VaccinePeptide])
      Output of ranked_vaccine_peptides(): sorted tuples of (Variant, list(VaccinePeptide))

    args : Namespace
      Parsed args from this run, relevant for output specs.
    """
    if not (args.output_csv or args.output_csv_report_dir):
        return

    make_csv_report(ranked_list_for_report,
        report_dir_path=args.output_csv_report_dir,
        combined_report_path=args.output_csv)


def output_templates(ranked_list_for_report, args, args_for_report):
    """
    Outputs the ranked list to templated reports, if requested by user. Some templated reports
    also print the args used to generate the data, so including those here.

    Parameters
    ----------
    ranked_list_for_report : list(Variant, [VaccinePeptide])
      Output of ranked_vaccine_peptides(): sorted tuples of (Variant, list(VaccinePeptide))

    args : Namespace
      Parsed args from this run, relevant for output specs and an input pickle file if present.

    args_for_report : Namespace
      Args relevant to generating the report templates, which will be displayed in the report.
      These are either from a previous run, if the data is read from a pickle, or from the
      current run if it's generated anew.
    """
    if not (args.output_ascii_report or args.output_html_report or args.output_pdf_report):
        return

    # create a dictionary of args that we want to display in the report - essentially
    # everything except output-related args
    args_to_display_in_report = {}
    for key, value in sorted(vars(args_for_report).items()):
        if not key.startswith('output'):
            args_to_display_in_report[key] = value

    # if necessary, propagate the input pickle data path to the output; by definition,
    # it didn't get stored when the pickle file was created previously
    if hasattr(args, 'input_pickle_file'):
        args_to_display_in_report['input_pickle_file'] = args.input_pickle_file

    output_values_for_template = {
        'args': sorted(args_to_display_in_report.items()),
        'patient_id': args.output_patient_id,
        'final_review': args.output_final_review,
        'reviewers': args.output_reviewed_by.split(','),
    }
    mhc_alleles = mhc_alleles_from_args(args_for_report)
    logger.info("MHC alleles: %s", mhc_alleles)
    variants = variant_collection_from_args(args_for_report)
    logger.info("Variants: %s", variants)
    template_data_creator = TemplateDataCreator(
        ranked_variants_with_vaccine_peptides=ranked_list_for_report,
        mhc_alleles=mhc_alleles,
        variants=variants,
        bam_path=args_for_report.bam,
        output_values=output_values_for_template)
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

    if "--input-pickle-file" in args_list:
        arg_parser = cached_run_arg_parser()
    else:
        arg_parser = new_run_arg_parser()

    args = arg_parser.parse_args(args_list)
    logger.info(args)
    check_args(args)

    data = ranked_variant_list_with_saved_args(args)
    output_data_frames(
        ranked_list_for_report=data['variants'],
        args=args)
    output_templates(
        ranked_list_for_report=data['variants'],
        args=args,
        args_for_report=data['saved_args'])
