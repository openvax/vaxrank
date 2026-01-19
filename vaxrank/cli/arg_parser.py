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



from argparse import ArgumentParser

from isovar.cli import make_isovar_arg_parser
from mhctools.cli import add_mhc_args


from .epitope_config_args import add_epitope_prediction_args
from .vaccine_config_args import add_vaccine_peptide_args
from ..version import __version__



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
    
    arg_parser.add_argument(
        "--config",
        default=None,
        help="Path to YAML file with options related to epitope prediction and vaccine design.")

    add_mhc_args(arg_parser)
    add_vaccine_peptide_args(arg_parser)
    add_epitope_prediction_args(arg_parser)
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



def add_supplemental_report_args(arg_parser):
    report_args_group = arg_parser.add_argument_group("Supplemental report options")
    # Primary option uses dashes; underscore version kept for backwards compatibility
    report_args_group.add_argument(
        "--cosmic-vcf-filename",
        "--cosmic_vcf_filename",  # legacy fallback
        dest="cosmic_vcf_filename",
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


def choose_arg_parser(args_list):
    # TODO: replace this with a command sub-parser
    if "--input-json-file" in args_list:
        return cached_run_arg_parser()
    else:
        return make_vaxrank_arg_parser()

def parse_vaxrank_args(args_list):
    arg_parser = choose_arg_parser(args_list)
    return arg_parser.parse_args(args_list)
