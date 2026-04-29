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



import sys
from argparse import SUPPRESS, Action, ArgumentParser, ArgumentTypeError
from importlib.resources import files

from isovar.cli import make_isovar_arg_parser
from mhctools.cli import add_mhc_args


from .epitope_config_args import add_epitope_prediction_args
from .vaccine_config_args import add_vaccine_peptide_args
from ..version import __version__


def _linker_arg(value):
    """argparse type=callable for --*-linker flags.

    Uppercases the input, then validates via vaccine_library.get_linker
    so the compositional grammar ((BASE)N / GnSm / AnY) resolves at
    parse time. Returns the canonical input string for downstream use.
    """
    from ..vaccine_library import get_linker
    upper = value.upper()
    try:
        get_linker(upper)
    except ValueError as e:
        raise ArgumentTypeError(str(e))
    return upper



class _PrintDefaultConfigAction(Action):
    """Dump the bundled default YAML config and exit.

    Useful as a starting point: ``vaxrank --print-default-config > my-config.yaml``
    produces a fully-commented config file the user can edit and then pass
    via ``--config my-config.yaml``.
    """
    def __init__(self, option_strings, dest=SUPPRESS, default=SUPPRESS,
                 help=None):
        super().__init__(
            option_strings=option_strings, dest=dest, default=default,
            nargs=0, help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        text = files("vaxrank.config").joinpath("default.yaml").read_text()
        sys.stdout.write(text)
        parser.exit()


def make_vaxrank_arg_parser():
    # create common parser with the --version flag
    parent_parser = ArgumentParser('parent', add_help=False)
    parent_parser.add_argument('--version', action='version', version='Vaxrank %s' % (__version__,))
    parent_parser.add_argument(
        '--print-default-config',
        action=_PrintDefaultConfigAction,
        help="Print the bundled default YAML config to stdout and exit. "
             "Pipe to a file (`> my-config.yaml`) to start from a fully "
             "documented config you can edit, then run with --config.")

    # inherit commandline options from Isovar
    arg_parser = make_isovar_arg_parser(
        prog="vaxrank",
        description=(
            "Rank personalized cancer neoantigens from somatic variants, "
            "tumor RNA, and patient HLA type, and emit them as analysis "
            "reports, peptide constructs, or mRNA constructs."),
        parents=[parent_parser],
    )
    
    arg_parser.add_argument(
        "--config",
        action="append",
        default=None,
        help="Path to YAML config. May be repeated; later files deep-merge "
             "over earlier ones. Each file may contain any of the top-level "
             "sections 'epitopes', 'vaccine_peptides', and 'manufacturability' "
             "— split them across files or keep them together.")
    add_config_override_args(arg_parser)

    arg_parser.add_argument(
        "--ensembl-release",
        default=None,
        type=int,
        help="Ensembl release number to use for gene annotations (e.g. 75, 93, 102). "
             "By default, pyensembl picks the most recent locally installed release "
             "for the reference assembly specified by --genome.")

    arg_parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        default=False,
        help="Print debug-level log messages to the console (by default only "
             "INFO and above are printed; DEBUG is always written to the log file).")

    arg_parser.add_argument(
        "--input-pvacseq",
        default=None,
        help="Path to a pVACseq aggregated TSV file (all_epitopes.aggregated.tsv). "
             "When provided without --vcf/--bam, vaxrank scores and ranks the "
             "pVACseq predictions and writes the requested output reports.")

    arg_parser.add_argument(
        "--input-lens",
        default=None,
        help="Path to a LENS report TSV file. "
             "When provided without --vcf/--bam, vaxrank scores and ranks the "
             "LENS predictions and writes the requested output reports.")

    add_mhc_args(arg_parser)
    add_vaccine_peptide_args(arg_parser)
    add_epitope_prediction_args(arg_parser)
    add_advanced_args(arg_parser)
    add_output_args(arg_parser)
    add_optional_output_args(arg_parser)
    add_supplemental_report_args(arg_parser)
    return arg_parser


def cached_run_arg_parser():
    arg_parser = ArgumentParser(
        prog="vaxrank",
        description=(
            "Rank personalized cancer neoantigens from somatic variants, "
            "tumor RNA, and patient HLA type, and emit them as analysis "
            "reports, peptide constructs, or mRNA constructs."),
    )
    arg_parser.add_argument(
        '--version',
        action='version',
        version='Vaxrank %s' % (__version__,))
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


def add_advanced_args(arg_parser):
    advanced_group = arg_parser.add_argument_group("Advanced options")
    advanced_group.add_argument(
        "--allow-dna-only-fallback",
        action="store_true",
        default=False,
        help="When a variant has no RNA support, attempt to construct vaccine "
             "peptides from DNA annotation alone using varcode's MutantTranscript. "
             "These peptides will have zero supporting RNA reads.")
    advanced_group.add_argument(
        "--prediction-cache",
        default=None,
        help="Path to pre-computed MHC binding predictions (topiary TSV/Parquet "
             "or CSV). Used as a lookup cache; peptides not in the cache fall "
             "through to the live MHC predictor.")
    advanced_group.add_argument(
        "--tumor-sample-name",
        default=None,
        help="Name of the tumor sample column in the VCF FORMAT fields, used "
             "to read DNA VAF (FORMAT/AF, falling back to AD). Required for "
             "multi-sample VCFs; auto-picked when only one sample is present.")


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
        "--pdf-backend",
        choices=["pdfkit", "weasyprint"],
        default="pdfkit",
        help="Library to use for PDF generation. 'weasyprint' is experimental "
             "and does not require wkhtmltopdf (default: pdfkit)")

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
        "--output-epitopes",
        default="",
        help="Path to save epitope predictions as a TSV/CSV file (format inferred "
             "from extension). This file can later be loaded with --input-epitopes.")

    output_args_group.add_argument(
        "--output-isovar-csv",
        default="",
        help="Path to CSV file containing raw RNA counts and filtering metadata "
             "for all variants (generated by Isovar)")

    add_mrna_output_args(output_args_group)
    add_peptide_output_args(output_args_group)

    output_args_group.add_argument(
        "--log-path",
        default="python.log",
        help="File path to write the vaxrank Python log to")

    output_args_group.add_argument(
        "--max-mutations-in-report",
        default=None,
        type=int,
        help="Number of mutations to report")



def add_mrna_output_args(group):
    """mRNA vaccine construct output (see vaxrank/mrna.py for assembly)."""
    from ..mrna_library import MITDS, SIGNAL_PEPTIDES, UTRS_3P, UTRS_5P
    group.add_argument(
        "--output-mrna",
        default="",
        help="Path to FASTA file of assembled mRNA vaccine constructs. "
             "When set, vaxrank emits one or more mRNA constructs containing "
             "the top vaccine peptides as concatenated antigens.")
    group.add_argument(
        "--output-mrna-manifest",
        default="",
        help="Optional path to a JSON manifest describing each mRNA construct "
             "(component names, length, contained antigens).")
    group.add_argument(
        "--mrna-signal-peptide",
        default="HLA_B",
        help="Signal peptide name from the mRNA library (one of: %s) or '' "
             "to omit. Default: HLA_B (BioNTech FixVac canonical, "
             "Sahin 2017)." % ", ".join(sorted(SIGNAL_PEPTIDES)))
    group.add_argument(
        "--mrna-linker",
        default="(G4S)2",
        type=_linker_arg,
        help="Linker name from the shared vocabulary (case-insensitive). "
             "Accepts named entries (G4S, AAY, EAAAK, P2A, ...) and "
             "compositional forms: (BASE)N / (BASE)xN / BASExN for repeats, "
             "GnSm for n-glycines + m-serines literal, AnY for n-alanines + Y. "
             "Examples: (G4S)2, G4Sx2, A3Y, G6S. "
             "Default: (G4S)2 (BioNTech FixVac canonical, Sahin 2017).")
    group.add_argument(
        "--mrna-include-mitd",
        action="store_true",
        default=True,
        help="Append the MHC-I trafficking domain (MITD) at the C-terminus to "
             "route antigens through the endolysosomal compartment. "
             "On by default to match BioNTech FixVac (Sahin 2017).")
    group.add_argument(
        "--mrna-no-mitd",
        action="store_false",
        dest="mrna_include_mitd",
        help="Disable MITD inclusion (rare for neoantigen vaccines).")
    group.add_argument(
        "--mrna-mitd",
        default="HLA_B",
        help="MITD name when --mrna-include-mitd is set (one of: %s). "
             "Default: HLA_B (BioNTech FixVac canonical)."
             % ", ".join(sorted(MITDS)))
    group.add_argument(
        "--mrna-5p-utr",
        default="HBB",
        help="5' UTR name (one of: %s). Default: HBB."
             % ", ".join(sorted(UTRS_5P)))
    group.add_argument(
        "--mrna-3p-utr",
        default="HBB_FI",
        help="3' UTR name (one of: %s). Default: HBB_FI (tandem 2× HBB / "
             "FI element, BioNTech FixVac canonical)."
             % ", ".join(sorted(UTRS_3P)))
    group.add_argument(
        "--mrna-codon-species",
        default="h_sapiens",
        help="Species key passed to DnaChisel CodonOptimize (e.g. h_sapiens, "
             "m_musculus). Default: h_sapiens.")
    group.add_argument(
        "--mrna-codon-method",
        default="use_best_codon",
        choices=["use_best_codon", "match_codon_usage", "harmonize_rca"],
        help="DnaChisel codon-optimization method. Default: use_best_codon.")
    group.add_argument(
        "--mrna-min-antigen-length-aa",
        default=15,
        type=int,
        help="Minimum amino-acid window per antigen. Default: 15.")
    group.add_argument(
        "--mrna-max-antigen-length-aa",
        default=25,
        type=int,
        help="Maximum amino-acid window per antigen. Default: 25 "
             "(BioNTech FixVac uses 27mers; 25 is the typical neoantigen "
             "vaccine SLP target).")
    group.add_argument(
        "--mrna-antigens-per-construct",
        default=5,
        type=int,
        help="Antigens packed into a single mRNA construct. Default: 5. "
             "Antigens beyond the cap spill into additional constructs.")
    group.add_argument(
        "--mrna-max-constructs",
        default=1,
        type=int,
        help="Maximum number of mRNA constructs in the vaccine. Default: 1 "
             "(one construct carrying all antigens). Higher values let "
             "additional constructs absorb spillover.")
    group.add_argument(
        "--mrna-candidates-per-slot",
        default=1,
        type=int,
        help="Alternative windows per variant slot in the mRNA construct. "
             "Default: 1. Higher emits alternative constructs with the "
             "same antigen lineup but different per-variant windows.")
    group.add_argument(
        "--mrna-max-length-nt",
        default=4000,
        type=int,
        help="Maximum total nucleotide length per construct (UTRs + CDS + stop). "
             "Antigens that would exceed this spill into additional constructs.")


def add_peptide_output_args(group):
    """Peptide vaccine construct output (see vaxrank/peptide.py)."""
    group.add_argument(
        "--output-peptide",
        default="",
        help="Path to FASTA file of peptide vaccine constructs. Default mode "
             "is one synthetic long peptide per ranked variant; see "
             "--peptide-mode.")
    group.add_argument(
        "--output-peptide-manifest",
        default="",
        help="Optional JSON manifest describing each peptide construct "
             "(matches the mRNA manifest schema for symmetric consumption).")
    group.add_argument(
        "--output-peptide-order-form",
        default="",
        help="Optional CSV order form ready to send to a peptide vendor.")
    group.add_argument(
        "--peptide-mode",
        default="slp",
        choices=["slp", "minimal_epitope", "multi_epitope"],
        help="slp (default): one synthetic long peptide per ranked variant. "
             "minimal_epitope: just the top mutant MHC ligand per variant. "
             "multi_epitope: concatenate antigens with --peptide-linker.")
    group.add_argument(
        "--peptide-linker",
        default="G4S3",
        type=_linker_arg,
        help="Linker used in --peptide-mode=multi_epitope (case-insensitive). "
             "Accepts named entries, compositional forms ((BASE)N / GnSm / "
             "AnY), and aliases. Shared with mRNA mode. Default: G4S3.")
    group.add_argument(
        "--peptide-min-antigen-length-aa",
        default=15,
        type=int,
        help="Minimum amino-acid window per antigen. Default: 15.")
    group.add_argument(
        "--peptide-max-antigen-length-aa",
        default=25,
        type=int,
        help="Maximum amino-acid window per antigen (SLP target). Default: 25.")
    group.add_argument(
        "--peptide-antigens-per-construct",
        default=1,
        type=int,
        help="Antigens per construct. Default: 1 (one antigen per peptide, "
             "the standard SLP layout). Higher values are only meaningful "
             "in --peptide-mode=multi_epitope.")
    group.add_argument(
        "--peptide-max-constructs",
        default=20,
        type=int,
        help="Maximum number of peptide constructs in the vaccine pool. "
             "Default: 20 (the typical PGV-001 personalized peptide vaccine "
             "pool size).")
    group.add_argument(
        "--peptide-candidates-per-slot",
        default=1,
        type=int,
        help="Alternative windows per variant slot. Default: 1. Set higher "
             "(e.g. 3) to emit alternates that the vaccine designer can "
             "choose between at synthesis time.")
    group.add_argument(
        "--peptide-n-terminal-acetyl",
        action="store_true",
        default=False,
        help="Mark constructs for N-terminal acetylation in the order form.")
    group.add_argument(
        "--peptide-c-terminal-amide",
        action="store_true",
        default=False,
        help="Mark constructs for C-terminal amidation in the order form.")
    group.add_argument(
        "--peptide-scale-mg",
        default=5.0,
        type=float,
        help="Synthesis scale in mg per peptide for the vendor order form. "
             "Default: 5.0.")
    group.add_argument(
        "--peptide-purity-percent",
        default=95.0,
        type=float,
        help="HPLC purity target (%%) for the vendor order form. Default: 95.0.")
    group.add_argument(
        "--peptide-counterion",
        default="TFA",
        help="Salt form / counterion for the vendor order form "
             "(e.g. TFA, acetate, HCl, free). Default: TFA.")


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
            args.output_isovar_csv or
            args.output_epitopes or
            args.output_mrna or
            args.output_peptide):
        raise ValueError(
            "Must specify at least one of: --output-csv, "
            "--output-xlsx-report, "
            "--output-ascii-report, "
            "--output-html-report, "
            "--output-pdf-report, "
            "--output-neoepitope-report, "
            "--output-json-file, "
            "--output-passing-variants-csv, "
            "--output-isovar-csv, "
            "--output-epitopes, "
            "--output-mrna, "
            "--output-peptide")


def external_input_arg_parser():
    """Parser for --input-pvacseq / --input-lens mode (no BAM/VCF required)."""
    arg_parser = ArgumentParser(
        prog="vaxrank",
        description=(
            "Score and rank pre-computed epitope predictions from "
            "pVACseq or LENS."),
    )
    arg_parser.add_argument(
        '--version',
        action='version',
        version='Vaxrank %s' % (__version__,))
    arg_parser.add_argument("--input-pvacseq", default=None)
    arg_parser.add_argument("--input-lens", default=None)
    arg_parser.add_argument(
        "--verbose", "-v", action="store_true", default=False)
    add_output_args(arg_parser)
    add_epitope_prediction_args(arg_parser)
    # config arg needed for epitope_config_from_args
    arg_parser.add_argument("--config", action="append", default=None)
    add_config_override_args(arg_parser)
    return arg_parser


def add_config_override_args(arg_parser):
    arg_parser.add_argument(
        "--config-value",
        dest="config_set_overrides",
        action=_ConfigOverrideAction,
        override_kind="set",
        default=None,
        metavar="PATH=VALUE",
        help=(
            "Override a config value using dotted.path=value. VALUE is parsed as "
            "YAML, so numbers, booleans, null, lists, and maps are supported."
        ),
    )
    arg_parser.add_argument(
        "--config-text",
        dest="config_expr_overrides",
        action=_ConfigOverrideAction,
        override_kind="expr",
        default=None,
        metavar="PATH=EXPR",
        help=(
            "Override a config value using dotted.path=text without YAML parsing. "
            "The right-hand side is stored as a raw string."
        ),
    )


def choose_arg_parser(args_list):
    # TODO: replace this with a command sub-parser
    if any(
            arg == "--input-json-file" or
            arg.startswith("--input-json-file=")
            for arg in args_list):
        return cached_run_arg_parser()
    elif any(
            arg in ("--input-pvacseq", "--input-lens") or
            arg.startswith("--input-pvacseq=") or
            arg.startswith("--input-lens=")
            for arg in args_list):
        return external_input_arg_parser()
    else:
        return make_vaxrank_arg_parser()

def parse_vaxrank_args(args_list):
    arg_parser = choose_arg_parser(args_list)
    return arg_parser.parse_args(args_list)
class _ConfigOverrideAction(Action):
    def __init__(self, option_strings, dest, **kwargs):
        self.override_kind = kwargs.pop("override_kind")
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        existing = getattr(namespace, self.dest, None)
        if existing is None:
            existing = []
        existing.append(values)
        setattr(namespace, self.dest, existing)

        ordered = getattr(namespace, "config_overrides", None)
        if ordered is None:
            ordered = []
        ordered.append((self.override_kind, values))
        setattr(namespace, "config_overrides", ordered)
