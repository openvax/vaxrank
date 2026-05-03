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
        "--vaccine-type",
        default="peptide",
        choices=["peptide", "mrna"],
        metavar="TYPE",
        help="Which vaccine type to design. One mode per run. The "
             "vaccine output destination goes in --vaccine-output; "
             "type-specific companion files in --vaccine-manifest, "
             "--vaccine-order-form (peptide), --vaccine-csv (mrna). "
             "Default: peptide.")

    # Unified vaccine output (replaces the old --output-peptide /
    # --output-mrna pair). The destination is interpreted by the
    # active --vaccine-type:
    #   peptide → FASTA file at PATH
    #   mrna    → directory at PATH (writes cds.fasta / no_polyA.fasta
    #             / full.fasta into it)
    output_args_group.add_argument(
        "--vaccine-output",
        default="",
        metavar="PATH",
        help="Destination for the assembled vaccine constructs. "
             "For --vaccine-type=peptide this is a FASTA file path; "
             "for --vaccine-type=mrna it is a directory that vaxrank "
             "writes cds.fasta, no_polyA.fasta, and full.fasta into. "
             "Required when designing a vaccine; omit for "
             "ranking-only / report-only runs.")
    output_args_group.add_argument(
        "--vaccine-manifest",
        default="",
        metavar="PATH",
        help="Optional JSON manifest path describing the assembled "
             "vaccine. Same schema for peptide and mRNA so downstream "
             "tools can consume either uniformly.")
    output_args_group.add_argument(
        "--vaccine-order-form",
        default="",
        metavar="PATH",
        help="Peptide-only. CSV order form ready to send to a peptide "
             "vendor. Errors fast if --vaccine-type is not 'peptide'.")
    output_args_group.add_argument(
        "--vaccine-csv",
        default="",
        metavar="PATH",
        help="mRNA-only. Long-format CSV: one row per (construct, "
             "element) with construct / element_kind / index / "
             "index_label / name / aa / nt / length_aa / length_nt / "
             "note. Errors fast if --vaccine-type is not 'mrna'.")
    output_args_group.add_argument(
        "--vaccine-csv-no-full-rows",
        dest="vaccine_csv_full_rows",
        action="store_false",
        default=True,
        help="mRNA-only. Suppress the per-construct cds / no_polyA / "
             "full summary rows in --vaccine-csv (the rows with the "
             "longest nt cells). Per-element rows are unaffected.")

    # Shared antigen-design axes — apply to BOTH peptide and mRNA
    # unless overridden by the per-type flag (--peptide-* / --mrna-*).
    # The two orthogonal axes:
    #   antigen_content: what each antigen *is*
    #   antigens_per_construct: how many to concatenate per construct
    output_args_group.add_argument(
        "--antigen-content",
        default=None,
        choices=["mutation_spanning", "minimal_epitope"],
        help="Default antigen content for both peptide and mRNA "
             "vaccines: 'mutation_spanning' (centered SLP-style window "
             "of --max-antigen-length-aa, default) or 'minimal_epitope' "
             "(top mutant MHC ligand(s), ~9 aa each). Per-type flags "
             "(--peptide-antigen-content / --mrna-antigen-content) "
             "override.")
    output_args_group.add_argument(
        "--epitopes-per-antigen",
        default=None,
        type=int,
        help="When --antigen-content=minimal_epitope, take this many "
             "top score-sorted ligands per ranked vaccine peptide as "
             "separate antigens. Default: 1 (the legacy single-top "
             "ligand semantics). Set >1 to pack multiple top ligands "
             "from the same variant.")

    # Issue #249: pepsickle proteasome-cleavage credibility annotations.
    # On by default; doesn't change vaccine ranking — purely surfaces
    # cleavage credibility on each predicted MHC ligand for clinical
    # review.
    output_args_group.add_argument(
        "--processing-aware-annotation",
        dest="processing_aware_annotation",
        action="store_true",
        default=True,
        help="Annotate each predicted MHC ligand with pepsickle's "
             "proteasome-cleavage credibility scores "
             "(pepsickle_c_term_cleavage_prob, "
             "pepsickle_max_internal_cut_prob, "
             "pepsickle_processing_score). On by default. "
             "Pepsickle runs in an isolated subprocess (issue #266) so "
             "torch's libomp doesn't clash with the parent's pandas / "
             "numpy / pyarrow OpenMP runtime. Doesn't change vaccine "
             "ranking — adds info for review.")
    output_args_group.add_argument(
        "--no-processing-aware-annotation",
        dest="processing_aware_annotation",
        action="store_false",
        help="Disable pepsickle credibility annotation (e.g. when "
             "pepsickle isn't installed or you want a faster run).")
    output_args_group.add_argument(
        "--pepsickle-human-only",
        action="store_true",
        default=False,
        help="Use pepsickle's human-only-trained model (default: "
             "off, all-mammal). Pass when running on human samples "
             "exclusively.")
    output_args_group.add_argument(
        "--pepsickle-threshold",
        type=float,
        default=0.5,
        help="Cleavage probability threshold used internally by "
             "pepsickle (default 0.5). Doesn't affect the continuous "
             "scores attached to each ligand, only pepsickle's "
             "internal binarization.")

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
    """mRNA-vaccine design knobs (see vaxrank/mrna.py for assembly).

    Output destinations are unified under ``--vaccine-output`` /
    ``--vaccine-manifest`` / ``--vaccine-csv`` (see
    ``add_output_args``); only mrna-specific design parameters live
    here.
    """
    from ..mrna_library import MITDS, SIGNAL_PEPTIDES, UTRS_3P, UTRS_5P
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
        "--mrna-poly-a-length",
        default=120,
        type=int,
        help="PolyA tail length in nt. Default: 120 (BioNTech FixVac / "
             "BNT122 per Sahin 2017). Set 0 to omit the polyA tail "
             "(full.fasta will then equal no_polyA.fasta).")
    group.add_argument(
        "--mrna-poly-a-segmented",
        action="store_true",
        default=False,
        help="Split the polyA into two segments separated by a 10-nt "
             "linker (BNT162b2 architecture: A30 + GCATATGACT + A70 per "
             "Xia 2021, PMC8310186). Increases mRNA stability by "
             "reducing template recombination during in vitro "
             "transcription.")
    group.add_argument(
        "--mrna-poly-a-first-segment",
        default=30,
        type=int,
        help="Length (nt) of the first polyA segment when "
             "--mrna-poly-a-segmented is set. Default: 30.")
    group.add_argument(
        "--mrna-poly-a-segment-linker",
        default="GCATATGACT",
        help="Linker sequence between segmented polyA segments. "
             "Default: GCATATGACT (BNT162b2 canonical).")
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
        "--mrna-optimize-linkers",
        dest="mrna_optimize_linkers",
        action="store_true",
        default=True,
        help="Enable per-junction linker optimization to minimize predicted "
             "MHC presentation of chimeric k-mers spanning antigen junctions "
             "(issue #247). On by default; requires --mhc-predictor + "
             "--mhc-alleles. The default linker is used as a fallback if no "
             "candidate outperforms it.")
    group.add_argument(
        "--mrna-no-optimize-linkers",
        dest="mrna_optimize_linkers",
        action="store_false",
        help="Disable per-junction linker optimization. The shared linker "
             "is used at every junction without optimization.")
    group.add_argument(
        "--mrna-junction-candidates",
        default="",
        help="Comma-separated linker names to try at each junction "
             "(e.g. 'G3S,G4S,(G3S)2,(G4S)2,AAA'). Empty = use the "
             "library default JUNCTION_SWAP_CANDIDATES.")
    group.add_argument(
        "--mrna-junction-rank-strong",
        default=0.5,
        type=float,
        help="Presentation rank below which a chimeric k-mer counts as "
             "a strong-binder hit (primary minimization target). "
             "Default: 0.5%%.")
    group.add_argument(
        "--mrna-junction-rank-mild",
        default=2.0,
        type=float,
        help="Presentation rank below which a chimeric k-mer counts as "
             "a mild-binder hit (secondary minimization target). "
             "Default: 2.0%%.")
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
        "--mrna-antigen-content",
        default=None,
        choices=["mutation_spanning", "minimal_epitope"],
        help="Override --antigen-content for mRNA. 'mutation_spanning' "
             "(default if --antigen-content unset) emits 25-aa SLP-style "
             "windows; 'minimal_epitope' emits top MHC ligands as "
             "antigens — concatenated minimal-epitope mRNA "
             "(\"string of beads\") is a first-class design.")
    group.add_argument(
        "--mrna-epitopes-per-antigen",
        default=None,
        type=int,
        help="Override --epitopes-per-antigen for mRNA. Only used when "
             "antigen content is 'minimal_epitope'.")
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
    """Peptide-vaccine design knobs (see vaxrank/peptide.py).

    Output destinations are unified under ``--vaccine-output`` /
    ``--vaccine-manifest`` / ``--vaccine-order-form`` (see
    ``add_output_args``); only peptide-specific design parameters
    live here.
    """
    group.add_argument(
        "--peptide-mode",
        default="slp",
        choices=["slp", "minimal_epitope", "multi_epitope"],
        help="Back-compat shorthand for the orthogonal axes "
             "(--antigen-content + --peptide-antigens-per-construct). "
             "slp (default): mutation_spanning, 1 per construct. "
             "minimal_epitope: minimal_epitope, 1 per construct. "
             "multi_epitope: mutation_spanning, N per construct (set N "
             "via --peptide-antigens-per-construct). Prefer the new "
             "axis flags for new designs.")
    group.add_argument(
        "--peptide-antigen-content",
        default=None,
        choices=["mutation_spanning", "minimal_epitope"],
        help="Override --antigen-content for peptide vaccines. When "
             "set, overrides whatever --peptide-mode would derive.")
    group.add_argument(
        "--peptide-epitopes-per-antigen",
        default=None,
        type=int,
        help="Override --epitopes-per-antigen for peptide vaccines.")
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


_PRIMARY_OUTPUT_FLAGS = (
    # (CLI flag, args attribute, one-line purpose)
    ("--output-csv",                    "output_csv",
        "ranked vaccine peptides as CSV"),
    ("--output-xlsx-report",            "output_xlsx_report",
        "ranked vaccine peptides as XLSX"),
    ("--output-ascii-report",           "output_ascii_report",
        "ASCII summary report"),
    ("--output-html-report",            "output_html_report",
        "HTML summary report"),
    ("--output-pdf-report",             "output_pdf_report",
        "PDF summary report"),
    ("--output-neoepitope-report",      "output_neoepitope_report",
        "per-(peptide, allele) neoepitope XLSX"),
    ("--output-json-file",              "output_json_file",
        "machine-readable JSON of ranked vaccine peptides"),
    ("--output-passing-variants-csv",   "output_passing_variants_csv",
        "filter-passing variants as CSV"),
    ("--output-isovar-csv",             "output_isovar_csv",
        "Isovar transcript-assembly intermediate CSV"),
    ("--output-epitopes",               "output_epitopes",
        "raw epitope predictions"),
    ("--vaccine-output",                "vaccine_output",
        "assembled vaccine constructs "
        "(FASTA for --vaccine-type=peptide, directory for mrna)"),
)


def check_args(args):
    """Fail fast when no primary --output-* flag is set, or when the
    user pairs LENS / pVACseq input with a pipeline-only output.

    Every output path is opt-in via its own flag; there is no default
    destination. Without this guard a user can run a long pipeline (or
    a LENS / pVACseq import) and end up with nothing on disk — exactly
    the surprise we want to prevent.

    Template reports (ASCII / HTML / PDF) need pyensembl ``Transcript``
    objects and the variant-counting metadata the pipeline path
    computes from VCF + BAM. LENS / pVACseq inputs carry transcript
    IDs as strings and no variant counts, so requesting those reports
    on the external path used to silently no-op. Reject the
    combination here with an explicit error pointing to the outputs
    that *are* reachable from external input. Tracked in
    openvax/vaxrank#268.
    """
    if not any(getattr(args, attr, '') for _, attr, _ in _PRIMARY_OUTPUT_FLAGS):
        flag_lines = "\n".join(
            "  %-32s %s" % (flag, purpose)
            for flag, _, purpose in _PRIMARY_OUTPUT_FLAGS)
        raise ValueError(
            "No output path specified. Pass at least one of the "
            "following --output-* flags so vaxrank knows where to "
            "write results:\n%s" % flag_lines)

    is_external = bool(
        getattr(args, 'input_lens', None) or
        getattr(args, 'input_pvacseq', None))
    if is_external:
        pipeline_only = [
            (flag, attr) for flag, attr in (
                ("--output-ascii-report", "output_ascii_report"),
                ("--output-html-report",  "output_html_report"),
                ("--output-pdf-report",   "output_pdf_report"),
            )
            if getattr(args, attr, '')]
        if pipeline_only:
            verb = "is" if len(pipeline_only) == 1 else "are"
            raise ValueError(
                "%s %s not reachable from --input-lens / "
                "--input-pvacseq input — those reports need full "
                "transcript objects and variant-call metadata that "
                "only the --vcf + --bam pipeline produces. Use "
                "--output-csv / --output-xlsx-report / "
                "--output-neoepitope-report / --output-json-file "
                "instead."
                % (", ".join(flag for flag, _ in pipeline_only), verb))

    # Reject companion flags meant for the other vaccine type — fail
    # fast rather than silently ignoring the file the user expected
    # to see written.
    vaccine_type = getattr(args, 'vaccine_type', 'peptide')
    if vaccine_type == 'mrna' and getattr(args, 'vaccine_order_form', ''):
        raise ValueError(
            "--vaccine-order-form is peptide-only (vendor synthesis "
            "form); --vaccine-type is 'mrna'. Drop --vaccine-order-form "
            "or switch --vaccine-type to 'peptide'.")
    if vaccine_type == 'peptide' and getattr(args, 'vaccine_csv', ''):
        raise ValueError(
            "--vaccine-csv is mrna-only (per-element layer table); "
            "--vaccine-type is 'peptide'. Drop --vaccine-csv or "
            "switch --vaccine-type to 'mrna'.")


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
