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
        help="Path to JSON file containing the cached vaxrank results "
             "(same shape as --output-json-file). Used to re-render "
             "reports without re-running the full pipeline.")
    add_output_args(arg_parser)
    add_optional_output_args(arg_parser)
    add_supplemental_report_args(arg_parser)
    return arg_parser



# Lets the user specify whether they want to see particular sections in the report.
def add_optional_output_args(arg_parser):
    # Manufacturability defaults to None so we can resolve it later
    # against ``--vaccine-type``: peptide-mode runs include the
    # GRAVY / Cys-content / N-terminal-Q manufacturability section
    # by default (it's relevant to peptide synthesis); mRNA-mode
    # runs exclude it (those features don't apply to mRNA
    # constructs). Users can still force either way explicitly.
    manufacturability_args = arg_parser.add_mutually_exclusive_group(required=False)
    manufacturability_args.add_argument(
        "--include-manufacturability-in-report",
        dest="manufacturability",
        action="store_true",
        help="Force-include the peptide manufacturability section "
             "(GRAVY scores, Cys content, N-terminal Q/E/C, etc.) "
             "in template reports. Default depends on --vaccine-type: "
             "on for peptide, off for mrna.")

    manufacturability_args.add_argument(
        "--no-manufacturability-in-report",
        dest="manufacturability",
        action="store_false",
        help="Force-exclude the peptide manufacturability section "
             "from template reports. See "
             "--include-manufacturability-in-report for the default.")
    arg_parser.set_defaults(manufacturability=None)

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
        nargs='+',
        default=["peptide", "mrna"],
        choices=["peptide", "mrna"],
        metavar="TYPE",
        help="Which vaccine type(s) to design. One or more of "
             "{peptide, mrna}. Default: 'peptide mrna' (both); pass "
             "a subset to design just one. Single-mode runs (one "
             "type) write directly into --output-dir. Multi-mode "
             "runs (≥2) write into per-modality subdirs "
             "(--output-dir/peptide/, --output-dir/mrna/, …) so the "
             "same flag set works regardless of how many modalities "
             "are active. Future modalities (DNA, …) plug in here.")

    # Unified output directory. ``--output-dir`` is always a
    # directory; vaxrank chooses canonical filenames inside it.
    # Single-mode: files go directly in DIR. Multi-mode: per-type
    # subdirs (DIR/peptide/, DIR/mrna/, …).
    #
    # Canonical filenames per type:
    #   peptide/  vaccine.fasta, manifest.json, order_form.csv
    #   mrna/     cds.fasta, no_polyA.fasta, full.fasta,
    #             manifest.json, mrna-sequence-parts.csv
    #
    # No file-extension-based mode switching: ``--output-dir`` is
    # never a file path. Tradeoff: single-file output ergonomics
    # gone; gain: one rule for every modality count and a clean
    # directory layout for hand-off.
    output_args_group.add_argument(
        "--output-dir",
        default="",
        metavar="DIR",
        help="Destination directory for the assembled vaccine "
             "constructs. Single-mode runs (one --vaccine-type) write "
             "files directly into DIR; multi-mode runs write to "
             "per-modality subdirs (DIR/peptide/, DIR/mrna/, …). "
             "Required when designing a vaccine; omit for "
             "ranking-only / report-only runs. Vaxrank picks "
             "canonical filenames inside (vaccine.fasta / cds.fasta / "
             "manifest.json / order_form.csv / mrna-sequence-parts.csv).")
    output_args_group.add_argument(
        "--mrna-csv-no-full-rows",
        dest="mrna_csv_full_rows",
        action="store_false",
        default=True,
        help="mRNA-only. Suppress the per-construct cds / no_polyA / "
             "full summary rows in the mRNA mrna-sequence-parts.csv (the rows "
             "with the longest nt cells). Per-element rows are "
             "unaffected.")

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
        help="Path to CSV with ranked antigens. Format depends on "
             "input source: pipeline (--vcf + --bam) emits "
             "per-variant rows (gene, mutation, antigen sequence, "
             "scores); LENS / pVACseq emit per-(peptide, allele) "
             "rows (the upstream tool's native granularity). "
             "Antigen-centric in both cases; content does not "
             "change with --vaccine-type.")

    output_args_group.add_argument(
        "--output-ascii-report",
        default="",
        help="Path to ASCII summary report — ranked antigens with "
             "their predicted MHC ligands, manufacturability metrics "
             "(peptide-mode only by default), and processing-credibility "
             "scores. Antigen-centric: content adapts to --vaccine-type "
             "but the same flag works for all modes.")

    output_args_group.add_argument(
        "--output-html-report",
        default="",
        help="Path to HTML summary report (same content as "
             "--output-ascii-report, rendered to HTML).")

    output_args_group.add_argument(
        "--output-pdf-report",
        default="",
        help="Path to PDF summary report (same content as "
             "--output-ascii-report, rendered to PDF).")

    output_args_group.add_argument(
        "--pdf-backend",
        choices=["pdfkit", "weasyprint"],
        default="pdfkit",
        help="Library to use for PDF generation. 'weasyprint' is experimental "
             "and does not require wkhtmltopdf (default: pdfkit)")

    output_args_group.add_argument(
        "--output-json-file",
        default="",
        help="Path to JSON dump of the ranked antigen set + "
             "command-line state. Same structure regardless of "
             "--vaccine-type.")

    output_args_group.add_argument(
        "--output-xlsx-report",
        default="",
        help="Path to XLSX summary report (one sheet per variant). "
             "Antigen-centric; intended for hand-off to a vaccine "
             "manufacturer / reviewer regardless of modality.")

    output_args_group.add_argument(
        "--output-neoepitope-report",
        default="",
        help="Path to per-(peptide, allele) neoepitope XLSX. "
             "Focused on the predicted minimal MHC ligands; available "
             "from any --vaccine-type and from the LENS / pVACseq "
             "external paths.")

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

    Output destinations are unified under ``--output-dir`` (see
    ``add_output_args``); the mRNA writer picks canonical filenames
    inside it (cds.fasta / no_polyA.fasta / full.fasta /
    manifest.json / mrna-sequence-parts.csv). Only mrna-specific design
    parameters live here.
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
        # ``G4Sx2`` is self-documenting; ``(G4S)2`` parses to the
        # same string but the parens require shell escaping. Avoid
        # the ambiguous bare ``GS2`` shorthand — see --peptide-linker.
        default="G4Sx2",
        type=_linker_arg,
        help="Linker name from the shared vocabulary (case-insensitive). "
             "Accepts named entries (G4S, AAY, EAAAK, P2A, ...) and "
             "compositional forms: (BASE)N / (BASE)xN / BASExN for repeats, "
             "GnSm for n-glycines + m-serines literal, AnY for n-alanines + Y. "
             "Examples: G4Sx2, (G4S)2, A3Y, G6S. "
             "Default: G4Sx2 = (G4S)2 = GGGGSGGGGS (BioNTech FixVac "
             "canonical, Sahin 2017).")
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

    Output destinations are unified under ``--output-dir`` (see
    ``add_output_args``); only peptide-specific design parameters
    live here.
    """
    group.add_argument(
        "--peptide-mode",
        default="slp",
        # ``type=str.lower`` so SLP / Slp / slp all reach the choice
        # comparison as the lowercase canonical form.
        type=str.lower,
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
        # ``G4Sx3`` is self-documenting: "G4S repeated 3x" =
        # GGGGSGGGGSGGGGS. The shorter alias ``GS3`` and the bare
        # form ``G4S3`` both also parse, but G4S3 means literally
        # ``GGGGSSS`` (4 G + 3 S, GnSm form) and GS3 reads as
        # ambiguous "G + S + 3"; G4Sx3 leaves no doubt about which
        # was meant.
        default="G4Sx3",
        type=_linker_arg,
        help="Linker used in --peptide-mode=multi_epitope (case-insensitive). "
             "Accepts named entries (P2A, AAY, EAAAK, …), compositional "
             "forms ((BASE)N or BASExN for repeats, GnSm for n-glycines + "
             "m-serines literal, AnY for n-alanines + Y), and aliases. "
             "Shared with mRNA mode. Default: G4Sx3 = (G4S)3 = "
             "GGGGSGGGGSGGGGS (15 aa, BioNTech FixVac canonical).")
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
    ("--output-dir",                    "output_dir",
        "directory for assembled vaccine constructs "
        "(per-modality subdirs when --vaccine-type has 2+ entries)"),
)


def check_args(args):
    """Fail fast when no primary --output-* flag is set, or when
    type-specific companion flags don't match --vaccine-type.

    Every output path is opt-in via its own flag; there is no default
    destination. Without this guard a user can run a long pipeline (or
    a LENS / pVACseq import) and end up with nothing on disk — exactly
    the surprise we want to prevent.
    """
    if not any(getattr(args, attr, '') for _, attr, _ in _PRIMARY_OUTPUT_FLAGS):
        flag_lines = "\n".join(
            "  %-32s %s" % (flag, purpose)
            for flag, _, purpose in _PRIMARY_OUTPUT_FLAGS)
        raise ValueError(
            "No output path specified. Pass at least one of the "
            "following flags so vaxrank knows where to write "
            "results:\n%s" % flag_lines)

    # ``--output-dir`` is contractually a directory; the file-path
    # form was removed when single-file outputs went away. Reject
    # paths that point at an existing file or look like a FASTA
    # (``out.fasta``) so the user sees the mismatch immediately
    # rather than burning ranking time and crashing in the writer.
    _reject_output_dir_misuse(args)

    # External input + template-report flag without ``--ensembl-release``
    # would render reports with empty effect annotations (transcript
    # IDs from LENS / pVACseq won't resolve). Reject pre-flight with a
    # build hint inferred from the file — empty-effect reports waste
    # the user's time and silently degrade quality, so make this a
    # hard error rather than a warning.
    _require_ensembl_release_for_template_reports(args)


# Suffixes that look like a single-file output but ``--output-dir``
# always wants a directory; same set ``mrna.py`` enforces locally for
# its own writer, hoisted here so both writers get the check.
#
# Scoped narrowly to FASTA-family extensions and archive containers:
# those are the formats this tree actually emits or that a user is
# reasonably likely to mistake for a vaccine output. Tabular suffixes
# (``.json`` / ``.csv`` / ``.tsv`` / ``.txt``) are deliberately *not*
# included — those happen to be common directory names too (e.g.
# ``runs/csv/`` for batch job outputs), and rejecting them
# unconditionally trips users who picked a sensible directory name.
# Mistakes in the tabular space surface as soon as the writer tries
# to ``os.makedirs`` an existing file (caught by the ``isfile`` check
# above).
_OUTPUT_DIR_REJECTED_SUFFIXES = (
    '.fasta', '.fa', '.fna', '.ffn', '.fas',
    '.zip', '.tar')


def _reject_output_dir_misuse(args):
    """Catch the "user passed a file path to --output-dir" case
    early. Three rejection categories:

    1. The path already exists and is a file (vaxrank would
       ``os.makedirs(path)`` and crash deep in the writer).
    2. The path doesn't exist but ends in a suffix that strongly
       suggests a single-file target — ``out.fasta``, ``out.zip``,
       etc. (See ``_OUTPUT_DIR_REJECTED_SUFFIXES`` for the exact
       list — kept narrow so common directory names like
       ``runs/csv/`` aren't rejected.) Better to refuse than create
       a directory called ``out.fasta/``.
    3. (Multi-mode-only) the path resolves to a sub-path of an
       existing modality subdir — handled by the writers.
    """
    import os
    output_dir = getattr(args, 'output_dir', '') or ''
    if not output_dir:
        return
    if os.path.isfile(output_dir):
        raise ValueError(
            "--output-dir points at an existing file (%r); the flag "
            "wants a directory. Pass a directory path or a path that "
            "doesn't exist yet." % output_dir)
    if any(output_dir.lower().endswith(s)
           for s in _OUTPUT_DIR_REJECTED_SUFFIXES):
        raise ValueError(
            "--output-dir is %r, which looks like a single-file "
            "target. Vaxrank picks canonical filenames inside the "
            "directory (vaccine.fasta / cds.fasta / manifest.json / "
            "etc.); pass a directory path instead (e.g. drop the "
            "suffix)." % output_dir)


def _require_ensembl_release_for_template_reports(args):
    needs_transcripts = any(
        getattr(args, attr, '') for attr in (
            'output_ascii_report',
            'output_html_report',
            'output_pdf_report',
        )
    )
    if not needs_transcripts:
        return
    if getattr(args, 'ensembl_release', None) is not None:
        return
    lens_path = getattr(args, 'input_lens', None)
    pvacseq_path = getattr(args, 'input_pvacseq', None)
    if not (lens_path or pvacseq_path):
        return
    # Try to infer the build from the LENS file so the hint can name
    # a plausible release. pVACseq aggregates don't carry a build
    # marker we can read here. Be honest about uncertainty:
    # ``origin_descriptor`` only pins the build (GRCh37 vs GRCh38),
    # not the Ensembl release — a build spans many releases. We
    # suggest something the user *already has installed* via the
    # pyensembl cache, falling back to the canonical build-final
    # release (GRCh37 → 75) when nothing is installed.
    hint = " Pass --ensembl-release N to populate them."
    if lens_path:
        from ..external_input import (
            infer_genome_build_from_lens,
            installed_ensembl_releases_for_build,
            _BUILD_TO_CANONICAL_RELEASE)
        try:
            build = infer_genome_build_from_lens(lens_path)
        except Exception:
            build = None
        if build:
            installed = installed_ensembl_releases_for_build(build)
            if installed:
                # Suggest the latest installed release for this build
                # — it's a real candidate the user can use without
                # downloading anything. Show the full list so they
                # know what else is available.
                hint = (
                    " The LENS file's origin_descriptor looks like "
                    "%s; you have these %s releases installed: %s. "
                    "Try --ensembl-release %d (or whichever matches "
                    "the build LENS used)."
                    % (build, build, installed, installed[-1]))
            elif build in _BUILD_TO_CANONICAL_RELEASE:
                # GRCh37 has a canonical "last release" (75); GRCh38
                # doesn't — we don't pick a number out of thin air.
                rel = _BUILD_TO_CANONICAL_RELEASE[build]
                hint = (
                    " The LENS file's origin_descriptor looks like "
                    "%s; the canonical final %s release is %d. "
                    "Install with `pyensembl install --release %d "
                    "--reference-name %s`."
                    % (build, build, rel, rel, build))
            else:
                hint = (
                    " The LENS file's origin_descriptor looks like "
                    "%s, but no %s releases are installed in your "
                    "pyensembl cache. Install one with `pyensembl "
                    "install --release N --reference-name %s` "
                    "(N must match the build LENS used)."
                    % (build, build, build))
    raise ValueError(
        "Template report(s) requested with --input-lens / "
        "--input-pvacseq but no --ensembl-release set. Transcript "
        "effect annotations would be empty (rendered as '—'), so "
        "vaxrank refuses to write a degraded report.%s" % hint)


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
    # ``--ensembl-release`` lets external-input runs resolve LENS /
    # pVACseq transcript_id strings to pyensembl ``Transcript`` objects
    # so template reports (ASCII / HTML / PDF) can render variant
    # effects. Without it, those reports fall back to a less detailed
    # rendering (no transcript-name / effect-description columns).
    arg_parser.add_argument(
        "--ensembl-release",
        default=None,
        type=int,
        help="Ensembl release for resolving transcript IDs (e.g. 75, 102). "
             "By default pyensembl picks the most recent locally installed "
             "release. Match the release LENS / pVACseq used.")
    add_output_args(arg_parser)
    # Template reports (ASCII / HTML / PDF) work on this path too once
    # transcripts are resolved (#268). The supplemental + optional groups
    # contribute the namespace shape ``TemplateDataCreator`` expects
    # (manufacturability / wt_epitopes / cosmic_vcf_filename).
    add_optional_output_args(arg_parser)
    add_supplemental_report_args(arg_parser)
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
    parsed = arg_parser.parse_args(args_list)
    # Snapshot of parser defaults so downstream code can detect
    # "user explicitly passed this flag" by ``args.X !=
    # args._parser_defaults['X']``. Used by the construct-kwargs
    # resolver to give CLI flags precedence over YAML config values.
    #
    # Apply ``action.type`` to string defaults the same way argparse
    # does at parse time — otherwise the snapshot stores the raw
    # default (e.g. ``'G4Sx3'``) while ``args.peptide_linker`` holds
    # the type-converted form (``'G4SX3'`` after _linker_arg
    # uppercases), and the equality check falsely flags the
    # unmodified default as user-supplied.
    defaults = {}
    for a in arg_parser._actions:
        if not a.dest or a.dest == 'help':
            continue
        val = a.default
        if isinstance(val, str) and a.type and callable(a.type):
            try:
                val = a.type(val)
            except Exception:
                pass
        defaults[a.dest] = val
    parsed._parser_defaults = defaults
    return parsed
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
