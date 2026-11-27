DEFAULT_LOGISTIC_EPITOPE_SCORE_MIDPOINT = 350.0
DEFAULT_LOGISTIC_EPITOPE_SCORE_WIDTH = 150.0
DEFAULT_MIN_EPITOPE_SCORE = 0.00001
DEFAULT_BINDING_AFFINITY_CUTOFF = 5000.0
DEFAULT_PERCENTILE_RANK_CUTOFF = 10.0

DEFAULT_PREFERRED_PEPTIDE_LENGTH = 25
DEFAULT_MIN_PEPTIDE_LENGTH = 25
DEFAULT_MAX_PEPTIDE_LENGTH = 25
DEFAULT_PADDING_AROUND_MUTATION = 5
DEFAULT_MAX_VACCINE_PEPTIDES_PER_VARIANT = 1
DEFAULT_NUM_MUTANT_EPITOPES_TO_KEEP = 1000
DEFAULT_VACCINE_PEPTIDE_SCORE_FRACTION_OF_BEST = 0.99

DEFAULT_MANUFACTURABILITY_MAX_C_TERMINAL_HYDROPATHY = 1.5
DEFAULT_MANUFACTURABILITY_MIN_KMER_HYDROPATHY = 0.0
DEFAULT_MANUFACTURABILITY_MAX_KMER_HYDROPATHY_LOW_PRIORITY = 1.5
DEFAULT_MANUFACTURABILITY_MAX_KMER_HYDROPATHY_HIGH_PRIORITY = 2.5

DEFAULT_MIN_KMER_LENGTH = 8
DEFAULT_MAX_KMER_LENGTH = 15

# Pure-stdlib CLI vocabularies shared between parser construction and domain modules
# so building argument parsers never triggers heavy runtime imports (#794).
ANTIGEN_KIND_MUTATION = "mutation"
ANTIGEN_KIND_FUSION = "fusion"
ANTIGEN_KIND_CTA = "CTA"
ANTIGEN_KIND_ERV = "ERV"
ANTIGEN_KIND_SPLICE = "splice"
ANTIGEN_KIND_VIRAL = "viral"
ANTIGEN_KIND_SV = "SV"
ANTIGEN_KIND_ISOFORM = "isoform"
ANTIGEN_KINDS = frozenset({
    ANTIGEN_KIND_MUTATION,
    ANTIGEN_KIND_FUSION,
    ANTIGEN_KIND_CTA,
    ANTIGEN_KIND_ERV,
    ANTIGEN_KIND_SPLICE,
    ANTIGEN_KIND_VIRAL,
    ANTIGEN_KIND_SV,
    ANTIGEN_KIND_ISOFORM,
})

ANTIGEN_SELECTION_TYPES = ANTIGEN_KINDS | {"ncORF"}
ANTIGEN_SELECTION_CHOICES = ANTIGEN_SELECTION_TYPES | {"variant"}

VARIANT_TYPES = frozenset({"SNV", "MNV", "INDEL", "SV"})
SV_TYPES = frozenset({"DEL", "DUP", "INV", "INS", "BND", "CNV"})

LENS_ANTIGEN_SOURCES = frozenset({
    "SNV", "INDEL", "FUSION", "SPLICE", "CTA/SELF", "ERV",
})

SHARED_SEQUENCE_RNA_POLICY_SOURCE = "source"
SHARED_SEQUENCE_RNA_POLICY_UNION = "within_library_union"
SHARED_SEQUENCE_RNA_POLICIES = frozenset({
    SHARED_SEQUENCE_RNA_POLICY_SOURCE,
    SHARED_SEQUENCE_RNA_POLICY_UNION,
})


def normalize_selection_type(value):
    """Canonicalize an antigen-selection type name case-insensitively."""
    return next(
        (kind for kind in ANTIGEN_SELECTION_CHOICES if kind.lower() == value.lower()),
        value,
    )


PRIMARY_OUTPUT_FLAGS = (
    # (CLI flag, args attribute, one-line purpose)
    ("--output-csv",                    "output_csv",
        "ranked candidate windows as CSV"),
    ("--output-tumor-board-dossier",    "output_tumor_board_dossier",
        "AI-ready tumor_board_dossier.json bundle"),
    ("--output-xlsx-report",            "output_xlsx_report",
        "ranked candidate windows as XLSX"),
    ("--output-ascii-report",           "output_ascii_report",
        "plain-text product and candidate report"),
    ("--output-html-report",            "output_html_report",
        "HTML summary report"),
    ("--output-pdf-report",             "output_pdf_report",
        "PDF summary report"),
    ("--output-neoepitope-report",      "output_neoepitope_report",
        "per-(peptide, allele) prediction table as CSV or XLSX"),
    ("--output-json-file",              "output_json_file",
        "legacy ranked-candidate JSON"),
    ("--output-passing-variants-csv",   "output_passing_variants_csv",
        "filter-passing variants as CSV"),
    ("--output-isovar-csv",             "output_isovar_csv",
        "Isovar transcript-assembly intermediate CSV"),
    ("--output-epitopes",               "output_epitopes",
        "native prediction dataset and replay evidence"),
    ("--output-antigens",               "output_antigens",
        "experimental resolved antigen sequences and admission evidence"),
    ("--output-input-predictions",      "output_input_predictions",
        "original external-table predictions"),
    ("--output-policy-dir",             "output_policy_dir",
        "resolved policy configuration and decision tables"),
    ("--output-dir",                    "output_dir",
        "directory for default analysis files and assembled products "
        "(per-modality subdirs when --vaccine-type has 2+ entries)"),
)
