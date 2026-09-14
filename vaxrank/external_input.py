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

"""Build vaccine candidates from normalized LENS and pVACseq reports.

All CLI spellings enter ``external_rescoring.load_unified_external`` for
source identity, prediction-mode validation and DSL scoring. Format-specific
adapters here supply context and antigen evidence to the shared
``(source, list[VaccinePeptide])`` output interface.

LENS supplies reported peptide windows; pVACseq tables generally supply only
the epitope itself. Neither path reconstructs missing protein sequence.
RNA counts and their derivations come from normalized source evidence;
expression values such as TPM are not substituted for read counts.
"""

import dataclasses
import logging
import os

import pandas as pd

from . import cells
from .amino_acids import has_only_standard_amino_acids
from .epitope_logic import slice_epitopes
from .candidate_epitope import (
    COMPARATOR_WT,
    SOURCE_CLASS_MUTATION,
    SOURCE_CLASS_SELF,
)
from .external_prediction import (
    _PVACSEQ_DNA_VAF_COLUMNS,
    _PVACSEQ_RNA_DEPTH_COLUMNS,
    _PVACSEQ_RNA_VAF_COLUMNS,
    external_text,
    external_values,
    lens_variant_id,
    pvacseq_variant_id,
)
from .external_report import GENOMIC_VARIANT_COLUMN, ExternalRecord
from .mutant_protein_fragment import MutantProteinFragment
from .ranking import DEFAULT_RANKING_RULES, rank_constructs
from .reference_proteome import self_reference_matches
from .vaccine_antigen import (
    ANTIGEN_KIND_CTA,
    ANTIGEN_KIND_ERV,
    ANTIGEN_KIND_FUSION,
    ANTIGEN_KIND_SPLICE,
    ATTESTATION_ADMITTED,
    ATTESTATION_OVERRIDDEN,
    AminoAcidInterval,
    TargetableMask,
    TumorSpecificityAttestation,
    TumorSpecificityEvidence,
    VaccineAntigen,
)
from .vaccine_config import (
    DEFAULT_COMBINED_SCORE_EXPR,
    DEFAULT_INCLUDED_ANTIGEN_SOURCES,
)
from .vaccine_library import truncate_at_stop_codon
from .vaccine_peptide import VaccinePeptide

logger = logging.getLogger(__name__)


# Marker appended to ``patient_info.mhc_alleles`` so reports + the
# linker-optimizer plumbing can tell user-supplied alleles from
# alleles inferred from a LENS / pVACseq report. Shared with
# ``vaxrank.cli.entry_point`` so producer + consumer can't drift.
LENS_PROVENANCE_MARKER = '(inferred from report)'


def _antigen_kind_sort_key(kind, count):
    """Stable antigen_source ordering — SNV / INDEL first, then by
    descending count, ``(missing)`` last. Shared by the input
    breakdown and the funnel's coord breakdown so both read the same.
    """
    priority = {'SNV': 0, 'INDEL': 1}
    if kind.upper() in priority:
        return (priority[kind.upper()], 0, kind)
    if kind == '(missing)':
        return (3, 0, kind)
    return (2, -count, kind)


def check_varcode_annotation(variant, transcript, provided_gene,
                             provided_is_frameshift):
    """Cross-check varcode against an external tool's own annotation.

    Shared by every external loader (LENS, pVACseq, …) — varcode is
    used here only to *validate* and interpret the columns the provider
    already gave us, never to supply information the provider lacks.

    Compares on the provider's own transcript (different isoforms
    renumber residues, so we check gene + effect *class*, never raw
    positions). Returns ``(gene_ok, effect_ok)``, or ``(None, None)``
    when varcode can't compute an effect (no resolved transcript) so
    the caller can skip the check rather than fabricate a verdict.
    """
    if transcript is None:
        return None, None
    try:
        effect = variant.effect_on_transcript(transcript)
    except Exception:
        # A sanity check must never crash the load; skip this variant.
        # Logged at DEBUG so a genuine effect-computation bug is still
        # diagnosable rather than silently swallowed.
        logger.debug(
            "varcode could not compute an effect for %s on %s; skipping "
            "the annotation cross-check.", variant, transcript,
            exc_info=True)
        return None, None
    gene_ok = (not provided_gene or not effect.gene_name
               or provided_gene == effect.gene_name)
    effect_ok = (
        bool(provided_is_frameshift) == (type(effect).__name__ == 'FrameShift'))
    return gene_ok, effect_ok


def log_varcode_agreement(results, source_name):
    """Summarize a varcode-vs-provider annotation cross-check in one log
    line. ``results`` is a list of ``(gene_ok, effect_ok, label)``
    tuples from :func:`check_varcode_annotation` (entries with
    ``gene_ok is None`` — varcode couldn't compute — are ignored).
    Shared by every external loader so the summary isn't re-implemented
    per modality.
    """
    checked = [(g, e, lbl) for g, e, lbl in results if g is not None]
    if not checked:
        return
    gene_bad = [lbl for g, e, lbl in checked if not g]
    effect_bad = [lbl for g, e, lbl in checked if not e]
    if gene_bad or effect_bad:
        # Distinct variants with any mismatch (a variant can be in both
        # lists), so the headline count is truthful.
        n_bad = len(set(gene_bad) | set(effect_bad))
        example = (gene_bad or effect_bad)[0]
        logger.warning(
            "varcode disagrees with %s on %d / %d checked variant(s): "
            "%d gene, %d effect-class mismatch(es) (e.g. %s). %s values "
            "are kept; check that the pyensembl release matches the build "
            "%s used.",
            source_name, n_bad, len(checked),
            len(gene_bad), len(effect_bad), example, source_name, source_name)
    else:
        logger.info(
            "varcode agrees with %s on all %d checked variant(s) "
            "(gene + effect class).", source_name, len(checked))


def variant_is_frameshift(variant):
    """True if the variant's genomic ref/alt imply a frameshift.

    Modality-agnostic — derived from the provider's own ref/alt
    nucleotides (LENS, pVACseq, VCF all build a varcode ``Variant`` the
    same way), not from a tool-specific annotation column. A frameshift
    indel changes coding length by a non-multiple of 3.
    """
    ref = variant.ref or ''
    alt = variant.alt or ''
    return abs(len(ref) - len(alt)) % 3 != 0


def maximal_mutant_span(rep_start, rep_end, peptides, context,
                        is_frameshift):
    """Mutant-residue span within ``context`` for an external antigen.

    Shared by every external loader so LENS and pVACseq mark mutations
    identically. For a **frameshift** the whole downstream sequence (to
    the new stop) is novel, but a single representative neoepitope only
    covers part of it — so we combine *every* neoepitope the provider
    listed for the variant: the union of their located spans tiles the
    novel tail and recovers the maximal mutant region. Non-frameshift
    variants keep the representative span (the mutation is local; a
    union would over-extend into wild-type flanks).

    Uses only provider-supplied peptides (no varcode-derived sequence).
    The start comes from combining rows (the earliest located neoepitope
    ≈ the frameshift onset); the end extends to the end of ``context``,
    since for a frameshift everything from the onset to the new stop
    (= the end of the translated context) is novel — even residues no
    single neoepitope happened to cover.

    Precondition: ``context`` is the mutant protein translation ending
    at the (new) stop codon, with no trailing wild-type sequence — true
    for LENS / pVACseq ``pep_context``. If a provider ever appended
    post-stop sequence, the extend-to-end step would over-mark it.
    """
    if not is_frameshift:
        return rep_start, rep_end
    start = rep_start
    for peptide in peptides:
        if not peptide:
            continue
        idx = context.find(peptide)
        if idx >= 0:
            start = min(start, idx)
    return start, len(context)


def log_transcript_resolution(n_with_ids, n_resolved, source_name,
                              id_label="transcript IDs"):
    """Warn once when some variants' transcript IDs didn't resolve
    against the configured pyensembl release (almost always a release
    mismatch). Shared by every external loader."""
    if n_with_ids and n_resolved < n_with_ids:
        logger.warning(
            "%d / %d variant(s) had %s that didn't resolve against the "
            "configured pyensembl release. Most often this is a release "
            "mismatch — pass --ensembl-release N to match the build %s "
            "used.", n_with_ids - n_resolved, n_with_ids, id_label,
            source_name)


@dataclasses.dataclass
class ExternalVariantEntry:
    """Represent one external variant and its translation status.

    LENS ingestion records the parsed variant and optional vaccine peptide
    alongside the transcript, annotation, VAF, and parsing outcomes that are
    aggregated into the final external-input report.
    """
    variant: object = None
    source: object = None
    vaccine_peptide: object = None
    had_transcript_ids: bool = False
    resolved_transcript: bool = False
    resolved_protein_context: bool = False
    annotation: object = None       # (gene_ok, effect_ok, label) or None
    dna_vaf: object = None          # float or None, from a DNA-qualified column
    # A variant fraction whose assay the source did not state. LENS names
    # its read columns rna_* and leaves `vaf` bare, so filing it as DNA
    # asserts something the file never said. Kept apart from dna_vaf so the
    # report can label it honestly rather than by guess.
    source_vaf: object = None       # float or None
    has_rna_support: bool = False
    unparseable: bool = False

    @property
    def ranking_source(self):
        """Opaque key carried beside this entry's vaccine peptides."""
        return self.source if self.source is not None else self.variant


@dataclasses.dataclass(frozen=True)
class ExternalInputSummary:
    """Filter-independent counts recovered from a parseable input report."""

    num_somatic_variants: int = 0
    num_coding_effect_variants: int = 0
    num_variants_with_rna_support: int = 0


@dataclasses.dataclass
class ExternalRankingResult:
    """Construct ranking plus metadata from the same external-file parse."""

    ranked: list = dataclasses.field(default_factory=list)
    entries: tuple = ()
    dna_vaf_by_variant: dict = dataclasses.field(default_factory=dict)
    source_vaf_by_variant: dict = dataclasses.field(default_factory=dict)
    input_summary: ExternalInputSummary = dataclasses.field(
        default_factory=ExternalInputSummary)


@dataclasses.dataclass(frozen=True)
class ExternalConstructSelection:
    """DSL-selected source window and the candidates it contains."""

    representative: ExternalRecord
    records: tuple[ExternalRecord, ...]

    @property
    def epitopes(self):
        """Candidates in the selected source window, without duplicates."""
        seen = set()
        result = []
        for record in self.records:
            identity = record.epitope.prediction_id
            if identity not in seen:
                seen.add(identity)
                result.append(record.epitope)
        return tuple(result)


def select_external_construct(records):
    """Select one source window by configured DSL score.

    The highest ``CandidateEpitope.epitope_score`` chooses the construct
    context. Other eligible epitopes from that exact context accompany it;
    evidence from another context can never inflate the selected construct.
    File order is the deterministic tie-breaker.
    """
    records = tuple(records)
    if not records:
        return None
    representative = records[0]
    for record in records[1:]:
        if record.epitope.epitope_score > representative.epitope.epitope_score:
            representative = record
    construct_id = representative.key.construct_identifier
    selected = tuple(
        record
        for record in records
        if record.key.construct_identifier == construct_id
    )
    return ExternalConstructSelection(
        representative=representative,
        records=selected,
    )


def parse_lens_variant_coordinates(coords):
    """Parse a LENS ``variant_coords`` string into ``(contig, pos)``
    or ``None`` if the cell is empty / malformed.

    Real LENS files emit several shapes:
      - ``chr1:26780312``         (2-part: chr + pos — LENS v1.9 dominant)
      - ``chr1:26780312:C``       (3-part: chr + pos + ref)
      - ``chr1:26780312:C:T``     (4-part: chr + pos + ref + alt)

    This helper *only* extracts the (chr, pos) tuple. Ref/alt are
    looked up from the dedicated per-antigen-source columns
    (``snv_ref_allele`` / ``snv_alt_allele`` for SNVs,
    ``indel_ref_allele`` / ``indel_alt_allele`` for indels) by
    :func:`variant_from_lens_row`, which builds the actual
    ``Variant``.

    NaN / empty / 'nan' returns ``None`` — non-SNV antigen rows
    (CTA / ERV / SPLICE / FUSION) genuinely don't carry genome
    coords here; LENS records their provenance in ``splice_coords``
    / ``fusion_*`` / ``erv_orf_id`` / etc.
    """
    if coords is None or not isinstance(coords, str):
        return None
    s = coords.strip()
    if not s or s.lower() == 'nan':
        return None
    parts = s.split(":")
    if len(parts) < 2:
        return None
    contig = parts[0]
    # Strip the ``chr`` prefix LENS emits — pyensembl's Ensembl-style
    # references use bare contigs (``3`` not ``chr3``), and varcode's
    # default contig normalization is bypassed downstream
    # (``normalize_contig_names=False``). Keeping the prefix means
    # ``variant.effect_on_transcript`` finds no transcripts and the
    # variant.short_description renders as ``chrchr3 …``.
    if contig.lower().startswith('chr'):
        contig = contig[3:]
    try:
        return contig, int(parts[1])
    except (ValueError, IndexError):
        return None


def normalize_lens_allele(value):
    """LENS records alt alleles as bracketed strings: ``'[T]'``,
    ``'[CA]'``. Strip the brackets and return the inner sequence;
    return None for missing / NaN / empty."""
    if value is None:
        return None
    if isinstance(value, float) and pd.isna(value):
        return None
    s = str(value).strip()
    if not s or s.lower() == 'nan':
        return None
    # Strip leading '[' and trailing ']' if present
    if s.startswith('['):
        s = s[1:]
    if s.endswith(']'):
        s = s[:-1]
    s = s.strip()
    return s or None


def resolve_external_transcripts(transcript_ids, genome):
    """Resolve a list of Ensembl transcript-ID strings to pyensembl
    ``Transcript`` objects.

    ``genome`` is a ``pyensembl.EnsemblRelease`` (or any object with
    ``transcript_by_id``). Versioned IDs (``ENST00000312960.4``) are
    stripped to bare IDs before lookup since pyensembl 2.x doesn't
    auto-strip them — passing the versioned form raises "Transcript
    not found." Unresolvable IDs are dropped (logged at DEBUG so the
    caller can summarize aggregate failures at INFO level rather
    than spamming per-row).

    Returns ``[]`` when ``genome`` is None — the "no genome plumbed
    through" case (e.g. unit tests bypassing the CLI, or external
    input paths run without ``--ensembl-release``). Downstream code
    tolerates an empty transcript list.

    Only catches the specific exceptions pyensembl raises for
    not-found IDs (``ValueError``, ``KeyError``); other exceptions
    propagate so genuine bugs aren't swallowed.
    """
    if genome is None or not transcript_ids:
        return []
    resolved = []
    for tid in transcript_ids:
        if not tid:
            continue
        # pyensembl 2.x doesn't strip Ensembl version suffixes
        # (``ENST00000312960.4`` errors as "not found"); strip
        # ourselves so LENS IDs that carry the suffix still resolve.
        bare = tid.split('.', 1)[0]
        try:
            resolved.append(genome.transcript_by_id(bare))
        except (ValueError, KeyError) as e:
            logger.debug(
                "Could not resolve transcript_id %r against the "
                "configured pyensembl release (%s); dropping.", tid, e)
    return resolved


def infer_genome_build_from_lens(lens_path, max_rows=500):
    """Peek at a LENS file's ``origin_descriptor`` column for build
    markers and return ``'GRCh38'`` / ``'GRCh37'`` / None.

    LENS encodes the genome build inside ERV-row ``origin_descriptor``
    values like ``Hsap38.chr2.156963765.156964472.-`` (build 38) or
    ``Hsap37.chr1.…`` (build 37). SNV / INDEL rows use the
    ``ENSGxxxxx.N:GENE`` form which doesn't encode the build directly.

    We sample the first ``max_rows`` and return the first build seen.
    Returns ``None`` when the column is missing or no Hsap37/Hsap38
    marker appears in the sample (often the case on synthetic test
    fixtures with only SNV rows).
    """
    try:
        df = pd.read_csv(
            lens_path, sep='\t', usecols=['origin_descriptor'],
            nrows=max_rows, low_memory=False)
    except (ValueError, FileNotFoundError):
        return None
    for desc in df['origin_descriptor'].dropna():
        s = str(desc)
        if 'Hsap38' in s:
            return 'GRCh38'
        if 'Hsap37' in s:
            return 'GRCh37'
    return None


# Map build → reasonable default Ensembl release. Conservative
# pinning: 75 is the canonical GRCh37 release; 102 is a stable
# mid-2020 GRCh38 release. Users on different installations should
# pass --ensembl-release explicitly; this mapping is only consulted
# for the pre-flight hint, not for silent auto-selection.
# Canonical "last release for this build" — only GRCh37 has one of
# these (release 75 was the final GRCh37 mainline release; everything
# after switched to GRCh38). For GRCh38 we don't pick a default
# number; instead we look at what the user actually has installed.
_BUILD_TO_CANONICAL_RELEASE = {
    'GRCh37': 75,
}


def installed_ensembl_releases_for_build(build):
    """Return the sorted list of Ensembl release numbers the user has
    *already downloaded* for ``build`` (e.g. 'GRCh38' / 'GRCh37'),
    derived from the pyensembl cache directory layout.

    pyensembl caches under ``<platformdirs>/pyensembl/<build>/ensembl<N>``
    so a directory listing is the authoritative answer to "what
    releases can the user pass to --ensembl-release without paying a
    download." Returns ``[]`` when the cache directory doesn't exist.

    Used to make ``--ensembl-release`` suggestions concrete:
    ``origin_descriptor`` only tells us the build (GRCh37 vs GRCh38),
    not the release — and a build spans many releases (GRCh38 covers
    Ensembl 76 through current ~113+). Suggesting an arbitrary number
    misleads the user; suggesting one they already have is honest.
    """
    import glob
    import re

    import platformdirs

    cache_root = platformdirs.user_cache_dir('pyensembl')
    pattern = os.path.join(cache_root, build, 'ensembl*')
    releases = []
    for path in glob.glob(pattern):
        m = re.search(r'ensembl(\d+)$', path)
        if m and os.path.isdir(path):
            releases.append(int(m.group(1)))
    return sorted(set(releases))


def variant_from_lens_row(row, genome=None):
    """Build a ``varcode.Variant`` from a LENS row using real ref/alt.

    LENS dedicates per-antigen-source columns for ref/alt:
      - SNV rows:   ``snv_ref_allele``, ``snv_alt_allele``
      - INDEL rows: ``indel_ref_allele``, ``indel_alt_allele``

    Both alt columns use bracket notation (``[T]``, ``[CA]``).
    ``variant_coords`` carries only ``chr:pos``; we glue the real
    alleles in here so downstream consumers (varcode-effect
    annotation, etc.) get a real biological genotype rather than a
    placeholder.

    Returns ``None`` when the row lacks parseable coords + alleles
    — non-SNV / non-INDEL rows (CTA / ERV / SPLICE / FUSION) have
    NaN ``variant_coords`` and are skipped upstream.
    """
    from varcode import Variant
    coords_parsed = parse_lens_variant_coordinates(row.get('variant'))
    if coords_parsed is None:
        return None
    contig, pos = coords_parsed
    antigen_source = (row.get('antigen_source') or '').upper()
    if antigen_source == 'SNV':
        ref = normalize_lens_allele(row.get('snv_ref_allele'))
        alt = normalize_lens_allele(row.get('snv_alt_allele'))
    elif antigen_source == 'INDEL':
        ref = normalize_lens_allele(row.get('indel_ref_allele'))
        alt = normalize_lens_allele(row.get('indel_alt_allele'))
    else:
        # SPLICE / FUSION / CTA-SELF / ERV rows have neither variant
        # coords (NaN handled above) nor SNV/INDEL alleles. Caller
        # shouldn't reach here for those, but defend.
        return None
    if not ref or not alt:
        logger.debug(
            "LENS row at %s:%d (%s) missing ref/alt alleles "
            "(ref=%r alt=%r); skipping.",
            contig, pos, antigen_source, ref, alt)
        return None
    try:
        return Variant(
            contig=contig, start=pos, ref=ref, alt=alt,
            genome=genome if genome is not None else cells.text(row.get('input_reference_assembly')) or None,
            normalize_contig_names=False)
    except Exception:
        logger.debug(
            "varcode rejected LENS row at %s:%d ref=%r alt=%r; skipping.",
            contig, pos, ref, alt, exc_info=True)
        return None


def peptide_offsets_in_context(peptide, peptide_context):
    """Locate the neoepitope inside its surrounding context.

    Returns ``(start, end)`` AA offsets of the peptide within
    ``peptide_context``. Returns ``(None, None)`` when the peptide can't
    be located — the caller should drop the row rather than fabricate
    a mutation span. Previously this defaulted to "the whole context
    is the mutation," which falsely told downstream code that every
    residue was mutated.
    """
    if not peptide_context or not peptide:
        return None, None
    idx = peptide_context.find(peptide)
    if idx < 0:
        return None, None
    return idx, idx + len(peptide)


def _read_counts_from_lens_row(row):
    """Read the counts topiary derived from a LENS row.

    Returns ``(n_overlapping_reads, n_alt_reads, n_ref_reads,
    n_alt_reads_supporting_protein_sequence)``.

    vaxrank used to derive these itself, and got the central one wrong.
    ``rna_reads_covering_genomic_origin_with_peptide_cds`` is a genuine
    count, but of reads overlapping the peptide's coding sequence — not of
    reads supporting the variant allele — and it was assigned straight to
    ``n_alt_reads``, with ``n_ref_reads`` derived from that. The two are
    different quantities and disagree in both directions:

      depth   cds_overlap   vaf      was n_alt   depth x vaf
       2337             2   0.022           2            51
        765            12   0.008          12             6
        291           166   0.258         166            75

    ``n_alt_reads`` fed the default combined score at the time, so this
    reordered vaccine peptides. It remains available for explicit read-based
    score expressions; the current default uses ``n_rna_alt``.

    topiary populates the fields with the derivation named per field
    (``rna_evidence_method``): depth x VAF for the alt/ref split, and the
    CDS-overlap count kept where it belongs, in
    ``n_alt_reads_supporting_protein_sequence``. A row whose source could
    not answer carries no estimate rather than a zero standing in for one;
    the zero appears here only because the fragment field is not nullable.
    """
    return (
        cells.integer(row.get('n_rna_overlapping'), default=0),
        cells.integer(row.get('n_rna_alt'), default=0),
        cells.integer(row.get('n_rna_ref'), default=0),
        # topiary dropped its supporting-count column as duplicative of
        # n_rna_* for single-unit sources, so this comes from the LENS
        # column directly. It is a pass-through of a number the file
        # states, not a derivation — and it counts reads overlapping the
        # peptide's CDS, which is why it lives in this field and not in
        # n_rna_alt.
        cells.integer(
            cells.first(
                row,
                # topiary 5.47.0 prefixes each tool's own columns with the
                # tool name; the bare spelling is the pre-5.47 one.
                'lens_rna_reads_covering_genomic_origin_with_peptide_cds',
                'rna_reads_covering_genomic_origin_with_peptide_cds'),
            default=0),
    )


def _dna_evidence_from_row(row):
    """DNA counts and their labels, or blanks when the source has none.

    Only pVACseq answers this, and only its all_epitopes flavour, which
    reports a DNA depth alongside the fraction. The aggregated flavour
    carries a fraction with no depth, so the counts stay absent rather than
    being invented from the fraction alone. LENS states no assay for its
    single fraction and gets no DNA evidence at all.
    """
    return {
        "n_dna_overlapping": cells.integer(
            row.get("n_dna_overlapping"), default=None),
        "n_dna_alt": cells.integer(row.get("n_dna_alt"), default=None),
        "n_dna_ref": cells.integer(row.get("n_dna_ref"), default=None),
        "dna_vaf": cells.number(row.get("dna_vaf")),
        "dna_evidence_method": cells.text(row.get("dna_evidence_method")),
        "dna_evidence_subject": cells.text(row.get("dna_evidence_subject")),
    }


def _read_provenance_from_lens_row(row):
    """``(rna_evidence_method, sequence_source, rna_evidence_subject)``.

    Read alongside the counts rather than derived: topiary labels each
    field with how it was obtained, and dropping the label leaves an
    estimate indistinguishable from a measurement (#375).
    """
    return (
        cells.text(row.get('rna_evidence_method')),
        cells.text(row.get('sequence_source')),
        cells.text(row.get('rna_evidence_subject')),
    )


@dataclasses.dataclass
class ExternalRankingAccumulator:
    """Fold per-variant entries into a ranking result.

    Both external formats walk their variant groups the same way: build the
    filter-independent metadata, build the construct, then tally the same six
    facts. Only *how* a group is turned into an entry differs per format, so
    that is the only thing the callers still supply. Keeping the tally here
    means the LENS and pVACseq ``PatientInfo`` counts cannot drift apart the
    way two hand-maintained counter blocks did.
    """

    require_target_epitopes: bool = True
    ranked: list = dataclasses.field(default_factory=list)
    entries: list = dataclasses.field(default_factory=list)
    dna_vaf_by_variant: dict = dataclasses.field(default_factory=dict)
    source_vaf_by_variant: dict = dataclasses.field(default_factory=dict)
    annotation_results: list = dataclasses.field(default_factory=list)
    n_unparseable: int = 0
    n_parseable: int = 0
    n_with_rna: int = 0
    n_with_transcript_ids: int = 0
    n_resolved_transcripts: int = 0
    n_resolved_protein_contexts: int = 0

    def add(self, entry):
        """Record one variant's metadata and construct outcome."""
        if entry is None or entry.unparseable:
            self.n_unparseable += 1
            return
        source = entry.ranking_source
        if source is None:
            return
        if (self.require_target_epitopes and entry.vaccine_peptide is not None
                and not entry.vaccine_peptide.contains_target_epitopes()):
            # Match direct-input admission before either final ranking or
            # repeated-source selection. Input evidence still enters the audit.
            entry = dataclasses.replace(entry, vaccine_peptide=None)
        self.entries.append(entry)
        self.n_parseable += 1
        if entry.annotation is not None:
            self.annotation_results.append(entry.annotation)
        if entry.had_transcript_ids:
            self.n_with_transcript_ids += 1
            if entry.resolved_transcript:
                self.n_resolved_transcripts += 1
        if entry.resolved_protein_context:
            self.n_resolved_protein_contexts += 1
        if entry.has_rna_support:
            self.n_with_rna += 1
        if entry.dna_vaf is not None and entry.variant is not None:
            self.dna_vaf_by_variant[entry.variant] = entry.dna_vaf
        if entry.source_vaf is not None and entry.variant is not None:
            self.source_vaf_by_variant[entry.variant] = entry.source_vaf
        if entry.vaccine_peptide is not None:
            self.ranked.append((source, [entry.vaccine_peptide]))

    def result(self, source_name, transcript_id_label="transcript IDs"):
        """Emit the summaries every external format owes its caller."""
        log_transcript_resolution(
            self.n_with_transcript_ids, self.n_resolved_transcripts,
            source_name, id_label=transcript_id_label)
        log_varcode_agreement(self.annotation_results, source_name)
        return ExternalRankingResult(
            entries=tuple(self.entries),
            ranked=rank_constructs(self.ranked),
            dna_vaf_by_variant=self.dna_vaf_by_variant,
            source_vaf_by_variant=self.source_vaf_by_variant,
            input_summary=ExternalInputSummary(
                num_somatic_variants=self.n_parseable,
                # Legacy field name. For source-agnostic inputs this counts
                # sources with a usable protein context, including fusions.
                num_coding_effect_variants=self.n_resolved_protein_contexts,
                num_variants_with_rna_support=self.n_with_rna,
            ),
        )


def lens_variant_metadata(variant_id, group_rows, genome=None):
    """Aggregate filter-independent facts from every raw LENS row."""
    parsed = [
        (row, variant_from_lens_row(row, genome=genome))
        for row in group_rows
    ]
    parsed = [(row, variant) for row, variant in parsed if variant is not None]
    if not parsed:
        logger.debug(
            "Could not build Variant for LENS identity %r; skipping.",
            variant_id)
        return ExternalVariantEntry(unparseable=True)
    representative_row, variant = parsed[0]

    transcript_ids = external_values(*(
        value
        for row in group_rows
        for value in (
            row.get("transcript_id"),
            row.get("all_transcript_ids_encoding_peptide"),
        )
    ))
    transcripts = resolve_external_transcripts(transcript_ids, genome)
    has_rna_support = any(
        counts[0] > 0 or counts[1] > 0
        for counts in map(_read_counts_from_lens_row, group_rows)
    )
    gene_names = external_values(*(
        value
        for row in group_rows
        for value in (
            row.get("gene"),
            row.get("all_gene_names_encoding_peptide"),
        )
    ))
    lens_is_frameshift = any(
        external_text(row.get("indel_type")).lower() == "frameshift"
        or "fs" in external_text(row.get("effect")).lower()
        for row in group_rows
    )
    entry = ExternalVariantEntry(
        variant=variant,
        had_transcript_ids=bool(transcript_ids),
        resolved_transcript=bool(transcripts),
        resolved_protein_context=bool(transcripts),
        has_rna_support=has_rna_support,
        annotation=(check_varcode_annotation(
            variant,
            transcripts[0] if transcripts else None,
            gene_names[0] if gene_names else "",
            lens_is_frameshift,
        ) + (str(variant_id),)),
    )
    for row in group_rows:
        # lens_vaf is topiary's source-prefixed original; `vaf` is the
        # pre-5.47 spelling. Deliberately not `rna_vaf`, which topiary also
        # publishes on LENS frames as a verbatim copy of this column —
        # reading that would inherit an assay assertion the file never made
        # and this function exists to avoid.
        value = cells.first(row, "lens_vaf", "vaf")
        if cells.missing(value):
            continue
        try:
            # LENS names its read columns rna_* explicitly and leaves this
            # one bare, so the assay is unstated. Recorded as the source's
            # own fraction rather than as a DNA VAF: calling it DNA asserts
            # an assay the file never claimed, and the report printed it
            # under "DNA VAF" on that basis. topiary reached the same
            # conclusion for its own use of this column, which is why the
            # LENS split is rna_depth_x_source_vaf.
            entry.source_vaf = float(value)
            break
        except (TypeError, ValueError):
            continue
    return entry


def lens_fusion_metadata(fusion_id, group_rows, genome=None):
    """Aggregate filter-independent facts for one caller-supplied fusion."""
    if not fusion_id:
        return ExternalVariantEntry(unparseable=True)
    transcript_ids = external_values(*(
        value
        for row in group_rows
        for value in (
            row.get("fusion_left_transcript"),
            row.get("fusion_right_transcript"),
            row.get("transcript_id"),
            row.get("all_transcript_ids_encoding_peptide"),
        )
    ))
    has_protein_context = any(
        bool(external_text(row.get("pep_context"))) for row in group_rows)
    transcripts = resolve_external_transcripts(transcript_ids, genome)
    # A LENS FUSION row is itself an RNA-derived antigen observation. Some
    # report versions do not populate the SNV-oriented rna_reads_* columns,
    # so absence there must not relabel the fusion as DNA-only.
    return ExternalVariantEntry(
        source=str(fusion_id),
        had_transcript_ids=bool(transcript_ids),
        resolved_transcript=bool(transcripts),
        resolved_protein_context=has_protein_context,
        has_rna_support=True,
    )


def lens_source_antigen_metadata(source_id, group_rows, genome=None):
    """Aggregate filter-independent facts for a non-coordinate antigen."""
    if not source_id:
        return ExternalVariantEntry(unparseable=True)
    transcript_ids = external_values(*(
        value
        for row in group_rows
        for value in (
            row.get("transcript_id"),
            row.get("all_transcript_ids_encoding_peptide"),
        )
    ))
    transcripts = resolve_external_transcripts(transcript_ids, genome)
    has_protein_context = any(
        bool(external_text(row.get("pep_context"))) for row in group_rows
    )
    has_rna_support = any(
        counts[0] > 0 or counts[1] > 0
        for counts in map(_read_counts_from_lens_row, group_rows)
    ) or any(
        not cells.missing(cells.first(
            row, "tpm", "erv_tumor_cpm", "snaf_exp"
        ))
        for row in group_rows
    )
    return ExternalVariantEntry(
        source=str(source_id),
        had_transcript_ids=bool(transcript_ids),
        resolved_transcript=bool(transcripts),
        resolved_protein_context=has_protein_context,
        has_rna_support=has_rna_support,
    )


@dataclasses.dataclass(frozen=True)
class ExternalConstructOptions:
    """Everything the vaccine-config layer contributes to one construct.

    External inputs used to build ``VaccinePeptide`` objects with nothing but
    a peptide length, so ``vaccine_peptides:`` config — combined-score
    expression, ranking rules, manufacturability thresholds, how many target
    epitopes to keep — applied to VCF runs and silently did nothing to LENS /
    pVACseq runs. Passing one object keeps the two paths honest.
    """

    vaccine_peptide_length: int = 25
    num_target_epitopes_to_keep: object = None
    combined_score_expr: object = None
    ranking_rules: object = None
    require_target_epitopes_in_variant: bool = True
    manufacturability_thresholds: dict = dataclasses.field(
        default_factory=dict)
    manufacturability_rules: object = None
    included_antigen_sources: tuple[str, ...] = (
        DEFAULT_INCLUDED_ANTIGEN_SOURCES
    )

    @classmethod
    def from_configs(cls, vaccine_config=None, manufacturability_config=None,
                     vaccine_peptide_length=None,
                     num_target_epitopes_to_keep=None):
        """Build options from the resolved config objects, if any."""
        if vaccine_config is None:
            length = vaccine_peptide_length or 25
            keep = num_target_epitopes_to_keep
            expr = None
            rules = None
            included_sources = DEFAULT_INCLUDED_ANTIGEN_SOURCES
        else:
            length = (
                vaccine_peptide_length
                or vaccine_config.preferred_peptide_length)
            keep = (
                num_target_epitopes_to_keep
                if num_target_epitopes_to_keep is not None
                else vaccine_config.num_target_epitopes_to_keep)
            expr = vaccine_config.combined_score_expr
            rules = vaccine_config.ranking_rules
            included_sources = vaccine_config.included_antigen_sources
        if manufacturability_config is None:
            thresholds, mfg_rules = {}, None
        else:
            thresholds = manufacturability_config.thresholds_dict()
            mfg_rules = manufacturability_config.rules
        return cls(
            vaccine_peptide_length=length,
            num_target_epitopes_to_keep=keep,
            combined_score_expr=expr,
            ranking_rules=rules,
            require_target_epitopes_in_variant=(
                vaccine_config.require_target_epitopes_in_variant
                if vaccine_config is not None else True),
            manufacturability_thresholds=thresholds,
            manufacturability_rules=mfg_rules,
            included_antigen_sources=tuple(included_sources),
        )


_SOURCE_AGNOSTIC_RANKING_RULES = (
    "target_epitope_score",
    "manufacturability",
    "self_epitope_score",
)


def source_agnostic_construct_options(options):
    """Use safe defaults for antigens with no mutation read-count fields.

    A source-agnostic expression supplied by the user is retained. The legacy
    mutation defaults are translated to their antigen equivalents; custom
    expressions that name mutation-only fields still fail explicitly in
    :class:`VaccinePeptide` instead of being silently rewritten.
    """
    expression = options.combined_score_expr
    if expression in (None, DEFAULT_COMBINED_SCORE_EXPR):
        expression = "target_epitope_score"
    rules = options.ranking_rules
    if rules is None or tuple(rules) == tuple(DEFAULT_RANKING_RULES):
        rules = _SOURCE_AGNOSTIC_RANKING_RULES
    return dataclasses.replace(
        options,
        combined_score_expr=expression,
        ranking_rules=tuple(rules),
    )


def _merged_epitope_intervals(epitopes):
    """Merge overlapping/adjacent candidate spans into a targetable mask."""
    spans = sorted(
        (epitope.offset, epitope.offset + len(epitope.sequence))
        for epitope in epitopes)
    merged = []
    for start, end in spans:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return tuple(AminoAcidInterval(start, end) for start, end in merged)


def _fusion_partner(value):
    """Split LENS's ``GENE^ENSG...`` partner spelling."""
    text = external_text(value)
    if not text:
        return "", ""
    gene, separator, gene_id = text.partition("^")
    return gene, gene_id if separator else ""


_LENS_ANTIGEN_KINDS = {
    "FUSION": ANTIGEN_KIND_FUSION,
    "SPLICE": ANTIGEN_KIND_SPLICE,
    "CTA/SELF": ANTIGEN_KIND_CTA,
    "ERV": ANTIGEN_KIND_ERV,
}


def _lens_source_metadata(row, antigen_source, *, sequence_source=""):
    """Retain the source fields needed to audit a LENS antigen."""
    names = {
        "antigen_source",
        "gene_detectable_normal_tissues",
        "gene_expression",
        "gene_id",
        "gene_main_subcellular_location",
        "gene_name",
        "gene_tpm",
        "gene_tpm_raw",
        "all_gene_ids_encoding_peptide",
        "all_gene_names_encoding_peptide",
        "all_transcript_ids_encoding_peptide",
        "lens_proportion_rna_reads_covering_genomic_origin_with_peptide_cds",
        "lens_rna_reads_covering_genomic_origin",
        "lens_rna_reads_covering_genomic_origin_with_peptide_cds",
        "mean_mtec_num_reads",
        "mean_mtec_tpm",
        "median_mtec_num_reads",
        "median_mtec_tpm",
        "mtec_p95_tpm",
        "norm_tissue_p95_tpm",
        "origin_descriptor",
        "primary_aln_rna_reads_covering_genomic_origin_with_peptide_cds",
        "snaf_exp",
        "source_sequence_name",
        "stdev_mtec_num_reads",
        "stdev_mtec_tpm",
        "tpm",
        "transcript_id",
    }
    if antigen_source == "FUSION":
        names.update({"fusion_annotation", "fusion_id", "fusion_type"})
    elif antigen_source == "SPLICE":
        names.update({
            "coding_sequence", "nt_context", "splice_coords",
            "splice_description",
        })
    elif antigen_source == "ERV":
        names.update(name for name in row if name.startswith("erv_"))
    metadata = {
        name: external_text(row.get(name))
        for name in names
        if external_text(row.get(name))
    }
    counts = _read_counts_from_lens_row(row)
    rna_method, stated_sequence_source, rna_subject = (
        _read_provenance_from_lens_row(row)
    )
    metadata.update({
        name: value
        for name, value in (
            ("rna_reads_covering_source", str(counts[0]) if counts[0] else ""),
            ("rna_reads_supporting_source", str(counts[1]) if counts[1] else ""),
            ("rna_reads_supporting_protein_sequence",
             str(counts[3]) if counts[3] else ""),
            ("rna_evidence_method", rna_method),
            ("rna_evidence_subject", rna_subject),
            ("sequence_source", stated_sequence_source or sequence_source),
        )
        if value
    })
    if antigen_source == "FUSION":
        left_gene, left_gene_id = _fusion_partner(
            row.get("fusion_left_gene")
        )
        right_gene, right_gene_id = _fusion_partner(
            row.get("fusion_right_gene")
        )
        metadata.update({
            name: value
            for name, value in (
                ("left_breakpoint", external_text(
                    row.get("fusion_left_breakpoint"))),
                ("left_gene", left_gene),
                ("left_gene_id", left_gene_id),
                ("left_transcript", external_text(
                    row.get("fusion_left_transcript"))),
                ("right_breakpoint", external_text(
                    row.get("fusion_right_breakpoint"))),
                ("right_gene", right_gene),
                ("right_gene_id", right_gene_id),
                ("right_transcript", external_text(
                    row.get("fusion_right_transcript"))),
                ("rna_reads_covering_breakpoint",
                 str(counts[0]) if counts[0] else ""),
                ("rna_reads_supporting_fusion",
                 str(counts[1]) if counts[1] else ""),
            )
            if value
        })
    metadata["antigen_source"] = antigen_source
    return tuple(sorted(metadata.items()))


def _lens_self_reference_matches(epitopes, antigen, genome):
    """Resolve exact self when Ensembl is installed; stay explicit otherwise."""
    peptides = tuple(dict.fromkeys(
        epitope.sequence for epitope in epitopes
    ))
    try:
        return self_reference_matches(peptides, antigen, genome)
    except ValueError as error:
        # External report mode can name an Ensembl release for annotation
        # provenance even when pyensembl's local GTF database is absent.
        # Transcript resolution already degrades to unresolved in that case;
        # exact-self must do the same rather than making an otherwise usable
        # LENS run fail.  The result remains provenance-incomplete, never a
        # claim that the reference was exhaustively searched.
        if "database needs to be created" not in str(error):
            raise
        logger.debug(
            "Could not build exact-self provenance for LENS antigen %s: %s",
            antigen.display_identifier,
            error,
        )
        return {
            peptide: antigen.self_reference_match(peptide, False)
            for peptide in peptides
        }


def lens_source_antigen_vaccine_entry(
        metadata, selection, genome=None, options=None):
    """Build a typed LENS fusion, splice, CTA, or ERV antigen.

    LENS identifies the caller-selected peptide but does not provide an exact
    amino-acid junction offset for fusion or splice rows.  Consequently only
    the supplied candidate intervals are marked targetable; this function
    never infers a junction or expands targetability to the full context.
    """
    if selection is None:
        return metadata
    options = source_agnostic_construct_options(
        options or ExternalConstructOptions())
    key = selection.representative.key
    antigen_source = key.antigen_source.upper()
    antigen_kind = _LENS_ANTIGEN_KINDS.get(antigen_source)
    if antigen_kind is None:
        return metadata
    peptide = truncate_at_stop_codon(key.peptide)
    context = truncate_at_stop_codon(key.source_sequence or key.peptide)
    if not peptide or not context:
        return metadata
    if not has_only_standard_amino_acids(context):
        logger.warning(
            "Dropped LENS %s construct %r: pep_context %r contains "
            "non-standard residues (allowed: 20 canonical AAs).",
            antigen_source, key.variant_id, context)
        return metadata
    start, end = peptide_offsets_in_context(peptide, context)
    if start is None:
        return metadata
    windowed, new_start, new_end = MutantProteinFragment.slp_window_around_mutation(
        context, start, end, options.vaccine_peptide_length)
    window_start = start - new_start
    window_end = window_start + len(windowed)
    if context[window_start:window_end] != windowed:
        raise ValueError("Could not locate the antigen SLP window in its context")
    epitopes = slice_epitopes(selection.epitopes, window_start, window_end)
    if not epitopes:
        return metadata

    row = selection.representative.row
    source_identifier = key.variant_id
    gene_name = key.primary_gene_name or (
        key.gene_names[0] if len(key.gene_names) == 1 else ""
    )
    gene_id = external_text(row.get("gene_id")) or (
        key.gene_ids[0] if len(key.gene_ids) == 1 else ""
    )
    transcript_ids = key.ordered_transcript_ids
    if antigen_source == "FUSION":
        left_gene, _left_gene_id = _fusion_partner(
            row.get("fusion_left_gene")
        )
        right_gene, _right_gene_id = _fusion_partner(
            row.get("fusion_right_gene")
        )
        gene_name = "::".join(
            value for value in (left_gene, right_gene) if value
        )
        transcript_ids = external_values(
            row.get("fusion_left_transcript"),
            row.get("fusion_right_transcript"),
            *transcript_ids,
        )
        source_identifier = (
            external_text(row.get("fusion_id")) or source_identifier
        )
    elif antigen_source == "SPLICE":
        source_identifier = (
            external_text(row.get("splice_description"))
            or external_text(row.get("splice_coords"))
            or source_identifier
        )
    elif antigen_source == "CTA/SELF":
        source_identifier = (
            external_text(row.get("origin_descriptor"))
            or gene_id
            or gene_name
            or source_identifier
        )
    elif antigen_source == "ERV":
        source_identifier = (
            external_text(row.get("erv_orf_id"))
            or external_text(row.get("origin_descriptor"))
            or source_identifier
        )
    source_metadata = _lens_source_metadata(
        row,
        antigen_source,
        sequence_source="caller-supplied LENS pep_context",
    )
    expression_derived = antigen_source in {"CTA/SELF", "ERV"}
    evidence_kind = {
        "FUSION": "caller_curated_fusion_neoantigen",
        "SPLICE": "caller_curated_aberrant_splice_neoantigen",
        "CTA/SELF": "caller_curated_CTA_self_expression",
        "ERV": "caller_curated_ERV_expression",
    }[antigen_source]
    antigen = VaccineAntigen(
        kind=antigen_kind,
        amino_acids=windowed,
        targetable_mask=TargetableMask(_merged_epitope_intervals(epitopes)),
        tumor_specificity=TumorSpecificityAttestation(
            status=(
                ATTESTATION_OVERRIDDEN
                if expression_derived
                else ATTESTATION_ADMITTED
            ),
            evidence_kind=evidence_kind,
            evidence_source="LENS report",
            patient_specific=True,
            rationale_code="lens_%s_antigen" % antigen_kind,
            requires_review=expression_derived,
            override_reason=(
                "explicit inclusion by antigen-source policy"
                if expression_derived
                else ""
            ),
            evidence_records=(TumorSpecificityEvidence(
                evidence_kind=evidence_kind,
                evidence_source="LENS report",
                subject_id=source_identifier,
                patient_specific=True,
                passed=None if expression_derived else True,
                details=source_metadata,
            ),),
        ),
        self_reference_excluded_gene_ids=(
            (gene_id,) if antigen_source == "CTA/SELF" and gene_id else ()
        ),
        gene_name=gene_name,
        gene_id=gene_id,
        transcript_ids=transcript_ids,
        species=key.species,
        source_identifier=source_identifier,
        source_metadata=source_metadata,
    )
    self_matches = _lens_self_reference_matches(epitopes, antigen, genome)
    source_class = (
        SOURCE_CLASS_SELF
        if antigen_source in {"CTA/SELF", "ERV"}
        else SOURCE_CLASS_MUTATION
    )
    epitopes = [
        dataclasses.replace(
            epitope,
            comparators={
                name: comparator
                for name, comparator in epitope.comparators.items()
                if name != COMPARATOR_WT
            },
            source_class=source_class,
            overlaps_mutation=False,
            overlaps_targetable=True,
            self_reference_match=self_matches[epitope.sequence],
        )
        for epitope in epitopes
    ]
    metadata.source = antigen
    metadata.vaccine_peptide = VaccinePeptide(
        antigen=antigen,
        epitopes=epitopes,
        num_target_epitopes_to_keep=options.num_target_epitopes_to_keep,
        manufacturability_thresholds=options.manufacturability_thresholds,
        manufacturability_rules=options.manufacturability_rules,
        combined_score_expr=options.combined_score_expr,
        ranking_rules=options.ranking_rules,
    )
    return metadata


def lens_fusion_vaccine_entry(metadata, selection, options=None, genome=None):
    """Backward-compatible fusion wrapper around the typed source path."""
    return lens_source_antigen_vaccine_entry(
        metadata, selection, genome=genome, options=options
    )


def external_vaccine_peptide(variant, selection, context, mutant_start,
                             mutant_end, gene_name, transcripts, counts,
                             options, provenance=("", "", ""),
                             dna_evidence=None):
    """Assemble one ``VaccinePeptide`` from an external construct window.

    Shared by every external format so a construct built from a LENS
    ``pep_context`` and one built from a pVACseq peptide obey the same
    invariant the VCF pipeline enforces in ``core_logic``: *every epitope
    attached to a vaccine peptide lies inside that vaccine peptide*.

    ``context`` may be longer than the configured peptide length (LENS
    ``pep_context`` is sometimes a 100+ aa protein prefix). Trimming it to an
    SLP window without also re-slicing the epitopes is what used to let a
    neoepitope 60 residues away from the mutation count toward the
    construct's ``target_epitope_score`` and render in its report table while
    being absent from the peptide that would actually be synthesized.

    ``counts`` is ``(n_overlapping, n_alt, n_ref, n_alt_supporting_protein)``.
    Returns ``None`` when nothing survives the window.
    """
    windowed, new_start, new_end = (
        MutantProteinFragment.slp_window_around_mutation(
            context, mutant_start, mutant_end,
            options.vaccine_peptide_length))
    # ``slp_window_around_mutation`` reports the mutation span rebased into
    # the window; recover the window's own start so the epitopes can be
    # rebased with it, and verify rather than trust the arithmetic.
    window_start = mutant_start - new_start
    window_end = window_start + len(windowed)
    if context[window_start:window_end] != windowed:
        raise ValueError(
            "Could not locate the SLP window inside its source context")

    epitopes = slice_epitopes(selection.epitopes, window_start, window_end)
    if not epitopes:
        logger.debug(
            "No candidate epitope from variant %r fits the %d-aa construct "
            "window; skipping.", variant, options.vaccine_peptide_length)
        return None
    epitopes = [
        epitope
        if epitope.overlaps_mutation
        else dataclasses.replace(epitope, overlaps_mutation=True)
        for epitope in epitopes
    ]
    n_total, n_alt, n_ref, n_alt_protein = counts
    # The subject comes from the reader rather than being assumed here:
    # both external sources happen to count reads, but that is a fact about
    # them, not something this function should assert on their behalf.
    rna_evidence_method, sequence_source, rna_evidence_subject = provenance
    dna_evidence = dna_evidence or {}
    fragment = MutantProteinFragment(
        variant=variant,
        gene_name=gene_name,
        amino_acids=windowed,
        mutant_amino_acid_start_offset=new_start,
        mutant_amino_acid_end_offset=new_end,
        supporting_reference_transcripts=transcripts,
        n_overlapping_reads=n_total,
        n_alt_reads=n_alt,
        n_ref_reads=n_ref,
        n_alt_reads_supporting_protein_sequence=n_alt_protein,
        placeholder_alleles=False,
        rna_evidence_method=rna_evidence_method,
        rna_evidence_subject=rna_evidence_subject,
        sequence_source=sequence_source,
        **dna_evidence,
    )
    return VaccinePeptide(
        mutant_protein_fragment=fragment,
        epitopes=epitopes,
        num_target_epitopes_to_keep=options.num_target_epitopes_to_keep,
        manufacturability_thresholds=options.manufacturability_thresholds,
        manufacturability_rules=options.manufacturability_rules,
        combined_score_expr=options.combined_score_expr,
        ranking_rules=options.ranking_rules,
    )


def lens_vaccine_entry(metadata, selection, genome=None, options=None):
    """Build one LENS vaccine peptide from a DSL-selected source window."""
    if selection is None:
        return metadata
    if selection.representative.key.antigen_source.upper() in _LENS_ANTIGEN_KINDS:
        return lens_source_antigen_vaccine_entry(
            metadata, selection, genome=genome, options=options
        )
    if metadata.variant is None:
        return metadata
    options = options or ExternalConstructOptions()
    key = selection.representative.key
    peptide = truncate_at_stop_codon(key.peptide)
    pep_context = truncate_at_stop_codon(key.source_sequence or key.peptide)
    if not peptide or not pep_context:
        return metadata
    if not has_only_standard_amino_acids(pep_context):
        logger.warning(
            "Dropped LENS construct for variant %r: pep_context %r contains "
            "non-standard residues (allowed: 20 canonical AAs).",
            key.variant_id, pep_context)
        return metadata
    start_off, end_off = peptide_offsets_in_context(peptide, pep_context)
    if start_off is None:
        return metadata
    start_off, end_off = maximal_mutant_span(
        start_off,
        end_off,
        [record.key.peptide for record in selection.records],
        pep_context,
        variant_is_frameshift(metadata.variant),
    )
    # One row's counts, not four independent maxima. Maximizing each field
    # separately can assemble a tuple no row ever reported — total from a
    # deep row, alt from a shallow one — and n_alt_reads feeds the
    # combined-score DSL, so an incoherent count reorders the ranking.
    # Best-supported row wins: most alt reads, then most coverage.
    # The provenance travels with the counts, from the same row. Reading
    # it separately would be the per-field maximum problem again, one field
    # over: a method label describing a row whose numbers were not used.
    observations = [
        (_read_counts_from_lens_row(record.row),
         _read_provenance_from_lens_row(record.row))
        for record in selection.records]
    (n_total, n_alt_reads, n_ref_reads, n_alt_protein), provenance = max(
        observations, key=lambda o: (o[0][1], o[0][0]),
        default=((0, 0, 0, 0), ("", "")))
    metadata.vaccine_peptide = external_vaccine_peptide(
        variant=metadata.variant,
        selection=selection,
        context=pep_context,
        mutant_start=start_off,
        mutant_end=end_off,
        gene_name=key.primary_gene_name,
        transcripts=resolve_external_transcripts(
            key.ordered_transcript_ids, genome),
        # Use the reference count topiary derived rather than recomputing
        # it as depth minus alt: that subtraction was only ever correct if
        # the alt count was variant support, which is exactly what was wrong.
        counts=(n_total, n_alt_reads, n_ref_reads, n_alt_protein),
        provenance=provenance,
        options=options,
    )
    return metadata


def lens_ranking_result(report, epitopes, genome=None, options=None):
    """Rank LENS constructs from an already-parsed report.

    Parameters
    ----------
    report : ExternalReport
        The single parse produced by ``read_lens_report``. Its ``records``
        already pair every source row with the identity derived from it, so
        no row is re-read and no identity is re-derived here.
    epitopes : list of CandidateEpitope
        Output of ``read_lens_report(path)``. Each CandidateEpitope groups all
        per-(allele, predictor) ``mhctools.Prediction`` records for
        one ``(peptide, source_sequence, offset)`` position.
    genome : varcode-compatible genome reference, optional
        Passed through to ``Variant`` construction. When None,
        Variants are constructed without genome resolution and
        downstream code that needs gene annotation may degrade.
    num_target_epitopes_to_keep : int, optional
        Forwarded to ``VaccinePeptide``. When None, all overlapping
        epitopes are kept.

    Returns
    -------
    ExternalRankingResult
    """
    if report is None or not report.rows:
        return ExternalRankingResult()
    options = options or ExternalConstructOptions()
    rows = list(report.rows)

    # Every ranking-eligible candidate, already bound to the source row whose
    # identity selected it. LENS uses lowercase snake_case columns.
    records_by_variant = {}
    for record in report.records_with_epitopes(epitopes):
        records_by_variant.setdefault(
            lens_variant_id(record.row), []).append(record)

    # Up-front antigen_source breakdown — logged BEFORE any filter
    # log lines so the operator sees the composition of the input
    # before they see what got dropped. Order: SNV / INDEL first
    # (the per-coord paths), then non-coord categories sorted by
    # count. Missing values surface as ``(missing)``.
    if rows:
        full_kinds: dict[str, int] = {}
        for r in rows:
            kind = r.get('antigen_source')
            kind_key = (
                str(kind).strip() if kind is not None and not (
                    isinstance(kind, float) and pd.isna(kind))
                else '(missing)')
            full_kinds[kind_key] = full_kinds.get(kind_key, 0) + 1
        ordered = sorted(
            full_kinds.items(),
            key=lambda kv: _antigen_kind_sort_key(kv[0], kv[1]))
        breakdown = ', '.join("%s=%d" % (k, v) for k, v in ordered)
    else:
        full_kinds = {}
        breakdown = ''

    groups = {}
    n_skipped_empty_coords = 0
    n_policy_excluded = 0
    # When a row has no variant_coords, the only sensible explanation
    # is a non-SNV / non-INDEL antigen kind (splice / fusion / ERV /
    # CTA-self / intron-retention). Verify that hypothesis instead of
    # asserting it — if any SNV / INDEL rows are missing coords, that's
    # an upstream bug worth surfacing distinctly.
    skipped_kinds = {}  # antigen_source value → count
    policy_excluded_kinds = {}
    for r in rows:
        coords = r.get('variant')
        kind = external_text(r.get('antigen_source')).upper()
        if kind and kind not in options.included_antigen_sources:
            n_policy_excluded += 1
            policy_excluded_kinds[kind] = (
                policy_excluded_kinds.get(kind, 0) + 1
            )
            continue
        if coords is None or (
                isinstance(coords, float) and pd.isna(coords)) or (
                isinstance(coords, str) and (
                    not coords.strip() or coords.strip().lower() == 'nan')):
            # Source identity lives outside variant_coords for these antigen
            # categories. They are first-class construct sources when the
            # explicit source-selection policy admits them.
            if kind in _LENS_ANTIGEN_KINDS and lens_variant_id(r):
                groups.setdefault(lens_variant_id(r), []).append(r)
                continue
            n_skipped_empty_coords += 1
            kind = r.get('antigen_source')
            kind_key = (
                str(kind).strip() if kind is not None and not (
                    isinstance(kind, float) and pd.isna(kind))
                else '(missing)')
            skipped_kinds[kind_key] = skipped_kinds.get(kind_key, 0) + 1
            continue
        groups.setdefault(lens_variant_id(r), []).append(r)
    skipped_breakdown = ', '.join(
        "%s=%d" % (k, v) for k, v in sorted(
            skipped_kinds.items(), key=lambda kv: -kv[1]))
    policy_excluded_breakdown = ', '.join(
        "%s=%d" % (k, v) for k, v in sorted(
            policy_excluded_kinds.items(), key=lambda kv: -kv[1]
        )
    )
    if n_skipped_empty_coords:
        # SNV / INDEL rows are *expected* to carry coords; flag them
        # separately so the user can chase upstream rather than assume
        # "typical". Remaining unsupported non-coordinate categories are
        # still scored in the neoepitope report but cannot yet enter construct
        # ranking. FUSION rows took the explicit source-antigen path above.
        unexpected = {k: v for k, v in skipped_kinds.items()
                      if k.upper() in ('SNV', 'INDEL')}
        if unexpected:
            logger.warning(
                "%d LENS row(s) declared antigen_source=%s but had "
                "empty variant_coords (SNV / INDEL antigens are "
                "expected to carry genome coords). This is likely an "
                "upstream LENS bug.",
                sum(unexpected.values()),
                '/'.join(sorted(unexpected)))

    tally = ExternalRankingAccumulator(
        require_target_epitopes=options.require_target_epitopes_in_variant)
    for variant_id, group_rows in groups.items():
        antigen_source = external_text(
            group_rows[0].get('antigen_source')).upper()
        metadata = (
            lens_fusion_metadata(variant_id, group_rows, genome=genome)
            if antigen_source == 'FUSION'
            else lens_source_antigen_metadata(
                variant_id, group_rows, genome=genome
            )
            if antigen_source in _LENS_ANTIGEN_KINDS
            else lens_variant_metadata(variant_id, group_rows, genome=genome)
        )
        tally.add(lens_vaccine_entry(
            metadata,
            select_external_construct(records_by_variant.get(variant_id, [])),
            genome=genome,
            options=options,
        ))
    n_unparseable = tally.n_unparseable
    ranked = tally.ranked

    if n_unparseable:
        logger.warning(
            "Skipped %d LENS variant(s): variant_coords couldn't be "
            "parsed as chr:pos OR the row's snv_*_allele / "
            "indel_*_allele columns were missing. See DEBUG log for "
            "the offenders.", n_unparseable)

    # Per-load summary of missing essential / important per-row data.
    # Essential = blocked the row (already counted as "skipped" above
    # or in the reader itself). Important-but-recoverable = warns so
    # users know what's degraded.
    n_no_pep_context = sum(
        1 for r in rows
        if not (r.get('pep_context') and not (
            isinstance(r.get('pep_context'), float) and pd.isna(r.get('pep_context')))))
    def _kind(r):
        k = r.get('antigen_source')
        return ('(missing)' if cells.missing(k) else str(k).strip())

    rows_no_gene_name = [r for r in rows if cells.missing(r.get('gene'))]
    rows_no_transcript = [r for r in rows if cells.missing(r.get('transcript_id'))]
    if n_no_pep_context:
        logger.warning(
            "%d / %d LENS row(s) lack pep_context — antigens for those "
            "rows degenerate to the bare neoepitope (no SLP context). "
            "Vaccine windows will be ~9 aa instead of ~25 aa.",
            n_no_pep_context, len(rows))
    # ``gene_name`` / ``transcript_id`` / ``rna_reads`` are *expected*
    # to be empty for ERV / CTA-SELF / SPLICE / FUSION antigens (no
    # canonical gene model / genome-coord origin). So we don't log the
    # full breakdown here (that just re-counts the same non-coord rows
    # the funnel already accounts for) — we only warn about SNV / INDEL
    # rows missing these, which would be a genuine upstream bug. The
    # counts feed the funnel's "note:" line.
    _SNV_OR_INDEL_KINDS = {'SNV', 'INDEL'}
    rows_no_reads = [
        r for r in rows
        if cells.missing(cells.first(
            r,
            # topiary 5.47.0 prefixes each tool's own columns with the tool
            # name. Reading only the bare spelling counted every row as
            # lacking RNA reads and warned about data that was present
            # under the other name — a rename reported as a defect in the
            # input file (#390).
            'lens_rna_reads_covering_genomic_origin',
            'rna_reads_covering_genomic_origin'))]
    n_snv_indel_no_gene = sum(
        1 for r in rows_no_gene_name if _kind(r).upper() in _SNV_OR_INDEL_KINDS)
    n_snv_indel_no_transcript = sum(
        1 for r in rows_no_transcript
        if _kind(r).upper() in _SNV_OR_INDEL_KINDS)
    n_snv_indel_no_reads = sum(
        1 for r in rows_no_reads if _kind(r).upper() in _SNV_OR_INDEL_KINDS)
    if n_snv_indel_no_gene:
        logger.warning(
            "%d SNV / INDEL row(s) lack gene_name — those antigen kinds "
            "are expected to carry one. Likely upstream LENS bug.",
            n_snv_indel_no_gene)
    if n_snv_indel_no_transcript:
        logger.warning(
            "%d SNV / INDEL row(s) lack transcript_id — those antigen "
            "kinds are expected to carry one. Likely upstream LENS bug.",
            n_snv_indel_no_transcript)
    if n_snv_indel_no_reads:
        logger.warning(
            "%d SNV / INDEL row(s) lack RNA-read counts "
            "(lens_rna_reads_covering_genomic_origin) — those kinds are "
            "expected to carry them.",
            n_snv_indel_no_reads)

    # ── Load funnel summary ──────────────────────────────────────────
    # One consolidated view of where the LENS rows went, replacing the
    # per-stage count lines that double-counted (the non-coord rows
    # "skipped" from variant ranking were the same rows re-counted as
    # "lack gene_name / transcript_id"). The rows feed two independent
    # uses: (1) every row minus mismatches is scored as a candidate
    # epitope for the neoepitope report; (2) rows with genomic coordinates
    # enter variant-based construct ranking, while fusion rows enter through
    # their paired-breakpoint antigen identity. Spelling that out here is what
    # makes the "skipped, then counted again" confusion go away.
    n_construct_source_rows = (
        len(rows) - n_skipped_empty_coords - n_policy_excluded
    )
    n_variants_ranked = sum(1 for _v, vps in ranked if vps)
    kept_kinds = {
        k: (
            full_kinds.get(k, 0)
            - skipped_kinds.get(k, 0)
            - policy_excluded_kinds.get(k.upper(), 0)
        )
        for k in full_kinds
    }
    coord_breakdown = ', '.join(
        "%s=%d" % (k, kept_kinds[k]) for k in sorted(
            kept_kinds, key=lambda k: _antigen_kind_sort_key(k, kept_kinds[k]))
        if kept_kinds[k] > 0) or '(none)'
    if n_snv_indel_no_gene or n_snv_indel_no_transcript or n_snv_indel_no_reads:
        note = (
            "%d missing gene_name, %d missing transcript_id, "
            "%d missing rna_reads" % (
                n_snv_indel_no_gene, n_snv_indel_no_transcript,
                n_snv_indel_no_reads))
    else:
        note = "all carry gene_name / transcript_id / rna_reads"
    funnel = [
        "LENS load funnel: %s" % os.path.basename(report.path),
        "  %d rows in: %s" % (len(rows), breakdown or '(none)'),
        "  → %d candidate epitopes eligible for construct ranking after "
        "epitope DSL filtering and the minimum-score gate" % len(epitopes),
        "  → %d construct-source row(s) (%s) → %d unique source(s) → "
        "%d ranked with vaccine peptide(s) for constructs" % (
            n_construct_source_rows, coord_breakdown, len(groups),
            n_variants_ranked),
    ]
    if n_skipped_empty_coords:
        funnel.append(
            "  → %d non-coord row(s) (%s) are report-only — no genome "
            "placement, so excluded from construct ranking" % (
                n_skipped_empty_coords, skipped_breakdown))
    if n_policy_excluded:
        funnel.append(
            "  → %d row(s) (%s) excluded from construct ranking by "
            "included_antigen_sources policy" % (
                n_policy_excluded, policy_excluded_breakdown
            )
        )
    funnel.append("  note (SNV/INDEL anomalies): %s" % note)
    logger.info("\n".join(funnel))

    # Top candidates first so they win greedy bin-packing in the construct
    # assemblers. The transcript-resolution and varcode-agreement summaries
    # every external format owes its caller are emitted by the accumulator.
    return tally.result("LENS")



def parse_pvacseq_variant(variant_id, genome=None):
    """Parse a pVACseq aggregate ``ID`` field into a ``varcode.Variant``.

    Common forms in the wild:
      - ``chr1-100000-100001-A-T`` (5-part dashed: contig-start-end-ref-alt)
      - ``chr1-100000-A-T`` (4-part dashed)
      - ``1.123.A.T`` (4-part dotted, legacy)

    All recognized forms supply real ref + alt nucleotides. Returns ``None``
    for unrecognized input.
    """
    from varcode import Variant
    if not variant_id:
        return None
    s = str(variant_id)
    contig = pos_s = ref = alt = None
    if '-' in s:
        parts = s.split('-')
        if len(parts) == 5:  # contig-start-end-ref-alt
            contig, pos_s, _, ref, alt = parts
        elif len(parts) == 4:  # contig-start-ref-alt
            contig, pos_s, ref, alt = parts
    if contig is None:
        # Dotted (legacy)
        parts = s.split('.')
        if len(parts) == 4:
            contig, pos_s, ref, alt = parts
    if contig is None:
        return None
    try:
        pos = int(pos_s)
    except (ValueError, TypeError):
        return None
    # Strip the ``chr`` prefix pVACseq emits, same rationale as LENS:
    # pyensembl uses bare contigs and ``normalize_contig_names=False``
    # below means varcode won't strip for us. Keeping the prefix
    # silently breaks downstream ``variant.effect_on_transcript``
    # lookups against the configured pyensembl release.
    if contig.lower().startswith('chr'):
        contig = contig[3:]
    try:
        v = Variant(
            contig=contig, start=pos, ref=ref, alt=alt,
            genome=genome, normalize_contig_names=False)
    except Exception:
        return None
    return v


def pvacseq_rna_depth(row):
    """RNA coverage reported by either pVACseq TSV flavor."""
    return cells.integer(
        cells.first(row, *_PVACSEQ_RNA_DEPTH_COLUMNS), default=0)


def pvacseq_rna_vaf(row):
    """RNA variant allele fraction reported by either pVACseq flavor.

    ``None`` when neither flavor's column carries one. Not 0.0: a fraction
    of zero says the variant was looked for and not seen, and an absent
    fraction says it was not looked for, and the two produce different
    read-count provenance.
    """
    return cells.first_number(row, *_PVACSEQ_RNA_VAF_COLUMNS)


def pvacseq_dna_vaf(row):
    """DNA variant allele fraction reported by either pVACseq flavor."""
    return cells.first_number(row, *_PVACSEQ_DNA_VAF_COLUMNS)


def pvacseq_genomic_variant(variant_id, group_rows, genome=None):
    """Resolve a pVACseq row group to a genomic ``varcode.Variant``.

    The join identity and the genomic variant are *different things* and must
    be derived differently. pVACseq's ``Index`` — which topiary prefers for
    the identity because it is the file's own stable row key — spells a
    variant as ``GENE.ENST….missense.806E/V``, which carries no coordinates
    at all. Deriving the variant from the identity therefore drops every
    all_epitopes file that ships an ``Index`` column, while deriving the
    identity from coordinates makes it disagree with topiary's. Both are
    needed, from their own sources: coordinates when the flavor supplies
    them, the ``ID`` / ``variant`` string (``chr1-100-A-T``) otherwise.
    """
    if genome is None:
        genome = next((cells.text(row.get('input_reference_assembly'))
                       for row in group_rows if cells.text(row.get('input_reference_assembly'))), None)
    for row in group_rows:
        genomic = cells.text(row.get(GENOMIC_VARIANT_COLUMN))
        if not genomic:
            continue
        variant = parse_pvacseq_variant(genomic, genome=genome)
        if variant is not None:
            return variant
    return parse_pvacseq_variant(variant_id, genome=genome)


def pvacseq_variant_metadata(variant_id, group_rows, genome=None):
    """Aggregate filter-independent facts from every raw pVACseq row."""
    variant = pvacseq_genomic_variant(variant_id, group_rows, genome=genome)
    if variant is None:
        # The aggregate warning tells operators to check DEBUG for the
        # offenders, so there has to be something here to find.
        logger.debug(
            "pVACseq group %r carried no genomic coordinates and its "
            "identifier could not be parsed as one; skipping.", variant_id)
        return ExternalVariantEntry(unparseable=True)
    transcript_ids = external_values(*(
        value
        for row in group_rows
        for value in (
            row.get("Best Transcript"),
            row.get("Transcript"),
            row.get("transcript"),
        )
    ))
    transcripts = resolve_external_transcripts(transcript_ids, genome)
    # depth * vaf can only exceed zero when depth already does (vaf is >= 0
    # by construction), so the product added nothing but a second pass.
    has_rna_support = any(
        pvacseq_rna_depth(row) > 0 for row in group_rows)
    gene_names = external_values(*(
        value
        for row in group_rows
        for value in (
            row.get("Gene"), row.get("Gene Name"), row.get("gene"))
    ))
    annotation = check_varcode_annotation(
        variant,
        transcripts[0] if transcripts else None,
        gene_names[0] if gene_names else "",
        variant_is_frameshift(variant),
    ) + (str(variant_id),)
    entry = ExternalVariantEntry(
        variant=variant,
        had_transcript_ids=bool(transcript_ids),
        resolved_transcript=bool(transcripts),
        resolved_protein_context=bool(transcripts),
        annotation=annotation,
        has_rna_support=has_rna_support,
    )
    for row in group_rows:
        dna_vaf = pvacseq_dna_vaf(row)
        if dna_vaf is not None:
            entry.dna_vaf = dna_vaf
            break
    return entry


def pvacseq_vaccine_entry(metadata, selection, genome=None, options=None):
    """Build one pVACseq vaccine peptide from a DSL-selected candidate."""
    if selection is None or metadata.variant is None:
        return metadata
    options = options or ExternalConstructOptions()
    key = selection.representative.key
    context = truncate_at_stop_codon(key.source_sequence or key.peptide)
    if not context or not has_only_standard_amino_acids(context):
        logger.warning(
            "Dropped pVACseq construct for variant %r: antigen sequence %r "
            "is empty or contains non-standard residues.",
            key.variant_id, key.source_sequence)
        return metadata
    # Read the counts topiary derived rather than deriving them here. The
    # LENS path has done this since #374; pVACseq kept a local depth x VAF
    # because topiary published no derived counts for it. topiary 5.47.0
    # does, so the last hand-rolled split goes (#383).
    #
    # The counts and the label describing them still come from one row: a
    # method describing a row whose numbers were not used is the per-field
    # maximum problem one field over.
    observations = [
        ((cells.integer(record.row.get("n_rna_overlapping"), default=0),
          cells.integer(record.row.get("n_rna_alt"), default=0),
          cells.integer(record.row.get("n_rna_ref"), default=0)),
         (cells.text(record.row.get("rna_evidence_method")),
          cells.text(record.row.get("sequence_source")),
          cells.text(record.row.get("rna_evidence_subject"))),
         _dna_evidence_from_row(record.row))
        for record in selection.records]
    (n_total, n_alt_reads, n_ref_reads), provenance, dna_evidence = max(
        observations, key=lambda o: (o[0][1], o[0][0]),
        default=((0, 0, 0), ("", "", ""), {}))
    metadata.vaccine_peptide = external_vaccine_peptide(
        variant=metadata.variant,
        selection=selection,
        # pVACseq ships no SLP context column, so the antigen window is the
        # peptide itself and every residue of it is treated as targetable.
        context=context,
        mutant_start=0,
        mutant_end=len(context),
        gene_name=key.primary_gene_name,
        transcripts=resolve_external_transcripts(
            key.ordered_transcript_ids, genome),
        # The fourth count is reads spanning the assembled protein
        # sequence. pVACseq does not report one, and setting it to
        # n_alt_reads asserted that every variant-supporting read also
        # spans the peptide — a measurement this source never made. Zero
        # says it was not measured; on the LENS path the field holds a
        # genuinely different count (cds_overlap_reads).
        counts=(n_total, n_alt_reads, n_ref_reads, 0),
        provenance=provenance,
        dna_evidence=dna_evidence,
        options=options,
    )
    return metadata


def pvacseq_ranking_result(report, epitopes, genome=None, options=None):
    """Rank pVACseq constructs from an already-parsed report.

    Consumes the topiary-normalized rows of ``report``: both pVACseq flavors
    have already been mapped onto one column vocabulary there, so the raw TSV
    is never re-read and the identities used for ranking are the identities
    the loader produced.

    Uses ``Best Peptide`` / ``MT Epitope Seq`` as the antigen sequence and
    treats the peptide itself as the mutation span (no SLP-context column).
    mRNA construct generation from pVACseq input therefore produces shorter
    antigen windows than the LENS path.
    """
    if report is None or not report.rows:
        return ExternalRankingResult()
    options = options or ExternalConstructOptions()

    rows = list(report.rows)
    groups = {}
    for r in rows:
        variant_id = pvacseq_variant_id(r)
        if not variant_id:
            continue
        groups.setdefault(variant_id, []).append(r)

    records_by_variant = {}
    for record in report.records_with_epitopes(epitopes):
        records_by_variant.setdefault(
            pvacseq_variant_id(record.row), []).append(record)

    tally = ExternalRankingAccumulator(
        require_target_epitopes=options.require_target_epitopes_in_variant)
    for variant_id, group_rows in groups.items():
        tally.add(pvacseq_vaccine_entry(
            pvacseq_variant_metadata(variant_id, group_rows, genome=genome),
            select_external_construct(records_by_variant.get(variant_id, [])),
            genome=genome,
            options=options,
        ))
    if tally.n_unparseable:
        logger.warning(
            "Skipped %d pVACseq row group(s) whose rows carried no genomic "
            "coordinates and whose identifier could not be parsed as one; "
            "see DEBUG log for details.", tally.n_unparseable)
    return tally.result("pVACseq", transcript_id_label="Best Transcript IDs")



def patient_info_from_external(ranked, source_path, patient_id,
                               input_summary,
                               input_label='External report',
                               predictions=None):
    """Build a :class:`PatientInfo` from external-input data.

    Input counts come from the filter-independent summary produced while
    parsing the external report:

      - ``num_somatic_variants`` = unique variants the input file
        produced antigens for (LENS / pVACseq are antigen-only files,
        so this is "variants that survived their pipeline"; silent
        / non-antigenic somatic calls aren't recoverable here)
      - ``num_coding_effect_variants`` = unique sources with a usable protein
        context (legacy field name): a resolved Transcript for point variants,
        or a caller-supplied translated context for fusion antigens
      - ``num_variants_with_rna_support`` = unique variants with at
        least one row carrying a non-zero RNA-read count

    Only the design-output count is derived from ``ranked``:

      - ``num_variants_with_vaccine_peptides`` = ``len(ranked)``,
        same definition as the pipeline path

    ``input_summary`` is required: deriving these counts from ``ranked``
    instead would report post-filter design output as if it were input
    composition, which is the opposite of what the header claims.
    """
    from .patient_info import PatientInfo
    n_with_peptides = sum(1 for _variant, vps in ranked if vps)
    # MHC alleles aren't carried as a separate header in LENS /
    # pVACseq files — they're implicit in the per-row predictions.
    # Infer the unique set from the predictions and mark "(inferred)"
    # so the user knows it came from the file content rather than a
    # header / explicit ``--mhc-alleles`` arg.
    mhc_alleles = []
    if predictions:
        seen = set()
        for ep in predictions:
            # Patient genotype is input provenance, not a side effect of
            # which peptide-allele groups survived the target filter.
            alleles = list(getattr(ep, 'patient_alleles', ()) or ())
            if not alleles:
                alleles = [
                    prediction.allele
                    for prediction in ep.predictions_flat()]
            for allele in alleles:
                if allele and allele not in seen:
                    seen.add(allele)
                    mhc_alleles.append(allele)
        if mhc_alleles:
            mhc_alleles = sorted(mhc_alleles) + [LENS_PROVENANCE_MARKER]
    return PatientInfo(
        patient_id=patient_id or '',
        # Leave the legacy vcf_paths empty — the external path's
        # input is *not* a VCF and labelling it as such was
        # actively misleading. The Inputs block in template reports
        # picks up ``inputs`` instead.
        vcf_paths=[],
        bam_path=None,
        mhc_alleles=mhc_alleles,
        num_somatic_variants=input_summary.num_somatic_variants,
        num_coding_effect_variants=input_summary.num_coding_effect_variants,
        num_variants_with_rna_support=(
            input_summary.num_variants_with_rna_support),
        num_variants_with_vaccine_peptides=n_with_peptides,
        inputs=([(input_label, source_path)] if source_path else []),
    )


def load_external_ranked(args, epitope_config=None, vaccine_config=None,
                         manufacturability_config=None):
    """Dispatch helper: load LENS / pVACseq based on args, return
    ``(ranked, report_df, predictions, patient_info, dna_vaf)`` or ``None``
    when no external input is supplied. Single-file aliases and repeatable
    ``--external-input`` use the same preparation and validation path.

    ``patient_info`` carries the variant-count metadata template reports
    (ASCII / HTML / PDF) need. Input counts are captured before filtering;
    only the vaccine-peptide count follows the filtered design output.

    ``epitope_config`` is evaluated by the selected loader before vaccine
    peptides are constructed, so ranking and template reports consume the
    configured DSL scores rather than the default affinity score.
    ``vaccine_config`` / ``manufacturability_config`` reach construct
    assembly through :class:`ExternalConstructOptions`. Occurrence/window
    selection uses target-epitope scores. The resulting constructs, including
    alternatives from repeated reports, use the common combined-score ordering
    in :func:`vaxrank.ranking.rank_constructs`.

    The returned ``predictions`` collection retains every loaded input group
    for audit reports and patient-genotype inference. ``ranked`` is built from
    separate copies narrowed to groups retained by the Topiary filter and
    meeting the configured minimum epitope score.
    """
    from .external_rescoring import external_inputs, load_unified_external

    if not external_inputs(args):
        return None
    options = ExternalConstructOptions.from_configs(
        vaccine_config=vaccine_config,
        manufacturability_config=manufacturability_config,
        vaccine_peptide_length=getattr(args, 'vaccine_peptide_length', None),
        num_target_epitopes_to_keep=getattr(
            args, 'num_epitopes_per_vaccine_peptide', None),
    )
    genome = getattr(args, 'genome', None)

    return load_unified_external(args, epitope_config, options, genome)
