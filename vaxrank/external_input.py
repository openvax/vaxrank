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

"""Synthesize ranked-vaccine-peptide structures from external predictions.

LENS and pVACseq files carry pre-computed (peptide, allele) MHC
predictions plus a peptide-context column (the SLP-style window
surrounding each neoepitope). Vaxrank's downstream code (reports,
peptide constructs, mRNA constructs) consumes a list of
``(varcode.Variant, list[VaccinePeptide])`` tuples — the same shape
the VCF/BAM pipeline emits.

This module bridges the two: it groups external predictions by
variant, picks one representative peptide per variant, and wraps the
context window in a ``MutantProteinFragment`` + ``VaccinePeptide`` so
**the same output dispatch** runs whether the input was a VCF or a
LENS report. Without this bridge the CLI hard-short-circuits external
inputs into a one-line CSV and never reaches the construct writers
(which is precisely the gap PR #253 closes).

Limitations:
- LENS / pVACseq don't carry RNA-supporting-fragment counts in a
  uniform way, so ``n_alt_reads`` is set from ``tpm`` (LENS) or
  defaults to 1. Combined-score ranking that depends on read counts
  collapses to epitope-only effectively.
- The "fragment" recovered from ``pep_context`` is just the SLP
  window; transcript context isn't preserved. Stop-loss / frameshift
  extensions aren't recoverable from external inputs.
"""

import dataclasses
import logging
import os

import pandas as pd

from .amino_acids import has_only_standard_amino_acids
from .epitope_io import detect_lens_predictors
from .mutant_protein_fragment import MutantProteinFragment
from .vaccine_library import truncate_at_stop_codon
from .vaccine_peptide import VaccinePeptide

# pVACseq scoring column priority. Sniffed against the row dict; the first
# column present + non-NaN wins. Covers both aggregate and all_epitopes TSVs.
_PVACSEQ_IC50_COLS = (
    'IC50 MT', 'Best IC50 Score MT',
    'Median MT IC50 Score', 'Best MT IC50 Score',
)
_PVACSEQ_RANK_COLS = (
    '%ile MT', 'Best Percentile MT',
    'Median MT Percentile', 'Best MT Percentile',
)

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


def ranked_sorted_by_target_score(ranked):
    """Order ``(variant, [VaccinePeptide])`` entries by descending
    target-epitope score, ties broken by variant string for
    determinism. Shared by every external loader so the ranking order
    is defined in one place."""
    return sorted(
        ranked,
        key=lambda pair: (
            -pair[1][0].target_epitope_score if pair[1] else 0.0,
            str(pair[0])))


@dataclasses.dataclass
class ExternalVariantEntry:
    """Represent one external variant and its translation status.

    LENS ingestion records the parsed variant and optional vaccine peptide
    alongside the transcript, annotation, VAF, and parsing outcomes that are
    aggregated into the final external-input report.
    """
    variant: object = None
    vaccine_peptide: object = None
    had_transcript_ids: bool = False
    resolved_transcript: bool = False
    annotation: object = None       # (gene_ok, effect_ok, label) or None
    dna_vaf: object = None          # float or None
    unparseable: bool = False


def _missing_cell(value):
    if value is None:
        return True
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def _first_cell(row, *names):
    for name in names:
        value = row.get(name)
        if not _missing_cell(value):
            return value
    return None


def _first_str(row, *names):
    value = _first_cell(row, *names)
    return str(value) if value is not None else ""


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
    import os
    import glob
    import re
    try:
        import platformdirs
        cache_root = platformdirs.user_cache_dir('pyensembl')
    except ImportError:
        # Fall back to the conventional macOS / XDG path so we still
        # work without platformdirs in the import graph.
        cache_root = os.path.expanduser('~/Library/Caches/pyensembl')
        if not os.path.isdir(cache_root):
            cache_root = os.path.expanduser('~/.cache/pyensembl')
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
    coords_parsed = parse_lens_variant_coordinates(row.get('variant_coords'))
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
            genome=genome, normalize_contig_names=False)
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


def _row_score(row, ic50_cols, rank_cols):
    """Numeric score for a row: smallest IC50 (strongest binder) wins,
    falling back to smallest percentile rank when IC50 is missing.

    ``ic50_cols`` / ``rank_cols`` are tuples of column names tried in
    priority order; the first non-NaN match is used. Returns a sortable
    ``(tier, value)`` tuple where tier 0 = real IC50, tier 1 = rank,
    tier 2 = no signal (rows tie at 0.0).
    """
    for col in ic50_cols:
        v = row.get(col)
        if v is not None and pd.notna(v):
            try:
                return (0, float(v))
            except (TypeError, ValueError):
                continue
    for col in rank_cols:
        v = row.get(col)
        if v is not None and pd.notna(v):
            try:
                return (1, float(v))
            except (TypeError, ValueError):
                continue
    return (2, 0.0)


def _pick_representative(group_rows, ic50_cols, rank_cols):
    """Pick the row whose peptide will represent this variant in the
    construct pipeline.

    Strategy: smallest IC50 (= strongest predicted binder). If no
    IC50 column is present, fall back to smallest percentile rank.

    Returns the strongest-binder row from ``group_rows``.
    """
    return min(group_rows, key=lambda r: _row_score(r, ic50_cols, rank_cols))


def _coerce_int(value, default=0):
    """Coerce a value (possibly NaN / None / str) to int.

    Returns ``default`` when the value is genuinely missing — pandas
    NaN, None, or a non-numeric string. Used to read RNA-read-count
    columns from LENS / pVACseq cleanly: when LENS omits the column
    or the row has no RNA support, we want 0 (honest "no signal"),
    not a fabricated stand-in.
    """
    if value is None:
        return default
    try:
        if pd.isna(value):
            return default
    except (TypeError, ValueError):
        pass
    try:
        return int(round(float(value)))
    except (TypeError, ValueError):
        return default


def _read_counts_from_lens_row(row):
    """Extract real RNA-read counts from a LENS row.

    LENS columns:
      - ``rna_reads_covering_genomic_origin``           → total reads at locus
      - ``rna_reads_covering_genomic_origin_with_peptide_cds`` → reads
        whose CDS produces this neoepitope (the alt-supporting count)
      - ``vaf``                                          → variant allele
        frequency (used only as a sanity check, not as a stand-in)

    Returns ``(n_overlapping_reads, n_alt_reads, n_ref_reads,
    n_alt_reads_supporting_protein_sequence)``. Missing columns yield
    0 — honest "no read signal" rather than a fabricated 1.
    """
    n_total = _coerce_int(
        row.get('rna_reads_covering_genomic_origin'), default=0)
    n_alt_cds = _coerce_int(
        row.get('rna_reads_covering_genomic_origin_with_peptide_cds'),
        default=0)
    n_alt_reads = n_alt_cds
    n_ref_reads = max(0, n_total - n_alt_reads)
    n_alt_supporting_protein = n_alt_cds
    return n_total, n_alt_reads, n_ref_reads, n_alt_supporting_protein


def _lens_ranked_entry(coords, group_rows, by_peptide, genome,
                       vaccine_peptide_length, affinity_cols, rank_cols,
                       num_target_epitopes_to_keep):
    """Build one variant's :class:`ExternalVariantEntry` from a LENS
    row group. Returns the entry (with ``unparseable`` / ``vaccine_peptide``
    None for rows that can't be turned into an antigen); the caller
    aggregates the stats and the ranked list."""
    rep = _pick_representative(group_rows, affinity_cols, rank_cols)
    # Build the Variant from REAL ref/alt columns (snv_*_allele /
    # indel_*_allele depending on antigen_source). variant_coords gives
    # only chr:pos in LENS v1.9; the alleles live in dedicated columns.
    variant = variant_from_lens_row(rep, genome=genome)
    if variant is None:
        logger.debug(
            "Could not build Variant from LENS row at coords=%r "
            "(antigen_source=%r); skipping.", coords,
            rep.get('antigen_source'))
        return ExternalVariantEntry(unparseable=True)

    peptide = rep.get('peptide') or ""
    pep_context = rep.get('pep_context') or ""
    if not pep_context:
        # No SLP window available — fall back to the neoepitope itself
        # as the antigen. The construct will be short but valid.
        pep_context = peptide
    # Translation halts at the first stop codon — anything after ``*``
    # doesn't exist as protein. Truncate so manufacturability /
    # hydropathy / codon-optimization see only real residues.
    pep_context = truncate_at_stop_codon(str(pep_context))
    peptide = truncate_at_stop_codon(str(peptide))
    if not pep_context or not peptide:
        logger.debug(
            "Dropped LENS row for variant %r: peptide / pep_context "
            "empty after stop-codon truncation.", coords)
        return ExternalVariantEntry()
    # Some LENS files emit non-standard residues (U / O / X / B / Z / J);
    # vaxrank's manufacturability / hydropathy code is keyed off the 20
    # canonical AAs, so drop the row rather than crash.
    if not has_only_standard_amino_acids(pep_context):
        logger.warning(
            "Dropped LENS row for variant %r: pep_context %r contains "
            "non-standard residues (allowed: 20 canonical AAs).",
            coords, pep_context)
        return ExternalVariantEntry()
    # Preserve LENS's own gene name; fall back to empty string (the
    # codebase convention for "not known").
    gene_name_raw = rep.get('gene_name')
    gene_name = (
        str(gene_name_raw) if gene_name_raw and not (
            isinstance(gene_name_raw, float) and pd.isna(gene_name_raw))
        else "")

    start_off, end_off = peptide_offsets_in_context(peptide, pep_context)
    if start_off is None:
        logger.debug(
            "Could not locate peptide %r in pep_context %r for variant "
            "%r; skipping (mutation span unknown).",
            peptide, pep_context, coords)
        return ExternalVariantEntry()

    # Frameshift → the whole downstream tail is novel; combine the
    # variant's neoepitope rows into the maximal mutant span. Frameshift
    # is derived from the variant's own ref/alt (general); LENS's own
    # frameshift annotation is kept separately only to sanity-check
    # varcode.
    start_off, end_off = maximal_mutant_span(
        start_off, end_off, [r.get('peptide') for r in group_rows],
        pep_context, variant_is_frameshift(variant))
    lens_is_frameshift = (
        str(rep.get('indel_type') or '').strip().lower() == 'frameshift'
        or 'fs' in str(rep.get('variant_effect') or '').lower())

    # Window pep_context to the canonical SLP fragment size; the same
    # operation the pipeline path's Isovar fragments satisfy by
    # construction (see MutantProteinFragment.slp_window_around_mutation).
    pep_context, start_off, end_off = (
        MutantProteinFragment.slp_window_around_mutation(
            pep_context, start_off, end_off, vaccine_peptide_length))

    # Real RNA-read counts from LENS columns. Missing → 0 (honest "no
    # signal"); never fabricated from TPM.
    (n_total, n_alt_reads, n_ref_reads,
     n_alt_supporting_protein) = _read_counts_from_lens_row(rep)

    # Real transcript IDs → pyensembl Transcript objects so template
    # reports have varcode context; shrinks to [] when unresolvable.
    tid = rep.get('transcript_id')
    all_tids = rep.get('all_transcript_ids_encoding_peptide')
    transcript_ids = []
    if tid and not (isinstance(tid, float) and pd.isna(tid)):
        transcript_ids.append(str(tid))
    if all_tids and isinstance(all_tids, str) and all_tids.lower() != 'nan':
        for t in all_tids.split(','):
            t = t.strip()
            if t and t not in transcript_ids:
                transcript_ids.append(t)
    transcripts = resolve_external_transcripts(transcript_ids, genome)

    entry = ExternalVariantEntry(
        variant=variant,
        had_transcript_ids=bool(transcript_ids),
        resolved_transcript=bool(transcripts),
        # Sanity-check LENS's annotation against varcode on LENS's own
        # transcript (validation only — LENS's values are kept).
        annotation=(check_varcode_annotation(
            variant, transcripts[0] if transcripts else None,
            gene_name, lens_is_frameshift) + (str(coords),)))

    fragment = MutantProteinFragment(
        variant=variant, gene_name=gene_name, amino_acids=str(pep_context),
        mutant_amino_acid_start_offset=start_off,
        mutant_amino_acid_end_offset=end_off,
        supporting_reference_transcripts=transcripts,
        n_overlapping_reads=n_total, n_alt_reads=n_alt_reads,
        n_ref_reads=n_ref_reads,
        n_alt_reads_supporting_protein_sequence=n_alt_supporting_protein,
        # Real ref/alt columns → real genotype, not a placeholder.
        placeholder_alleles=False)

    # Collect the variant's CandidateEpitopes (deduped on identity); LENS
    # rows are mutant-derived, so every one overlaps the mutation.
    seen_peptides = set()
    seen_epitope_ids = set()
    variant_epitopes = []
    for r in group_rows:
        pep = r.get('peptide') or ""
        if not pep or pep in seen_peptides:
            continue
        seen_peptides.add(pep)
        for e in by_peptide.get(pep, []):
            if id(e) in seen_epitope_ids:
                continue
            seen_epitope_ids.add(id(e))
            if not e.overlaps_mutation:
                e = dataclasses.replace(e, overlaps_mutation=True)
            variant_epitopes.append(e)
    if not variant_epitopes:
        # Transcript stats + annotation already recorded above; just no
        # vaccine peptide for this variant.
        return entry

    entry.vaccine_peptide = VaccinePeptide(
        mutant_protein_fragment=fragment,
        epitopes=variant_epitopes,
        num_target_epitopes_to_keep=num_target_epitopes_to_keep)

    # DNA VAF for the report (representative row; LENS dedupes by variant).
    vaf_raw = rep.get('vaf')
    if vaf_raw is not None and not (
            isinstance(vaf_raw, float) and pd.isna(vaf_raw)):
        try:
            entry.dna_vaf = float(vaf_raw)
        except (TypeError, ValueError):
            pass
    return entry


def ranked_from_lens_predictions(epitopes, lens_tsv_path, genome=None,
                                 num_target_epitopes_to_keep=None,
                                 vaccine_peptide_length=25):
    """Group LENS epitopes by variant; emit ranked vaccine peptides.

    Parameters
    ----------
    epitopes : list of CandidateEpitope
        Output of ``load_lens(path)``. Each CandidateEpitope groups all
        per-(allele, predictor) ``mhctools.Prediction`` records for
        one ``(peptide, source_sequence, offset)`` position.
    lens_tsv_path : str
        Path to the original LENS TSV. We re-read the raw rows here
        because the per-(peptide, allele) ``report_df`` returned by
        ``load_lens`` doesn't carry ``pep_context`` / ``variant_coords``
        in their original lower-snake form — they get rendered into
        display columns. The raw TSV preserves them.
    genome : varcode-compatible genome reference, optional
        Passed through to ``Variant`` construction. When None,
        Variants are constructed without genome resolution and
        downstream code that needs gene annotation may degrade.
    num_target_epitopes_to_keep : int, optional
        Forwarded to ``VaccinePeptide``. When None, all overlapping
        epitopes are kept.

    Returns
    -------
    list[(varcode.Variant, list[VaccinePeptide])]
    """
    if not lens_tsv_path:
        return [], {}
    df = pd.read_csv(lens_tsv_path, sep="\t", low_memory=False)
    if df.empty:
        return [], {}

    # Auto-detect which LENS predictor columns are present so the
    # representative-row pick reads real values. Without this, every
    # row scored 0.0 and stable-sort returned the file-order first
    # row (typically the alphabetically-first allele) instead of the
    # strongest binder.
    detected = detect_lens_predictors(df.columns)
    affinity_cols = tuple(
        d.value_col for d in detected if d.kind == 'pMHC_affinity')
    rank_cols = tuple(
        d.percentile_col for d in detected
        if d.percentile_col and d.kind == 'pMHC_affinity')
    if not affinity_cols and not rank_cols:
        # No pMHC_affinity predictor detected — typically a
        # presentation-only LENS file. _row_score will tie every row
        # at (2, 0.0) and the file-order first row "wins". Warn so
        # the user knows the representative pick is arbitrary.
        kinds = sorted({d.kind for d in detected}) or ['(none)']
        logger.warning(
            "LENS file has no pMHC_affinity scoring columns "
            "(detected kinds: %s); per-variant representative-peptide "
            "pick will fall back to file order. Consider running with "
            "a predictor that emits pMHC affinity (e.g. mhcflurry, "
            "netmhcpan-ba).", kinds)

    # Index epitopes by their mutant peptide sequence so we can attach
    # the right CandidateEpitope(s) to each candidate vaccine peptide. One
    # peptide can map to multiple Epitopes when the same sequence
    # appears with different source contexts / offsets across
    # variants — we'll filter by (peptide, allele) presence below.
    by_peptide: dict[str, list] = {}
    for e in epitopes:
        by_peptide.setdefault(e.sequence, []).append(e)

    # Group rows by variant. LENS uses lowercase snake_case columns.
    rows = df.to_dict('records')

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
    # When a row has no variant_coords, the only sensible explanation
    # is a non-SNV / non-INDEL antigen kind (splice / fusion / ERV /
    # CTA-self / intron-retention). Verify that hypothesis instead of
    # asserting it — if any SNV / INDEL rows are missing coords, that's
    # an upstream bug worth surfacing distinctly.
    skipped_kinds = {}  # antigen_source value → count
    for r in rows:
        coords = r.get('variant_coords')
        if coords is None or (
                isinstance(coords, float) and pd.isna(coords)) or (
                isinstance(coords, str) and (
                    not coords.strip() or coords.strip().lower() == 'nan')):
            n_skipped_empty_coords += 1
            kind = r.get('antigen_source')
            kind_key = (
                str(kind).strip() if kind is not None and not (
                    isinstance(kind, float) and pd.isna(kind))
                else '(missing)')
            skipped_kinds[kind_key] = skipped_kinds.get(kind_key, 0) + 1
            continue
        groups.setdefault(coords, []).append(r)
    skipped_breakdown = ', '.join(
        "%s=%d" % (k, v) for k, v in sorted(
            skipped_kinds.items(), key=lambda kv: -kv[1]))
    if n_skipped_empty_coords:
        # SNV / INDEL rows are *expected* to carry coords; flag them
        # separately so the user can chase upstream rather than assume
        # "typical". (The non-coord rows aren't dropped from the run —
        # they're still scored as candidate epitopes for the neoepitope
        # report; they just can't be placed on the genome for the
        # variant-based vaccine-construct ranking. The funnel below
        # spells this out.)
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

    ranked = []
    n_unparseable = 0
    # Aggregate transcript-resolution stats so we summarize once at
    # the end instead of spamming a log line per row. ``n_with_ids`` is
    # variants that *had* at least one transcript_id in the LENS file;
    # ``n_resolved`` is the subset where pyensembl actually returned a
    # Transcript. The gap is almost always a release mismatch.
    n_with_ids = 0
    n_resolved = 0
    # Sanity-check varcode against LENS's own annotation. Collected
    # per variant, summarized once by log_varcode_agreement.
    annotation_results = []
    # ``vaf`` column carries DNA VAF in LENS v1.9 (sits between
    # ``mhcflurry_agretopicity`` and ``totcopynum`` / ``multiplicity``
    # / ``ccf`` — all DNA-clonal annotations; values track DNA VAF
    # rather than reads-derived RNA VAF on real Pt02 data, e.g.
    # 0.141 vs computed 0.138). Plumb to ``dna_vaf_by_variant`` so
    # the report's "DNA VAF" field gets populated instead of "n/a".
    dna_vaf_by_variant = {}
    for coords, group_rows in groups.items():
        entry = _lens_ranked_entry(
            coords, group_rows, by_peptide, genome, vaccine_peptide_length,
            affinity_cols, rank_cols, num_target_epitopes_to_keep)
        if entry.unparseable:
            n_unparseable += 1
            continue
        if entry.annotation is not None:
            annotation_results.append(entry.annotation)
        if entry.had_transcript_ids:
            n_with_ids += 1
            if entry.resolved_transcript:
                n_resolved += 1
        if entry.dna_vaf is not None:
            dna_vaf_by_variant[entry.variant] = entry.dna_vaf
        if entry.vaccine_peptide is not None:
            ranked.append((entry.variant, [entry.vaccine_peptide]))

    if n_unparseable:
        logger.warning(
            "Skipped %d LENS variant(s): variant_coords couldn't be "
            "parsed as chr:pos OR the row's snv_*_allele / "
            "indel_*_allele columns were missing. See DEBUG log for "
            "the offenders.", n_unparseable)

    # Per-load summary of missing essential / important per-row data.
    # Essential = blocked the row (already counted as "skipped" above
    # or in load_lens itself). Important-but-recoverable = warns so
    # users know what's degraded.
    n_no_pep_context = sum(
        1 for r in rows
        if not (r.get('pep_context') and not (
            isinstance(r.get('pep_context'), float) and pd.isna(r.get('pep_context')))))
    def _is_missing(v):
        return v is None or v == '' or (
            isinstance(v, float) and pd.isna(v)) or (
            isinstance(v, str) and v.strip().lower() == 'nan')

    def _kind(r):
        k = r.get('antigen_source')
        return ('(missing)' if _is_missing(k) else str(k).strip())

    rows_no_gene_name = [r for r in rows if _is_missing(r.get('gene_name'))]
    rows_no_transcript = [r for r in rows if _is_missing(r.get('transcript_id'))]
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
        if _is_missing(r.get('rna_reads_covering_genomic_origin'))]
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
            "%d SNV / INDEL row(s) lack rna_reads_covering_genomic_origin "
            "— those kinds are expected to carry RNA-read counts.",
            n_snv_indel_no_reads)

    # Transcript-resolution summary. The "no genome configured" case
    # is now caught pre-flight in ``arg_parser.check_args`` for any
    # run that requests template reports; if we still hit it here
    # we're in a "loaded but reports will be empty" path that's
    # acceptable (e.g., a CSV-only run — no transcript effects
    # rendered anywhere). Only the release-mismatch case is worth
    # logging here.
    log_transcript_resolution(n_with_ids, n_resolved, "LENS")
    log_varcode_agreement(annotation_results, "LENS")

    # ── Load funnel summary ──────────────────────────────────────────
    # One consolidated view of where the LENS rows went, replacing the
    # per-stage count lines that double-counted (the non-coord rows
    # "skipped" from variant ranking were the same rows re-counted as
    # "lack gene_name / transcript_id"). The rows feed two independent
    # uses: (1) every row minus mismatches is scored as a candidate
    # epitope for the neoepitope report; (2) only rows with genomic
    # variant_coords can be placed on the genome for variant-based
    # vaccine-construct ranking. Spelling that out here is what makes
    # the "skipped, then counted again" confusion go away.
    n_coord_rows = len(rows) - n_skipped_empty_coords
    n_variants_ranked = sum(1 for _v, vps in ranked if vps)
    kept_kinds = {
        k: full_kinds.get(k, 0) - skipped_kinds.get(k, 0) for k in full_kinds}
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
        "LENS load funnel: %s" % os.path.basename(lens_tsv_path),
        "  %d rows in: %s" % (len(rows), breakdown or '(none)'),
        "  → %d candidate epitopes scored for the neoepitope report "
        "(all antigen sources)" % len(epitopes),
        "  → %d row(s) with genomic variant_coords (%s) → %d unique "
        "variant(s) → %d ranked with vaccine peptide(s) for constructs" % (
            n_coord_rows, coord_breakdown, len(groups), n_variants_ranked),
    ]
    if n_skipped_empty_coords:
        funnel.append(
            "  → %d non-coord row(s) (%s) are report-only — no genome "
            "placement, so excluded from construct ranking" % (
                n_skipped_empty_coords, skipped_breakdown))
    funnel.append("  note (SNV/INDEL anomalies): %s" % note)
    logger.info("\n".join(funnel))

    # Top candidates first so they win greedy bin-packing in the
    # construct assemblers (shared ordering helper).
    return ranked_sorted_by_target_score(ranked), dna_vaf_by_variant


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


def _pvacseq_variant_key(row):
    """Variant key for aggregated or all_epitopes pVACseq rows."""
    existing = _first_cell(row, 'ID')
    if existing is not None:
        return existing
    chrom = _first_cell(row, 'Chromosome')
    start = _first_cell(row, 'Start')
    ref = _first_cell(row, 'Reference')
    alt = _first_cell(row, 'Variant')
    if None in (chrom, start, ref, alt):
        index = _first_cell(row, 'Index')
        return index
    stop = _first_cell(row, 'Stop')
    if stop is not None:
        return f"{chrom}-{start}-{stop}-{ref}-{alt}"
    return f"{chrom}-{start}-{ref}-{alt}"


def _pvacseq_source_key(row):
    """Source key matching load_pvacseq's CandidateEpitope source."""
    existing = _first_cell(row, 'ID')
    if existing is not None:
        return str(existing)
    index = _first_cell(row, 'Index')
    if index is not None:
        return str(index)
    chrom = _first_cell(row, 'Chromosome')
    start = _first_cell(row, 'Start')
    ref = _first_cell(row, 'Reference')
    alt = _first_cell(row, 'Variant')
    if None in (chrom, start, ref, alt):
        return ""
    return f"{chrom}-{start}-{ref}-{alt}"


def _pvacseq_peptide(row):
    return _first_str(row, 'Best Peptide', 'MT Epitope Seq', 'peptide')


def _pvacseq_gene(row):
    return _first_str(row, 'Gene', 'Gene Name', 'gene')


def _pvacseq_transcript_ids(row):
    transcript = _first_str(row, 'Best Transcript', 'Transcript', 'transcript')
    return [transcript] if transcript else []


def _pvacseq_rna_depth(row):
    return _coerce_int(_first_cell(row, 'RNA Depth', 'Tumor RNA Depth'), default=0)


def _pvacseq_rna_vaf(row):
    raw = _first_cell(row, 'RNA VAF', 'Tumor RNA VAF')
    try:
        return float(raw) if raw is not None else 0.0
    except (TypeError, ValueError):
        return 0.0


def _pvacseq_dna_vaf(row):
    raw = _first_cell(row, 'DNA VAF', 'Tumor DNA VAF')
    try:
        return float(raw) if raw is not None else None
    except (TypeError, ValueError):
        return None


def ranked_from_pvacseq_predictions(epitopes, pvacseq_tsv_path,
                                    genome=None,
                                    num_target_epitopes_to_keep=None):
    """pVACseq variant of :func:`ranked_from_lens_predictions`.

    Re-reads the raw pVACseq TSV (the per-(peptide, allele)
    ``report_df`` returned by ``load_pvacseq`` carries display columns
    rather than every original aggregate / all_epitopes column).

    Uses ``Best Peptide`` / ``MT Epitope Seq`` as the antigen sequence
    and treats the peptide itself as the mutation span (no SLP-context column).
    mRNA construct generation from pVACseq input therefore produces
    shorter antigen windows than the LENS path.
    """
    if not pvacseq_tsv_path:
        return [], {}
    df = pd.read_csv(pvacseq_tsv_path, sep="\t", low_memory=False)
    if df.empty:
        return [], {}

    by_source_peptide: dict[tuple[str, str], list] = {}
    for e in epitopes:
        by_source_peptide.setdefault((e.source_sequence, e.sequence), []).append(e)

    rows = df.to_dict('records')
    groups = {}
    for r in rows:
        key = _pvacseq_variant_key(r)
        if _missing_cell(key):
            continue
        groups.setdefault(key, []).append(r)

    ranked = []
    n_skipped = 0
    # Aggregate transcript-resolution stats (see LENS path for rationale).
    n_with_ids = 0
    n_resolved = 0
    annotation_results = []
    # pVACseq carries an explicit ``DNA VAF`` column (separate from
    # the ``RNA VAF`` we already consume for read counts). Plumb it
    # through to the report's DNA-VAF field.
    dna_vaf_by_variant = {}
    for vid, group_rows in groups.items():
        rep = _pick_representative(
            group_rows, _PVACSEQ_IC50_COLS, _PVACSEQ_RANK_COLS)
        best_pep = _pvacseq_peptide(rep)
        # Preserve pVACseq's own gene name; empty string when missing
        # (codebase convention for "not known"). No 'unknown' invention.
        gene = _pvacseq_gene(rep)

        variant = parse_pvacseq_variant(vid, genome=genome)
        if variant is None:
            n_skipped += 1
            logger.debug(
                "Could not parse pVACseq ID %r as a Variant; skipping.",
                vid)
            continue

        # Real RNA-read counts from pVACseq columns.
        # ``RNA Depth`` / ``Tumor RNA Depth`` = total coverage at the
        # variant position. ``RNA VAF`` / ``Tumor RNA VAF`` = variant
        # allele frequency (alt / total).
        # n_alt_reads = round(depth × vaf); ref = depth − alt. Missing
        # values yield 0 — honest "no read signal," never fabricated.
        n_total = _pvacseq_rna_depth(rep)
        vaf = _pvacseq_rna_vaf(rep)
        n_alt_reads = int(round(n_total * vaf)) if n_total > 0 else 0
        n_ref_reads = max(0, n_total - n_alt_reads)

        # pVACseq aggregate's "Best Peptide" is the minimal-epitope
        # neoantigen, not an SLP context. The whole peptide IS the
        # mutation span (the aggregate doesn't split mutant vs.
        # flanking residues), so claiming offsets [0, len(best_pep))
        # is structurally honest — a transcript-aware loader could
        # narrow this further but the aggregate doesn't carry the
        # data to do so.
        # pVACseq's aggregate carries one ``Best Transcript`` column
        # per row (the transcript backing the chosen Best Peptide).
        # Resolve to a pyensembl ``Transcript`` so template reports
        # have real effect context; falls back to [] when the genome
        # isn't plumbed through or the ID can't be resolved.
        transcript_ids = _pvacseq_transcript_ids(rep)
        transcripts = resolve_external_transcripts(transcript_ids, genome)
        if transcript_ids:
            n_with_ids += 1
            if transcripts:
                n_resolved += 1
        # Sanity-check varcode against the variant. pVACseq's aggregate
        # carries no explicit effect annotation, so the "provided"
        # frameshift comes from the variant's own ref/alt (provider
        # data) — this still cross-checks that varcode's effect class
        # is consistent with the genomic alleles, and that the gene
        # agrees. Shared with the LENS path via check_varcode_annotation
        # + log_varcode_agreement.
        gene_ok, effect_ok = check_varcode_annotation(
            variant, transcripts[0] if transcripts else None,
            gene, variant_is_frameshift(variant))
        annotation_results.append((gene_ok, effect_ok, str(vid)))
        fragment = MutantProteinFragment(
            variant=variant, gene_name=gene,
            amino_acids=str(best_pep),
            mutant_amino_acid_start_offset=0,
            mutant_amino_acid_end_offset=len(best_pep),
            supporting_reference_transcripts=transcripts,
            n_overlapping_reads=n_total, n_alt_reads=n_alt_reads,
            n_ref_reads=n_ref_reads,
            n_alt_reads_supporting_protein_sequence=n_alt_reads,
        )
        seen_peptides = set()
        seen_epitope_ids = set()
        variant_epitopes = []
        for r in group_rows:
            pep = _pvacseq_peptide(r)
            if not pep or pep in seen_peptides:
                continue
            seen_peptides.add(pep)
            source_key = _pvacseq_source_key(r)
            for e in by_source_peptide.get((source_key, pep), []):
                if id(e) in seen_epitope_ids:
                    continue
                seen_epitope_ids.add(id(e))
                if not e.overlaps_mutation:
                    e = dataclasses.replace(e, overlaps_mutation=True)
                variant_epitopes.append(e)
        if not variant_epitopes:
            continue
        ranked.append((variant, [VaccinePeptide(
            mutant_protein_fragment=fragment,
            epitopes=variant_epitopes,
            num_target_epitopes_to_keep=num_target_epitopes_to_keep,
        )]))

        dna_vaf = _pvacseq_dna_vaf(rep)
        if dna_vaf is not None:
            dna_vaf_by_variant[variant] = dna_vaf

    if n_skipped:
        logger.warning(
            "Skipped %d pVACseq row group(s) with unparseable IDs; "
            "see DEBUG log for details.", n_skipped)

    log_transcript_resolution(
        n_with_ids, n_resolved, "pVACseq", id_label="Best Transcript IDs")
    log_varcode_agreement(annotation_results, "pVACseq")
    return ranked_sorted_by_target_score(ranked), dna_vaf_by_variant


def patient_info_from_external(ranked, source_path, patient_id,
                               input_label='External report',
                               predictions=None):
    """Build a :class:`PatientInfo` from external-input data.

    Counts are derived from the ranked output:
      - ``num_somatic_variants`` = unique variants the input file
        produced antigens for (LENS / pVACseq are antigen-only files,
        so this is "variants that survived their pipeline"; silent
        / non-antigenic somatic calls aren't recoverable here)
      - ``num_coding_effect_variants`` = unique variants whose
        representative fragment resolved at least one Transcript
        (the rest are ERV / non-genic / unresolvable IDs)
      - ``num_variants_with_rna_support`` = unique variants with at
        least one row carrying a non-zero RNA-read count
      - ``num_variants_with_vaccine_peptides`` = ``len(ranked)``,
        same definition as the pipeline path
    """
    from .patient_info import PatientInfo
    n_with_transcript = 0
    n_with_rna = 0
    n_with_peptides = 0
    for _variant, vps in ranked:
        if not vps:
            continue
        n_with_peptides += 1
        # ``vps[0]`` is sufficient on the LENS / pVACseq paths because
        # both loaders emit exactly one VaccinePeptide per variant
        # (the ``MutantProteinFragment`` is built from the
        # representative row's pep_context). If a future loader emits
        # multiple VPs per variant, this should iterate ``vps``.
        frag = vps[0].mutant_protein_fragment
        if frag.supporting_reference_transcripts:
            n_with_transcript += 1
        if (frag.n_alt_reads or 0) > 0 or (frag.n_overlapping_reads or 0) > 0:
            n_with_rna += 1
    # MHC alleles aren't carried as a separate header in LENS /
    # pVACseq files — they're implicit in the per-row predictions.
    # Infer the unique set from the predictions and mark "(inferred)"
    # so the user knows it came from the file content rather than a
    # header / explicit ``--mhc-alleles`` arg.
    mhc_alleles = []
    if predictions:
        seen = set()
        for ep in predictions:
            # ``predictions`` here is a list of CandidateEpitope, which
            # has no ``.allele`` — its alleles are the keys of
            # ``per_allele_scores`` (populated by the DSL at load).
            # Fall back to the raw Prediction.allele list for any loader
            # that doesn't populate per_allele_scores.
            alleles = list((getattr(ep, 'per_allele_scores', None) or {}))
            if not alleles:
                alleles = [
                    getattr(p, 'allele', '')
                    for p in (getattr(ep, 'predictions', None) or [])]
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
        num_somatic_variants=len(ranked),
        num_coding_effect_variants=n_with_transcript,
        num_variants_with_rna_support=n_with_rna,
        num_variants_with_vaccine_peptides=n_with_peptides,
        inputs=([(input_label, source_path)] if source_path else []),
    )


def load_external_ranked(args):
    """Dispatch helper: load LENS / pVACseq based on args, return
    ``(ranked, report_df, predictions, patient_info)`` or ``None`` if
    neither flag is set.

    ``patient_info`` carries the variant-count metadata template
    reports (ASCII / HTML / PDF) need; counts are proxies derived
    from the ranked output (see :func:`patient_info_from_external`).
    """
    from .epitope_io import load_lens, load_pvacseq
    patient_id = getattr(args, 'output_patient_id', '') or ''
    vaccine_peptide_length = getattr(args, 'vaccine_peptide_length', None) or 25
    if getattr(args, 'input_lens', None):
        report_df, predictions = load_lens(args.input_lens)
        ranked, dna_vaf_by_variant = ranked_from_lens_predictions(
            predictions, args.input_lens,
            genome=getattr(args, 'genome', None),
            vaccine_peptide_length=vaccine_peptide_length)
        patient_info = patient_info_from_external(
            ranked, args.input_lens, patient_id,
            input_label='LENS report',
            predictions=predictions)
        return (ranked, report_df, predictions, patient_info,
                dna_vaf_by_variant)
    if getattr(args, 'input_pvacseq', None):
        report_df, predictions = load_pvacseq(args.input_pvacseq)
        ranked, dna_vaf_by_variant = ranked_from_pvacseq_predictions(
            predictions, args.input_pvacseq,
            genome=getattr(args, 'genome', None))
        patient_info = patient_info_from_external(
            ranked, args.input_pvacseq, patient_id,
            input_label='pVACseq report',
            predictions=predictions)
        return (ranked, report_df, predictions, patient_info,
                dna_vaf_by_variant)
    return None
