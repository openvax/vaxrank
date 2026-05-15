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

import pandas as pd

from .epitope_io import detect_lens_predictors
from .mutant_protein_fragment import MutantProteinFragment
from .vaccine_library import (
    has_only_standard_amino_acids,
    truncate_at_stop_codon,
)
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


def _parse_variant_coords(coords):
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
    :func:`_variant_from_lens_row`, which builds the actual
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


def _strip_lens_allele(value):
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


def _resolve_transcripts(transcript_ids, genome):
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


def _variant_from_lens_row(row, genome=None):
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
    coords_parsed = _parse_variant_coords(row.get('variant_coords'))
    if coords_parsed is None:
        return None
    contig, pos = coords_parsed
    antigen_source = (row.get('antigen_source') or '').upper()
    if antigen_source == 'SNV':
        ref = _strip_lens_allele(row.get('snv_ref_allele'))
        alt = _strip_lens_allele(row.get('snv_alt_allele'))
    elif antigen_source == 'INDEL':
        ref = _strip_lens_allele(row.get('indel_ref_allele'))
        alt = _strip_lens_allele(row.get('indel_alt_allele'))
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


def _mut_offsets_in_context(peptide, pep_context):
    """Locate the neoepitope inside its surrounding context.

    Returns ``(start, end)`` AA offsets of the peptide within
    ``pep_context``. Returns ``(None, None)`` when the peptide can't
    be located — the caller should drop the row rather than fabricate
    a mutation span. Previously this defaulted to "the whole context
    is the mutation," which falsely told downstream code that every
    residue was mutated.
    """
    if not pep_context or not peptide:
        return None, None
    idx = pep_context.find(peptide)
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
        return []
    df = pd.read_csv(lens_tsv_path, sep="\t", low_memory=False)
    if df.empty:
        return []

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
        # Stable order: SNV / INDEL up-front (the variant-coord paths),
        # then everything else by descending count, with (missing) last.
        priority = {'SNV': 0, 'INDEL': 1}
        def _bucket_order(item):
            k, count = item
            if k.upper() in priority:
                return (priority[k.upper()], 0, k)
            if k == '(missing)':
                return (3, 0, k)
            return (2, -count, k)
        ordered = sorted(full_kinds.items(), key=_bucket_order)
        breakdown = ', '.join("%s=%d" % (k, v) for k, v in ordered)
        logger.info(
            "LENS report contains %d row(s); antigen_source "
            "breakdown: %s.", len(rows), breakdown)

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
    if n_skipped_empty_coords:
        breakdown = ', '.join(
            "%s=%d" % (k, v) for k, v in sorted(
                skipped_kinds.items(), key=lambda kv: -kv[1]))
        logger.info(
            "Skipped %d LENS row(s) with empty variant_coords; "
            "antigen_source breakdown: %s.",
            n_skipped_empty_coords, breakdown)
        # SNV / INDEL rows are *expected* to carry coords; flag them
        # separately so the user can chase upstream rather than assume
        # "typical".
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
    # ``vaf`` column carries DNA VAF in LENS v1.9 (sits between
    # ``mhcflurry_agretopicity`` and ``totcopynum`` / ``multiplicity``
    # / ``ccf`` — all DNA-clonal annotations; values track DNA VAF
    # rather than reads-derived RNA VAF on real Pt02 data, e.g.
    # 0.141 vs computed 0.138). Plumb to ``dna_vaf_by_variant`` so
    # the report's "DNA VAF" field gets populated instead of "n/a".
    dna_vaf_by_variant = {}
    for coords, group_rows in groups.items():
        rep = _pick_representative(group_rows, affinity_cols, rank_cols)
        # Build the Variant from REAL ref/alt columns
        # (snv_ref_allele/snv_alt_allele or indel_ref_allele/indel_alt_allele
        # depending on antigen_source). variant_coords gives only chr:pos in
        # LENS v1.9; the alleles live in dedicated columns.
        variant = _variant_from_lens_row(rep, genome=genome)
        if variant is None:
            n_unparseable += 1
            logger.debug(
                "Could not build Variant from LENS row at coords=%r "
                "(antigen_source=%r); skipping.",
                coords, rep.get('antigen_source'))
            continue
        peptide = rep.get('peptide') or ""
        pep_context = rep.get('pep_context') or ""
        if not pep_context:
            # No SLP window available — fall back to using the
            # neoepitope itself as the antigen. The construct will be
            # short but valid.
            pep_context = peptide
        # Translation halts at the first stop codon — anything after
        # ``*`` doesn't exist as protein. Truncate so manufacturability
        # / hydropathy / codon-optimization see only real residues.
        pep_context = truncate_at_stop_codon(str(pep_context))
        peptide = truncate_at_stop_codon(str(peptide))
        # If the neoepitope itself was past the stop, we can't make a
        # vaccine peptide out of it — drop the row.
        if not pep_context or not peptide:
            n_skipped_post_stop = locals().get('n_skipped_post_stop', 0) + 1
            logger.debug(
                "Dropped LENS row for variant %r: peptide / pep_context "
                "empty after stop-codon truncation.", coords)
            continue
        # Some LENS files emit non-standard residues (selenocysteine
        # 'U', pyrrolysine 'O', ambiguous 'X' / 'B' / 'Z' / 'J').
        # Vaxrank's manufacturability / hydropathy code is keyed off
        # the 20 canonical AAs, so drop the row rather than crash.
        if not has_only_standard_amino_acids(pep_context):
            logger.warning(
                "Dropped LENS row for variant %r: pep_context %r "
                "contains non-standard residues (allowed: 20 canonical "
                "AAs).", coords, pep_context)
            continue
        # Preserve LENS's own gene name; fall back to empty string
        # (the codebase convention for "not known") rather than
        # invent "unknown". Empty propagates through
        # iter_named_antigens which already handles it.
        gene_name_raw = rep.get('gene_name')
        gene_name = (
            str(gene_name_raw) if gene_name_raw and not (
                isinstance(gene_name_raw, float) and pd.isna(gene_name_raw))
            else "")

        start_off, end_off = _mut_offsets_in_context(peptide, pep_context)
        if start_off is None:
            logger.debug(
                "Could not locate peptide %r in pep_context %r for "
                "variant %r; skipping (mutation span unknown).",
                peptide, pep_context, coords)
            continue

        # Window pep_context down to the canonical SLP fragment size
        # (--vaccine-peptide-length, default 25). LENS sometimes
        # emits the entire protein prefix here, which would render
        # as a 100+aa "vaccine peptide" — not what
        # ``MutantProteinFragment`` represents elsewhere in vaxrank.
        # The pipeline path gets correctly-sized fragments from
        # Isovar by construction; the same operation lives on
        # ``MutantProteinFragment.slp_window_around_mutation`` so
        # both paths share one definition of "an SLP-windowed
        # protein fragment around a mutation."
        pep_context, start_off, end_off = (
            MutantProteinFragment.slp_window_around_mutation(
                pep_context, start_off, end_off, vaccine_peptide_length))

        # Real RNA-read counts from LENS columns. Missing → 0
        # (honest "no signal"); never fabricate from TPM.
        (n_total, n_alt_reads, n_ref_reads,
         n_alt_supporting_protein) = _read_counts_from_lens_row(rep)

        # Pull real transcript IDs when available, then resolve to
        # pyensembl ``Transcript`` objects so downstream
        # ``predicted_effect()`` (used by template reports) has real
        # varcode context. When ``genome`` is None or an ID can't be
        # resolved, the resolved list shrinks accordingly — the
        # report renderer tolerates the empty case.
        tid = rep.get('transcript_id')
        all_tids = rep.get('all_transcript_ids_encoding_peptide')
        transcript_ids = []
        if tid and not (isinstance(tid, float) and pd.isna(tid)):
            transcript_ids.append(str(tid))
        if (all_tids and isinstance(all_tids, str)
                and all_tids.lower() != 'nan'):
            for t in all_tids.split(','):
                t = t.strip()
                if t and t not in transcript_ids:
                    transcript_ids.append(t)
        transcripts = _resolve_transcripts(transcript_ids, genome)
        if transcript_ids:
            n_with_ids += 1
            if transcripts:
                n_resolved += 1

        fragment = MutantProteinFragment(
            variant=variant,
            gene_name=gene_name,
            amino_acids=str(pep_context),
            mutant_amino_acid_start_offset=start_off,
            mutant_amino_acid_end_offset=end_off,
            supporting_reference_transcripts=transcripts,
            n_overlapping_reads=n_total,
            n_alt_reads=n_alt_reads,
            n_ref_reads=n_ref_reads,
            n_alt_reads_supporting_protein_sequence=n_alt_supporting_protein,
            # Real LENS ref/alt columns are consulted directly, so
            # the synthesized Variant carries a real biological
            # genotype — placeholder_alleles is False.
            placeholder_alleles=False,
        )

        # Collect all Epitopes whose mutant peptide appears in any of
        # this variant's rows. The load_lens adapter has already
        # grouped all per-(allele, predictor) Prediction records into
        # one CandidateEpitope per (peptide, source, offset), so dedup is on
        # CandidateEpitope identity — typically yields one CandidateEpitope per unique
        # peptide in the variant's row group.
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
                # CandidateEpitope is frozen — flag mutation-overlap via
                # dataclasses.replace. LENS rows are mutant-derived
                # by construction, so every CandidateEpitope reached this way
                # overlaps the mutation.
                if not e.overlaps_mutation:
                    e = dataclasses.replace(e, overlaps_mutation=True)
                variant_epitopes.append(e)

        if not variant_epitopes:
            continue

        vp = VaccinePeptide(
            mutant_protein_fragment=fragment,
            epitopes=variant_epitopes,
            num_target_epitopes_to_keep=num_target_epitopes_to_keep,
        )
        ranked.append((variant, [vp]))

        # Capture DNA VAF for the report's "DNA VAF" field. Take from
        # the representative row; LENS dedupes by variant so the
        # value is consistent across the group.
        vaf_raw = rep.get('vaf')
        if vaf_raw is not None and not (
                isinstance(vaf_raw, float) and pd.isna(vaf_raw)):
            try:
                dna_vaf_by_variant[variant] = float(vaf_raw)
            except (TypeError, ValueError):
                pass

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

    def _kind_breakdown(rows_subset):
        b = {}
        for r in rows_subset:
            k = _kind(r)
            b[k] = b.get(k, 0) + 1
        return ', '.join(
            "%s=%d" % (k, v) for k, v in sorted(
                b.items(), key=lambda kv: -kv[1]))

    rows_no_gene_name = [r for r in rows if _is_missing(r.get('gene_name'))]
    rows_no_transcript = [r for r in rows if _is_missing(r.get('transcript_id'))]
    n_no_gene_name = len(rows_no_gene_name)
    n_no_transcript = len(rows_no_transcript)
    n_no_reads = sum(
        1 for r in rows
        if _is_missing(r.get('rna_reads_covering_genomic_origin')))
    if n_no_pep_context:
        logger.warning(
            "%d / %d LENS row(s) lack pep_context — antigens for those "
            "rows degenerate to the bare neoepitope (no SLP context). "
            "Vaccine windows will be ~9 aa instead of ~25 aa.",
            n_no_pep_context, len(rows))
    # ``gene_name`` and ``transcript_id`` are expected to be empty for
    # ERV / CTA-SELF / SPLICE / FUSION antigens (no canonical gene model
    # to point at). Show the actual antigen_source breakdown rather
    # than asserting "typical for X" — and warn separately if any SNV
    # / INDEL rows lack these (which would be an upstream bug).
    _SNV_OR_INDEL_KINDS = {'SNV', 'INDEL'}
    if n_no_gene_name:
        logger.info(
            "%d / %d LENS row(s) lack gene_name; antigen_source "
            "breakdown: %s.",
            n_no_gene_name, len(rows), _kind_breakdown(rows_no_gene_name))
        anomalous = [r for r in rows_no_gene_name
                     if _kind(r).upper() in _SNV_OR_INDEL_KINDS]
        if anomalous:
            logger.warning(
                "%d SNV / INDEL row(s) lack gene_name — those antigen "
                "kinds are expected to carry one. Likely upstream "
                "LENS bug.", len(anomalous))
    if n_no_transcript:
        logger.info(
            "%d / %d LENS row(s) lack transcript_id; antigen_source "
            "breakdown: %s. Downstream varcode-effect annotation won't "
            "have a transcript context for these.",
            n_no_transcript, len(rows), _kind_breakdown(rows_no_transcript))
        anomalous = [r for r in rows_no_transcript
                     if _kind(r).upper() in _SNV_OR_INDEL_KINDS]
        if anomalous:
            logger.warning(
                "%d SNV / INDEL row(s) lack transcript_id — those "
                "antigen kinds are expected to carry one. Likely "
                "upstream LENS bug.", len(anomalous))
    if n_no_reads:
        rows_no_reads = [
            r for r in rows
            if _is_missing(r.get('rna_reads_covering_genomic_origin'))
        ]
        # Antigen kinds where ``rna_reads_covering_genomic_origin`` is
        # legitimately missing: ERV / SPLICE / FUSION don't have a
        # genome-coord origin in the same sense as SNV / INDEL —
        # those rows are RNA-evidenced via different LENS columns
        # (e.g. ERV expression). Surface SNV / INDEL gaps as the
        # noteworthy anomaly and treat the rest as informational.
        anomalous_no_reads = [
            r for r in rows_no_reads
            if _kind(r).upper() in _SNV_OR_INDEL_KINDS
        ]
        logger.info(
            "%d / %d LENS row(s) lack rna_reads_covering_genomic_origin; "
            "antigen_source breakdown: %s. Those rows score with "
            "n_alt_reads=0 — combined-score ranking falls back to the "
            "epitope-only branch for them.",
            n_no_reads, len(rows), _kind_breakdown(rows_no_reads))
        if anomalous_no_reads:
            logger.warning(
                "%d SNV / INDEL row(s) lack rna_reads_covering_genomic_origin "
                "— those kinds are expected to carry RNA-read counts.",
                len(anomalous_no_reads))

    # Transcript-resolution summary. The "no genome configured" case
    # is now caught pre-flight in ``arg_parser.check_args`` for any
    # run that requests template reports; if we still hit it here
    # we're in a "loaded but reports will be empty" path that's
    # acceptable (e.g., a CSV-only run — no transcript effects
    # rendered anywhere). Only the release-mismatch case is worth
    # logging here.
    if n_with_ids and n_resolved < n_with_ids:
        n_unresolved = n_with_ids - n_resolved
        logger.warning(
            "%d / %d variant(s) had transcript IDs that didn't resolve "
            "against the configured pyensembl release. Most often this "
            "is a release mismatch — pass --ensembl-release N to match "
            "the build LENS used.", n_unresolved, n_with_ids)

    # Order by target_epitope_score descending so the top candidates
    # win greedy bin-packing in the construct assemblers. Ties broken
    # alphabetically by variant coordinates for determinism.
    ranked.sort(
        key=lambda pair: (
            -pair[1][0].target_epitope_score if pair[1] else 0.0,
            str(pair[0]),
        ))
    return ranked, dna_vaf_by_variant


def _parse_pvacseq_id(vid, genome=None):
    """Parse a pVACseq aggregate ``ID`` field into a
    ``(varcode.Variant, alleles_real)`` pair, matching the
    :func:`_parse_variant_coords` shape used on the LENS path.

    Common forms in the wild:
      - ``chr1-100000-100001-A-T`` (5-part dashed: contig-start-end-ref-alt)
      - ``chr1-100000-A-T`` (4-part dashed)
      - ``1.123.A.T`` (4-part dotted, legacy)

    All recognized forms supply real ref + alt nucleotides, so
    ``alleles_real`` is always ``True`` on success. Returns
    ``(None, False)`` for unrecognized input.

    The shape symmetry with :func:`_parse_variant_coords` lets future
    code that handles both LENS and pVACseq paths uniformly read a
    single ``alleles_real`` flag without case analysis.
    """
    from varcode import Variant
    if not vid:
        return None, False
    s = str(vid)
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
        return None, False
    try:
        pos = int(pos_s)
    except (ValueError, TypeError):
        return None, False
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
        return None, False
    return v, True


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
        return []
    df = pd.read_csv(pvacseq_tsv_path, sep="\t", low_memory=False)
    if df.empty:
        return []

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

        variant, alleles_real = _parse_pvacseq_id(vid, genome=genome)
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
        transcripts = _resolve_transcripts(transcript_ids, genome)
        if transcript_ids:
            n_with_ids += 1
            if transcripts:
                n_resolved += 1
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

    # The "no genome configured" case is now caught pre-flight in
    # ``arg_parser.check_args`` for runs that request template
    # reports; only the release-mismatch case is worth logging here.
    if n_with_ids and n_resolved < n_with_ids:
        n_unresolved = n_with_ids - n_resolved
        logger.warning(
            "%d / %d variant(s) had Best Transcript IDs that didn't "
            "resolve against the configured pyensembl release. Most "
            "often this is a release mismatch — pass --ensembl-release "
            "N to match the build pVACseq used.",
            n_unresolved, n_with_ids)

    ranked.sort(
        key=lambda pair: (
            -pair[1][0].target_epitope_score if pair[1] else 0.0,
            str(pair[0]),
        ))
    return ranked, dna_vaf_by_variant


def _patient_info_from_external(ranked, source_path, patient_id,
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
        for p in predictions:
            allele = getattr(p, 'allele', '') or ''
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
    from the ranked output (see ``_patient_info_from_external``).
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
        patient_info = _patient_info_from_external(
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
        patient_info = _patient_info_from_external(
            ranked, args.input_pvacseq, patient_id,
            input_label='pVACseq report',
            predictions=predictions)
        return (ranked, report_df, predictions, patient_info,
                dna_vaf_by_variant)
    return None
