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

import logging

import pandas as pd

from .epitope_io import detect_lens_predictors, normalize_hla_allele
from .mutant_protein_fragment import MutantProteinFragment
from .vaccine_library import (
    has_only_standard_amino_acids,
    truncate_at_stop_codon,
)
from .vaccine_peptide import VaccinePeptide

# pVACseq aggregate-report scoring column priority. Sniffed against
# the row dict; the first column present + non-NaN wins.
_PVACSEQ_IC50_COLS = ('IC50 MT', 'Best IC50 Score MT')
_PVACSEQ_RANK_COLS = ('%ile MT', 'Best Percentile MT')

logger = logging.getLogger(__name__)


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


def ranked_from_lens_predictions(predictions, lens_tsv_path, genome=None,
                                 num_mutant_epitopes_to_keep=None):
    """Group LENS predictions by variant; emit ranked vaccine peptides.

    Parameters
    ----------
    predictions : list of EpitopePrediction
        Output of ``load_lens(path)``. Already one prediction per
        (peptide, allele, predictor).
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
    num_mutant_epitopes_to_keep : int, optional
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

    # Index predictions by (peptide, allele) so we can attach all
    # per-predictor predictions to each candidate vaccine peptide.
    by_key = {}
    for p in predictions:
        by_key.setdefault((p.peptide_sequence, p.allele), []).append(p)

    # Group rows by variant. LENS uses lowercase snake_case columns.
    rows = df.to_dict('records')
    groups = {}
    n_skipped_empty_coords = 0
    for r in rows:
        coords = r.get('variant_coords')
        if coords is None or (
                isinstance(coords, float) and pd.isna(coords)) or (
                isinstance(coords, str) and (
                    not coords.strip() or coords.strip().lower() == 'nan')):
            # Non-SNV antigen rows (splice, fusion, intron retention)
            # genuinely lack genome coords — silent skip.
            n_skipped_empty_coords += 1
            continue
        groups.setdefault(coords, []).append(r)
    if n_skipped_empty_coords:
        logger.info(
            "Skipped %d LENS row(s) with empty variant_coords (typical "
            "for non-SNV antigen sources: splice / fusion / intron "
            "retention).", n_skipped_empty_coords)

    ranked = []
    n_unparseable = 0
    # Aggregate transcript-resolution stats so we summarize once at
    # the end instead of spamming a log line per row. ``n_with_ids`` is
    # variants that *had* at least one transcript_id in the LENS file;
    # ``n_resolved`` is the subset where pyensembl actually returned a
    # Transcript. The gap is almost always a release mismatch.
    n_with_ids = 0
    n_resolved = 0
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

        # Collect all predictions for any peptide associated with this
        # variant (not just the representative) so reports / scoring
        # see the full landscape.
        epitope_preds = []
        seen_keys = set()
        for r in group_rows:
            pep = r.get('peptide') or ""
            # load_lens normalizes 'HLA-A02:01' → 'HLA-A*02:01';
            # match the same form here so the (peptide, allele)
            # lookup hits.
            allele = normalize_hla_allele(str(r.get('allele') or ""))
            key = (pep, allele)
            if key in seen_keys:
                continue
            seen_keys.add(key)
            epitope_preds.extend(by_key.get(key, []))

        # Make sure at least one prediction overlaps the mutation so
        # the VaccinePeptide isn't pruned to zero mutant epitopes.
        # LENS predictions are mutant-derived — flag them all.
        for p in epitope_preds:
            if not p.overlaps_mutation:
                p.overlaps_mutation = True

        if not epitope_preds:
            continue

        vp = VaccinePeptide(
            mutant_protein_fragment=fragment,
            epitope_predictions=epitope_preds,
            num_mutant_epitopes_to_keep=num_mutant_epitopes_to_keep,
        )
        ranked.append((variant, [vp]))

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
    n_no_gene_name = sum(
        1 for r in rows
        if not (r.get('gene_name') and not (
            isinstance(r.get('gene_name'), float) and pd.isna(r.get('gene_name')))))
    n_no_transcript = sum(
        1 for r in rows
        if not (r.get('transcript_id') and not (
            isinstance(r.get('transcript_id'), float)
            and pd.isna(r.get('transcript_id')))))
    n_no_reads = sum(
        1 for r in rows
        if (r.get('rna_reads_covering_genomic_origin') is None
            or (isinstance(r.get('rna_reads_covering_genomic_origin'), float)
                and pd.isna(r.get('rna_reads_covering_genomic_origin')))))
    if n_no_pep_context:
        logger.warning(
            "%d / %d LENS row(s) lack pep_context — antigens for those "
            "rows degenerate to the bare neoepitope (no SLP context). "
            "Vaccine windows will be ~9 aa instead of ~25 aa.",
            n_no_pep_context, len(rows))
    if n_no_gene_name:
        logger.info(
            "%d / %d LENS row(s) lack gene_name (typical for ERV / "
            "non-genic antigens).", n_no_gene_name, len(rows))
    if n_no_transcript:
        logger.info(
            "%d / %d LENS row(s) lack transcript_id; downstream "
            "varcode-effect annotation won't have a transcript context "
            "for these.", n_no_transcript, len(rows))
    if n_no_reads:
        logger.warning(
            "%d / %d LENS row(s) lack RNA-read counts "
            "(rna_reads_covering_genomic_origin). Combined-score "
            "ranking will collapse to epitope-only for those rows; "
            "consider --combined-score-mode=epitope_only if RNA data "
            "is consistently absent.", n_no_reads, len(rows))

    # Transcript-resolution summary. Three failure modes worth surfacing
    # at INFO/WARN level — the per-row debug logs in
    # :func:`_resolve_transcripts` capture the individual IDs.
    if n_with_ids and genome is None:
        logger.warning(
            "%d variant(s) had transcript_id(s) in LENS but no "
            "pyensembl genome was configured; ASCII / HTML / PDF "
            "report effect annotations will be empty. Pass "
            "--ensembl-release to populate them.", n_with_ids)
    elif n_with_ids and n_resolved < n_with_ids:
        n_unresolved = n_with_ids - n_resolved
        logger.warning(
            "%d / %d variant(s) had transcript IDs that didn't resolve "
            "against the configured pyensembl release. Most often this "
            "is a release mismatch — pass --ensembl-release N to match "
            "the build LENS used.", n_unresolved, n_with_ids)

    # Order by mutant_epitope_score descending so the top candidates
    # win greedy bin-packing in the construct assemblers. Ties broken
    # alphabetically by variant coordinates for determinism.
    ranked.sort(
        key=lambda pair: (
            -pair[1][0].mutant_epitope_score if pair[1] else 0.0,
            str(pair[0]),
        ))
    return ranked


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


def ranked_from_pvacseq_predictions(predictions, pvacseq_tsv_path,
                                    genome=None,
                                    num_mutant_epitopes_to_keep=None):
    """pVACseq variant of :func:`ranked_from_lens_predictions`.

    Re-reads the raw aggregate TSV (the per-(peptide, allele)
    ``report_df`` returned by ``load_pvacseq`` carries display columns
    rather than the original ``ID`` / ``Best Peptide`` / ``IC50 MT``).

    Uses ``Best Peptide`` as the antigen sequence and treats the
    peptide itself as the mutation span (no SLP-context column).
    mRNA construct generation from pVACseq input therefore produces
    shorter antigen windows than the LENS path.
    """
    if not pvacseq_tsv_path:
        return []
    df = pd.read_csv(pvacseq_tsv_path, sep="\t", low_memory=False)
    if df.empty:
        return []

    by_key = {}
    for p in predictions:
        by_key.setdefault((p.peptide_sequence, p.allele), []).append(p)

    rows = df.to_dict('records')
    groups = {}
    for r in rows:
        key = r.get('ID') or r.get('Index')
        if key is None or pd.isna(key):
            continue
        groups.setdefault(key, []).append(r)

    ranked = []
    n_skipped = 0
    # Aggregate transcript-resolution stats (see LENS path for rationale).
    n_with_ids = 0
    n_resolved = 0
    for vid, group_rows in groups.items():
        rep = _pick_representative(
            group_rows, _PVACSEQ_IC50_COLS, _PVACSEQ_RANK_COLS)
        best_pep = rep.get('Best Peptide') or rep.get('peptide') or ""
        # Preserve pVACseq's own gene name; empty string when missing
        # (codebase convention for "not known"). No 'unknown' invention.
        gene_raw = rep.get('Gene')
        gene = (
            str(gene_raw) if gene_raw and not (
                isinstance(gene_raw, float) and pd.isna(gene_raw))
            else "")

        variant, alleles_real = _parse_pvacseq_id(vid, genome=genome)
        if variant is None:
            n_skipped += 1
            logger.debug(
                "Could not parse pVACseq ID %r as a Variant; skipping.",
                vid)
            continue

        # Real RNA-read counts from pVACseq aggregate columns.
        # ``RNA Depth`` = total coverage at the variant position.
        # ``RNA VAF`` = variant allele frequency (alt / total).
        # n_alt_reads = round(depth × vaf); ref = depth − alt. Missing
        # values yield 0 — honest "no read signal," never fabricated.
        n_total = _coerce_int(rep.get('RNA Depth'), default=0)
        rna_vaf = rep.get('RNA VAF')
        try:
            vaf = float(rna_vaf) if rna_vaf is not None and not (
                isinstance(rna_vaf, float) and pd.isna(rna_vaf)) else 0.0
        except (TypeError, ValueError):
            vaf = 0.0
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
        best_tid = rep.get('Best Transcript')
        transcript_ids = (
            [str(best_tid)]
            if best_tid and not (isinstance(best_tid, float) and pd.isna(best_tid))
            else [])
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
        epitope_preds = []
        seen_keys = set()
        for r in group_rows:
            pep = r.get('Best Peptide') or r.get('peptide') or ""
            allele_raw = r.get('Allele') or r.get('allele') or ""
            allele = normalize_hla_allele(str(allele_raw))
            key = (pep, allele)
            if key in seen_keys:
                continue
            seen_keys.add(key)
            epitope_preds.extend(by_key.get(key, []))
        for p in epitope_preds:
            if not p.overlaps_mutation:
                p.overlaps_mutation = True
        if not epitope_preds:
            continue
        ranked.append((variant, [VaccinePeptide(
            mutant_protein_fragment=fragment,
            epitope_predictions=epitope_preds,
            num_mutant_epitopes_to_keep=num_mutant_epitopes_to_keep,
        )]))

    if n_skipped:
        logger.warning(
            "Skipped %d pVACseq row group(s) with unparseable IDs; "
            "see DEBUG log for details.", n_skipped)

    if n_with_ids and genome is None:
        logger.warning(
            "%d variant(s) had Best Transcript IDs in pVACseq but no "
            "pyensembl genome was configured; ASCII / HTML / PDF "
            "report effect annotations will be empty. Pass "
            "--ensembl-release to populate them.", n_with_ids)
    elif n_with_ids and n_resolved < n_with_ids:
        n_unresolved = n_with_ids - n_resolved
        logger.warning(
            "%d / %d variant(s) had Best Transcript IDs that didn't "
            "resolve against the configured pyensembl release. Most "
            "often this is a release mismatch — pass --ensembl-release "
            "N to match the build pVACseq used.",
            n_unresolved, n_with_ids)

    ranked.sort(
        key=lambda pair: (
            -pair[1][0].mutant_epitope_score if pair[1] else 0.0,
            str(pair[0]),
        ))
    return ranked


def _patient_info_from_external(ranked, source_path, patient_id):
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
    return PatientInfo(
        patient_id=patient_id or '',
        vcf_paths=[source_path] if source_path else [],
        bam_path=None,
        num_somatic_variants=len(ranked),
        num_coding_effect_variants=n_with_transcript,
        num_variants_with_rna_support=n_with_rna,
        num_variants_with_vaccine_peptides=n_with_peptides,
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
    if getattr(args, 'input_lens', None):
        report_df, predictions = load_lens(args.input_lens)
        ranked = ranked_from_lens_predictions(
            predictions, args.input_lens,
            genome=getattr(args, 'genome', None))
        patient_info = _patient_info_from_external(
            ranked, args.input_lens, patient_id)
        return ranked, report_df, predictions, patient_info
    if getattr(args, 'input_pvacseq', None):
        report_df, predictions = load_pvacseq(args.input_pvacseq)
        ranked = ranked_from_pvacseq_predictions(
            predictions, args.input_pvacseq,
            genome=getattr(args, 'genome', None))
        patient_info = _patient_info_from_external(
            ranked, args.input_pvacseq, patient_id)
        return ranked, report_df, predictions, patient_info
    return None
