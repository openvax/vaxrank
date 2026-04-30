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
from .vaccine_peptide import VaccinePeptide

# pVACseq aggregate-report scoring column priority. Sniffed against
# the row dict; the first column present + non-NaN wins.
_PVACSEQ_IC50_COLS = ('IC50 MT', 'Best IC50 Score MT')
_PVACSEQ_RANK_COLS = ('%ile MT', 'Best Percentile MT')

logger = logging.getLogger(__name__)


def _parse_variant_coords(coords, genome=None, ref_alt_fallback=("A", "C")):
    """Parse a LENS ``variant_coords`` cell into a ``(Variant, alleles_real)``
    pair, where ``alleles_real`` is True iff ref + alt came from the
    LENS row (vs. placeholder fallback).

    Real LENS files emit several forms; we accept whichever is present:

      - ``chr1:26780312:C:T``  (4-part: chr, pos, ref, alt — test fixtures
                                 + occasional real rows. ``alleles_real=True``)
      - ``chr1:26780312:C``    (3-part: chr, pos, ref. alt is placeholder.
                                 ``alleles_real=False``)
      - ``chr1:26780312``      (2-part: chr, pos only — dominant in real
                                 LENS v1.9 reports; both ref and alt are
                                 placeholders. ``alleles_real=False``)

    Placeholder nucleotides (``A``/``C`` by default) satisfy varcode's
    validator without claiming a specific genotype. Downstream construct
    assembly keys off the (chr, pos) tuple, not the alleles. **Important**:
    when ``alleles_real=False`` the resulting ``Variant`` should NOT be
    fed to varcode effect annotation or any code path that interprets
    ref/alt as biology — the genotype is fictional. Callers that need
    real ref/alt (e.g. reannotating against a transcript) should drop
    rows where ``alleles_real`` is False and ask the upstream LENS run
    for a more complete report.

    Returns ``(None, False)`` for empty / NaN / malformed input. NaN
    is a normal artifact of non-SNV antigen rows (splice, fusion,
    intron retention) where ``variant_coords`` is genuinely empty —
    callers should treat None as "skip silently," not as a warning
    condition.
    """
    from varcode import Variant
    if coords is None or not isinstance(coords, str):
        return None, False
    s = coords.strip()
    if not s or s.lower() == 'nan':
        return None, False
    parts = s.split(":")
    if len(parts) < 2:
        return None, False
    contig = parts[0]
    try:
        start = int(parts[1])
    except (ValueError, IndexError):
        return None, False
    # varcode rejects 'N' and empty strings as nucleotides. Track
    # whether ref AND alt came from the LENS row so downstream code
    # can detect the placeholder case.
    ref_real = (
        len(parts) >= 3 and parts[2] and parts[2] != 'N')
    alt_real = (
        len(parts) >= 4 and parts[3] and parts[3] != 'N')
    ref = parts[2] if ref_real else ref_alt_fallback[0]
    alt = parts[3] if alt_real else ref_alt_fallback[1]
    if ref == alt:  # varcode rejects no-op variants
        alt = 'C' if ref == 'A' else 'A'
        alt_real = False
    try:
        v = Variant(
            contig=contig, start=start, ref=ref, alt=alt,
            genome=genome, normalize_contig_names=False)
    except Exception:
        return None, False
    return v, (ref_real and alt_real)


def _mut_offsets_in_context(peptide, pep_context):
    """Locate the neoepitope inside its surrounding context.

    Returns ``(start, end)`` AA offsets of the peptide within
    ``pep_context``. If the peptide can't be located, the whole
    context is treated as the mutation span (conservative — keeps the
    centering helper from cropping out the mutation).
    """
    if not pep_context or not peptide:
        return 0, len(pep_context or "")
    idx = pep_context.find(peptide)
    if idx < 0:
        return 0, len(pep_context)
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


def _coerce_n_alt_reads(tpm):
    """LENS doesn't carry alt-read count directly. Use ``ceil(tpm)`` as
    a stand-in so combined-score ranking has a non-degenerate
    expression term. Falls back to 1 when tpm is missing/NaN."""
    if tpm is None or pd.isna(tpm):
        return 1
    try:
        return max(1, int(round(float(tpm))))
    except (TypeError, ValueError):
        return 1


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
    n_placeholder_alleles = 0
    for coords, group_rows in groups.items():
        variant, alleles_real = _parse_variant_coords(coords, genome=genome)
        if variant is None:
            n_unparseable += 1
            logger.debug(
                "Could not parse LENS variant_coords %r; skipping.", coords)
            continue
        if not alleles_real:
            n_placeholder_alleles += 1
        rep = _pick_representative(group_rows, affinity_cols, rank_cols)
        peptide = rep.get('peptide') or ""
        pep_context = rep.get('pep_context') or ""
        if not pep_context:
            # No SLP window available — fall back to using the
            # neoepitope itself as the antigen. The construct will be
            # short but valid.
            pep_context = peptide
        gene_name = rep.get('gene_name') or 'unknown'
        tpm = rep.get('tpm')

        start_off, end_off = _mut_offsets_in_context(peptide, pep_context)
        n_alt = _coerce_n_alt_reads(tpm)

        fragment = MutantProteinFragment(
            variant=variant,
            gene_name=str(gene_name),
            amino_acids=str(pep_context),
            mutant_amino_acid_start_offset=start_off,
            mutant_amino_acid_end_offset=end_off,
            supporting_reference_transcripts=[],
            n_overlapping_reads=n_alt,
            n_alt_reads=n_alt,
            n_ref_reads=0,
            n_alt_reads_supporting_protein_sequence=n_alt,
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
            "Skipped %d LENS variant_coords value(s) that didn't parse "
            "as chr:pos[:ref:alt]; see DEBUG log for the offenders.",
            n_unparseable)
    if n_placeholder_alleles:
        logger.info(
            "%d LENS variant(s) used placeholder ref/alt nucleotides "
            "because the source row only carried chr:pos. The synthesized "
            "MutantProteinFragment.variant.ref / .alt are placeholders "
            "(not real biology) and must NOT be fed to varcode-effect "
            "annotation. Construct assembly only uses (chr, pos) and is "
            "unaffected.", n_placeholder_alleles)

    # Order by mutant_epitope_score descending so the top candidates
    # win greedy bin-packing in the construct assemblers. Ties broken
    # alphabetically by variant coordinates for determinism.
    ranked.sort(
        key=lambda pair: (
            -pair[1][0].mutant_epitope_score if pair[1] else 0.0,
            str(pair[0]),
        ))
    return ranked


def _parse_pvacseq_id(vid):
    """Parse a pVACseq aggregate ``ID`` field into a varcode.Variant.

    Common forms in the wild:
      - ``chr1-100000-100001-A-T`` (5-part dashed: contig-start-end-ref-alt)
      - ``chr1-100000-A-T`` (4-part dashed)
      - ``1.123.A.T`` (4-part dotted, legacy)

    Returns ``(contig, start, ref, alt)`` tuple or ``None`` if the
    string doesn't match any recognized form.
    """
    if not vid:
        return None
    s = str(vid)
    # Try dashed first (modern pVACseq aggregate output).
    if '-' in s:
        parts = s.split('-')
        if len(parts) == 5:  # contig-start-end-ref-alt
            try:
                return parts[0], int(parts[1]), parts[3], parts[4]
            except ValueError:
                return None
        if len(parts) == 4:  # contig-start-ref-alt
            try:
                return parts[0], int(parts[1]), parts[2], parts[3]
            except ValueError:
                return None
    # Dotted (legacy)
    parts = s.split('.')
    if len(parts) == 4:
        try:
            return parts[0], int(parts[1]), parts[2], parts[3]
        except ValueError:
            return None
    return None


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
    from varcode import Variant
    n_skipped = 0
    for vid, group_rows in groups.items():
        rep = _pick_representative(
            group_rows, _PVACSEQ_IC50_COLS, _PVACSEQ_RANK_COLS)
        best_pep = rep.get('Best Peptide') or rep.get('peptide') or ""
        gene = rep.get('Gene') or 'unknown'

        parsed = _parse_pvacseq_id(vid)
        if parsed is None:
            n_skipped += 1
            logger.debug(
                "Could not parse pVACseq ID %r as a Variant; skipping.",
                vid)
            continue
        contig, pos, ref, alt = parsed
        try:
            variant = Variant(
                contig=contig, start=pos, ref=ref, alt=alt,
                genome=genome, normalize_contig_names=False)
        except Exception:
            n_skipped += 1
            logger.debug(
                "Variant construction failed for pVACseq ID %r; skipping.",
                vid, exc_info=True)
            continue

        fragment = MutantProteinFragment(
            variant=variant, gene_name=str(gene),
            amino_acids=str(best_pep),
            mutant_amino_acid_start_offset=0,
            mutant_amino_acid_end_offset=len(best_pep),
            supporting_reference_transcripts=[],
            n_overlapping_reads=1, n_alt_reads=1,
            n_ref_reads=0, n_alt_reads_supporting_protein_sequence=1,
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

    ranked.sort(
        key=lambda pair: (
            -pair[1][0].mutant_epitope_score if pair[1] else 0.0,
            str(pair[0]),
        ))
    return ranked


def load_external_ranked(args):
    """Dispatch helper: load LENS / pVACseq based on args, return
    ``(ranked, report_df, predictions)`` or ``None`` if neither flag
    is set.
    """
    from .epitope_io import load_lens, load_pvacseq
    if getattr(args, 'input_lens', None):
        report_df, predictions = load_lens(args.input_lens)
        ranked = ranked_from_lens_predictions(
            predictions, args.input_lens,
            genome=getattr(args, 'genome', None))
        return ranked, report_df, predictions
    if getattr(args, 'input_pvacseq', None):
        report_df, predictions = load_pvacseq(args.input_pvacseq)
        ranked = ranked_from_pvacseq_predictions(
            predictions, args.input_pvacseq,
            genome=getattr(args, 'genome', None))
        return ranked, report_df, predictions
    return None
