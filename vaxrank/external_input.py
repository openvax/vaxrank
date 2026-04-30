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

from .mutant_protein_fragment import MutantProteinFragment
from .vaccine_peptide import VaccinePeptide

logger = logging.getLogger(__name__)


def _parse_variant_coords(coords, genome=None):
    """Parse ``chr:pos:ref:alt`` (LENS) into a ``varcode.Variant``.

    Returns None if the string is malformed.
    """
    from varcode import Variant
    if not isinstance(coords, str) or coords.count(":") < 3:
        return None
    parts = coords.split(":")
    contig, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
    try:
        start = int(pos)
    except ValueError:
        return None
    return Variant(contig=contig, start=start, ref=ref, alt=alt,
                   genome=genome, normalize_contig_names=False)


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


def _pick_representative(group_rows):
    """Pick the row whose peptide will represent this variant in the
    construct pipeline.

    Strategy: smallest IC50 (= strongest predicted binder). If the
    LENS file doesn't carry IC50 (mhcflurry-presentation only), fall
    back to lowest percentile rank.
    """
    def _key(r):
        ic50 = r.get('ic50')
        if ic50 is not None and pd.notna(ic50):
            return (0, float(ic50))
        rank = r.get('percentile_rank')
        if rank is not None and pd.notna(rank):
            return (1, float(rank))
        return (2, 0.0)
    return sorted(group_rows, key=_key)[0]


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

    # Index predictions by (peptide, allele) so we can attach all
    # per-predictor predictions to each candidate vaccine peptide.
    by_key = {}
    for p in predictions:
        by_key.setdefault((p.peptide_sequence, p.allele), []).append(p)

    # Group rows by variant. LENS uses lowercase snake_case columns.
    rows = df.to_dict('records')
    groups = {}
    for r in rows:
        coords = r.get('variant_coords')
        if not coords:
            continue
        groups.setdefault(coords, []).append(r)

    ranked = []
    for coords, group_rows in groups.items():
        variant = _parse_variant_coords(coords, genome=genome)
        if variant is None:
            logger.warning(
                "Could not parse LENS variant_coords %r; skipping.", coords)
            continue
        rep = _pick_representative(group_rows)
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
            allele_raw = r.get('allele') or ""
            # load_lens normalizes 'HLA-A02:01' → 'HLA-A*02:01';
            # match the same form here so the (peptide, allele)
            # lookup hits.
            from .epitope_io import _normalize_hla_allele
            allele = _normalize_hla_allele(str(allele_raw))
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

    # Order by mutant_epitope_score descending so the top candidates
    # win greedy bin-packing in the construct assemblers. Ties broken
    # alphabetically by variant coordinates for determinism.
    ranked.sort(
        key=lambda pair: (
            -pair[1][0].mutant_epitope_score if pair[1] else 0.0,
            str(pair[0]),
        ))
    return ranked


def ranked_from_pvacseq_predictions(predictions, report_df, genome=None,
                                    num_mutant_epitopes_to_keep=None):
    """pVACseq variant of :func:`ranked_from_lens_predictions`.

    pVACseq aggregate reports (the format ``load_pvacseq`` parses)
    don't carry an SLP-context column in the same shape; we use
    ``Best Peptide`` as the antigen sequence and treat the peptide
    itself as the mutation span. mRNA construct generation from
    pVACseq input therefore produces shorter antigen windows than
    the LENS path.
    """
    if report_df is None or report_df.empty:
        return []

    by_key = {}
    for p in predictions:
        by_key.setdefault((p.peptide_sequence, p.allele), []).append(p)

    rows = report_df.to_dict('records')
    groups = {}
    for r in rows:
        # pVACseq aggregate uses 'ID' as a stable per-variant key.
        key = r.get('ID') or r.get('Index')
        if key is None:
            continue
        groups.setdefault(key, []).append(r)

    ranked = []
    from varcode import Variant
    for vid, group_rows in groups.items():
        rep = _pick_representative(group_rows)
        best_pep = rep.get('Best Peptide') or rep.get('peptide') or ""
        gene = rep.get('Gene') or 'unknown'
        # pVACseq IDs look like "1.123.A.T" or similar; try to parse
        # for variant location, else fall back to a placeholder.
        contig, pos, ref, alt = '?', 0, 'N', 'N'
        try:
            parts = str(vid).split('.')
            if len(parts) == 4:
                contig, pos_s, ref, alt = parts
                pos = int(pos_s)
        except (TypeError, ValueError):
            pass
        try:
            variant = Variant(
                contig=contig, start=pos, ref=ref, alt=alt,
                genome=genome, normalize_contig_names=False)
        except Exception:
            logger.warning("Could not parse pVACseq ID %r as a Variant; "
                           "skipping.", vid)
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
            allele = r.get('Allele') or r.get('allele') or ""
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
            predictions, report_df,
            genome=getattr(args, 'genome', None))
        return ranked, report_df, predictions
    return None
