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

"""
Serialization, deserialization, and import of epitope predictions.

Supports:
  - vaxrank native (CSV/TSV roundtrip of one row per leaf prediction)
  - pVACseq aggregated TSV (all_epitopes.aggregated.tsv)
  - LENS report TSV

All loader functions emit ``list[vaxrank.candidate_epitope.CandidateEpitope]``;
the flat per-(peptide, allele, predictor) row shape lives only inside
the CSV / DataFrame layer for round-trip fidelity. Loaders feed each
row to :func:`~vaxrank.candidate_epitope.candidate_epitopes_from_rows`,
which groups by ``(peptide, source, offset)`` into ``CandidateEpitope``
objects at the end.
"""

import logging
import os
import re
from dataclasses import dataclass

import pandas as pd
from mhctools.pred import Prediction

from .epitope_dsl import _kind_for_method
from .candidate_epitope import (
    SOURCE_CLASS_MUTATION, candidate_epitopes_from_rows,
)

logger = logging.getLogger(__name__)


# ── Vaxrank native format ────────────────────────────────────────────────────

VAXRANK_COLUMNS = [
    "allele",
    "peptide_sequence",
    "wt_peptide_sequence",
    "ic50",
    "wt_ic50",
    "percentile_rank",
    "prediction_method_name",
    "source_class",
    "overlaps_mutation",
    "source_sequence",
    "offset",
    "occurs_in_reference",
    "occurs_in_non_CTA_reference",
    "predictor_version",
]


def _epitope_to_rows(epitope):
    """Flatten one ``CandidateEpitope`` to per-(predictor, allele) row dicts
    matching :data:`VAXRANK_COLUMNS`. Drives ``predictions_to_dataframe``
    and ``save_predictions``."""
    wt = epitope.wt
    wt_peptide_sequence = wt.sequence if wt is not None else ""
    # Build (allele, predictor_name) -> wt_ic50 lookup so each per-allele
    # candidate row gets paired with its matching WT IC50, if any.
    wt_ic50_by_key = {}
    if wt is not None:
        for p in wt.predictions_flat():
            if p.kind == 'pMHC_affinity' and p.value is not None:
                wt_ic50_by_key.setdefault((p.allele, p.predictor_name), p.value)
    rows = []
    for p in epitope.predictions_flat():
        if p.kind != 'pMHC_affinity':
            # The vaxrank-native CSV format only carries affinity rows;
            # presentation / stability / processing live elsewhere.
            continue
        wt_ic50 = wt_ic50_by_key.get((p.allele, p.predictor_name))
        rows.append({
            "allele": p.allele,
            "peptide_sequence": epitope.sequence,
            "wt_peptide_sequence": wt_peptide_sequence,
            "ic50": p.value,
            "wt_ic50": wt_ic50,
            "percentile_rank": p.percentile_rank,
            "prediction_method_name": p.predictor_name,
            "source_class": epitope.source_class,
            "overlaps_mutation": epitope.overlaps_mutation,
            "source_sequence": epitope.source_sequence,
            "offset": epitope.offset,
            "occurs_in_reference": epitope.occurs_in_reference,
            "occurs_in_non_CTA_reference": epitope.occurs_in_non_CTA_reference,
            "predictor_version": p.predictor_version,
        })
    return rows


def predictions_to_dataframe(epitopes):
    """Convert a collection of ``CandidateEpitope`` objects to a flat DataFrame
    matching :data:`VAXRANK_COLUMNS` — one row per leaf
    ``pMHC_affinity`` prediction.

    Parameters
    ----------
    epitopes : list of CandidateEpitope, dict, or iterable
        Most callers pass a plain ``list[CandidateEpitope]``. A ``dict`` is
        accepted for back-compat with pre-3.0 call sites that keyed
        an ordered map by ``(peptide, allele)`` — the values are used.
    """
    if isinstance(epitopes, dict):
        epitopes = list(epitopes.values())
    rows = []
    for e in epitopes:
        rows.extend(_epitope_to_rows(e))
    return pd.DataFrame(rows, columns=VAXRANK_COLUMNS)


def save_predictions(epitopes, path):
    """Save ``CandidateEpitope`` objects to a CSV/TSV file. One row per
    ``(epitope, allele, predictor)`` leaf prediction.

    Format is inferred from the file extension (.tsv -> tab, else comma).
    """
    df = predictions_to_dataframe(epitopes)
    sep = "\t" if str(path).endswith(".tsv") else ","
    df.to_csv(path, sep=sep, index=False)
    logger.info("Saved %d prediction row(s) to %s", len(df), path)


def load_predictions(path):
    """Load a vaxrank-native CSV/TSV file and return ``list[CandidateEpitope]``."""
    sep = "\t" if str(path).endswith(".tsv") else ","
    df = pd.read_csv(path, sep=sep)

    def _str_or_empty(val):
        if val is None or (isinstance(val, float) and pd.isna(val)):
            return ""
        return str(val)

    rows = []
    for _, row in df.iterrows():
        wt_ic50 = row.get("wt_ic50")
        if pd.isna(wt_ic50):
            wt_ic50 = None
        percentile_rank = row.get("percentile_rank")
        if pd.isna(percentile_rank):
            percentile_rank = None
        peptide = row["peptide_sequence"]
        wt_peptide = _str_or_empty(row.get("wt_peptide_sequence", ""))
        method = _str_or_empty(row.get("prediction_method_name", ""))
        version = _str_or_empty(row.get("predictor_version", ""))
        allele = row["allele"]
        kind = _kind_for_method(method)
        mutant = Prediction(
            kind=kind, predictor_name=method, predictor_version=version,
            allele=allele, peptide=peptide, value=float(row["ic50"]),
            score=0.0,
            percentile_rank=(
                float(percentile_rank) if percentile_rank is not None
                else None),
        )
        wt = None
        if wt_ic50 is not None:
            wt = Prediction(
                kind=kind, predictor_name=method,
                predictor_version=version, allele=allele,
                peptide=wt_peptide, value=float(wt_ic50), score=0.0,
                percentile_rank=None,
            )
        occurs_in_ref = bool(row.get("occurs_in_reference", False))
        rows.append({
            'peptide': peptide,
            'source': _str_or_empty(row.get("source_sequence", "")),
            'offset': int(row.get("offset", 0)),
            'mutant': mutant,
            'wt': wt,
            'source_class': _str_or_empty(row.get("source_class", "")) or None,
            'overlaps_mutation': bool(row.get("overlaps_mutation", True)),
            'occurs_in_reference': occurs_in_ref,
            # Fall back to raw occurs_in_reference when the CSV is
            # from a pre-3.0 producer that didn't emit this column.
            'occurs_in_non_CTA_reference': bool(
                row.get("occurs_in_non_CTA_reference", occurs_in_ref)),
        })
    return candidate_epitopes_from_rows(rows)


# ── pVACseq import ───────────────────────────────────────────────────────────

def _is_missing(val):
    """True for scalar missing values from pandas / pVACseq TSVs."""
    if val is None:
        return True
    try:
        return bool(pd.isna(val))
    except (TypeError, ValueError):
        return False


def _format_nm(value):
    # Missing affinity → blank, consistent with every other missing
    # value in the neoepitope CSVs (see _build_lens_report_row). A
    # blank cell in a "Predicted … affinity" column reads as "not
    # predicted" without mixing a string sentinel into numeric data.
    value = _safe_float(value)
    return '%.2f nM' % value if value is not None else ''


# Canonical core of the per-(peptide, allele) neoepitope table. All
# three input paths (VCF pipeline, LENS, pVACseq) build their report
# rows on top of this so the shared columns — and the missing-value
# convention — are produced by one piece of code and can't drift
# (this is where the old "No prediction" vs blank inconsistency came
# from: three hand-rolled row dicts). Affinity args are raw nM floats
# (or None → blank). Callers append their source-specific columns.
NEOEPITOPE_CORE_COLUMNS = (
    'Allele', 'Mutant peptide sequence', 'Score',
    'Predicted mutant pMHC affinity', 'Wildtype sequence',
    'Predicted wildtype pMHC affinity', 'Gene name', 'Genomic variant')


def neoepitope_core_row(allele, mutant_peptide, mutant_affinity,
                        wt_peptide, wt_affinity, gene_name, variant,
                        score=None):
    """Build the shared core columns of one neoepitope-report row.

    ``score`` is included only when provided (the pipeline path
    surfaces it inline; the import paths add ``vaxrank_score`` in
    :func:`write_neoepitope_report`). Affinities are raw nM floats or
    None; formatting + the missing-value rule live in ``_format_nm``.
    """
    row = {
        'Allele': allele or '',
        'Mutant peptide sequence': mutant_peptide or '',
    }
    if score is not None:
        row['Score'] = score
    row['Predicted mutant pMHC affinity'] = _format_nm(mutant_affinity)
    row['Wildtype sequence'] = wt_peptide or ''
    row['Predicted wildtype pMHC affinity'] = _format_nm(wt_affinity)
    row['Gene name'] = gene_name or ''
    row['Genomic variant'] = variant or ''
    return row


def _first_float(*values):
    for value in values:
        coerced = _safe_float(value)
        if coerced is not None:
            return coerced
    return None


def _topiary_group_source(row):
    """Source key matching topiary.ranking.EvalContext grouping."""
    return (
        _safe_str(row.get("variant"))
        or _safe_str(row.get("source_sequence_name"))
        or _safe_str(row.get("peptide"))
    )


def _topiary_pvacseq_to_epitope_rows(df):
    """Convert topiary.read_pvacseq long-form rows to CandidateEpitope rows."""
    epitope_rows = []
    n_rows = 0
    for _, row in df.iterrows():
        peptide = _safe_str(row.get("peptide"))
        allele = _safe_str(row.get("allele"))
        value = _safe_float(row.get("value"))
        if not peptide or not allele or value is None:
            continue

        method = _safe_str(row.get("prediction_method_name")) or "pvacseq"
        version = _safe_str(row.get("predictor_version"))
        kind = _safe_str(row.get("kind")) or _kind_for_method(method)
        percentile_rank = _safe_float(row.get("percentile_rank"))
        score = _safe_float(row.get("score"))
        mutant = Prediction(
            kind=kind, predictor_name=method, predictor_version=version,
            allele=allele, peptide=peptide, value=value,
            score=score if score is not None else 0.0,
            percentile_rank=percentile_rank,
        )

        wt = None
        wt_value = _safe_float(row.get("wt_value"))
        if wt_value is not None:
            wt_method = _safe_str(row.get("wt_prediction_method_name")) or method
            wt_version = _safe_str(row.get("wt_predictor_version")) or version
            wt = Prediction(
                kind=kind, predictor_name=wt_method,
                predictor_version=wt_version, allele=allele,
                peptide=_safe_str(row.get("wt_peptide")),
                value=wt_value, score=0.0,
                percentile_rank=_safe_float(row.get("wt_percentile_rank")),
            )

        epitope_rows.append({
            'peptide': peptide,
            'source': _topiary_group_source(row),
            'offset': int(_safe_float(row.get("peptide_offset")) or 0),
            'mutant': mutant,
            'wt': wt,
            'source_class': SOURCE_CLASS_MUTATION,
            'overlaps_mutation': _safe_bool(
                row.get("contains_mutant_residues"), default=True),
            'occurs_in_reference': _safe_bool(row.get("pvacseq_ref_match")),
            'occurs_in_non_CTA_reference': _safe_bool(
                row.get("pvacseq_ref_match")),
        })
        n_rows += 1
    return epitope_rows, n_rows


def _build_pvacseq_report_row(row):
    """Build one user-facing report row from a topiary pVACseq row."""
    rna_expr = _first_float(
        row.get("rna_transcript_expression"),
        row.get("transcript_expression"),
        row.get("gene_expression"))
    rna_vaf = _first_float(row.get("rna_vaf"), row.get("tumor_rna_vaf"))
    dna_vaf = _first_float(row.get("dna_vaf"), row.get("tumor_dna_vaf"))
    out = neoepitope_core_row(
        allele=_safe_str(row.get("allele")),
        mutant_peptide=_safe_str(row.get("peptide")),
        mutant_affinity=row.get("value"),
        wt_peptide=_safe_str(row.get("wt_peptide")),
        wt_affinity=row.get("wt_value"),
        gene_name=_safe_str(row.get("gene")),
        variant=_safe_str(row.get("variant")))
    out.update({
        'Tier': _safe_str(row.get("pvacseq_tier")),
        'Ref Match': _safe_bool(row.get("pvacseq_ref_match")),
        'RNA Expr': rna_expr,
        'RNA VAF': rna_vaf,
        'DNA VAF': dna_vaf,
        '%ile MT': _safe_float(row.get("percentile_rank")),
        'Source sequence name': _topiary_group_source(row),
        'Peptide offset': int(_safe_float(row.get("peptide_offset")) or 0),
        'MHC class': _safe_str(row.get("mhc_class")),
        'Contains mutant residues': _safe_bool(
            row.get("contains_mutant_residues"), default=True),
    })
    for col, value in row.items():
        if col.startswith("pvacseq_") and col not in {
                "pvacseq_ref_match", "pvacseq_tier"}:
            out[col] = _safe_float(value)
    return out


def load_pvacseq(path):
    """
    Import a pVACseq TSV and return a neoepitope report DataFrame ready
    for output.

    Both pVACseq output flavors are accepted:
    ``*all_epitopes.aggregated.tsv`` and ``*all_epitopes.tsv``.

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    pandas.DataFrame
        Columns: Allele, Mutant peptide sequence, Predicted mutant pMHC
        affinity, Wildtype sequence, Predicted wildtype pMHC affinity,
        Gene name, Genomic variant, Tier, Ref Match, RNA Expr, DNA VAF,
        vaxrank_score
    list of CandidateEpitope
        One ``vaxrank.candidate_epitope.CandidateEpitope`` per unique
        ``(peptide, source_sequence, offset)`` group.
    """
    from topiary import read_pvacseq

    from .epitope_dsl import drop_empty_sample_name
    result = read_pvacseq(path)
    topiary_df = drop_empty_sample_name(result.df)
    epitope_rows, n_rows = _topiary_pvacseq_to_epitope_rows(topiary_df)
    report_rows = [
        _build_pvacseq_report_row(row)
        for _, row in topiary_df.iterrows()
        if _safe_str(row.get("peptide")) and _safe_str(row.get("allele"))
    ]

    report_df = pd.DataFrame(report_rows) if report_rows else pd.DataFrame()
    report_df.attrs["topiary_df"] = topiary_df
    epitopes = candidate_epitopes_from_rows(epitope_rows)
    from .epitope_dsl import attach_per_allele_scores
    epitopes = attach_per_allele_scores(epitopes, topiary_df=topiary_df)
    logger.info(
        "Loaded %d epitope(s) (%d row(s), %s flavor) from pVACseq file %s",
        len(epitopes), n_rows, result.extra.get("pvacseq_format"), path)
    return report_df, epitopes


# ── Shared cell-coercion helpers ─────────────────────────────────────────────

def _safe_str(val):
    """Coerce a cell to str, mapping NaN/None/'NA' to empty string."""
    if _is_missing(val):
        return ""
    s = str(val)
    if s == "NA":
        return ""
    return s


def _safe_float(val):
    """Convert to float, returning None for NaN/missing."""
    if _is_missing(val):
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def _safe_bool(val, default=False):
    """Coerce common TSV boolean spellings; missing values use *default*."""
    if _is_missing(val):
        return default
    if isinstance(val, str):
        lowered = val.strip().lower()
        if lowered in {"true", "t", "yes", "y", "1"}:
            return True
        if lowered in {"false", "f", "no", "n", "0"}:
            return False
        return default
    return bool(val)


# ── LENS import ──────────────────────────────────────────────────────────────

# LENS columns follow the convention "<tool>_<version>.<metric>", e.g.
# "mhcflurry_2.1.1.aff", "netmhcpan_4.1b.perc_rank_el",
# "netmhcstabpan_1.0.halflife_hours". We sniff which tools are present by
# regex, bucket columns by (tool, version), and consult the registry below
# to know which metric is the "value" and which are percentile ranks.
_LENS_COLUMN_RE = re.compile(r"^([A-Za-z][A-Za-z0-9]*)_(\d[\w.]*)\.(.+)$")

# Registry describing how each LENS-emitted predictor maps into the Topiary
# DSL. ``kind`` is the topiary prediction Kind (pMHC_affinity etc.).
# ``value_cols`` and ``percentile_cols`` are LENS metric suffixes in
# preference order — the first present metric wins. MHCflurry's integrated
# model contributes two additional kinds: allele-specific presentation and
# allele-independent antigen processing. Their columns stay attached to the
# same detected tool/version so the loader can emit all three canonical
# mhctools prediction kinds.
_LENS_PREDICTOR_REGISTRY = {
    "mhcflurry": {
        "kind": "pMHC_affinity",
        "value_cols": ("aff",),
        "percentile_cols": ("aff_perc",),
        "presentation_score_cols": ("pres_score",),
        "presentation_percentile_cols": ("pres_perc",),
        "processing_score_cols": ("proc_score",),
    },
    "netmhcpan": {
        "kind": "pMHC_affinity",
        "value_cols": ("aff_nm",),
        "percentile_cols": ("perc_rank_el", "perc_rank_ba"),
    },
    "netmhcstabpan": {
        "kind": "pMHC_stability",
        "value_cols": ("halflife_hours",),
        "percentile_cols": ("perc_rank_stab",),
    },
}

@dataclass(slots=True)
class DetectedPredictor:
    """A predictor whose columns were found in a LENS TSV.

    ``value_col`` is the LENS column carrying the predicted value
    (IC50 for affinity, half-life hours for stability). The other
    ``*_col`` attributes are ``None`` when that signal isn't emitted
    for the detected ``(tool, version)``. MHCflurry presentation and
    processing columns are modeled separately rather than overloading
    its affinity percentile.
    """
    tool: str
    version: str
    kind: str
    value_col: str
    percentile_col: str | None = None
    agretopicity_col: str | None = None
    presentation_score_col: str | None = None
    presentation_percentile_col: str | None = None
    processing_score_col: str | None = None


def detect_lens_predictors(columns):
    """Return a list of :class:`DetectedPredictor` for every known tool present.

    Unknown tools (not in ``_LENS_PREDICTOR_REGISTRY``) are ignored; we'd have
    no way to know which of their columns is the affinity / percentile / kind.
    """
    # Group columns by (tool, version)
    buckets = {}
    for col in columns:
        match = _LENS_COLUMN_RE.match(col)
        if not match:
            continue
        tool, version, metric = match.group(1), match.group(2), match.group(3)
        buckets.setdefault((tool, version), {})[metric] = col

    detected = []
    for (tool, version), metrics in buckets.items():
        spec = _LENS_PREDICTOR_REGISTRY.get(tool)
        if spec is None:
            logger.debug(
                "Skipping unknown LENS predictor %r (version %s)", tool, version)
            continue
        value_col = next(
            (metrics[m] for m in spec["value_cols"] if m in metrics), None)
        if value_col is None:
            # A tool prefix was present but without its canonical value metric
            # — treat as not emitted rather than raising.
            continue
        percentile_col = next(
            (metrics[m] for m in spec["percentile_cols"] if m in metrics), None)
        presentation_score_col = next(
            (metrics[m] for m in spec.get("presentation_score_cols", ())
             if m in metrics), None)
        presentation_percentile_col = next(
            (metrics[m] for m in spec.get(
                "presentation_percentile_cols", ()) if m in metrics), None)
        processing_score_col = next(
            (metrics[m] for m in spec.get("processing_score_cols", ())
             if m in metrics), None)
        # Agretopicity columns aren't versioned (just "<tool>_agretopicity"),
        # so we look them up directly on the full column set.
        agretopicity_col = f"{tool}_agretopicity"
        if agretopicity_col not in columns:
            agretopicity_col = None
        detected.append(DetectedPredictor(
            tool=tool,
            version=version,
            kind=spec["kind"],
            value_col=value_col,
            percentile_col=percentile_col,
            agretopicity_col=agretopicity_col,
            presentation_score_col=presentation_score_col,
            presentation_percentile_col=presentation_percentile_col,
            processing_score_col=processing_score_col,
        ))
    # Deterministic order: tool name alphabetical
    detected.sort(key=lambda d: (d.tool, d.version))
    return detected


def _pick_canonical_predictor(affinity_preds):
    """Pick the canonical pMHC_affinity :class:`DetectedPredictor` from a list.

    Used to populate the report's display columns when multiple affinity
    predictors are present. Preference order: ``mhcflurry`` (present
    across all LENS versions), then ``netmhcpan``, then alphabetical on
    tool name for any remainder.
    """
    if not affinity_preds:
        return None
    for preferred in ("mhcflurry", "netmhcpan"):
        for d in affinity_preds:
            if d.tool == preferred:
                return d
    return sorted(affinity_preds, key=lambda d: d.tool)[0]


def normalize_hla_allele(allele):
    """LENS emits alleles as 'HLA-A01:01'; vaxrank output uses 'HLA-A*01:01'."""
    if not allele:
        return allele
    return re.sub(r"^(HLA-[A-Z]{1,3})(\d)", r"\1*\2", allele)


def load_lens(path):
    """
    Import a LENS report TSV and return a neoepitope report DataFrame
    plus a list of ``vaxrank.candidate_epitope.CandidateEpitope`` objects.

    LENS column schemas vary across versions (v1.4, v1.5, v1.9+); we sniff
    which predictors a given file actually emits rather than requiring a
    fixed schema. Internally we build one flat per-(peptide, allele,
    detected predictor, prediction kind) record per row, then group them
    into ``CandidateEpitope`` objects keyed by
    ``(peptide, source_sequence, offset)`` so Topiary
    DSL expressions can combine multi-predictor predictions via
    ``affinity['mhcflurry']`` / ``affinity['netmhcpan']`` or use an
    unqualified ``affinity`` fallback via the epitope config's
    ``default_methods`` setting.

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    pandas.DataFrame
        Report-ready DataFrame (one row per peptide × allele).
    list of CandidateEpitope
        One ``CandidateEpitope`` per ``(peptide, source_sequence, offset)``
        group, each carrying its per-(allele, predictor, kind)
        ``mhctools.Prediction`` records inside ``epitope``. MHCflurry
        inputs may carry distinct ``pMHC_affinity``, ``pMHC_presentation``,
        and allele-independent ``antigen_processing`` leaves.
    """
    # low_memory=False avoids mixed-dtype warnings — LENS reports have many
    # optional columns that are mostly NA for a given antigen type.
    df = pd.read_csv(path, sep="\t", low_memory=False)

    required = {"peptide", "allele"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"LENS file {path} missing required columns: {missing}")

    detected = detect_lens_predictors(df.columns)
    if not detected:
        raise ValueError(
            f"LENS file {path} has no recognized predictor columns "
            f"(expected columns like 'mhcflurry_<version>.aff' or "
            f"'netmhcpan_<version>.aff_nm')")

    logger.info(
        "LENS file %s: detected predictors=%s",
        path, [(d.tool, d.version) for d in detected])

    # Report columns like 'Predicted mutant pMHC affinity' and '%ile rank'
    # render a single predictor's value. Pick the canonical pMHC_affinity
    # predictor for that display role (every per-tool raw value is also
    # exposed in its own column further down). Stability / presentation
    # predictors are never chosen for display.
    affinity_preds = [d for d in detected if d.kind == "pMHC_affinity"]
    display_pred = (
        _pick_canonical_predictor(affinity_preds) if affinity_preds else None)
    chosen = list(detected)

    epitope_rows = []
    report_rows = []
    processing_scores_by_position = {}
    # Rows whose peptide isn't a substring of a non-empty pep_context —
    # peptide and pep_context came from different isoforms / annotation
    # snapshots upstream. Dropped here (not carried into the report,
    # constructs, or pepsickle); summarized once after the loop.
    n_dropped_peptide_context_mismatch = 0
    mismatch_examples = []
    for _, row in df.iterrows():
        peptide = row.get("peptide", "")
        if not peptide or pd.isna(peptide):
            continue

        allele_raw = row.get("allele", "")
        if pd.isna(allele_raw) or not allele_raw:
            continue
        allele = normalize_hla_allele(str(allele_raw))

        # mhcflurry_agretopicity = MT_IC50 / WT_IC50 in LENS's
        # convention (small values ≪ 1 = mutation strengthened
        # binding). WT_IC50 = MT_IC50 / agretopicity is the only
        # way to recover wildtype affinity from a LENS row, since
        # LENS doesn't emit a WT IC50 column directly. Computed
        # once per row (shared across detected predictors).
        agretopicity = _safe_float(row.get("mhcflurry_agretopicity"))
        # Locate the neoepitope inside its surrounding context once per
        # row (pep_context is a row-level column, not per-predictor).
        # LENS centers the peptide within pep_context but doesn't emit
        # the offset directly; substring search recovers it. When
        # pep_context is non-empty but doesn't contain the peptide, the
        # two were built from different isoforms / annotation snapshots
        # (an upstream LENS issue) — drop the row here rather than carry
        # a fabricated offset downstream into the report, constructs,
        # and pepsickle. An empty pep_context is a different case (no
        # SLP window at all): keep it, offset 0, bare-neoepitope fallback.
        pep_context = _safe_str(row.get("pep_context"))
        if pep_context and str(peptide) not in pep_context:
            n_dropped_peptide_context_mismatch += 1
            if len(mismatch_examples) < 1:
                mismatch_examples.append((str(peptide), pep_context))
            continue
        offset = pep_context.find(str(peptide)) if pep_context else 0
        row_added = False
        shared_epitope_fields = {
            'peptide': str(peptide),
            'source': pep_context,
            'offset': offset,
            'source_class': SOURCE_CLASS_MUTATION,
            # LENS pre-curates rows as neoepitopes and pre-filters
            # against the patient reference proteome — both flags
            # are structurally true for any surviving row.
            'overlaps_mutation': True,
            'occurs_in_reference': False,
            'occurs_in_non_CTA_reference': False,
        }
        for d in chosen:
            value = _safe_float(row.get(d.value_col))
            percentile_rank = (
                _safe_float(row.get(d.percentile_col))
                if d.percentile_col else None)
            # Derive WT IC50 from agretopicity, only for the
            # mhcflurry-affinity predictor (the value matches that
            # tool's IC50 scale). Other predictors leave wt_ic50=None.
            if value is not None:
                wt_ic50 = None
                if (d.tool == "mhcflurry" and d.kind == "pMHC_affinity"
                        and agretopicity is not None and agretopicity > 0):
                    wt_ic50 = value / agretopicity
                mutant = Prediction(
                    kind=d.kind, predictor_name=d.tool,
                    predictor_version=d.version, allele=allele,
                    peptide=str(peptide), value=value, score=0.0,
                    percentile_rank=percentile_rank,
                )
                wt = None
                if wt_ic50 is not None:
                    # LENS doesn't emit a WT peptide column — pass an
                    # empty string. ``candidate_epitopes_from_rows`` keeps
                    # the WT context (anonymous-WT signal) because wt_ic50
                    # is non-None even though the sequence is unknown.
                    wt = Prediction(
                        kind=d.kind, predictor_name=d.tool,
                        predictor_version=d.version, allele=allele,
                        peptide="", value=wt_ic50, score=0.0,
                        percentile_rank=None,
                    )
                epitope_rows.append({
                    **shared_epitope_fields,
                    'mutant': mutant,
                    'wt': wt,
                })
                row_added = True

            presentation_score = (
                _safe_float(row.get(d.presentation_score_col))
                if d.presentation_score_col else None)
            presentation_percentile = (
                _safe_float(row.get(d.presentation_percentile_col))
                if d.presentation_percentile_col else None)
            if (presentation_score is not None
                    or presentation_percentile is not None):
                # mhctools represents MHCflurry presentation as its own
                # allele-specific kind. ``Prediction.score`` is required;
                # retain a percentile-only legacy row with a neutral 0.0
                # score rather than discarding the percentile evidence.
                epitope_rows.append({
                    **shared_epitope_fields,
                    'mutant': Prediction(
                        kind="pMHC_presentation",
                        predictor_name=d.tool,
                        predictor_version=d.version,
                        allele=allele,
                        peptide=str(peptide),
                        value=None,
                        score=(
                            presentation_score
                            if presentation_score is not None else 0.0),
                        percentile_rank=presentation_percentile,
                    ),
                    'wt': None,
                })
                row_added = True

            processing_score = (
                _safe_float(row.get(d.processing_score_col))
                if d.processing_score_col else None)
            if processing_score is not None:
                # MHCflurry processing depends on peptide + flanks, not MHC.
                # LENS repeats it on every allele row, while mhctools emits
                # one allele-less prediction per peptide position. Preserve
                # that canonical shape and reject inconsistent repeats.
                processing_key = (
                    str(peptide), pep_context, offset, d.tool, d.version)
                if processing_key in processing_scores_by_position:
                    if processing_scores_by_position[
                            processing_key] != processing_score:
                        raise ValueError(
                            "Conflicting LENS processing scores for one "
                            "peptide position and predictor")
                else:
                    processing_scores_by_position[
                        processing_key] = processing_score
                    epitope_rows.append({
                        **shared_epitope_fields,
                        'mutant': Prediction(
                            kind="antigen_processing",
                            predictor_name=d.tool,
                            predictor_version=d.version,
                            allele="",
                            peptide=str(peptide),
                            value=None,
                            score=processing_score,
                            percentile_rank=None,
                        ),
                        'wt': None,
                    })
                row_added = True

        if not row_added:
            continue

        report_rows.append(_build_lens_report_row(
            row=row,
            allele=allele,
            peptide=peptide,
            source_sequence_name=pep_context or str(peptide),
            peptide_offset=offset,
            detected=detected,
            display_pred=display_pred,
            chosen_tools=[d.tool for d in chosen],
        ))

    if n_dropped_peptide_context_mismatch:
        ex_peptide, ex_context = mismatch_examples[0]
        logger.warning(
            "Dropped %d LENS row(s) whose peptide isn't a substring of "
            "its (non-empty) pep_context — peptide and pep_context were "
            "built from different isoforms / annotation snapshots "
            "upstream. These are excluded from the report and "
            "constructs. Example: peptide=%r not in pep_context=%r.",
            n_dropped_peptide_context_mismatch, ex_peptide, ex_context)

    report_df = pd.DataFrame(report_rows) if report_rows else pd.DataFrame()
    epitopes = candidate_epitopes_from_rows(epitope_rows)
    # See load_pvacseq for the parallel rationale: the 3.1 refactor moved
    # per-(peptide, allele) scoring to the DSL and stored the result on
    # CandidateEpitope.per_allele_scores. Loaders must populate it or
    # downstream ranking sees zero-scored epitopes and drops everything.
    from .epitope_dsl import attach_per_allele_scores
    epitopes = attach_per_allele_scores(epitopes)
    logger.info(
        "Loaded %d epitope(s) (%d row(s) × %d predictor(s)) from %s",
        len(epitopes), len(report_df), len(chosen), path)
    return report_df, epitopes


def _build_lens_report_row(row, allele, peptide, source_sequence_name,
                           peptide_offset, detected, display_pred, chosen_tools):
    """Assemble one row of the LENS neoepitope report DataFrame.

    ``display_pred`` (if not None) drives the legacy Affinity / %ile
    columns; every *detected* predictor's raw value is also exposed as
    '<Tool> value (nM/hours)' so downstream DSL / scripts can see both.
    """
    antigen_source = _safe_str(row.get("antigen_source"))
    gene = _safe_str(row.get("gene_name"))
    variant_pos = _safe_str(row.get("variant_coords"))
    if not variant_pos:
        variant_pos = _safe_str(row.get("mut_aa_pos"))
    if not variant_pos:
        variant_pos = _safe_str(row.get("origin_descriptor"))

    display_value = (
        _safe_float(row.get(display_pred.value_col)) if display_pred else None)
    display_percentile = (
        _safe_float(row.get(display_pred.percentile_col))
        if display_pred and display_pred.percentile_col else None)
    display_agretopicity = (
        _safe_float(row.get(display_pred.agretopicity_col))
        if display_pred and display_pred.agretopicity_col else None)

    # Shared core columns (WT affinity isn't predicted on the import
    # path → None → blank); then LENS-specific extras.
    out = neoepitope_core_row(
        allele=allele, mutant_peptide=str(peptide),
        mutant_affinity=display_value, wt_peptide='', wt_affinity=None,
        gene_name=gene, variant=variant_pos)
    out.update({
        'Antigen source': antigen_source,
        'Predictors used': ','.join(chosen_tools),
        '%ile rank': display_percentile,
        'Agretopicity': display_agretopicity,
        'TPM': _safe_float(row.get("tpm")),
        'Source sequence name': source_sequence_name,
        'Peptide offset': peptide_offset,
    })
    # Every detected predictor gets a raw-value column so users can see
    # per-tool signals side-by-side in the report, even if not chosen
    # for DSL scoring.
    unit_by_kind = {"pMHC_affinity": "nM", "pMHC_stability": "hours"}
    for d in detected:
        unit = unit_by_kind.get(d.kind, "")
        label = f"{d.tool} value ({unit})" if unit else f"{d.tool} value"
        out[label] = _safe_float(row.get(d.value_col))
        if d.kind == "pMHC_affinity" and d.percentile_col:
            out[f"{d.tool} affinity percentile rank"] = _safe_float(
                row.get(d.percentile_col))
        if d.presentation_score_col:
            out[f"{d.tool} presentation score"] = _safe_float(
                row.get(d.presentation_score_col))
        if d.presentation_percentile_col:
            out[f"{d.tool} presentation percentile rank"] = _safe_float(
                row.get(d.presentation_percentile_col))
        if d.processing_score_col:
            out[f"{d.tool} processing score"] = _safe_float(
                row.get(d.processing_score_col))
    return out


# ── Shared report writer ─────────────────────────────────────────────────────

def _ensure_parent_dir(path):
    """Create the parent directory of an output file if it's missing.

    pandas' ``to_csv`` / ``to_excel`` refuse to write into a
    non-existent directory; mirror the vaccine writers, which already
    ``os.makedirs`` their target dir, so ``--output-csv foo/bar.csv``
    works without the user pre-creating ``foo/``.
    """
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def write_neoepitope_report(report_df, epitopes, excel_report_path=None,
                            csv_report_path=None, epitope_config=None,
                            topiary_df=None):
    """
    Score epitopes and write a neoepitope report from imported data.

    Scoring runs through the Topiary DSL (``filter_expr`` / ``score_expr``
    on the :class:`EpitopeConfig`), so external inputs (LENS, pVACseq)
    pick up the same formulas as the main VCF/BAM pipeline. When both
    expressions are unset the default node exactly reproduces the legacy
    per-prediction logistic score byte-for-byte.

    When multiple Predictions share a (peptide, allele) — as emitted by
    LENS files exposing both MHCflurry and netMHCpan — the DSL receives
    all of them and collapses to one score per group; that single score
    lands on the corresponding row of ``report_df``.

    Parameters
    ----------
    report_df : pandas.DataFrame
        Report-ready DataFrame from load_pvacseq or load_lens
        (one row per (peptide, allele) pair).
    epitopes : list of CandidateEpitope
        Output of ``load_lens`` / ``load_pvacseq`` — one CandidateEpitope per
        unique ``(peptide, source, offset)``.
    excel_report_path : str, optional
    csv_report_path : str, optional
    epitope_config : EpitopeConfig, optional
    topiary_df : pandas.DataFrame, optional
        Pre-built topiary long-form rows to use for DSL validation and
        scoring. ``load_pvacseq`` stores this on ``report_df.attrs`` so
        pVACseq annotation columns remain available to custom DSL
        expressions.
    """
    from .epitope_config import EpitopeConfig
    from .epitope_dsl import (
        epitopes_to_topiary_df,
        score_predictions,
        validate_dsl_against_predictions,
    )

    if epitope_config is None:
        epitope_config = EpitopeConfig()

    peptide_col = 'Mutant peptide sequence'
    allele_col = 'Allele'

    # The score-merge below keys by (peptide, allele) and broadcasts the
    # score back across every matching report_df row. Real-world LENS
    # files emit the same (peptide, allele) pair from multiple sources
    # (alternative transcripts, homologous regions, multiple variants)
    # — that's expected, not an error. Each row retains its own source
    # provenance (variant_coords, gene_name, transcript) and gets the
    # same score, which is the correct behavior for a per-row report.
    # The duplicate count is only useful when chasing a downstream
    # row-count mismatch — log at DEBUG, not INFO.
    if len(report_df) > 0 and logger.isEnabledFor(logging.DEBUG):
        dup_count = int(report_df.duplicated(
            subset=[peptide_col, allele_col]).sum())
        if dup_count:
            logger.debug(
                "report_df has %d duplicate (peptide, allele) row(s) "
                "from multi-source input; the same score broadcasts "
                "across each duplicate. Per-source provenance is "
                "preserved by the other columns.", dup_count)

    # Build the topiary DataFrame once and share it between validator and
    # scorer. pVACseq import keeps the loader-produced frame in
    # ``report_df.attrs`` so passthrough columns such as
    # ``pvacseq_mhcflurry_ic50_mt`` remain DSL-addressable; LENS and
    # vaxrank-native callers fall back to rebuilding from CandidateEpitope.
    if topiary_df is None:
        topiary_df = report_df.attrs.get("topiary_df")
    if topiary_df is None:
        topiary_df = epitopes_to_topiary_df(epitopes)

    validate_dsl_against_predictions(
        epitope_config, epitopes, topiary_df=topiary_df)

    score_series = score_predictions(
        epitopes, epitope_config, topiary_df=topiary_df)

    # score_series is indexed by topiary's group key
    # (source/variant, peptide, peptide_offset, allele). pVACseq reports
    # carry that source key through for precise alignment. If an exact
    # source-keyed row is absent, it was filtered out and must stay at 0.
    # Older LENS-shaped reports lack source/offset columns and keep the
    # legacy (peptide, allele) broadcast.
    scores_by_key = {}
    scores_by_pair = {}
    for idx_tuple, score in score_series.items():
        source_name, peptide, offset, allele = idx_tuple
        scores_by_key[(source_name, peptide, int(offset), allele)] = score
        scores_by_pair[(peptide, allele)] = score

    report_df = report_df.copy()
    source_col = 'Source sequence name'
    offset_col = 'Peptide offset'
    has_source_keys = (
        source_col in report_df.columns and offset_col in report_df.columns)
    scores = []
    for _, row in report_df.iterrows():
        score = None
        if has_source_keys:
            offset = _safe_float(row.get(offset_col))
            if offset is not None:
                score = scores_by_key.get((
                    row.get(source_col), row[peptide_col], int(offset),
                    row[allele_col]))
            if score is None:
                score = 0.0
        else:
            score = scores_by_pair.get((row[peptide_col], row[allele_col]), 0.0)
        scores.append(round(float(score), 6))
    report_df.insert(2, 'vaxrank_score', scores)

    report_df = report_df.sort_values('vaxrank_score', ascending=False)
    report_df.insert(0, 'rank', range(1, len(report_df) + 1))

    if csv_report_path:
        _ensure_parent_dir(csv_report_path)
        report_df.to_csv(csv_report_path, index=False)
        logger.info('Wrote CSV neoepitope report to %s', csv_report_path)

    if excel_report_path:
        _ensure_parent_dir(excel_report_path)
        writer = pd.ExcelWriter(excel_report_path, engine='openpyxl')
        report_df.to_excel(writer, sheet_name='Neoepitopes', index=False)
        worksheet = writer.sheets['Neoepitopes']
        worksheet.column_dimensions['C'].width = 23
        worksheet.column_dimensions['E'].width = 27
        worksheet.column_dimensions['F'].width = 17
        worksheet.column_dimensions['G'].width = 30
        worksheet.column_dimensions['H'].width = 18
        writer.close()
        logger.info('Wrote XLSX neoepitope report to %s', excel_report_path)
