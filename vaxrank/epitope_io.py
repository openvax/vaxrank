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
  - vaxrank native (CSV/TSV roundtrip of EpitopePrediction fields)
  - pVACseq aggregated TSV (all_epitopes.aggregated.tsv)
  - LENS report TSV
"""

import logging
import re

import pandas as pd

from .epitope_prediction import EpitopePrediction

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
    "overlaps_mutation",
    "source_sequence",
    "offset",
    "occurs_in_reference",
]


def predictions_to_dataframe(predictions):
    """
    Convert a collection of EpitopePrediction objects to a DataFrame.

    Parameters
    ----------
    predictions : dict or list
        Either an OrderedDict keyed by (peptide, allele) -> EpitopePrediction,
        or a plain list of EpitopePrediction objects.

    Returns
    -------
    pandas.DataFrame
    """
    if isinstance(predictions, dict):
        predictions = list(predictions.values())
    rows = []
    for p in predictions:
        rows.append({
            "allele": p.allele,
            "peptide_sequence": p.peptide_sequence,
            "wt_peptide_sequence": p.wt_peptide_sequence,
            "ic50": p.ic50,
            "wt_ic50": p.wt_ic50,
            "percentile_rank": p.percentile_rank,
            "prediction_method_name": p.prediction_method_name,
            "overlaps_mutation": p.overlaps_mutation,
            "source_sequence": p.source_sequence,
            "offset": p.offset,
            "occurs_in_reference": p.occurs_in_reference,
        })
    return pd.DataFrame(rows, columns=VAXRANK_COLUMNS)


def save_predictions(predictions, path):
    """
    Save EpitopePrediction objects to a CSV/TSV file.

    Format is inferred from the file extension (.tsv -> tab, else comma).
    """
    df = predictions_to_dataframe(predictions)
    sep = "\t" if str(path).endswith(".tsv") else ","
    df.to_csv(path, sep=sep, index=False)
    logger.info("Saved %d predictions to %s", len(df), path)


def load_predictions(path):
    """
    Load EpitopePrediction objects from a vaxrank-native CSV/TSV file.

    Returns
    -------
    list of EpitopePrediction
    """
    sep = "\t" if str(path).endswith(".tsv") else ","
    df = pd.read_csv(path, sep=sep)
    return _dataframe_to_predictions(df)


def _dataframe_to_predictions(df):
    """Convert a vaxrank-format DataFrame to a list of EpitopePrediction."""
    predictions = []
    for _, row in df.iterrows():
        wt_ic50 = row.get("wt_ic50")
        if pd.isna(wt_ic50):
            wt_ic50 = None
        percentile_rank = row.get("percentile_rank")
        if pd.isna(percentile_rank):
            percentile_rank = None
        predictions.append(EpitopePrediction(
            allele=row["allele"],
            peptide_sequence=row["peptide_sequence"],
            wt_peptide_sequence=row.get("wt_peptide_sequence", ""),
            ic50=float(row["ic50"]),
            wt_ic50=float(wt_ic50) if wt_ic50 is not None else None,
            percentile_rank=float(percentile_rank) if percentile_rank is not None else None,
            prediction_method_name=row.get("prediction_method_name", ""),
            overlaps_mutation=bool(row.get("overlaps_mutation", True)),
            source_sequence=row.get("source_sequence", ""),
            offset=int(row.get("offset", 0)),
            occurs_in_reference=bool(row.get("occurs_in_reference", False)),
        ))
    return predictions


# ── pVACseq import ───────────────────────────────────────────────────────────

def load_pvacseq(path):
    """
    Import a pVACseq aggregated TSV (``*all_epitopes.aggregated.tsv``)
    and return a neoepitope report DataFrame ready for output.

    Each row is one variant's best epitope.

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
    list of EpitopePrediction
        Corresponding EpitopePrediction objects
    """
    df = pd.read_csv(path, sep="\t")
    required = {"Best Peptide", "Allele", "IC50 MT"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"pVACseq file {path} missing required columns: {missing}")

    predictions = []
    report_rows = []
    for _, row in df.iterrows():
        peptide = row.get("Best Peptide", "")
        if not peptide or pd.isna(peptide):
            continue

        allele = row.get("Allele", "")
        if pd.isna(allele):
            continue

        ic50_mt = row.get("IC50 MT")
        if pd.isna(ic50_mt):
            continue
        ic50_mt = float(ic50_mt)

        ic50_wt = row.get("IC50 WT")
        wt_ic50 = float(ic50_wt) if not pd.isna(ic50_wt) else None

        percentile = row.get("%ile MT")
        percentile_rank = float(percentile) if not pd.isna(percentile) else None

        wt_peptide = row.get("WT Epitope Seq", "")
        if pd.isna(wt_peptide):
            wt_peptide = ""

        ref_match = row.get("Ref Match", False)
        if isinstance(ref_match, str):
            occurs_in_reference = ref_match.strip().lower() == "true"
        else:
            occurs_in_reference = bool(ref_match) if not pd.isna(ref_match) else False

        method = row.get("Best MT IC50 Score Method", "pvacseq")
        if pd.isna(method):
            method = "pvacseq"

        gene = row.get("Gene", "")
        if pd.isna(gene):
            gene = ""
        variant_id = row.get("ID", "")
        if pd.isna(variant_id):
            variant_id = ""
        tier = row.get("Tier", "")
        if pd.isna(tier):
            tier = ""

        pred = EpitopePrediction(
            allele=str(allele),
            peptide_sequence=str(peptide),
            wt_peptide_sequence=str(wt_peptide),
            ic50=ic50_mt,
            wt_ic50=wt_ic50,
            percentile_rank=percentile_rank,
            prediction_method_name=str(method),
            overlaps_mutation=True,
            source_sequence="",
            offset=0,
            occurs_in_reference=occurs_in_reference,
        )
        predictions.append(pred)

        # Build report row preserving pVACseq metadata
        wt_ic50_str = '%.2f nM' % wt_ic50 if wt_ic50 is not None else 'No prediction'
        report_rows.append({
            'Allele': str(allele),
            'Mutant peptide sequence': str(peptide),
            'Predicted mutant pMHC affinity': '%.2f nM' % ic50_mt,
            'Wildtype sequence': str(wt_peptide),
            'Predicted wildtype pMHC affinity': wt_ic50_str,
            'Gene name': str(gene),
            'Genomic variant': str(variant_id),
            'Tier': str(tier),
            'Ref Match': occurs_in_reference,
            'RNA Expr': _safe_float(row.get("RNA Expr")),
            'RNA VAF': _safe_float(row.get("RNA VAF")),
            'DNA VAF': _safe_float(row.get("DNA VAF")),
            '%ile MT': percentile_rank,
        })

    report_df = pd.DataFrame(report_rows) if report_rows else pd.DataFrame()
    logger.info("Loaded %d epitope predictions from pVACseq file %s", len(predictions), path)
    return report_df, predictions


# LENS columns follow the convention "<tool>_<version>.<metric>", e.g.
# "mhcflurry_2.1.1.aff", "netmhcpan_4.1b.perc_rank_el",
# "netmhcstabpan_1.0.halflife_hours". We sniff which tools are present by
# regex, bucket columns by (tool, version), and consult the registry below
# to know which metric is the "value" and which are percentile ranks.
_LENS_COLUMN_RE = re.compile(r"^([A-Za-z][A-Za-z0-9]*)_(\d[\w.]*)\.(.+)$")

# Registry describing how each LENS-emitted predictor maps into the Topiary
# DSL. ``kind`` is the topiary prediction Kind (pMHC_affinity etc.).
# ``value_cols`` and ``percentile_cols`` are LENS metric suffixes in
# preference order — the first present metric wins.
_LENS_PREDICTOR_REGISTRY = {
    "mhcflurry": {
        "kind": "pMHC_affinity",
        "value_cols": ("aff",),
        "percentile_cols": ("pres_perc", "aff_perc"),
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

# Free-standing agretopicity columns (not versioned): "<tool>_agretopicity"
_LENS_AGRETOPICITY_COL = "{tool}_agretopicity"


class DetectedPredictor:
    """A predictor detected in a LENS TSV.

    Attributes
    ----------
    tool : str            - registry key, e.g. "mhcflurry"
    version : str         - version string parsed from the column name
    kind : str            - topiary prediction Kind (e.g. "pMHC_affinity")
    value_col : str       - LENS column for the predicted value (IC50 / hours)
    percentile_col : str or None
    agretopicity_col : str or None
    """
    __slots__ = (
        "tool", "version", "kind", "value_col",
        "percentile_col", "agretopicity_col",
    )

    def __init__(self, tool, version, kind, value_col,
                 percentile_col=None, agretopicity_col=None):
        self.tool = tool
        self.version = version
        self.kind = kind
        self.value_col = value_col
        self.percentile_col = percentile_col
        self.agretopicity_col = agretopicity_col

    def __repr__(self):
        return (f"DetectedPredictor(tool={self.tool!r}, version={self.version!r}, "
                f"kind={self.kind!r})")


def _detect_lens_predictors(columns):
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
        agretopicity_col = _LENS_AGRETOPICITY_COL.format(tool=tool)
        if agretopicity_col not in columns:
            agretopicity_col = None
        detected.append(DetectedPredictor(
            tool=tool,
            version=version,
            kind=spec["kind"],
            value_col=value_col,
            percentile_col=percentile_col,
            agretopicity_col=agretopicity_col,
        ))
    # Deterministic order: tool name alphabetical
    detected.sort(key=lambda d: (d.tool, d.version))
    return detected


def _pick_canonical_predictor(detected):
    """For --lens-predictor=auto, pick a single pMHC_affinity predictor.

    Preference order: mhcflurry (present across all LENS versions), then
    netmhcpan, then any remaining affinity predictor alphabetically.
    """
    affinity_preds = [d for d in detected if d.kind == "pMHC_affinity"]
    if not affinity_preds:
        return None
    for preferred in ("mhcflurry", "netmhcpan"):
        for d in affinity_preds:
            if d.tool == preferred:
                return d
    return affinity_preds[0]


def _normalize_hla_allele(allele):
    """LENS emits alleles as 'HLA-A01:01'; vaxrank output uses 'HLA-A*01:01'."""
    if not allele:
        return allele
    return re.sub(r"^(HLA-[A-Z]{1,3})(\d)", r"\1*\2", allele)


def load_lens(path, predictor="auto"):
    """
    Import a LENS report TSV and return a neoepitope report DataFrame
    plus EpitopePrediction objects.

    LENS column schemas vary across versions (v1.4, v1.5, v1.9+); we sniff
    which predictors a given file actually emits rather than requiring a
    fixed schema.

    Parameters
    ----------
    path : str or Path
    predictor : str
        ``"auto"`` picks a single canonical predictor (mhcflurry > netmhcpan).
        ``"all"`` emits one :class:`EpitopePrediction` per (peptide, allele,
        detected predictor) so that Topiary DSL expressions can combine them
        via ``affinity['mhcflurry']`` / ``affinity['netmhcpan']`` etc.
        A specific tool name (e.g. ``"mhcflurry"``) uses only that predictor
        and errors if absent.

    Returns
    -------
    pandas.DataFrame
        Report-ready DataFrame (one row per peptide × allele).
    list of EpitopePrediction
        Per-predictor predictions. In ``all`` mode there are N predictions
        per report row (N = number of detected predictors).
    """
    # low_memory=False avoids mixed-dtype warnings — LENS reports have many
    # optional columns that are mostly NA for a given antigen type.
    df = pd.read_csv(path, sep="\t", low_memory=False)

    required = {"peptide", "allele"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"LENS file {path} missing required columns: {missing}")

    detected = _detect_lens_predictors(df.columns)
    if not detected:
        raise ValueError(
            f"LENS file {path} has no recognized predictor columns "
            f"(expected columns like 'mhcflurry_<version>.aff' or "
            f"'netmhcpan_<version>.aff_nm')")

    available_tools = sorted({d.tool for d in detected})

    if predictor == "auto":
        canonical = _pick_canonical_predictor(detected)
        if canonical is None:
            raise ValueError(
                f"LENS file {path} has no pMHC affinity predictor "
                f"(detected: {[d.tool for d in detected]})")
        chosen = [canonical]
    elif predictor == "all":
        chosen = list(detected)
    else:
        chosen = [d for d in detected if d.tool == predictor]
        if not chosen:
            raise ValueError(
                f"LENS file {path} does not contain predictor {predictor!r} "
                f"(available: {available_tools})")

    logger.info(
        "LENS file %s: detected=%s, emitting for %s",
        path, [d.tool for d in detected], [d.tool for d in chosen])

    # Use the first chosen pMHC_affinity predictor as the "display" one
    # whose value lands in the report's Affinity / %ile columns. Stability
    # predictors are never chosen for display; they ride along in the
    # per-predictor extra columns.
    display_pred = next(
        (d for d in chosen if d.kind == "pMHC_affinity"), None)

    predictions = []
    report_rows = []
    for _, row in df.iterrows():
        peptide = row.get("peptide", "")
        if not peptide or pd.isna(peptide):
            continue

        allele_raw = row.get("allele", "")
        if pd.isna(allele_raw) or not allele_raw:
            continue
        allele = _normalize_hla_allele(str(allele_raw))

        # Build one EpitopePrediction per chosen predictor for this row.
        # Predictions with no usable value for this row are skipped.
        row_preds = []
        for d in chosen:
            value = _safe_float(row.get(d.value_col))
            if value is None:
                continue
            percentile_rank = (
                _safe_float(row.get(d.percentile_col))
                if d.percentile_col else None)
            pep_context = _safe_str(row.get("pep_context"))
            row_preds.append(EpitopePrediction(
                allele=allele,
                peptide_sequence=str(peptide),
                wt_peptide_sequence="",
                ic50=value,
                wt_ic50=None,
                percentile_rank=percentile_rank,
                prediction_method_name=d.tool,
                overlaps_mutation=True,
                source_sequence=pep_context,
                offset=0,
                occurs_in_reference=False,
            ))

        if not row_preds:
            continue
        predictions.extend(row_preds)

        report_rows.append(_build_lens_report_row(
            row=row,
            allele=allele,
            peptide=peptide,
            detected=detected,
            display_pred=display_pred,
            chosen_tools=[d.tool for d in chosen],
        ))

    report_df = pd.DataFrame(report_rows) if report_rows else pd.DataFrame()
    logger.info(
        "Loaded %d predictions (%d rows × %d predictors) from %s",
        len(predictions), len(report_df), len(chosen), path)
    return report_df, predictions


def _build_lens_report_row(row, allele, peptide, detected, display_pred,
                           chosen_tools):
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

    out = {
        'Allele': allele,
        'Mutant peptide sequence': str(peptide),
        'Predicted mutant pMHC affinity': (
            '%.2f nM' % display_value if display_value is not None
            else 'No prediction'),
        'Wildtype sequence': '',
        'Predicted wildtype pMHC affinity': 'No prediction',
        'Gene name': gene,
        'Genomic variant': variant_pos,
        'Antigen source': antigen_source,
        'Predictors used': ','.join(chosen_tools),
        '%ile rank': display_percentile,
        'Agretopicity': display_agretopicity,
        'TPM': _safe_float(row.get("tpm")),
    }
    # Every detected predictor gets a raw-value column so users can see
    # per-tool signals side-by-side in the report, even if not chosen
    # for DSL scoring.
    unit_by_kind = {"pMHC_affinity": "nM", "pMHC_stability": "hours"}
    for d in detected:
        unit = unit_by_kind.get(d.kind, "")
        label = f"{d.tool} value ({unit})" if unit else f"{d.tool} value"
        out[label] = _safe_float(row.get(d.value_col))
    return out


def _safe_str(val):
    """Coerce a cell to str, mapping NaN/None/'NA' to empty string."""
    if val is None:
        return ""
    if isinstance(val, float) and pd.isna(val):
        return ""
    s = str(val)
    if s == "NA":
        return ""
    return s


def _safe_float(val):
    """Convert to float, returning None for NaN/missing."""
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def write_neoepitope_report(report_df, predictions, excel_report_path=None,
                            csv_report_path=None, epitope_config=None):
    """
    Score predictions and write a neoepitope report from imported data.

    Scoring runs through the Topiary DSL (``filter_expr`` / ``score_expr``
    on the :class:`EpitopeConfig`), so external inputs (LENS, pVACseq)
    pick up the same formulas as the main VCF/BAM pipeline. When both
    expressions are unset the default node exactly reproduces the legacy
    :meth:`EpitopePrediction.logistic_epitope_score` output.

    When multiple EpitopePrediction objects share a (peptide, allele) —
    as emitted by ``load_lens(..., predictor='all')`` — the DSL receives
    all of them and collapses to one score per group; that single score
    lands on the corresponding row of ``report_df``.

    Parameters
    ----------
    report_df : pandas.DataFrame
        Report-ready DataFrame from load_pvacseq or load_lens
        (one row per (peptide, allele) pair).
    predictions : list of EpitopePrediction
    excel_report_path : str, optional
    csv_report_path : str, optional
    epitope_config : EpitopeConfig, optional
    """
    from .epitope_config import EpitopeConfig
    from .epitope_dsl import (
        score_predictions,
        validate_dsl_against_predictions,
    )

    if epitope_config is None:
        epitope_config = EpitopeConfig()

    # Fail loud if filter_expr / score_expr references a predictor / column
    # the loaded data doesn't expose.
    validate_dsl_against_predictions(epitope_config, predictions)

    score_series = score_predictions(predictions, epitope_config)

    # score_series is indexed by (source_sequence_name, peptide,
    # peptide_offset, allele). We built the topiary DF with
    # source_sequence_name = peptide, peptide_offset = 0, so the effective
    # group key is (peptide, allele) — align back onto report_df using that.
    scores_by_key = {}
    for idx_tuple, score in score_series.items():
        # idx_tuple = (source_seq_name, peptide, peptide_offset, allele)
        _, peptide, _, allele = idx_tuple
        scores_by_key[(peptide, allele)] = score

    report_df = report_df.copy()
    peptide_col = 'Mutant peptide sequence'
    allele_col = 'Allele'
    scores = [
        round(float(scores_by_key.get((row[peptide_col], row[allele_col]), 0.0)), 6)
        for _, row in report_df.iterrows()
    ]
    report_df.insert(2, 'vaxrank_score', scores)

    report_df = report_df.sort_values('vaxrank_score', ascending=False)
    report_df.insert(0, 'rank', range(1, len(report_df) + 1))

    if csv_report_path:
        report_df.to_csv(csv_report_path, index=False)
        logger.info('Wrote CSV neoepitope report to %s', csv_report_path)

    if excel_report_path:
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
