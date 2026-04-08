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


def load_lens(path):
    """
    Import a LENS report TSV and return a neoepitope report DataFrame
    ready for output.

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    pandas.DataFrame
        Report-ready DataFrame
    list of EpitopePrediction
        Corresponding EpitopePrediction objects
    """
    df = pd.read_csv(path, sep="\t")
    required = {"peptide", "mhc_allele", "binding_affinity"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"LENS file {path} missing required columns: {missing}")

    predictions = []
    report_rows = []
    for _, row in df.iterrows():
        peptide = row.get("peptide", "")
        if not peptide or pd.isna(peptide):
            continue

        allele = row.get("mhc_allele", "")
        if pd.isna(allele):
            continue

        binding_affinity = row.get("binding_affinity")
        if pd.isna(binding_affinity):
            continue
        ic50 = float(binding_affinity)

        wt_ic50 = None
        ref_affinity = row.get("reference_binding_affinity")
        if ref_affinity is not None and not pd.isna(ref_affinity):
            wt_ic50 = float(ref_affinity)

        wt_peptide = row.get("reference_peptide", "")
        if pd.isna(wt_peptide):
            wt_peptide = ""

        percentile_ba = row.get("percent_rank_ba")
        percentile_el = row.get("percent_rank_el")
        if percentile_el is not None and not pd.isna(percentile_el):
            percentile_rank = float(percentile_el)
        elif percentile_ba is not None and not pd.isna(percentile_ba):
            percentile_rank = float(percentile_ba)
        else:
            percentile_rank = None

        has_wt = bool(wt_peptide) and wt_peptide != peptide

        gene = row.get("gene_name", "")
        if pd.isna(gene):
            gene = ""
        variant_pos = row.get("variant_position", "")
        if pd.isna(variant_pos):
            variant_pos = ""
        antigen_source = row.get("antigen_source", "")
        if pd.isna(antigen_source):
            antigen_source = ""

        pred = EpitopePrediction(
            allele=str(allele),
            peptide_sequence=str(peptide),
            wt_peptide_sequence=str(wt_peptide),
            ic50=ic50,
            wt_ic50=wt_ic50,
            percentile_rank=percentile_rank,
            prediction_method_name="lens",
            overlaps_mutation=has_wt,
            source_sequence="",
            offset=0,
            occurs_in_reference=False,
        )
        predictions.append(pred)

        wt_ic50_str = '%.2f nM' % wt_ic50 if wt_ic50 is not None else 'No prediction'
        report_rows.append({
            'Allele': str(allele),
            'Mutant peptide sequence': str(peptide),
            'Predicted mutant pMHC affinity': '%.2f nM' % ic50,
            'Wildtype sequence': str(wt_peptide),
            'Predicted wildtype pMHC affinity': wt_ic50_str,
            'Gene name': str(gene),
            'Genomic variant': str(variant_pos),
            'Antigen source': str(antigen_source),
            '%ile EL': _safe_float(row.get("percent_rank_el")),
            '%ile BA': _safe_float(row.get("percent_rank_ba")),
            'Agretopicity': _safe_float(row.get("agretopicity")),
            'TPM': _safe_float(row.get("tpm")),
        })

    report_df = pd.DataFrame(report_rows) if report_rows else pd.DataFrame()
    logger.info("Loaded %d epitope predictions from LENS file %s", len(predictions), path)
    return report_df, predictions


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

    Parameters
    ----------
    report_df : pandas.DataFrame
        Report-ready DataFrame from load_pvacseq or load_lens

    predictions : list of EpitopePrediction

    excel_report_path : str, optional
        Path for XLSX output

    csv_report_path : str, optional
        Path for CSV output

    epitope_config : EpitopeConfig, optional
    """
    from .epitope_config import EpitopeConfig
    if epitope_config is None:
        epitope_config = EpitopeConfig()

    # Add vaxrank score column
    scores = []
    for p in predictions:
        scores.append(p.logistic_epitope_score(
            midpoint=epitope_config.logistic_epitope_score_midpoint,
            width=epitope_config.logistic_epitope_score_width,
            ic50_cutoff=epitope_config.binding_affinity_cutoff,
            scoring_mode=epitope_config.scoring_mode,
        ))
    report_df = report_df.copy()
    report_df.insert(2, 'vaxrank_score', [round(s, 6) for s in scores])

    # Sort by score descending
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
