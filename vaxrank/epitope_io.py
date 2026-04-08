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

Supports three formats:
  - vaxrank native (CSV/TSV roundtrip of EpitopePrediction fields)
  - pVACseq aggregated TSV (all_epitopes.aggregated.tsv)
  - LENS report TSV

The common goal: produce a list of EpitopePrediction objects that can
be fed into vaxrank's vaccine peptide ranking pipeline.
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

    Format is inferred from the file extension (.tsv → tab, else comma).
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

# Regex for HLA allele columns in pVACseq aggregated files (e.g. A*29:02)
_PVAC_ALLELE_COL_RE = re.compile(r"^[ABC]\*\d+:\d+$|^D[PQR]")


def load_pvacseq(path):
    """
    Import epitope predictions from a pVACseq aggregated TSV file
    (``*all_epitopes.aggregated.tsv``).

    Each row in the pVACseq file is one variant; we extract the best
    epitope per variant and convert it to an EpitopePrediction.

    Parameters
    ----------
    path : str or Path
        Path to the pVACseq aggregated TSV file.

    Returns
    -------
    list of EpitopePrediction
    """
    df = pd.read_csv(path, sep="\t")
    required = {"Best Peptide", "Allele", "IC50 MT"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"pVACseq file {path} missing required columns: {missing}")

    predictions = []
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

        # pVACseq doesn't always include the WT peptide sequence in the
        # aggregated file, but the per-epitope file does. Use empty string
        # as fallback.
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

        predictions.append(EpitopePrediction(
            allele=str(allele),
            peptide_sequence=str(peptide),
            wt_peptide_sequence=str(wt_peptide),
            ic50=ic50_mt,
            wt_ic50=wt_ic50,
            percentile_rank=percentile_rank,
            prediction_method_name=str(method),
            overlaps_mutation=True,  # pVACseq only reports mutant epitopes
            source_sequence="",  # not available in aggregated format
            offset=0,
            occurs_in_reference=occurs_in_reference,
        ))

    logger.info("Loaded %d epitope predictions from pVACseq file %s", len(predictions), path)
    return predictions


# ── LENS import ──────────────────────────────────────────────────────────────

def load_lens(path):
    """
    Import epitope predictions from a LENS report TSV file.

    Parameters
    ----------
    path : str or Path
        Path to the LENS report TSV file.

    Returns
    -------
    list of EpitopePrediction
    """
    df = pd.read_csv(path, sep="\t")
    required = {"peptide", "mhc_allele", "binding_affinity"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"LENS file {path} missing required columns: {missing}")

    predictions = []
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
        # Prefer EL percentile rank if available
        if percentile_el is not None and not pd.isna(percentile_el):
            percentile_rank = float(percentile_el)
        elif percentile_ba is not None and not pd.isna(percentile_ba):
            percentile_rank = float(percentile_ba)
        else:
            percentile_rank = None

        has_wt = bool(wt_peptide) and wt_peptide != peptide

        predictions.append(EpitopePrediction(
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
        ))

    logger.info("Loaded %d epitope predictions from LENS file %s", len(predictions), path)
    return predictions


# ── Format auto-detection ────────────────────────────────────────────────────

def detect_format(path):
    """
    Auto-detect epitope file format by inspecting column headers.

    Returns
    -------
    str
        One of "vaxrank", "pvacseq", "lens"

    Raises
    ------
    ValueError
        If format cannot be determined.
    """
    sep = "\t" if str(path).endswith(".tsv") else ","
    df = pd.read_csv(path, sep=sep, nrows=0)
    cols = set(df.columns)

    if "peptide_sequence" in cols and "allele" in cols:
        return "vaxrank"
    if "Best Peptide" in cols and "IC50 MT" in cols:
        return "pvacseq"
    if "peptide" in cols and "mhc_allele" in cols and "binding_affinity" in cols:
        return "lens"

    raise ValueError(
        f"Cannot detect epitope file format from columns: {sorted(cols)[:10]}... "
        "Expected vaxrank (peptide_sequence, allele), "
        "pVACseq (Best Peptide, IC50 MT), or "
        "LENS (peptide, mhc_allele, binding_affinity).")


def load_epitopes(path):
    """
    Load epitope predictions from any supported format (auto-detected).

    Parameters
    ----------
    path : str or Path
        Path to epitope predictions file.

    Returns
    -------
    list of EpitopePrediction
    """
    fmt = detect_format(path)
    logger.info("Detected epitope file format: %s", fmt)
    if fmt == "vaxrank":
        return load_predictions(path)
    elif fmt == "pvacseq":
        return load_pvacseq(path)
    elif fmt == "lens":
        return load_lens(path)
    else:
        raise ValueError(f"Unknown epitope file format: {fmt}")
