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
which groups native rows by position and external rows by their complete
``ExternalPredictionKey`` identity into ``CandidateEpitope`` objects.
"""

import json
import logging
import os
import re
from dataclasses import dataclass

import pandas as pd
from mhctools.pred import Prediction

from . import cells
from .epitope_dsl import prediction_kind_for_method
from .external_prediction import (
    _PVACSEQ_DNA_VAF_COLUMNS,
    _PVACSEQ_RNA_VAF_COLUMNS,
    ExternalPredictionKey,
)
from .external_report import (
    GENOMIC_VARIANT_COLUMN, ExternalRecord, ExternalReport,
)
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
    # Targetability and the antigen-aware exact-self result are *decisions*,
    # not raw inputs: a reload that dropped them silently re-admitted
    # held-out epitopes as vaccine targets, because both fields fall back to
    # the permissive legacy behavior when absent.
    "overlaps_targetable",
    "self_reference_match",
    "source_sequence",
    "source_name",
    # Flanking residues drive processing/cleavage predictors, so a reload
    # without them cannot reproduce the scores of the file it read.
    "n_flank",
    "c_flank",
    "prediction_id",
    "offset",
    "occurs_in_reference",
    "occurs_in_non_CTA_reference",
    "predictor_version",
    # Generic prediction fields. The legacy ``ic50`` column remains above
    # for compatibility, but cannot represent presentation, processing, or
    # any future non-affinity kind.
    "prediction_kind",
    "prediction_peptide",
    "prediction_value",
    "prediction_score",
    "prediction_tcr",
    "prediction_n_flank",
    "prediction_c_flank",
    "prediction_source_sequence_name",
    "prediction_offset",
    # Explicit input-allele membership is necessary for canonical
    # allele-independent processing predictions.
    "patient_alleles",
    # Preserve the exact configured DSL result so a native reload does not
    # silently revert to an unrelated default scoring formula.
    "per_allele_scores",
    # Generic WT fields parallel the mutant prediction fields. The legacy
    # ``wt_ic50`` column remains populated for affinity consumers.
    "wt_prediction_kind",
    "wt_prediction_peptide",
    "wt_prediction_value",
    "wt_prediction_score",
    "wt_percentile_rank",
    "wt_prediction_tcr",
    "wt_prediction_n_flank",
    "wt_prediction_c_flank",
    "wt_prediction_source_sequence_name",
    "wt_prediction_offset",
]


def epitope_prediction_rows(epitope):
    """Serialize one candidate to one native row per mutant prediction.

    Every prediction kind is emitted. Affinity values are duplicated into
    the legacy ``ic50`` column for compatibility; the generic prediction
    columns retain the kind, value, normalized score, percentile, and
    provenance needed to reconstruct the original ``Prediction`` records.
    """
    wt = epitope.wt
    wt_peptide_sequence = wt.sequence if wt is not None else ""
    wt_by_key = {}
    if wt is not None:
        for p in wt.predictions_flat():
            key = (
                p.kind, p.predictor_name, p.predictor_version, p.allele)
            wt_by_key.setdefault(key, []).append(p)
    rows = []
    for p in epitope.predictions_flat():
        key = (p.kind, p.predictor_name, p.predictor_version, p.allele)
        matching_wt = wt_by_key.get(key, [])
        wt_prediction = matching_wt.pop(0) if matching_wt else None
        rows.append({
            "allele": p.allele,
            "peptide_sequence": epitope.sequence,
            "wt_peptide_sequence": wt_peptide_sequence,
            "ic50": p.value if p.kind == "pMHC_affinity" else None,
            "wt_ic50": (
                wt_prediction.value
                if (
                    wt_prediction is not None
                    and wt_prediction.kind == "pMHC_affinity"
                )
                else None),
            "percentile_rank": p.percentile_rank,
            "prediction_method_name": p.predictor_name,
            "source_class": epitope.source_class,
            "overlaps_mutation": epitope.overlaps_mutation,
            "overlaps_targetable": epitope.overlaps_targetable,
            "self_reference_match": (
                epitope.self_reference_match.to_json()
                if epitope.self_reference_match is not None else ""),
            "source_sequence": epitope.source_sequence,
            "source_name": epitope.source_name,
            "n_flank": epitope.n_flank,
            "c_flank": epitope.c_flank,
            "prediction_id": epitope.prediction_id,
            "offset": epitope.offset,
            "occurs_in_reference": epitope.occurs_in_reference,
            "occurs_in_non_CTA_reference": epitope.occurs_in_non_CTA_reference,
            "predictor_version": p.predictor_version,
            "prediction_kind": p.kind,
            "prediction_peptide": p.peptide,
            "prediction_value": p.value,
            "prediction_score": p.score,
            "prediction_tcr": p.tcr,
            "prediction_n_flank": p.n_flank,
            "prediction_c_flank": p.c_flank,
            "prediction_source_sequence_name": p.source_sequence_name,
            "prediction_offset": p.offset,
            "patient_alleles": json.dumps(
                list(epitope.patient_alleles), separators=(",", ":")),
            "per_allele_scores": json.dumps(
                epitope.per_allele_scores,
                sort_keys=True,
                separators=(",", ":"),
            ),
            "wt_prediction_kind": (
                wt_prediction.kind if wt_prediction is not None else ""),
            "wt_prediction_peptide": (
                wt_prediction.peptide if wt_prediction is not None else ""),
            "wt_prediction_value": (
                wt_prediction.value if wt_prediction is not None else None),
            "wt_prediction_score": (
                wt_prediction.score if wt_prediction is not None else None),
            "wt_percentile_rank": (
                wt_prediction.percentile_rank
                if wt_prediction is not None else None),
            "wt_prediction_tcr": (
                wt_prediction.tcr if wt_prediction is not None else ""),
            "wt_prediction_n_flank": (
                wt_prediction.n_flank if wt_prediction is not None else ""),
            "wt_prediction_c_flank": (
                wt_prediction.c_flank if wt_prediction is not None else ""),
            "wt_prediction_source_sequence_name": (
                wt_prediction.source_sequence_name
                if wt_prediction is not None else None),
            "wt_prediction_offset": (
                wt_prediction.offset if wt_prediction is not None else None),
        })
    unpaired_wt = [
        prediction
        for predictions in wt_by_key.values()
        for prediction in predictions
    ]
    if unpaired_wt:
        raise ValueError(
            "Cannot serialize WT predictions without matching mutant "
            "prediction leaves; refusing to lose comparator evidence")
    return rows


def predictions_to_dataframe(epitopes):
    """Convert a collection of ``CandidateEpitope`` objects to a flat DataFrame
    matching :data:`VAXRANK_COLUMNS` — one row per mutant prediction leaf,
    regardless of prediction kind.

    Parameters
    ----------
    epitopes : iterable of CandidateEpitope
    """
    rows = []
    for e in epitopes:
        rows.extend(epitope_prediction_rows(e))
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
    """Load a vaxrank-native CSV/TSV file.

    Files written before generic prediction columns were introduced remain
    supported as affinity-only inputs. New-format rows must carry an explicit
    prediction kind and score; missing scientific values are not converted to
    implicit zeros.
    """
    sep = "\t" if str(path).endswith(".tsv") else ","
    df = pd.read_csv(path, sep=sep)

    def string_or_empty(val):
        if val is None or (isinstance(val, float) and pd.isna(val)):
            return ""
        return str(val)

    def optional_float(val):
        if val is None or pd.isna(val):
            return None
        return float(val)

    has_generic_predictions = "prediction_kind" in df.columns

    rows = []
    for _, row in df.iterrows():
        wt_ic50 = optional_float(row.get("wt_ic50"))
        percentile_rank = optional_float(row.get("percentile_rank"))
        peptide = row["peptide_sequence"]
        wt_peptide = string_or_empty(row.get("wt_peptide_sequence", ""))
        method = string_or_empty(row.get("prediction_method_name", ""))
        version = string_or_empty(row.get("predictor_version", ""))
        allele = string_or_empty(row.get("allele", ""))
        if has_generic_predictions:
            kind = string_or_empty(row.get("prediction_kind"))
            if not kind:
                raise ValueError(
                    "Native epitope row is missing prediction_kind")
            score = optional_float(row.get("prediction_score"))
            if score is None:
                raise ValueError(
                    "Native epitope row is missing prediction_score")
            value = optional_float(row.get("prediction_value"))
            prediction_peptide = string_or_empty(
                row.get("prediction_peptide")) or peptide
        else:
            # Legacy files only represented affinity IC50 and did not retain
            # normalized scores. Keep that historical 0.0 compatibility
            # behavior only when the generic schema is entirely absent.
            kind = prediction_kind_for_method(method)
            score = 0.0
            value = float(row["ic50"])
            prediction_peptide = peptide
        mutant = Prediction(
            kind=kind, predictor_name=method, predictor_version=version,
            allele=allele, peptide=prediction_peptide, value=value,
            score=score,
            percentile_rank=(
                float(percentile_rank) if percentile_rank is not None
                else None),
            tcr=string_or_empty(row.get("prediction_tcr", "")),
            n_flank=string_or_empty(row.get("prediction_n_flank", "")),
            c_flank=string_or_empty(row.get("prediction_c_flank", "")),
            source_sequence_name=(
                string_or_empty(row.get("prediction_source_sequence_name"))
                or None),
            offset=int(optional_float(row.get("prediction_offset")) or 0),
        )
        wt = None
        wt_kind = string_or_empty(row.get("wt_prediction_kind", ""))
        if wt_kind:
            wt_score = optional_float(row.get("wt_prediction_score"))
            if wt_score is None:
                raise ValueError(
                    "Native epitope row is missing wt_prediction_score")
            wt = Prediction(
                kind=wt_kind, predictor_name=method,
                predictor_version=version, allele=allele,
                peptide=(
                    string_or_empty(row.get("wt_prediction_peptide"))
                    or wt_peptide),
                value=optional_float(row.get("wt_prediction_value")),
                score=wt_score,
                percentile_rank=optional_float(row.get("wt_percentile_rank")),
                tcr=string_or_empty(row.get("wt_prediction_tcr", "")),
                n_flank=string_or_empty(
                    row.get("wt_prediction_n_flank", "")),
                c_flank=string_or_empty(
                    row.get("wt_prediction_c_flank", "")),
                source_sequence_name=(
                    string_or_empty(
                        row.get("wt_prediction_source_sequence_name"))
                    or None),
                offset=int(optional_float(
                    row.get("wt_prediction_offset")) or 0),
            )
        elif wt_ic50 is not None:
            wt = Prediction(
                kind=kind, predictor_name=method,
                predictor_version=version, allele=allele,
                peptide=wt_peptide, value=wt_ic50, score=0.0,
                percentile_rank=None,
            )
        patient_alleles_raw = string_or_empty(row.get("patient_alleles", ""))
        patient_alleles = (
            tuple(json.loads(patient_alleles_raw))
            if patient_alleles_raw else ())
        per_allele_scores_raw = string_or_empty(
            row.get("per_allele_scores", ""))
        per_allele_scores = (
            json.loads(per_allele_scores_raw)
            if per_allele_scores_raw else {})
        occurs_in_ref = bool(row.get("occurs_in_reference", False))
        # ``overlaps_targetable`` is tri-state: True / False / "producer did
        # not decide". A blank cell must reload as None, not as False.
        targetable_raw = row.get("overlaps_targetable")
        overlaps_targetable = (
            None if cells.missing(targetable_raw)
            else cells.boolean(targetable_raw))
        self_reference_raw = string_or_empty(
            row.get("self_reference_match", ""))
        self_reference_match = None
        if self_reference_raw:
            from .vaccine_antigen import SelfReferenceMatch
            self_reference_match = SelfReferenceMatch.from_json(
                self_reference_raw)
        rows.append({
            'peptide': peptide,
            'source': string_or_empty(row.get("source_sequence", "")),
            'source_name': string_or_empty(row.get("source_name", "")),
            'n_flank': string_or_empty(row.get("n_flank", "")),
            'c_flank': string_or_empty(row.get("c_flank", "")),
            'prediction_id': string_or_empty(row.get("prediction_id", "")),
            'offset': int(row.get("offset", 0)),
            'overlaps_targetable': overlaps_targetable,
            'self_reference_match': self_reference_match,
            'mutant': mutant,
            'wt': wt,
            'source_class': string_or_empty(
                row.get("source_class", "")) or None,
            'overlaps_mutation': bool(row.get("overlaps_mutation", True)),
            'occurs_in_reference': occurs_in_ref,
            # Fall back to raw occurs_in_reference when the CSV is
            # from a pre-3.0 producer that didn't emit this column.
            'occurs_in_non_CTA_reference': bool(
                row.get("occurs_in_non_CTA_reference", occurs_in_ref)),
            'patient_alleles': patient_alleles,
            'per_allele_scores': per_allele_scores,
        })
    return candidate_epitopes_from_rows(rows)


# ── pVACseq import ───────────────────────────────────────────────────────────

def _format_nm(value):
    # Missing affinity → blank, consistent with every other missing
    # value in the neoepitope CSVs (see _build_lens_report_row). A
    # blank cell in a "Predicted … affinity" column reads as "not
    # predicted" without mixing a string sentinel into numeric data.
    value = cells.number(value)
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


def _topiary_group_source(row):
    """Source key matching topiary.ranking.EvalContext grouping."""
    return (
        cells.text(row.get("variant"))
        or cells.text(row.get("source_sequence_name"))
        or cells.text(row.get("peptide"))
    )


def _topiary_pvacseq_to_epitope_rows(rows):
    """Convert topiary.read_pvacseq long-form rows to CandidateEpitope rows."""
    epitope_rows = []
    n_rows = 0
    for row in rows:
        peptide = cells.text(row.get("peptide"))
        allele = cells.text(row.get("allele"))
        value = cells.number(row.get("value"))
        if not peptide or not allele or value is None:
            continue

        method = cells.text(row.get("prediction_method_name")) or "pvacseq"
        version = cells.text(row.get("predictor_version"))
        kind = (
            cells.text(row.get("kind"))
            or prediction_kind_for_method(method))
        percentile_rank = cells.number(row.get("percentile_rank"))
        score = cells.number(row.get("score"))
        mutant = Prediction(
            kind=kind, predictor_name=method, predictor_version=version,
            allele=allele, peptide=peptide, value=value,
            score=score if score is not None else 0.0,
            percentile_rank=percentile_rank,
        )

        wt = None
        wt_value = cells.number(row.get("wt_value"))
        if wt_value is not None:
            wt_method = cells.text(row.get("wt_prediction_method_name")) or method
            wt_version = cells.text(row.get("wt_predictor_version")) or version
            wt = Prediction(
                kind=kind, predictor_name=wt_method,
                predictor_version=wt_version, allele=allele,
                peptide=cells.text(row.get("wt_peptide")),
                value=wt_value, score=0.0,
                percentile_rank=cells.number(row.get("wt_percentile_rank")),
            )

        epitope_rows.append({
            'peptide': peptide,
            # pVACseq ships no source-protein context, so the peptide *is*
            # its own source window — the same thing the prediction identity
            # records. The variant string is provenance and belongs in
            # ``source_name``; storing it in ``source_sequence`` made a
            # non-sequence look like a sequence, and any consumer that
            # sliced that window (construct assembly does) got a fragment of
            # the variant's *name* back.
            'source': peptide,
            'source_name': _topiary_group_source(row),
            'prediction_id': cells.text(row.get("prediction_id")),
            # 0, not the frame's peptide_offset, and for the same reason
            # ExternalPredictionKey.from_pvacseq_row pins it: the peptide is
            # its own source window, so it sits at offset 0 of that window.
            # Reading a source-protein offset here would let the epitope and
            # the identity that names it disagree about the same position.
            'offset': 0,
            'mutant': mutant,
            'wt': wt,
            'source_class': SOURCE_CLASS_MUTATION,
            'overlaps_mutation': cells.boolean(
                row.get("contains_mutant_residues"), default=True),
            'occurs_in_reference': cells.boolean(row.get("pvacseq_ref_match")),
            'occurs_in_non_CTA_reference': cells.boolean(
                row.get("pvacseq_ref_match")),
        })
        n_rows += 1
    return epitope_rows, n_rows


def _build_pvacseq_report_row(row):
    """Build one user-facing report row from a topiary pVACseq row."""
    # Same helper, same precedence, as the ranking path uses for these
    # fields — so the report and the ranking cannot disagree about a row.
    rna_expr = cells.first_number(
        row, "rna_transcript_expression", "transcript_expression",
        "gene_expression")
    rna_vaf = cells.first_number(row, *_PVACSEQ_RNA_VAF_COLUMNS)
    dna_vaf = cells.first_number(row, *_PVACSEQ_DNA_VAF_COLUMNS)
    out = neoepitope_core_row(
        allele=cells.text(row.get("allele")),
        mutant_peptide=cells.text(row.get("peptide")),
        mutant_affinity=row.get("value"),
        wt_peptide=cells.text(row.get("wt_peptide")),
        wt_affinity=row.get("wt_value"),
        gene_name=cells.text(row.get("gene")),
        variant=cells.text(row.get("variant")))
    out.update({
        'Tier': cells.text(row.get("pvacseq_tier")),
        'Ref Match': cells.boolean(row.get("pvacseq_ref_match")),
        'RNA Expr': rna_expr,
        'RNA VAF': rna_vaf,
        'DNA VAF': dna_vaf,
        '%ile MT': cells.number(row.get("percentile_rank")),
        'Prediction identity': cells.text(row.get("prediction_id")),
        'Source sequence name': _topiary_group_source(row),
        'Peptide offset': int(cells.number(row.get("peptide_offset")) or 0),
        'MHC class': cells.text(row.get("mhc_class")),
        'Contains mutant residues': cells.boolean(
            row.get("contains_mutant_residues"), default=True),
    })
    for col, value in row.items():
        if col.startswith("pvacseq_") and col not in {
                "pvacseq_ref_match", "pvacseq_tier"}:
            out[col] = cells.number(value)
    return out


# pVACseq's all_epitopes flavor carries both a row key (``Index``) and genomic
# coordinates. ``topiary.read_pvacseq`` prefers ``Index`` for its ``variant``
# column and does not carry the coordinate columns through — so the identity
# survives the normalization and the genome position does not. Construct
# ranking needs the position, so we recover it here, at load, keyed by the
# same identity topiary produced. This is the only place the raw TSV is
# consulted, and what it recovers is stated rather than assumed.
_PVACSEQ_COORDINATE_COLUMNS = ("Chromosome", "Start", "Reference", "Variant")


def pvacseq_genomic_coordinates(path):
    """Map each pVACseq row identity to its ``chrom-start-ref-alt`` string.

    Returns an empty map for flavors that ship no coordinate columns (the
    aggregated flavor encodes the position in its ``ID`` instead).
    """
    from .external_prediction import pvacseq_variant_id
    try:
        header = pd.read_csv(path, sep="\t", nrows=0)
    except (ValueError, OSError):
        return {}
    if not set(_PVACSEQ_COORDINATE_COLUMNS) <= set(header.columns):
        return {}
    id_columns = [
        column for column in ("Index", "ID")
        if column in header.columns]
    frame = pd.read_csv(
        path, sep="\t", low_memory=False,
        usecols=list(_PVACSEQ_COORDINATE_COLUMNS) + id_columns)
    coordinates = {}
    for _, raw in frame.iterrows():
        row = dict(raw)
        variant_id = pvacseq_variant_id(row)
        if not variant_id:
            continue
        parts = [cells.text(row.get(c)) for c in _PVACSEQ_COORDINATE_COLUMNS]
        if all(parts):
            coordinates.setdefault(variant_id, "-".join(parts))
    return coordinates


def read_pvacseq_report(path, epitope_config=None):
    """
    Import a pVACseq TSV and return a neoepitope report DataFrame ready
    for output.

    Both pVACseq output flavors are accepted:
    ``*all_epitopes.aggregated.tsv`` and ``*all_epitopes.tsv``.

    Parameters
    ----------
    path : str or Path
    epitope_config : EpitopeConfig, optional
        Effective DSL configuration to evaluate before returning candidates.

    Returns
    -------
    pandas.DataFrame
        Columns: Allele, Mutant peptide sequence, Predicted mutant pMHC
        affinity, Wildtype sequence, Predicted wildtype pMHC affinity,
        Gene name, Genomic variant, Tier, Ref Match, RNA Expr, DNA VAF,
        vaxrank_score
    list of CandidateEpitope
        One ``vaxrank.candidate_epitope.CandidateEpitope`` per unique
        external prediction identity.
    """
    from topiary import read_pvacseq

    result = read_pvacseq(path)
    topiary_df = result.df.copy()

    # Topiary has already normalized both pVACseq flavors onto one column
    # vocabulary (``variant``, ``gene``, ``transcript``, ``rna_depth`` /
    # ``tumor_rna_depth``, …). Every later stage reads these rows, never the
    # raw TSV — so an identity built here and an identity used for construct
    # ranking are the same object rather than two derivations that have to be
    # kept in agreement.
    rows = [dict(row) for _, row in topiary_df.iterrows()]
    coordinates = pvacseq_genomic_coordinates(path)
    records = []
    prediction_ids = []
    for row in rows:
        genomic = coordinates.get(cells.text(row.get("variant")))
        if genomic:
            row[GENOMIC_VARIANT_COLUMN] = genomic
        key = ExternalPredictionKey.from_pvacseq_row(row)
        prediction_ids.append(key.identifier if key is not None else "")
        if key is not None:
            records.append(ExternalRecord(key=key, row=row))
    topiary_df["prediction_id"] = prediction_ids
    for row, prediction_id in zip(rows, prediction_ids):
        row["prediction_id"] = prediction_id

    epitope_rows, n_rows = _topiary_pvacseq_to_epitope_rows(rows)
    report_rows = [
        _build_pvacseq_report_row(row)
        for row in rows
        if cells.text(row.get("peptide")) and cells.text(row.get("allele"))
    ]

    report_df = pd.DataFrame(report_rows) if report_rows else pd.DataFrame()
    report_df.attrs["topiary_df"] = topiary_df
    epitopes = candidate_epitopes_from_rows(epitope_rows)
    from .epitope_dsl import attach_per_allele_scores
    epitopes = attach_per_allele_scores(
        epitopes, epitope_config, topiary_df=topiary_df)
    logger.info(
        "Loaded %d epitope(s) (%d row(s), %s flavor) from pVACseq file %s",
        len(epitopes), n_rows, result.extra.get("pvacseq_format"), path)
    return ExternalReport(
        source_format="pvacseq",
        path=str(path),
        report_df=report_df,
        epitopes=tuple(epitopes),
        records=tuple(records),
        rows=tuple(rows),
    )


# ── LENS import ──────────────────────────────────────────────────────────────

# ── LENS import ──────────────────────────────────────────────────────────────
#
# LENS column names vary across versions — v1.4 writes ``mhcflurry_2.1.1.aff``,
# v1.9 writes something else, and each tool spells its metrics differently.
# vaxrank used to absorb that itself with a ``<tool>_<version>.<metric>`` regex
# and a per-tool registry of metric-name preferences.
#
# ``topiary.read_lens`` already does it, for every LENS version topiary knows,
# and normalizes onto a regular ``<tool>_<kind>_<metric>`` scheme. Reading that
# instead of re-deriving it means new LENS versions arrive with a topiary
# upgrade rather than a vaxrank patch — and it surfaces signals the local
# registry never listed, notably netmhcpan presentation and netmhcstabpan
# score.

# Two printings of one processing score are "the same value" below this.
# Wide enough to absorb a file that writes 0.8 on one row and 0.8000 on the
# next, tight enough that a genuinely different score still gets reported.
PROCESSING_SCORE_TOLERANCE = 1e-9

# (topiary metric suffix, mhctools Kind, which Prediction field it fills).
# ``value`` is the predicted quantity (IC50 nM, stability half-life);
# ``score`` is a normalized 0-1 signal; ``rank`` is a percentile.
_LENS_KIND_METRICS = (
    ("affinity", "pMHC_affinity"),
    ("presentation", "pMHC_presentation"),
    ("stability", "pMHC_stability"),
    ("antigen_processing", "antigen_processing"),
)


@dataclass(slots=True)
class DetectedPredictor:
    """One ``(tool, kind)`` signal found in a topiary-normalized LENS frame.

    Built from the columns ``topiary.read_lens`` produced rather than from
    raw LENS spellings. ``value_col`` / ``score_col`` / ``rank_col`` are
    ``None`` where that metric is absent for this tool and kind.
    """
    tool: str
    kind: str
    value_col: str | None = None
    score_col: str | None = None
    rank_col: str | None = None
    agretopicity_col: str | None = None
    # Retained so report columns and log lines keep naming a version.
    version: str = ""


# Bridge for openvax/topiary#206: read_lens normalizes a binding column only
# when it recognizes the version string, and passes an unrecognized one
# through verbatim — so a LENS file written by netMHCpan "4.1" rather than
# "4.1b" would lose its entire affinity axis with nothing in the logs. Map
# any residual predictor-shaped column by its metric, which is what actually
# determines meaning, and say so. Delete once topiary matches structurally.
_RESIDUAL_LENS_COLUMN_RE = re.compile(
    r"^([A-Za-z][A-Za-z0-9]*)_(\d[\w.]*)\.(.+)$")

# LENS metric suffix -> (topiary kind segment, metric field).
_RESIDUAL_METRIC_MAP = {
    "aff": ("affinity", "value"),
    "aff_nm": ("affinity", "value"),
    "aff_perc": ("affinity", "rank"),
    "perc_rank_ba": ("affinity", "rank"),
    "score_ba": ("affinity", "score"),
    "perc_rank_el": ("presentation", "rank"),
    "score_el": ("presentation", "score"),
    "pres_score": ("presentation", "score"),
    "pres_perc": ("presentation", "rank"),
    "proc_score": ("antigen_processing", "score"),
    "halflife_hours": ("stability", "value"),
    "perc_rank_stab": ("stability", "rank"),
}


def _map_residual_lens_columns(df):
    """Normalize predictor columns topiary left unmapped, and report them."""
    renames = {}
    unmapped = []
    for column in df.columns:
        match = _RESIDUAL_LENS_COLUMN_RE.match(str(column))
        if not match:
            continue
        tool, version, metric = match.groups()
        mapped = _RESIDUAL_METRIC_MAP.get(metric)
        if mapped is None:
            unmapped.append(str(column))
            continue
        kind, field = mapped
        target = "%s_%s_%s" % (tool, kind, field)
        if target not in df.columns:
            renames[column] = target
    if renames:
        logger.warning(
            "topiary did not normalize %d LENS binding column(s), most "
            "likely an unrecognized predictor version (openvax/topiary#206). "
            "Mapping them by metric so the signal is not lost: %s",
            len(renames),
            ", ".join("%s -> %s" % kv for kv in sorted(renames.items())))
        df = df.rename(columns=renames)
    if unmapped:
        logger.warning(
            "LENS column(s) %s look like predictor output but their metric "
            "is unknown to vaxrank and to topiary; they are ignored.",
            ", ".join(sorted(unmapped)))
    return df


def detect_lens_predictors(columns, versions=None):
    """Return a :class:`DetectedPredictor` per ``(tool, kind)`` present.

    ``columns`` are the columns of a ``topiary.read_lens`` frame; ``versions``
    is topiary's ``{method: version}`` map when available. Tools are whatever
    topiary emitted — vaxrank no longer keeps its own list, so a predictor
    topiary learns to read needs no change here.
    """
    columns = set(columns)
    versions = versions or {}
    detected = []
    tools = sorted({
        column.split("_", 1)[0]
        for column in columns
        for metric, _kind in _LENS_KIND_METRICS
        if column.startswith(("%s_" % metric,)) is False
        and ("_%s_" % metric) in column
    })
    for tool in tools:
        for metric, kind in _LENS_KIND_METRICS:
            names = {
                field: "%s_%s_%s" % (tool, metric, field)
                for field in ("value", "score", "rank")
            }
            present = {
                field: name for field, name in names.items()
                if name in columns}
            if not present:
                continue
            agretopicity = "%s_agretopicity" % tool
            detected.append(DetectedPredictor(
                tool=tool,
                kind=kind,
                value_col=present.get("value"),
                score_col=present.get("score"),
                rank_col=present.get("rank"),
                agretopicity_col=(
                    agretopicity if agretopicity in columns else None),
                version=str(versions.get(tool, "") or ""),
            ))
    detected.sort(key=lambda d: (d.tool, d.kind))
    return detected


def _pick_canonical_predictor(affinity_preds):
    """Pick the :class:`DetectedPredictor` whose value drives display columns.

    Uses topiary's ``CANONICAL_METHOD_PREFERENCE`` so the report's headline
    affinity names the same model the DSL resolves an unqualified reference
    to. vaxrank used to keep its own preference order here, and it had
    drifted from topiary's.
    """
    if not affinity_preds:
        return None
    from topiary.ranking import CANONICAL_METHOD_PREFERENCE

    def rank(detected):
        try:
            return (CANONICAL_METHOD_PREFERENCE.index(detected.tool),
                    detected.tool)
        except ValueError:
            return (len(CANONICAL_METHOD_PREFERENCE), detected.tool)

    return sorted(affinity_preds, key=rank)[0]


def normalize_hla_allele(allele):
    """LENS emits alleles as 'HLA-A01:01'; vaxrank output uses 'HLA-A*01:01'."""
    if not allele:
        return allele
    return re.sub(r"^(HLA-[A-Z]{1,3})(\d)", r"\1*\2", allele)


def lens_epitope_position(peptide, peptide_context):
    """Return the position identity shared by LENS import and ranking.

    The identity is ``(peptide, source context, offset)``. Missing context is
    represented by an empty source and offset zero; a peptide that does not
    occur in a non-empty context has no valid position and returns ``None``.
    """
    peptide = cells.text(peptide)
    peptide_context = cells.text(peptide_context)
    if not peptide:
        return None
    if peptide_context:
        offset = peptide_context.find(peptide)
        if offset < 0:
            return None
    else:
        offset = 0
    return peptide, peptide_context, offset


def read_lens_report(path, epitope_config=None):
    """
    Import a LENS report TSV and return a neoepitope report DataFrame
    plus a list of ``vaxrank.candidate_epitope.CandidateEpitope`` objects.

    LENS column schemas vary across versions (v1.4, v1.5, v1.9+); we sniff
    which predictors a given file actually emits rather than requiring a
    fixed schema. Internally we build one flat per-(peptide, allele,
    detected predictor, prediction kind) record per row, then group them
    into ``CandidateEpitope`` objects keyed by complete variant, gene,
    transcript, and peptide-position provenance so Topiary
    DSL expressions can combine multi-predictor predictions via
    ``affinity['mhcflurry']`` / ``affinity['netmhcpan']`` or use an
    unqualified ``affinity`` fallback via the epitope config's
    ``default_methods`` setting.

    Parameters
    ----------
    path : str or Path
    epitope_config : EpitopeConfig, optional
        Effective DSL configuration to evaluate before returning candidates.
        Supplying it here ensures vaccine ranking and template reports consume
        the same scores as the neoepitope CSV/XLSX report.

    Returns
    -------
    pandas.DataFrame
        Report-ready DataFrame (one row per peptide × allele).
    list of CandidateEpitope
        One ``CandidateEpitope`` per external prediction identity, each
        carrying its per-(allele, predictor, kind)
        ``mhctools.Prediction`` records inside ``epitope``. MHCflurry
        inputs may carry distinct ``pMHC_affinity``, ``pMHC_presentation``,
        and allele-independent ``antigen_processing`` leaves.
    """
    # low_memory=False avoids mixed-dtype warnings — LENS reports have many
    # optional columns that are mostly NA for a given antigen type.
    from topiary import read_lens

    # topiary owns LENS parsing: it detects the version, absorbs the
    # per-version column spellings, and normalizes every predictor onto
    # ``<tool>_<kind>_<metric>``. Metadata columns arrive under topiary's
    # names too — ``gene``, ``variant``, ``effect``, ``gene_tpm`` — which is
    # the vocabulary the rest of this module now reads.
    result = read_lens(path)
    df = _map_residual_lens_columns(result.df)
    # topiary detects each binding model's version while parsing and records
    # it on the metadata. vaxrank's version-qualified DSL references
    # (``affinity['mhcflurry', '2.1.1']``) depend on those reaching the
    # Prediction leaves, so read them from there rather than re-deriving.
    lens_versions = dict(getattr(result.metadata, "models", None) or {})

    required = {"peptide", "allele"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"LENS file {path} missing required columns: {missing}")

    detected = detect_lens_predictors(df.columns, lens_versions)
    if not detected:
        raise ValueError(
            f"LENS file {path} has no recognized predictor columns. "
            f"topiary normalizes LENS binding columns to "
            f"'<tool>_<kind>_<metric>' (e.g. 'mhcflurry_affinity_value'); "
            f"none were present after reading.")

    logger.info(
        "LENS file %s (version %s): detected %s",
        path, (result.extra or {}).get("lens_version", "unknown"),
        sorted({(d.tool, d.kind) for d in detected}))

    # Report columns like 'Predicted mutant pMHC affinity' and '%ile rank'
    # render a single predictor's value. Pick the canonical pMHC_affinity
    # predictor for that display role (every per-tool raw value is also
    # exposed in its own column further down). Stability / presentation
    # predictors are never chosen for display.
    affinity_preds = [d for d in detected if d.kind == "pMHC_affinity"]
    display_pred = (
        _pick_canonical_predictor(affinity_preds) if affinity_preds else None)
    chosen = list(detected)

    # One materialization of the source rows. Construct ranking consumes
    # these same dicts, so LENS evidence columns (read counts, VAF, alleles)
    # never require a second read of the file.
    rows = [dict(row) for _, row in df.iterrows()]
    records = []
    epitope_rows = []
    report_rows = []
    processing_rows_by_position = {}
    # Rows whose peptide isn't a substring of a non-empty pep_context —
    # peptide and pep_context came from different isoforms / annotation
    # snapshots upstream. Dropped here (not carried into the report,
    # constructs, or pepsickle); summarized once after the loop.
    n_dropped_peptide_context_mismatch = 0
    mismatch_examples = []
    # (peptide, tool, kept score, conflicting score) for allele rows that
    # disagreed about an allele-independent processing score.
    processing_conflicts = []
    for row in rows:
        peptide = cells.text(row.get("peptide"))
        if not peptide:
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
        agretopicities = {
            d.tool: cells.number(row.get(d.agretopicity_col))
            for d in detected if d.agretopicity_col}
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
        pep_context = cells.text(row.get("pep_context"))
        position = lens_epitope_position(peptide, pep_context)
        if position is None:
            n_dropped_peptide_context_mismatch += 1
            if len(mismatch_examples) < 1:
                mismatch_examples.append((peptide, pep_context))
            continue
        peptide, pep_context, offset = position
        prediction_key = ExternalPredictionKey.from_lens_row(row)
        if prediction_key is None:
            continue
        prediction_id = prediction_key.identifier
        row["prediction_id"] = prediction_id
        records.append(ExternalRecord(key=prediction_key, row=row))
        row_added = False
        shared_epitope_fields = {
            'peptide': str(peptide),
            'source': pep_context,
            'prediction_id': prediction_id,
            'offset': offset,
            'source_class': SOURCE_CLASS_MUTATION,
            # LENS pre-curates rows as neoepitopes and pre-filters
            # against the patient reference proteome — both flags
            # are structurally true for any surviving row.
            'overlaps_mutation': True,
            'occurs_in_reference': False,
            'occurs_in_non_CTA_reference': False,
        }
        # One detected entry per (tool, kind); each fills whichever of
        # value / score / rank topiary's normalized frame carries for it.
        for d in chosen:
            value = cells.number(row.get(d.value_col)) if d.value_col else None
            score = cells.number(row.get(d.score_col)) if d.score_col else None
            rank = cells.number(row.get(d.rank_col)) if d.rank_col else None
            if value is None and score is None and rank is None:
                continue

            if d.kind == "antigen_processing":
                if score is None:
                    continue
                # Processing depends on peptide + flanks, not MHC. LENS
                # repeats it on every allele row, while mhctools emits one
                # allele-less prediction per peptide position; preserve that
                # canonical shape.
                processing_key = (prediction_id, d.tool, d.version)
                existing = processing_rows_by_position.get(processing_key)
                if existing is not None:
                    kept = existing['mutant'].score
                    if abs(kept - score) > PROCESSING_SCORE_TOLERANCE:
                        processing_conflicts.append(
                            (peptide, d.tool, kept, score))
                    existing['patient_alleles'].add(allele)
                else:
                    processing_row = {
                        **shared_epitope_fields,
                        'mutant': Prediction(
                            kind=d.kind, predictor_name=d.tool,
                            predictor_version=d.version, allele="",
                            peptide=str(peptide), value=None, score=score,
                            percentile_rank=None,
                        ),
                        'wt': None,
                        'patient_alleles': {allele},
                    }
                    processing_rows_by_position[processing_key] = (
                        processing_row)
                    epitope_rows.append(processing_row)
                row_added = True
                continue

            # Allele-scoped kinds. ``Prediction.score`` is required, so a
            # signal that only carries a percentile keeps a neutral 0.0
            # rather than discarding the percentile evidence.
            mutant = Prediction(
                kind=d.kind, predictor_name=d.tool,
                predictor_version=d.version, allele=allele,
                peptide=str(peptide), value=value,
                score=score if score is not None else 0.0,
                percentile_rank=rank,
            )

            # Agretopicity is MT/WT on this tool's own value scale, so it
            # recovers a WT comparator only for that tool's affinity — LENS
            # emits no WT column of its own.
            wt = None
            agretopicity = agretopicities.get(d.tool)
            if (d.kind == "pMHC_affinity" and value is not None
                    and agretopicity is not None and agretopicity > 0):
                # No WT peptide sequence exists upstream; an empty sequence
                # is the anonymous-WT signal candidate_epitopes_from_rows
                # keeps because the value is non-None.
                wt = Prediction(
                    kind=d.kind, predictor_name=d.tool,
                    predictor_version=d.version, allele=allele,
                    peptide="", value=value / agretopicity, score=0.0,
                    percentile_rank=None,
                )

            epitope_rows.append({
                **shared_epitope_fields, 'mutant': mutant, 'wt': wt})
            row_added = True

        if not row_added:
            continue

        report_rows.append(_build_lens_report_row(
            row=row,
            allele=allele,
            peptide=peptide,
            prediction_id=prediction_id,
            source_sequence_name=pep_context or str(peptide),
            peptide_offset=offset,
            detected=detected,
            display_pred=display_pred,
            chosen_tools=sorted({d.tool for d in chosen}),
        ))

    if processing_conflicts:
        ex_peptide, ex_tool, ex_kept, ex_other = processing_conflicts[0]
        logger.warning(
            "%d peptide/predictor pair(s) reported different %s processing "
            "scores across their allele rows. That score depends on peptide "
            "and flanks, not on MHC, so the repeats should agree; most often "
            "the file printed one value at two precisions. Kept the "
            "first-seen value. Example: peptide %s, %s scored %r on one row "
            "and %r on another.",
            len(processing_conflicts), ex_tool, ex_peptide, ex_tool,
            ex_kept, ex_other)

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
    # See read_pvacseq_report for the parallel rationale: the 3.1 refactor moved
    # per-(peptide, allele) scoring to the DSL and stored the result on
    # CandidateEpitope.per_allele_scores. Loaders must populate it or
    # downstream ranking sees zero-scored epitopes and drops everything.
    from .epitope_dsl import attach_per_allele_scores
    epitopes = attach_per_allele_scores(epitopes, epitope_config)
    logger.info(
        "Loaded %d epitope(s) (%d row(s) × %d predictor(s)) from %s",
        len(epitopes), len(report_df), len(chosen), path)
    return ExternalReport(
        source_format="lens",
        path=str(path),
        report_df=report_df,
        epitopes=tuple(epitopes),
        records=tuple(records),
        rows=tuple(rows),
    )


def _build_lens_report_row(row, allele, peptide, prediction_id,
                           source_sequence_name, peptide_offset, detected,
                           display_pred, chosen_tools):
    """Assemble one row of the LENS neoepitope report DataFrame.

    ``display_pred`` (if not None) drives the legacy Affinity / %ile
    columns; every *detected* predictor's raw value is also exposed as
    '<Tool> value (nM/hours)' so downstream DSL / scripts can see both.
    """
    antigen_source = cells.text(row.get("antigen_source"))
    gene = cells.text(row.get("gene"))
    variant_pos = cells.text(row.get("variant"))
    if not variant_pos:
        variant_pos = cells.text(row.get("mut_aa_pos"))
    if not variant_pos:
        variant_pos = cells.text(row.get("origin_descriptor"))

    display_value = (
        cells.number(row.get(display_pred.value_col))
        if display_pred and display_pred.value_col else None)
    display_percentile = (
        cells.number(row.get(display_pred.rank_col))
        if display_pred and display_pred.rank_col else None)
    display_agretopicity = (
        cells.number(row.get(display_pred.agretopicity_col))
        if display_pred and display_pred.agretopicity_col else None)

    # Shared core columns (WT affinity isn't predicted on the import
    # path → None → blank); then LENS-specific extras.
    out = neoepitope_core_row(
        allele=allele, mutant_peptide=str(peptide),
        mutant_affinity=display_value, wt_peptide='', wt_affinity=None,
        gene_name=gene, variant=variant_pos)
    out.update({
        'Prediction identity': prediction_id,
        'Antigen source': antigen_source,
        'Predictors used': ','.join(chosen_tools),
        '%ile rank': display_percentile,
        'Agretopicity': display_agretopicity,
        'TPM': cells.number(row.get("gene_tpm")),
        'Source sequence name': source_sequence_name,
        'Peptide offset': peptide_offset,
    })
    # Every detected predictor gets a raw-value column so users can see
    # per-tool signals side-by-side in the report, even if not chosen
    # for DSL scoring.
    unit_by_kind = {"pMHC_affinity": "nM", "pMHC_stability": "hours"}
    for d in detected:
        unit = unit_by_kind.get(d.kind, "")
        if d.value_col:
            label = f"{d.tool} value ({unit})" if unit else f"{d.tool} value"
            out[label] = cells.number(row.get(d.value_col))
        if d.kind == "pMHC_affinity" and d.rank_col:
            out[f"{d.tool} affinity percentile rank"] = cells.number(
                row.get(d.rank_col))
        if d.kind == "pMHC_presentation":
            if d.score_col:
                out[f"{d.tool} presentation score"] = cells.number(
                    row.get(d.score_col))
            if d.rank_col:
                out[f"{d.tool} presentation percentile rank"] = cells.number(
                    row.get(d.rank_col))
        if d.kind == "antigen_processing" and d.score_col:
            out[f"{d.tool} processing score"] = cells.number(
                row.get(d.score_col))
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
        Report-ready DataFrame from an external report reader
        (one row per (peptide, allele) pair).
    epitopes : list of CandidateEpitope
        Output of ``read_lens_report`` / ``read_pvacseq_report`` — one per
        complete external prediction identity.
    excel_report_path : str, optional
    csv_report_path : str, optional
    epitope_config : EpitopeConfig, optional
    topiary_df : pandas.DataFrame, optional
        Pre-built topiary long-form rows to use for DSL validation and
        scoring. ``read_pvacseq_report`` stores this on ``report_df.attrs`` so
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

    # Duplicate peptide/allele pairs are expected across variants and
    # transcripts; the identity join keeps their scores distinct.
    if len(report_df) > 0 and logger.isEnabledFor(logging.DEBUG):
        dup_count = int(report_df.duplicated(
            subset=[peptide_col, allele_col]).sum())
        if dup_count:
            logger.debug(
                "report_df has %d duplicate (peptide, allele) row(s) "
                "from multi-source input; explicit prediction identities "
                "keep their scores source-specific.", dup_count)

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

    # score_series is indexed by the stable prediction identity, peptide,
    # offset, and allele. If an exact row is absent, the DSL filtered it out.
    # There is no (peptide, allele) fallback: it broadcast one source's score
    # onto every row sharing a sequence, which is the merge the identity
    # exists to prevent. A frame without the identity columns is refused
    # below rather than silently broadcast.
    scores_by_key = {}
    for idx_tuple, score in score_series.items():
        prediction_id, peptide, offset, allele = idx_tuple
        scores_by_key[(prediction_id, peptide, int(offset), allele)] = score

    report_df = report_df.copy()
    identity_col = 'Prediction identity'
    offset_col = 'Peptide offset'
    missing = [
        column for column in (identity_col, offset_col)
        if column not in report_df.columns]
    if missing:
        raise ValueError(
            "Neoepitope report frame is missing %s. Scores join on the exact "
            "prediction identity; a (peptide, allele) fallback would "
            "broadcast one source's score onto another source that happens "
            "to share a sequence." % " and ".join(missing))
    scores = []
    filter_passed = []
    rank_eligible = []
    exclusion_reasons = []
    for _, row in report_df.iterrows():
        offset = cells.number(row.get(offset_col))
        score_key = (
            (row.get(identity_col), row[peptide_col], int(offset),
             row[allele_col])
            if offset is not None else None)
        passed = score_key in scores_by_key
        score = scores_by_key.get(score_key)
        eligible = (
            passed
            and score is not None
            and float(score) >= epitope_config.min_epitope_score
        )
        scores.append(round(float(score), 6) if passed else None)
        filter_passed.append(passed)
        rank_eligible.append(eligible)
        exclusion_reasons.append(
            "" if eligible else (
                "dsl_filter" if not passed else "min_epitope_score"))
    report_df.insert(2, 'vaxrank_score', scores)
    report_df.insert(3, 'vaxrank_filter_passed', filter_passed)
    report_df.insert(4, 'vaxrank_rank_eligible', rank_eligible)
    report_df.insert(5, 'vaxrank_exclusion_reason', exclusion_reasons)

    report_df = report_df.sort_values(
        'vaxrank_score', ascending=False, na_position='last')
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
