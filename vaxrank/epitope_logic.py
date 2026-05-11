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


import traceback
import logging
from typing import Optional

import numpy as np
from pyensembl import Genome
from topiary import TopiaryPredictor
from topiary.ranking import EvalContext, apply_filter

from .config.defaults import DEFAULT_MIN_KMER_LENGTH
from .epitope_config import EpitopeConfig
from .epitope_dsl import (
    _legacy_predictions_to_epitopes,
    build_filter_node,
    build_score_node,
)
from .epitope_prediction import EpitopePrediction
from .mutant_protein_fragment import MutantProteinFragment
from .peptide_context import Epitope
from .reference_proteome import ReferenceProteome

logger = logging.getLogger(__name__)


def slice_epitopes(epitopes, start_offset, end_offset):
    """Return subset of ``Epitope`` objects whose mutant peptide lies
    fully within ``[start_offset, end_offset)``, with each one's source
    window narrowed to that range and offset rebased.

    Drop-in replacement for the legacy ``slice_epitope_predictions``
    that operated on flat ``EpitopePrediction`` records.
    """
    sliced = (e.sliced(start_offset, end_offset) for e in epitopes)
    return [e for e in sliced if e is not None]

def predict_epitopes(
        mhc_predictor,
        protein_fragment : MutantProteinFragment,
        epitope_config : Optional[EpitopeConfig] = None,
        genome : Optional[Genome] = None) -> list[Epitope]:
    """
    Parameters
    ----------
    mhc_predictor
        Either a topiary.TopiaryPredictor (which can wrap multiple
        models) or a mhctools BasePredictor. A bare BasePredictor is
        wrapped in a single-model TopiaryPredictor; pass a multi-model
        TopiaryPredictor (or use ``--mhc-predictor`` with multiple
        names on the CLI) to get all alleles' / predictors' results
        rolled into each Epitope's mutant context.

    protein_fragment
        Protein sub-sequence to run MHC binding predictor over

    epitope_config
        Configuration object with parameters for scoring epitopes, if
        missing then default values are used

    genome
        Genome whose proteome to use for reference peptide filtering

    Returns
    -------
    list[Epitope]
        One ``Epitope`` per (peptide, peptide_offset) — each carrying
        all per-(allele, predictor) ``mhctools.Prediction`` records
        in its mutant ``PeptideContext``, plus the WT comparator
        context when the peptide overlaps the mutation.

    Uses the input genome to evaluate whether the epitope occurs in reference.
    """
    if epitope_config is None:
        epitope_config = EpitopeConfig()

    # Internal scratchpad of flat EpitopePrediction records; converted
    # to ``list[Epitope]`` at the end via the migration adapter so this
    # producer's body doesn't need a structural rewrite. (Issue #284 —
    # adapter drops out once the body emits Predictions / contexts
    # natively.)
    results: list[EpitopePrediction] = []
    reference_proteome = ReferenceProteome(genome)

    # Wrap bare mhctools predictors in a TopiaryPredictor
    if not isinstance(mhc_predictor, TopiaryPredictor):
        topiary_predictor = TopiaryPredictor(models=[mhc_predictor])
    else:
        topiary_predictor = mhc_predictor

    # Run predictions via Topiary
    try:
        predictions_df = topiary_predictor.predict_from_named_sequences(
            {protein_fragment.gene_name: protein_fragment.amino_acids})
    except Exception:
        logger.error(
            'MHC prediction errored for protein fragment %s, with traceback: %s',
            protein_fragment, traceback.format_exc())
        return results

    if predictions_df.empty:
        return results

    # ``default_methods`` (per-kind ``prediction_method_name`` defaults
    # for unqualified DSL refs) is required when multi-predictor data
    # is in play — without it, ``affinity.value`` is ambiguous and
    # topiary raises. Single-predictor runs leave ``default_methods``
    # unset and the DSL resolves the only available method.
    default_methods = getattr(epitope_config, 'default_methods', None) or None

    # Apply the configured filter (default: affinity or percentile-rank cutoff)
    # so that rows dropped by the filter are never scored or WT-predicted.
    filter_node = build_filter_node(epitope_config)
    if filter_node is not None:
        predictions_df = apply_filter(
            predictions_df, filter_node, default_methods=default_methods)
        if predictions_df.empty:
            return results

    # Evaluate the score expression once; indexed by
    # (source_sequence_name, peptide, peptide_offset, allele) group tuple.
    score_node = build_score_node(epitope_config)
    score_ctx = EvalContext(predictions_df, default_methods=default_methods)
    score_series = (
        score_node.eval(score_ctx).reindex(score_ctx.group_index).fillna(0.0)
    )

    # Compute WT epitopes for peptides that overlap the mutation
    wt_peptides = {}
    for _, row in predictions_df.iterrows():
        peptide = row["peptide"]
        peptide_start_offset = row["peptide_offset"]
        peptide_length = row["peptide_length"]
        peptide_end_offset = peptide_start_offset + peptide_length

        overlaps_mutation = protein_fragment.interval_overlaps_mutation(
            start_offset=peptide_start_offset,
            end_offset=peptide_end_offset)

        if overlaps_mutation and peptide not in wt_peptides:
            full_reference_protein_sequence = (
                protein_fragment.predicted_effect().original_protein_sequence
            )
            global_epitope_start_pos = (
                protein_fragment.global_start_pos() + peptide_start_offset
            )
            wt_peptide = full_reference_protein_sequence[
                global_epitope_start_pos:global_epitope_start_pos + peptide_length]
            wt_peptides[peptide] = wt_peptide

    # Predict binding for WT peptides via Topiary
    wt_predictions_grouped = {}
    min_peptide_length = min(predictions_df["peptide_length"]) if len(predictions_df) > 0 else DEFAULT_MIN_KMER_LENGTH
    valid_wt_peptides = {
        f"wt_{i}": seq
        for i, seq in enumerate(set(wt_peptides.values()))
        if len(seq) >= min_peptide_length
    }
    if valid_wt_peptides:
        try:
            wt_df = topiary_predictor.predict_from_named_sequences(valid_wt_peptides)
            # Index WT predictions by (peptide, allele)
            for _, row in wt_df.iterrows():
                key = (row["peptide"], row["allele"])
                wt_predictions_grouped[key] = row
        except Exception:
            logger.error(
                'MHC prediction for WT peptides errored, with traceback: %s',
                traceback.format_exc())

    # Convert Topiary DataFrame rows to EpitopePrediction objects
    num_total = 0
    num_occurs_in_reference = 0
    num_low_scoring = 0
    for _, row in predictions_df.iterrows():
        num_total += 1
        peptide = row["peptide"]
        peptide_length = row["peptide_length"]
        peptide_start_offset = row["peptide_offset"]
        peptide_end_offset = peptide_start_offset + peptide_length

        overlaps_mutation = protein_fragment.interval_overlaps_mutation(
            start_offset=peptide_start_offset,
            end_offset=peptide_end_offset)

        occurs_in_reference = reference_proteome.contains(peptide)
        if occurs_in_reference:
            logger.debug('Peptide %s occurs in reference', peptide)
            num_occurs_in_reference += 1

        # IC50 value: use the "affinity" column if available, otherwise "value"
        ic50 = row.get("affinity")
        if ic50 is None or (isinstance(ic50, float) and np.isnan(ic50)):
            ic50 = row["value"]

        percentile_rank = row.get("percentile_rank")
        if percentile_rank is not None and isinstance(percentile_rank, float) and np.isnan(percentile_rank):
            percentile_rank = None

        # Compute WT epitope sequence and binding
        if overlaps_mutation:
            wt_peptide = wt_peptides[peptide]
            wt_row = wt_predictions_grouped.get(
                (wt_peptide, row["allele"]))
            wt_ic50 = None
            if wt_row is None:
                if len(wt_peptide) < min_peptide_length:
                    logger.info(
                        'No prediction for too-short WT epitope %s: possible stop-loss variant',
                        wt_peptide)
            else:
                wt_affinity = wt_row.get("affinity")
                if wt_affinity is not None and not (isinstance(wt_affinity, float) and np.isnan(wt_affinity)):
                    wt_ic50 = wt_affinity
                else:
                    wt_ic50 = wt_row["value"]
        else:
            wt_peptide = peptide
            wt_ic50 = ic50

        epitope_prediction = EpitopePrediction(
            allele=row["allele"],
            peptide_sequence=peptide,
            wt_peptide_sequence=wt_peptide,
            ic50=ic50,
            wt_ic50=wt_ic50,
            percentile_rank=percentile_rank,
            prediction_method_name=row.get("prediction_method_name", ""),
            overlaps_mutation=overlaps_mutation,
            source_sequence=protein_fragment.amino_acids,
            offset=peptide_start_offset,
            occurs_in_reference=occurs_in_reference)
        group_key = (
            row["source_sequence_name"],
            peptide,
            peptide_start_offset,
            row["allele"],
        )
        # reindex + fillna above guarantee every group tuple is in the series.
        epitope_score = float(score_series[group_key])

        if epitope_score >= epitope_config.min_epitope_score:
            results.append(epitope_prediction)
        else:
            num_low_scoring += 1

    logger.info(
        "%d total peptides: %d occur in reference, %d failed score threshold",
        num_total,
        num_occurs_in_reference,
        num_low_scoring)
    return _legacy_predictions_to_epitopes(results)
