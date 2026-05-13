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


import math
import traceback
import logging
from typing import Optional

from mhctools.pred import Prediction
from pyensembl import Genome
from topiary import TopiaryPredictor
from topiary.ranking import EvalContext, apply_filter

from .config.defaults import DEFAULT_MIN_KMER_LENGTH
from .epitope_config import EpitopeConfig
from .epitope_dsl import build_filter_node, build_score_node
from .mutant_protein_fragment import MutantProteinFragment
from .candidate_epitope import (
    CandidateEpitope, SOURCE_CLASS_MUTATION,
    candidate_epitopes_from_rows,
)
from .reference_proteome import ReferenceProteome

logger = logging.getLogger(__name__)


def _finite_or_none(x):
    """Coerce ``None`` / NaN / Inf floats to ``None``.

    Topiary frames legitimately carry NaN in the ``value`` column for
    non-affinity kinds (``pMHC_presentation`` carries its probability
    in ``score``, not ``value``). Pass-through writes NaN into
    ``Prediction.value`` which then poisons every downstream consumer
    — sorting, scoring, and especially JSON serialization (simplejson
    rejects NaN / Inf for strict compliance). Producers should emit
    ``None`` to mean "no IC50 for this kind", not ``NaN``.
    """
    if x is None:
        return None
    if isinstance(x, float) and (math.isnan(x) or math.isinf(x)):
        return None
    return x


def slice_epitopes(epitopes, start_offset, end_offset):
    """Return subset of ``CandidateEpitope`` objects whose mutant peptide lies
    fully within ``[start_offset, end_offset)``, with each one's source
    window narrowed to that range and offset rebased.

    Drop-in replacement for the legacy ``slice_epitope_predictions``
    that operated on pre-3.0 flat per-(peptide, allele) records.
    """
    sliced = (e.sliced(start_offset, end_offset) for e in epitopes)
    return [e for e in sliced if e is not None]

def predict_epitopes(
        mhc_predictor,
        protein_fragment : MutantProteinFragment,
        epitope_config : Optional[EpitopeConfig] = None,
        genome : Optional[Genome] = None) -> list[CandidateEpitope]:
    """
    Parameters
    ----------
    mhc_predictor
        Either a topiary.TopiaryPredictor (which can wrap multiple
        models) or a mhctools BasePredictor. A bare BasePredictor is
        wrapped in a single-model TopiaryPredictor; pass a multi-model
        TopiaryPredictor (or use ``--mhc-predictor`` with multiple
        names on the CLI) to get all alleles' / predictors' results
        rolled into each CandidateEpitope's mutant context.

    protein_fragment
        Protein sub-sequence to run MHC binding predictor over

    epitope_config
        Configuration object with parameters for scoring epitopes, if
        missing then default values are used

    genome
        Genome whose proteome to use for reference peptide filtering

    Returns
    -------
    list[CandidateEpitope]
        One ``CandidateEpitope`` per (peptide, peptide_offset) — each carrying
        all per-(allele, predictor) ``mhctools.Prediction`` records
        in its mutant ``Peptide``, plus the WT comparator
        context when the peptide overlaps the mutation.

    Uses the input genome to evaluate whether the epitope occurs in reference.
    """
    if epitope_config is None:
        epitope_config = EpitopeConfig()

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
        return []

    if predictions_df.empty:
        return []

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
            return []

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

    # Walk topiary frame rows; build per-(allele, predictor) leaf
    # ``Prediction`` records into a flat list of row dicts.
    # ``candidate_epitopes_from_rows`` groups them by
    # (peptide, source, offset) into one CandidateEpitope per
    # peptide position at the end.
    rows = []
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

        # reindex + fillna above guarantee every group tuple is in the series.
        group_key = (
            row["source_sequence_name"],
            peptide,
            peptide_start_offset,
            row["allele"],
        )
        epitope_score = float(score_series[group_key])
        if epitope_score < epitope_config.min_epitope_score:
            num_low_scoring += 1
            continue

        # IC50 value: use the "affinity" column when present (affinity
        # rows), otherwise fall back to "value". For non-affinity kinds
        # (e.g. pMHC_presentation), both columns are legitimately NaN
        # in the topiary frame — coerce to None so the Prediction
        # carries an honest "no IC50" instead of a NaN poison pill.
        ic50 = _finite_or_none(row.get("affinity"))
        if ic50 is None:
            ic50 = _finite_or_none(row.get("value"))

        percentile_rank = _finite_or_none(row.get("percentile_rank"))

        # Resolve WT comparator only when the peptide overlaps the
        # mutation — non-overlapping peptides aren't neoepitopes and
        # don't get a meaningful WT pair.
        wt_pred = None
        if overlaps_mutation:
            wt_peptide = wt_peptides[peptide]
            wt_row = wt_predictions_grouped.get((wt_peptide, row["allele"]))
            if wt_row is None:
                if len(wt_peptide) < min_peptide_length:
                    logger.info(
                        'No prediction for too-short WT epitope %s: possible stop-loss variant',
                        wt_peptide)
            else:
                wt_ic50 = _finite_or_none(wt_row.get("affinity"))
                if wt_ic50 is None:
                    wt_ic50 = _finite_or_none(wt_row.get("value"))
                tool = row.get("prediction_method_name", "")
                wt_pred = Prediction(
                    kind=row.get("kind") or "pMHC_affinity",
                    predictor_name=tool,
                    predictor_version=row.get("predictor_version", "") or "",
                    allele=row["allele"],
                    peptide=wt_peptide,
                    value=wt_ic50,
                    score=0.0,
                    percentile_rank=None,
                )

        tool = row.get("prediction_method_name", "")
        mutant_pred = Prediction(
            kind=row.get("kind") or "pMHC_affinity",
            predictor_name=tool,
            predictor_version=row.get("predictor_version", "") or "",
            allele=row["allele"],
            peptide=peptide,
            value=ic50,
            score=0.0,
            percentile_rank=percentile_rank,
        )
        rows.append({
            'peptide': peptide,
            'source': protein_fragment.amino_acids,
            'offset': peptide_start_offset,
            'mutant': mutant_pred,
            'wt': wt_pred,
            'source_class': SOURCE_CLASS_MUTATION,
            'overlaps_mutation': overlaps_mutation,
            'occurs_in_reference': occurs_in_reference,
            # ReferenceProteome doesn't yet know about CTAs, so the
            # CTA-aware flag mirrors the raw one. When the CTA set
            # is populated (via pirlygenes), this branch will diverge
            # for CTA-matching peptides.
            'occurs_in_non_CTA_reference': occurs_in_reference,
            # Per-(peptide, allele) DSL score; threaded onto the
            # CandidateEpitope as ``per_allele_scores[allele]`` by
            # ``candidate_epitopes_from_rows``. Single source of
            # truth — VaccinePeptide reads it directly and never
            # recomputes from ic50 / percentile_rank.
            'allele_score': epitope_score,
        })

    logger.info(
        "%d total peptides: %d occur in reference, %d failed score threshold",
        num_total,
        num_occurs_in_reference,
        num_low_scoring)
    return candidate_epitopes_from_rows(rows)
