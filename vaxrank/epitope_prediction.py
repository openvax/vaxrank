# Copyright (c) 2016-2018. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import absolute_import, print_function, division
from collections import namedtuple, OrderedDict
import traceback
import logging

import numpy as np

from .reference_proteome import ReferenceProteome

logger = logging.getLogger(__name__)

EpitopePredictionBase = namedtuple("EpitopePrediction", [
    "allele",
    "peptide_sequence",
    "wt_peptide_sequence",
    "length",
    "ic50",
    "wt_ic50",
    "percentile_rank",
    "prediction_method_name",
    "overlaps_mutation",
    "source_sequence",
    "offset",
    "occurs_in_reference",
])


class EpitopePrediction(EpitopePredictionBase):

    def logistic_epitope_score(
            self,
            midpoint=350.0,
            width=150.0,
            ic50_cutoff=5000.0):  # TODO: add these default values into CLI as arguments
        """
        Map from IC50 values to score where 1.0 = strong binder, 0.0 = weak binder
        Default midpoint and width for logistic determined by max likelihood fit
        for data from Alessandro Sette's 1994 paper:

           "The relationship between class I binding affinity
            and immunogenicity of potential cytotoxic T cell epitopes.

        TODO: Use a large dataset to find MHC binding range predicted to #
        correlate with immunogenicity
        """
        if self.ic50 >= ic50_cutoff:
            return 0.0

        rescaled = (float(self.ic50) - midpoint) / width
        # simplification of 1.0 - logistic(x) = logistic(-x)
        logistic = 1.0 / (1.0 + np.exp(rescaled))

        # since we're scoring IC50 values, let's normalize the output
        # so IC50 near 0.0 always returns a score of 1.0
        normalizer = 1.0 / (1.0 + np.exp(-midpoint / width))

        return logistic / normalizer


def predict_epitopes(
        mhc_predictor,
        protein_fragment,
        min_epitope_score=0.0,
        genome=None):
    """
    Parameters
    ----------
    mhc_predictor : mhctools.BasePredictor
        Object with predict_peptides method

    protein_fragment : MutantProteinFragment

    peptide_length : list of int
        Lengths of peptides to make pMHC binding predictions for

    min_epitope_score : float
        Ignore peptides with binding predictions whose normalized score is less
        than this.

    genome : pyensembl.Genome
        Genome whose proteome to use for reference peptide filtering

    Returns an OrderedDict of EpitopePrediction objects, keyed by a
    (peptide sequence, allele) tuple, that have a normalized score greater
    than min_epitope_score.

    Uses the input genome to evaluate whether the epitope occurs in reference.
    """
    results = OrderedDict()
    reference_proteome = ReferenceProteome(genome)

    # sometimes the predictors will fail, and we don't want to crash vaxrank in that situation
    # TODO: make more specific or remove when we fix error handling in mhctools
    try:
        mhctools_binding_predictions = mhc_predictor.predict_subsequences(
            {protein_fragment.gene_name: protein_fragment.amino_acids})
    except:
        logger.error(
            'MHC prediction errored for protein fragment %s, with traceback: %s',
            protein_fragment, traceback.format_exc())
        return results

    # compute the WT epitopes for each mutant fragment's epitopes; mutant -> WT
    wt_peptides = {}
    for binding_prediction in mhctools_binding_predictions:
        peptide = binding_prediction.peptide
        peptide_length = binding_prediction.length
        peptide_start_offset = binding_prediction.offset
        peptide_end_offset = peptide_start_offset + peptide_length

        overlaps_mutation = protein_fragment.interval_overlaps_mutation(
            start_offset=peptide_start_offset,
            end_offset=peptide_end_offset)

        if overlaps_mutation:
            full_reference_protein_sequence = (
                protein_fragment.predicted_effect().original_protein_sequence
            )
            global_epitope_start_pos = (
                protein_fragment.global_start_pos() + peptide_start_offset
            )
            wt_peptide = full_reference_protein_sequence[
                global_epitope_start_pos:global_epitope_start_pos + peptide_length]
            wt_peptides[peptide] = wt_peptide

    wt_predictions = []
    try:
        # filter to minimum peptide lengths
        valid_wt_peptides = [
            x for x in wt_peptides.values() if len(x) >= mhc_predictor.min_peptide_length
        ]
        if len(valid_wt_peptides) > 0:
            wt_predictions = mhc_predictor.predict_peptides(valid_wt_peptides)
    except:
        logger.error(
            'MHC prediction for WT peptides errored, with traceback: %s',
            traceback.format_exc())

    wt_predictions_grouped = {}
    # break it out: (peptide, allele) -> prediction
    for wt_prediction in wt_predictions:
        wt_predictions_grouped[(wt_prediction.peptide, wt_prediction.allele)] = wt_prediction

    # convert from mhctools.BindingPrediction objects to EpitopePrediction
    # which differs primarily by also having a boolean field
    # 'overlaps_mutation' that indicates whether the epitope overlaps
    # mutant amino acids or both sides of a deletion
    num_total = 0
    num_occurs_in_reference = 0
    num_low_scoring = 0
    for binding_prediction in mhctools_binding_predictions:
        num_total += 1
        peptide = binding_prediction.peptide
        peptide_length = binding_prediction.length
        peptide_start_offset = binding_prediction.offset
        peptide_end_offset = peptide_start_offset + peptide_length

        overlaps_mutation = protein_fragment.interval_overlaps_mutation(
            start_offset=peptide_start_offset,
            end_offset=peptide_end_offset)

        occurs_in_reference = reference_proteome.contains(peptide)
        if occurs_in_reference:
            logger.debug('Peptide %s occurs in reference', peptide)
            num_occurs_in_reference += 1

        # compute WT epitope sequence, if this epitope overlaps the mutation
        if overlaps_mutation:
            wt_peptide = wt_peptides[peptide]
            wt_prediction = wt_predictions_grouped.get((wt_peptide, binding_prediction.allele))
            wt_ic50 = None
            if wt_prediction is None:
                # this can happen in a stop-loss variant: do we want to check that here?
                if len(wt_peptide) < mhc_predictor.min_peptide_length:
                    logger.info(
                        'No prediction for too-short WT epitope %s: possible stop-loss variant',
                        wt_peptide)
            else:
                wt_ic50 = wt_prediction.value

        else:
            wt_peptide = peptide
            wt_ic50 = binding_prediction.value

        epitope_prediction = EpitopePrediction(
            allele=binding_prediction.allele,
            peptide_sequence=peptide,
            wt_peptide_sequence=wt_peptide,
            length=len(peptide),
            ic50=binding_prediction.value,
            wt_ic50=wt_ic50,
            percentile_rank=binding_prediction.percentile_rank,
            prediction_method_name=binding_prediction.prediction_method_name,
            overlaps_mutation=overlaps_mutation,
            source_sequence=protein_fragment.amino_acids,
            offset=peptide_start_offset,
            occurs_in_reference=occurs_in_reference)
        if epitope_prediction.logistic_epitope_score() >= min_epitope_score:
            key = (epitope_prediction.peptide_sequence, epitope_prediction.allele)
            results[key] = epitope_prediction
        else:
            num_low_scoring += 1

    logger.info(
        "%d total peptides: %d occur in reference, %d failed score threshold",
        num_total,
        num_occurs_in_reference,
        num_low_scoring)
    return results


def slice_epitope_predictions(
        epitope_predictions,
        start_offset,
        end_offset):
    """
    Return subset of EpitopePrediction objects which overlap the given interval
    and slice through their source sequences and adjust their offset.
    """
    return [
        EpitopePrediction(
            allele=p.allele,
            peptide_sequence=p.peptide_sequence,
            wt_peptide_sequence=p.wt_peptide_sequence,
            length=p.length,
            ic50=p.ic50,
            wt_ic50=p.wt_ic50,
            percentile_rank=p.percentile_rank,
            prediction_method_name=p.prediction_method_name,
            overlaps_mutation=p.overlaps_mutation,
            source_sequence=p.source_sequence[start_offset:end_offset],
            offset=p.offset - start_offset,
            occurs_in_reference=p.occurs_in_reference)
        for p in epitope_predictions
        if p.offset >= start_offset and p.offset + p.length <= end_offset
    ]
