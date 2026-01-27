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


from collections import OrderedDict
import traceback
import logging


from pyensembl import Genome
from mhctools.base_predictor import BasePredictor


from .epitope_config import EpitopeConfig
from .epitope_prediction import EpitopePrediction
from .mutant_protein_fragment import MutantProteinFragment
from .reference_proteome import ReferenceProteome

logger = logging.getLogger(__name__)


def slice_epitope_predictions(
        epitope_predictions,
        start_offset,
        end_offset):
    """
    Return subset of EpitopePrediction objects which overlap the given interval
    and slice through their source sequences and adjust their offset.
    """
    return [
        p.slice_source_sequence(start_offset, end_offset)
        for p in epitope_predictions
        if p.offset >= start_offset and p.offset + p.length <= end_offset
    ]

Peptide = str
Allele = str

def predict_epitopes(
        mhc_predictor : BasePredictor,
        protein_fragment : MutantProteinFragment,
        epitope_config : EpitopeConfig = None,
        genome :  Genome = None) -> dict[tuple[Peptide, Allele], EpitopePrediction]:
    """
    Parameters
    ----------
    mhc_predictor
        Object with predict_peptides method

    protein_fragment 
        Protein sub-sequence to run MHC binding predictor over 

    epitope_config
        Configuration object with parameters for scoring epitopes, if 
        missing then default values are used

    genome
        Genome whose proteome to use for reference peptide filtering

    Returns an OrderedDict of EpitopePrediction objects, keyed by a
    (peptide sequence, allele) tuple, that have a normalized score greater
    than min_epitope_score.

    Uses the input genome to evaluate whether the epitope occurs in reference.
    """
    if epitope_config is None:
        epitope_config = EpitopeConfig()

    results = OrderedDict()
    reference_proteome = ReferenceProteome(genome)

    # sometimes the predictors will fail, and we don't want to crash vaxrank
    # in that situation
    # TODO: make more specific or remove when we fix error handling in mhctools
    try:
        mhctools_binding_predictions = mhc_predictor.predict_subsequences(
            {protein_fragment.gene_name: protein_fragment.amino_acids})
    except Exception:
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
    # filter to minimum peptide lengths
    valid_wt_peptides = [
        x for x in wt_peptides.values() if len(x) >= mhc_predictor.min_peptide_length
    ]
    if len(valid_wt_peptides) > 0:

        try:
            wt_predictions = mhc_predictor.predict_peptides(valid_wt_peptides)
        except Exception:
            logger.error(
                'MHC prediction for WT peptides errored, with traceback: %s',
                traceback.format_exc())

    # break it out: (peptide, allele) -> prediction
    wt_predictions_grouped = {
        (wt_prediction.peptide, wt_prediction.allele): wt_prediction
        for wt_prediction in wt_predictions
    }

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
            wt_prediction = wt_predictions_grouped.get(
                (wt_peptide, binding_prediction.allele))
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
            ic50=binding_prediction.value,
            wt_ic50=wt_ic50,
            percentile_rank=binding_prediction.percentile_rank,
            prediction_method_name=binding_prediction.prediction_method_name,
            overlaps_mutation=overlaps_mutation,
            source_sequence=protein_fragment.amino_acids,
            offset=peptide_start_offset,
            occurs_in_reference=occurs_in_reference)
        epitope_score = epitope_prediction.logistic_epitope_score(
            midpoint=epitope_config.logistic_epitope_score_midpoint,
            width=epitope_config.logistic_epitope_score_width,
            ic50_cutoff=epitope_config.binding_affinity_cutoff)

        if epitope_score >= epitope_config.min_epitope_score:
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
