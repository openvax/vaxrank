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

from __future__ import absolute_import, print_function, division

import logging

from numpy import isclose

from .mutant_protein_fragment import MutantProteinFragment
from .epitope_prediction import predict_epitopes, slice_epitope_predictions
from .vaccine_peptide import VaccinePeptide
from .vaxrank_results import VaxrankResults

logger = logging.getLogger(__name__)

def run_vaxrank(
        isovar_results,
        mhc_predictor,
        vaccine_peptide_length=25,
        max_vaccine_peptides_per_variant=1,
        num_mutant_epitopes_to_keep=10000,
        min_epitope_score=0.0):
    """
    Parameters
    ----------
    isovar_results : list of isovar.IsovarResult
        Each IsovarResult corresponds to one somatic variant and its collection
         of protein sequences determined from RNA.

    mhc_predictor : mhctools.BasePredictor
        Object with predict_peptides method, used for making pMHC binding
        predictions

    vaccine_peptide_length : int
        Length of vaccine SLP to construct

    max_vaccine_peptides_per_variant : int
        Number of vaccine peptides to generate for each mutation.

    num_mutant_epitopes_to_keep : int, optional
        Number of top-ranking epitopes for each vaccine peptide to include in
        computation.

    min_epitope_score : float, optional
        Ignore peptides with binding predictions whose normalized score is less
        than this.
    """
    variant_to_vaccine_peptides_dict = create_vaccine_peptides_dict(
        isovar_results=isovar_results,
        mhc_predictor=mhc_predictor,
        vaccine_peptide_length=vaccine_peptide_length,
        max_vaccine_peptides_per_variant=max_vaccine_peptides_per_variant,
        num_mutant_epitopes_to_keep=num_mutant_epitopes_to_keep,
        min_epitope_score=min_epitope_score)
    ranked_list = ranked_vaccine_peptides(variant_to_vaccine_peptides_dict)

    return VaxrankResults(
        isovar_results=isovar_results,
        variant_to_vaccine_peptides_dict=variant_to_vaccine_peptides_dict,
        ranked_vaccine_peptides=ranked_list)


def create_vaccine_peptides_dict(
        isovar_results,
        mhc_predictor,
        vaccine_peptide_length=25,
        max_vaccine_peptides_per_variant=1,
        num_mutant_epitopes_to_keep=10 ** 5,
        min_epitope_score=0.0):
    """
    Parameters
    ----------
    isovar_results : list of isovar.IsovarResult
        List with one object per variant optionally containing protein sequences

    mhc_predictor : mhctools.BasePredictor
        Object with predict_peptides method, used for making pMHC binding
        predictions

    vaccine_peptide_length : int
        Length of vaccine SLP to construct

    max_vaccine_peptides_per_variant : int
        Number of vaccine peptides to generate for each mutation.

    num_mutant_epitopes_to_keep : int, optional
        Number of top-ranking epitopes for each vaccine peptide to include in
        computation.

    min_epitope_score : float, optional
        Ignore peptides with binding predictions whose normalized score is less
        than this.

    Returns
    -------
    Returns a dictionary of varcode.Variant objects to a list of
    VaccinePeptides.
    """
    vaccine_peptides_dict = {}
    for isovar_result in isovar_results:
        variant = isovar_result.variant
        vaccine_peptides = vaccine_peptides_for_variant(
            isovar_result=isovar_result,
            mhc_predictor=mhc_predictor,
            vaccine_peptide_length=vaccine_peptide_length,
            max_vaccine_peptides_per_variant=max_vaccine_peptides_per_variant,
            num_mutant_epitopes_to_keep=num_mutant_epitopes_to_keep,
            min_epitope_score=min_epitope_score)

        if any(x.contains_mutant_epitopes() for x in vaccine_peptides):
            vaccine_peptides_dict[variant] = vaccine_peptides

    return vaccine_peptides_dict

def vaccine_peptides_for_variant(
        isovar_result,
        mhc_predictor,
        vaccine_peptide_length,
        max_vaccine_peptides_per_variant,
        num_mutant_epitopes_to_keep=None,
        min_epitope_score=0.0):
    """
    Parameters
    ----------
    isovar_result : isovar.IsovarResult

    mhc_predictor : mhctools.BasePredictor
        Object with predict_peptides method, used for making pMHC binding
        predictions

    vaccine_peptide_length : int
        Length of vaccine SLP to construct

    max_vaccine_peptides_per_variant : int
        Number of vaccine peptides to generate for each mutation.

    num_mutant_epitopes_to_keep : int, optional
        Number of top-ranking epitopes for each vaccine peptide to include in
        computation.

    min_epitope_score : float, optional
        Ignore peptides with binding predictions whose normalized score is less
        than this.

    Returns
    -------
    Sorted list of VaccinePeptide objects. If there are no suitable vaccine
    peptides (no strong MHC binder subsequences), returns an empty list.
    """
    if not isovar_result.passes_all_filters:
        # don't consider candidate vaccine peptides from variants which either
        # failed their filters or don't have an RNA-derived protein sequence
        return []

    variant = isovar_result.variant
    long_protein_fragment = MutantProteinFragment.from_isovar_result(isovar_result)

    logger.info(
        "Mutant protein fragment for %s: %s",
        variant,
        long_protein_fragment)

    epitope_predictions = predict_epitopes(
        mhc_predictor=mhc_predictor,
        protein_fragment=long_protein_fragment,
        min_epitope_score=min_epitope_score,
        genome=variant.ensembl).values()

    candidate_vaccine_peptides = []

    for offset, candidate_fragment in long_protein_fragment.sorted_subsequences(
            subsequence_length=vaccine_peptide_length):

        subsequence_epitope_predictions = slice_epitope_predictions(
            epitope_predictions,
            start_offset=offset,
            end_offset=offset + len(candidate_fragment))
        # filter out peptides that have no epitopes
        if not subsequence_epitope_predictions:
            logger.info(
                "No epitope predictions for mutant protein fragment %s",
                candidate_fragment)
            continue

        assert all(
            p.source_sequence == candidate_fragment.amino_acids
            for p in subsequence_epitope_predictions)

        candidate_vaccine_peptide = VaccinePeptide(
            mutant_protein_fragment=candidate_fragment,
            epitope_predictions=subsequence_epitope_predictions,
            num_mutant_epitopes_to_keep=num_mutant_epitopes_to_keep)

        logger.debug(
            "%s, combined score: %0.4f",
            candidate_vaccine_peptide,
            candidate_vaccine_peptide.combined_score)
        candidate_vaccine_peptides.append(candidate_vaccine_peptide)

    n_total_candidates = len(candidate_vaccine_peptides)
    if n_total_candidates == 0:
        logger.info("No candidate peptides for variant %s", variant.short_description)
        return []

    max_score = max(vp.combined_score for vp in candidate_vaccine_peptides)
    if isclose(max_score, 0.0):
        filtered_candidate_vaccine_peptides = candidate_vaccine_peptides
    else:
        # only keep candidate vaccines that are within 1% of the maximum
        # combined score
        filtered_candidate_vaccine_peptides = [
            vp
            for vp in candidate_vaccine_peptides
            if vp.combined_score / max_score > 0.99
        ]
    n_filtered = len(filtered_candidate_vaccine_peptides)
    logger.info(
        "Keeping %d/%d vaccine peptides for %s",
        n_filtered,
        n_total_candidates,
        variant)

    if n_filtered == 0:
        return []

    filtered_candidate_vaccine_peptides.sort(key=VaccinePeptide.lexicographic_sort_key)

    logger.debug("Top vaccine peptides for %s:", variant)
    for i, vaccine_peptide in enumerate(filtered_candidate_vaccine_peptides):
        logger.debug(
            "%d) %s (combined score = %0.4f)",
            i + 1,
            vaccine_peptide,
            vaccine_peptide.combined_score)

    return filtered_candidate_vaccine_peptides[:max_vaccine_peptides_per_variant]


def ranked_vaccine_peptides(variant_to_vaccine_peptides_dict):
    """
    This function returns a sorted list whose first element is a Variant and whose second
    element is a list of VaccinePeptide objects.

    Parameters
    ----------
    variant_to_vaccine_peptides_dict : dict
        Dictionary from varcode.Variant to list of VaccinePeptide

    Returns list of (varcode.Variant, VaccinePeptide list) tuples
    """
    result_list = list(variant_to_vaccine_peptides_dict.items())

    def sort_key(variant_and_vaccine_peptides_pair):
        vaccine_peptides = variant_and_vaccine_peptides_pair[1]
        if len(vaccine_peptides) == 0:
            return 0.0
        else:
            top_vaccine_peptide = vaccine_peptides[0]
            return top_vaccine_peptide.combined_score

    # sort in descending order of combined (expression * mhc binding) scores
    result_list.sort(key=sort_key, reverse=True)
    return result_list
