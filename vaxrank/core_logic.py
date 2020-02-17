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
from .protein_sequences import create_variant_to_protein_sequence_dict

logger = logging.getLogger(__name__)

def run_vaxrank(
        variants,
        reads_generator,
        mhc_predictor,
        vaccine_peptide_length=25,
        padding_around_mutation=5,
        max_vaccine_peptides_per_variant=1,
        min_alt_rna_reads=1,
        min_variant_sequence_coverage=1,
        variant_sequence_assembly=True,
        num_mutant_epitopes_to_keep=10000,
        min_epitope_score=0.0):
    """
    Parameters
    ----------
    variants : VariantCollection
        Variant objects to evaluate for vaccine inclusion

    reads_generator : generator
        Yields sequence of varcode.Variant objects paired with sequences of
        AlleleRead objects that support that variant.

    mhc_predictor : mhctools.BasePredictor
        Object with predict_peptides method, used for making pMHC binding
        predictions

    vaccine_peptide_length : int
        Length of vaccine SLP to construct

    padding_around_mutation : int
        Number of off-center windows around the mutation to consider as vaccine
        peptides.

    max_vaccine_peptides_per_variant : int
        Number of vaccine peptides to generate for each mutation.

    min_alt_rna_reads : int
        Drop variant sequences at loci with fewer than this number of reads
        supporting the alt allele.

    min_variant_sequence_coverage : int
        Trim variant sequences to positions supported by at least this number
        of RNA reads.

    variant_sequence_assembly : bool
        If True, then assemble variant cDNA sequences based on overlap of RNA
        reads. If False, then variant cDNA sequences must be fully spanned and
        contained within RNA reads.

    num_mutant_epitopes_to_keep : int, optional
        Number of top-ranking epitopes for each vaccine peptide to include in
        computation.

    min_epitope_score : float, optional
        Ignore peptides with binding predictions whose normalized score is less
        than this.
    """
    # total number of amino acids is the vaccine peptide length plus the
    # number of off-center windows around the mutation
    protein_fragment_sequence_length = (
            vaccine_peptide_length + 2 * padding_around_mutation)
    variant_to_protein_sequences_dict = create_variant_to_protein_sequence_dict(
        reads_generator=reads_generator,
        sequence_length=protein_fragment_sequence_length,
        min_alt_rna_reads=min_alt_rna_reads,
        min_variant_sequence_coverage=min_variant_sequence_coverage,
        variant_sequence_assembly=variant_sequence_assembly)
    variant_to_vaccine_peptides_dict = create_vaccine_peptides_dict(
        protein_sequence_dict=variant_to_protein_sequences_dict,
        mhc_predictor=mhc_predictor,
        vaccine_peptide_length=vaccine_peptide_length,
        padding_around_mutation=padding_around_mutation,
        max_vaccine_peptides_per_variant=max_vaccine_peptides_per_variant,
        num_mutant_epitopes_to_keep=num_mutant_epitopes_to_keep,
        min_epitope_score=min_epitope_score)
    ranked_list = ranked_vaccine_peptides(variant_to_vaccine_peptides_dict)
    return VaxrankResults(
        variants=variants,
        variant_to_protein_sequences_dict=variant_to_protein_sequences_dict,
        variant_to_vaccine_peptides_dict=variant_to_vaccine_peptides_dict,
        ranked_vaccine_peptides=ranked_list)


def create_vaccine_peptides_dict(
        protein_sequence_dict,
        mhc_predictor,
        vaccine_peptide_length=25,
        padding_around_mutation=10,
        max_vaccine_peptides_per_variant=1,
        num_mutant_epitopes_to_keep=10 ** 5,
        min_epitope_score=0.0):
    """
    Parameters
    ----------
    protein_sequence_dict : dict
        Dictionary from varcode.Variant to isovar ProteinSequence

    mhc_predictor : mhctools.BasePredictor
        Object with predict_peptides method, used for making pMHC binding
        predictions

    vaccine_peptide_length : int
        Length of vaccine SLP to construct

    padding_around_mutation : int
        Number of off-center windows around the mutation to consider as vaccine
        peptides.

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
    for variant, isovar_protein_sequence in protein_sequence_dict.items():
        vaccine_peptides = vaccine_peptides_for_variant(
            variant=variant,
            isovar_protein_sequence=isovar_protein_sequence,
            mhc_predictor=mhc_predictor,
            vaccine_peptide_length=vaccine_peptide_length,
            padding_around_mutation=padding_around_mutation,
            max_vaccine_peptides_per_variant=max_vaccine_peptides_per_variant,
            num_mutant_epitopes_to_keep=num_mutant_epitopes_to_keep,
            min_epitope_score=min_epitope_score)

        if any(x.contains_mutant_epitopes() for x in vaccine_peptides):
            vaccine_peptides_dict[variant] = vaccine_peptides

    return vaccine_peptides_dict

def vaccine_peptides_for_variant(
        variant,
        isovar_protein_sequence,
        mhc_predictor,
        vaccine_peptide_length,
        padding_around_mutation,
        max_vaccine_peptides_per_variant,
        num_mutant_epitopes_to_keep=None,
        min_epitope_score=0.0):
    """
    Parameters
    ----------
    variant : varcode.Variant

    isovar_protein_sequence : isovar. ProteinSequence

    mhc_predictor : mhctools.BasePredictor
        Object with predict_peptides method, used for making pMHC binding
        predictions

    vaccine_peptide_length : int
        Length of vaccine SLP to construct

    padding_around_mutation : int
        Number of off-center windows around the mutation to consider as vaccine
        peptides.

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

    At this point, we know the variant has RNA support, as per the
    isovar_protein_sequence.
    """
    protein_fragment = MutantProteinFragment.from_isovar_protein_sequence(
        variant=variant,
        protein_sequence=isovar_protein_sequence)

    logger.info(
        "Mutant protein fragment for %s: %s",
        variant,
        protein_fragment)

    epitope_predictions = predict_epitopes(
        mhc_predictor=mhc_predictor,
        protein_fragment=protein_fragment,
        min_epitope_score=min_epitope_score,
        genome=variant.ensembl).values()

    candidate_vaccine_peptides = []

    for offset, candidate_fragment in protein_fragment.top_k_subsequences(
            subsequence_length=vaccine_peptide_length,
            k=2 * padding_around_mutation + 1):

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
