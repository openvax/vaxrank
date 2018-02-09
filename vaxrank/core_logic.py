# Copyright (c) 2016. Mount Sinai School of Medicine
#
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
from collections import defaultdict
import logging

from numpy import isclose

from isovar.protein_sequences import (
    reads_generator_to_protein_sequences_generator,
)

from .mutant_protein_fragment import MutantProteinFragment
from .epitope_prediction import predict_epitopes, slice_epitope_predictions
from .vaccine_peptide import VaccinePeptide


logger = logging.getLogger(__name__)


def vaccine_peptides_for_variant(
        variant,
        isovar_protein_sequence,
        mhc_predictor,
        vaccine_peptide_length,
        padding_around_mutation,
        max_vaccine_peptides_per_variant,
        min_epitope_score,
        num_mutant_epitopes_to_keep):
    """
    Returns sorted list of VaccinePeptide objects.
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

    logger.info(
        "Keeping %d/%d vaccine peptides for %s",
        len(filtered_candidate_vaccine_peptides),
        n_total_candidates,
        variant)
    filtered_candidate_vaccine_peptides.sort(key=VaccinePeptide.lexicographic_sort_key)

    if len(filtered_candidate_vaccine_peptides) > 0:
        logger.debug("Top vaccine peptides for %s:", variant)
        for i, vaccine_peptide in enumerate(filtered_candidate_vaccine_peptides):
            logger.debug(
                "%d) %s (combined score = %0.4f)",
                i + 1,
                vaccine_peptide,
                vaccine_peptide.combined_score)
    return filtered_candidate_vaccine_peptides[:max_vaccine_peptides_per_variant]

def generate_vaccine_peptides(
        reads_generator,
        mhc_predictor,
        vaccine_peptide_length,
        padding_around_mutation,
        max_vaccine_peptides_per_variant,
        min_alt_rna_reads,
        min_variant_sequence_coverage,
        variant_sequence_assembly,
        num_mutant_epitopes_to_keep=10000,
        min_epitope_score=0):
    """
    Returns a tuple of two values:
    - dictionary mapping each variant to list of VaccinePeptide objects
    - dictionary containing some variant counts for report display
    """

    # total number of amino acids is the vaccine peptide length plus the
    # number of off-center windows around the mutation
    protein_fragment_sequence_length = (
        vaccine_peptide_length + 2 * padding_around_mutation)

    protein_sequences_generator = reads_generator_to_protein_sequences_generator(
        reads_generator,
        transcript_id_whitelist=None,
        protein_sequence_length=protein_fragment_sequence_length,
        min_alt_rna_reads=min_alt_rna_reads,
        min_variant_sequence_coverage=min_variant_sequence_coverage,
        variant_sequence_assembly=variant_sequence_assembly,
        max_protein_sequences_per_variant=1)

    result_dict = {}
    counts_dict = defaultdict(int)
    for variant, isovar_protein_sequences in protein_sequences_generator:
        if len(variant.effects().drop_silent_and_noncoding()) > 0:
            counts_dict['num_coding_effect_variants'] += 1
        isovar_protein_sequences = list(isovar_protein_sequences)
        if len(isovar_protein_sequences) == 0:
            # this means the variant RNA support is below threshold
            logger.info("No protein sequences for %s", variant)
            continue

        # use the first protein sequence - why?
        counts_dict['num_variants_with_rna_support'] += 1
        isovar_protein_sequence = isovar_protein_sequences[0]
        vaccine_peptides = vaccine_peptides_for_variant(
            variant=variant,
            isovar_protein_sequence=isovar_protein_sequence,
            mhc_predictor=mhc_predictor,
            vaccine_peptide_length=vaccine_peptide_length,
            padding_around_mutation=padding_around_mutation,
            max_vaccine_peptides_per_variant=max_vaccine_peptides_per_variant,
            num_mutant_epitopes_to_keep=num_mutant_epitopes_to_keep,
            min_epitope_score=min_epitope_score)

        # do any of this variant's vaccine peptides contain mutant epitopes?
        any_mutant_epitopes = False
        for vaccine_peptide in vaccine_peptides:
            if vaccine_peptide.contains_mutant_epitopes():
                any_mutant_epitopes = True
                break
        if any_mutant_epitopes:
            counts_dict['num_variants_with_vaccine_peptides'] += 1
        result_dict[variant] = vaccine_peptides

    for key, value in counts_dict.items():
        logger.info('%s: %d', key, value)

    return result_dict, counts_dict

def ranked_vaccine_peptides(
        reads_generator,
        mhc_predictor,
        vaccine_peptide_length,
        padding_around_mutation,
        max_vaccine_peptides_per_variant,
        min_alt_rna_reads,
        min_variant_sequence_coverage,
        variant_sequence_assembly,
        num_mutant_epitopes_to_keep=10000,
        min_epitope_score=0):
    """
    Returns a tuple of two values:
    - sorted list whose first element is a Variant and whose second
    element is a list of VaccinePeptide objects
    - dictionary containing "funnel" variant counts for report display. Keys are e.g.
    "num_variants_with_rna_support", values are integer counts.
    """
    variants_to_vaccine_peptides_dict, variant_counts_dict = generate_vaccine_peptides(
        reads_generator=reads_generator,
        mhc_predictor=mhc_predictor,
        vaccine_peptide_length=vaccine_peptide_length,
        padding_around_mutation=padding_around_mutation,
        max_vaccine_peptides_per_variant=max_vaccine_peptides_per_variant,
        min_alt_rna_reads=min_alt_rna_reads,
        min_variant_sequence_coverage=min_variant_sequence_coverage,
        variant_sequence_assembly=variant_sequence_assembly,
        num_mutant_epitopes_to_keep=num_mutant_epitopes_to_keep,
        min_epitope_score=min_epitope_score)
    result_list = list(variants_to_vaccine_peptides_dict.items())
    # TODO: move this sort key into its own function, also make less nuts
    result_list.sort(
        key=lambda x: x[1][0].combined_score if len(x[1]) > 0 else 0.0,
        reverse=True)
    return result_list, variant_counts_dict
