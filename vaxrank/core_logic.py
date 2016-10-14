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
import logging
from collections import OrderedDict

import pandas as pd
from isovar.protein_sequences import (
    reads_generator_to_protein_sequences_generator,
)

from .mutant_protein_fragment import MutantProteinFragment
from .epitope_prediction import predict_epitopes, slice_epitope_predictions
from .vaccine_peptide import VaccinePeptide


logger = logging.getLogger(__name__)


def vaccine_peptides_for_variant(
        variant,
        isovar_protein_sequences,
        mhc_predictor,
        vaccine_peptide_length,
        padding_around_mutation,
        max_vaccine_peptides_per_variant):
    """
    Returns list containing (MutantProteinFragment, VaccinePeptideMetrics) pairs
    """
    isovar_protein_sequences = list(isovar_protein_sequences)
    if len(isovar_protein_sequences) == 0:
        logger.info("No protein sequences for %s", variant)
        return []

    protein_fragment = MutantProteinFragment.from_isovar_protein_sequence(
        variant=variant,
        protein_sequence=isovar_protein_sequences[0])

    logger.info("Mutant protein fragment for %s: %s",
        variant,
        protein_fragment)

    epitope_predictions = predict_epitopes(
        mhc_predictor=mhc_predictor,
        protein_fragment=protein_fragment)

    candidate_vaccine_peptides = []

    for offset, candidate_fragment in protein_fragment.top_k_subsequences(
            subsequence_length=vaccine_peptide_length,
            k=2 * padding_around_mutation + 1):
        subsequence_epitope_predictions = slice_epitope_predictions(
            epitope_predictions,
            start_offset=offset,
            end_offset=offset + len(candidate_fragment))

        assert all(
            p.source_sequence == candidate_fragment.amino_acids
            for p in subsequence_epitope_predictions)

        candidate_vaccine_peptide = VaccinePeptide(
            mutant_protein_fragment=candidate_fragment,
            epitope_predictions=subsequence_epitope_predictions)

        logger.debug("%s, combined score: %0.4f",
            candidate_vaccine_peptide,
            candidate_vaccine_peptide.combined_score)
        candidate_vaccine_peptides.append(candidate_vaccine_peptide)

    max_score = max(vp.combined_score for vp in candidate_vaccine_peptides)
    n_total_candidates = len(candidate_vaccine_peptides)

    # only keep candidate vaccines that are within 1% of the maximum
    # combined score
    filtered_candidate_vaccine_peptides = [
        vp
        for vp in candidate_vaccine_peptides
        if vp.combined_score / max_score > 0.99
    ]

    logger.info("Keeping %d/%d vaccine peptides for %s",
        len(filtered_candidate_vaccine_peptides),
        n_total_candidates,
        variant)
    filtered_candidate_vaccine_peptides.sort(key=VaccinePeptide.lexicographic_sort_key)

    if len(filtered_candidate_vaccine_peptides) > 0:
        logger.debug("Top vaccine peptides for %s:", variant)
        for i, vaccine_peptide in enumerate(filtered_candidate_vaccine_peptides):
            logger.debug("%d) %s (combined score = %0.4f)",
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
        min_reads_supporting_cdna_sequence):
    """
    Returns dictionary mapping each variant to list of VaccinePeptide objects.
    """

    # total number of amino acids is the vaccine peptide length plus the
    # number of off-center windows around the mutation
    protein_fragment_sequence_length = (
        vaccine_peptide_length + 2 * padding_around_mutation)

    protein_sequences_generator = reads_generator_to_protein_sequences_generator(
        reads_generator,
        transcript_id_whitelist=None,
        protein_sequence_length=protein_fragment_sequence_length,
        min_reads_supporting_cdna_sequence=min_reads_supporting_cdna_sequence,
        max_protein_sequences_per_variant=1)

    result_dict = {}
    for variant, isovar_protein_sequences in protein_sequences_generator:
        vaccine_peptides = vaccine_peptides_for_variant(
            variant=variant,
            isovar_protein_sequences=isovar_protein_sequences,
            mhc_predictor=mhc_predictor,
            vaccine_peptide_length=vaccine_peptide_length,
            padding_around_mutation=padding_around_mutation,
            max_vaccine_peptides_per_variant=max_vaccine_peptides_per_variant)
        result_dict[variant] = vaccine_peptides
    return result_dict

def ranked_vaccine_peptides(
        reads_generator,
        mhc_predictor,
        vaccine_peptide_length,
        padding_around_mutation,
        max_vaccine_peptides_per_variant,
        min_reads_supporting_cdna_sequence):
    """
    Returns sorted list whose first element is a Variant and whose second
    element is a list of VaccinePeptide objects.
    """
    variants_to_vaccine_peptides_dict = generate_vaccine_peptides(
        reads_generator=reads_generator,
        mhc_predictor=mhc_predictor,
        vaccine_peptide_length=vaccine_peptide_length,
        padding_around_mutation=padding_around_mutation,
        max_vaccine_peptides_per_variant=max_vaccine_peptides_per_variant,
        min_reads_supporting_cdna_sequence=min_reads_supporting_cdna_sequence)
    result_list = list(variants_to_vaccine_peptides_dict.items())
    result_list.sort(
        key=lambda x: x[1][0].combined_score if len(x[1]) > 0 else 0.0,
        reverse=True)
    return result_list


def dataframe_from_ranked_list(ranked_list_of_variants_with_vaccine_peptides):
    """
    Converts list of (Variant, [VaccinePeptide]) into a DataFrame, dropping
    any variants without vaccine peptide predictions.
    """
    columns = OrderedDict([
        ("genome", []),
        ("chr", []),
        ("pos", []),
        ("ref", []),
        ("alt", []),

    ])
    for field in MutantProteinFragment._fields:
        columns[field] = []
    columns["variant_rank"] = []
    columns["peptide_secondary_rank"] = []
    columns["mutant_epitope_score"] = []
    columns["wildtype_epitope_score"] = []
    columns["expression_score"] = []
    columns["combined_score"] = []

    for i, (variant, vaccine_peptides) in enumerate(
            ranked_list_of_variants_with_vaccine_peptides):
        for j, vaccine_peptide in enumerate(vaccine_peptides):
            columns["genome"].append(variant.reference_name)
            columns["chr"].append(variant.contig)
            columns["pos"].append(variant.original_start)
            columns["ref"].append(variant.original_ref)
            columns["alt"].append(variant.original_alt)

            for field in MutantProteinFragment._fields:
                columns[field].append(getattr(
                    vaccine_peptide.mutant_protein_fragment, field))
            columns["variant_rank"].append(i + 1)
            columns["peptide_secondary_rank"].append(j + 1)
            columns["mutant_epitope_score"].append(vaccine_peptide.mutant_epitope_score)
            columns["wildtype_epitope_score"].append(vaccine_peptide.wildtype_epitope_score)
            columns["expression_score"].append(vaccine_peptide.expression_score)
            columns["combined_score"].append(vaccine_peptide.combined_score)
    return pd.DataFrame(columns, columns=columns.keys())
