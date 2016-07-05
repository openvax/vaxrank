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
import heapq
from collections import OrderedDict

import pandas as pd
from isovar.protein_sequences import (
    reads_generator_to_protein_sequences_generator,
)

from .mutant_protein_fragment import MutantProteinFragment
from .epitope_prediction import predict_epitopes, slice_epitope_predictions
from .vaccine_peptide_metrics import VaccinePeptideMetrics

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
        logging.info("No protein sequences for %s" % (variant,))
        return []

    protein_fragment = MutantProteinFragment.from_isovar_protein_sequence(
        variant=variant,
        protein_sequence=isovar_protein_sequences[0])

    epitope_predictions = predict_epitopes(
        mhc_predictor=mhc_predictor,
        protein_fragment=protein_fragment)

    candidate_vaccine_peptides_to_metrics = {}
    for offset, candidate_vaccine_peptide in protein_fragment.top_k_subsequences(
            subsequence_length=vaccine_peptide_length,
            k=2 * padding_around_mutation + 1):
        subsequence_epitope_predictions = slice_epitope_predictions(
            epitope_predictions,
            start_offset=offset,
            end_offset=offset + len(candidate_vaccine_peptide))
        assert all(
            p.source_sequence == candidate_vaccine_peptide.amino_acids
            for p in subsequence_epitope_predictions)

        vaccine_peptide_metrics = VaccinePeptideMetrics.from_epitope_predictions(
            epitope_predictions=subsequence_epitope_predictions,
            mutant_protein_fragment=candidate_vaccine_peptide)

        logging.info("%s: %s, combined score: %0.4f" % (
            candidate_vaccine_peptide,
            vaccine_peptide_metrics,
            vaccine_peptide_metrics.combined_score()))

        candidate_vaccine_peptides_to_metrics[
            candidate_vaccine_peptide] = vaccine_peptide_metrics

    max_score = max(
        metrics.combined_score()
        for metrics in candidate_vaccine_peptides_to_metrics.values())
    n_total_candidates = len(candidate_vaccine_peptides_to_metrics)
    # only keep candidate vaccines that are within 1% of the maximum
    # combined score
    filtered_candidate_vaccine_peptides_with_metrics = [
        (vaccine_peptide, metrics)
        for (vaccine_peptide, metrics)
        in candidate_vaccine_peptides_to_metrics.items()
        if metrics.combined_score() / max_score > 0.99
    ]
    print("Keeping %d/%d vaccine peptides for %s" % (
        len(filtered_candidate_vaccine_peptides_with_metrics),
        n_total_candidates,
        variant))
    filtered_candidate_vaccine_peptides_with_metrics.sort(
        key=lambda x: x[1].lexicographic_sort_key())
    if len(filtered_candidate_vaccine_peptides_with_metrics) > 0:
        print("\n\nTop vaccine peptides for %s:" % variant)
        for i, (vaccine_peptide, metrics) in enumerate(
                filtered_candidate_vaccine_peptides_with_metrics):
            print("%d) %s: %s (combined score = %0.4f)\n" % (
                i + 1,
                vaccine_peptide,
                metrics,
                metrics.combined_score()))
    return filtered_candidate_vaccine_peptides_with_metrics[
        :max_vaccine_peptides_per_variant]

def generate_vaccine_peptides(
        reads_generator,
        mhc_predictor,
        vaccine_peptide_length,
        padding_around_mutation,
        max_vaccine_peptides_per_variant,
        min_reads_supporting_cdna_sequence):
    """
    Returns dictionary mapping each variant to list containing
    (MutantProteinFragment, VaccinePeptideMetrics) pairs
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
        vaccine_peptides_and_metrics = vaccine_peptides_for_variant(
            variant=variant,
            isovar_protein_sequences=isovar_protein_sequences,
            mhc_predictor=mhc_predictor,
            vaccine_peptide_length=vaccine_peptide_length,
            padding_around_mutation=padding_around_mutation,
            max_vaccine_peptides_per_variant=max_vaccine_peptides_per_variant)
        result_dict[variant] = vaccine_peptides_and_metrics
    return result_dict

def select_vaccine_peptides(
        reads_generator,
        mhc_predictor,
        vaccine_peptide_length,
        padding_around_mutation,
        max_vaccine_peptides_per_variant,
        min_reads_supporting_cdna_sequence,
        max_variants_selected):
    """
    Returns sorted list of at most `max_variants_selected` tuples whose
    first element is a Variant and whose second element is a list of
    (MutantProteinFragment, VaccinePeptideMetrics) pairs.
    """
    variants_to_vaccine_peptides_and_metrics_dict = generate_vaccine_peptides(
        reads_generator=reads_generator,
        mhc_predictor=mhc_predictor,
        vaccine_peptide_length=vaccine_peptide_length,
        padding_around_mutation=padding_around_mutation,
        max_vaccine_peptides_per_variant=max_vaccine_peptides_per_variant,
        min_reads_supporting_cdna_sequence=min_reads_supporting_cdna_sequence)

    variants_to_combined_scores_dict = {
        variant: vaccine_peptides_and_metrics[0][1].combined_score()
        for (variant, vaccine_peptides_and_metrics) in
        variants_to_vaccine_peptides_and_metrics_dict.items()
        if len(vaccine_peptides_and_metrics) > 0
    }
    return [
        (v, variants_to_vaccine_peptides_and_metrics_dict[v])
        for v in heapq.nlargest(
            max_variants_selected,
            variants_to_combined_scores_dict,
            key=variants_to_combined_scores_dict.get)
    ]

def select_vaccine_peptides_dataframe(**kwargs):
    """
    Takes the same arguments as `select_vaccine_peptides` but turns result
    type into a DataFrame.
    """
    columns = OrderedDict([
        ("genome", []),
        ("chr", []),
        ("pos", []),
        ("ref", []),
        ("alt", []),

    ])
    for field in MutantProteinFragment._fields + VaccinePeptideMetrics._fields:
        columns[field] = []
    columns["variant_rank"] = []
    columns["peptide_secondary_rank"] = []
    columns["expression_score"] = []
    columns["combined_score"] = []

    for i, (variant, vaccine_peptide_and_metrics_list) in enumerate(
            select_vaccine_peptides(**kwargs)):
        for j, (vaccine_peptide, metrics) in enumerate(
                vaccine_peptide_and_metrics_list):
            columns["genome"].append(variant.reference_name)
            columns["chr"].append(variant.contig)
            columns["pos"].append(variant.original_start)
            columns["ref"].append(variant.original_ref)
            columns["alt"].append(variant.original_alt)

            for field in MutantProteinFragment._fields:
                columns[field].append(getattr(vaccine_peptide, field))
            for field in VaccinePeptideMetrics._fields:
                if field in MutantProteinFragment._fields:
                    continue
                columns[field].append(getattr(metrics, field))
            columns["variant_rank"].append(i + 1)
            columns["peptide_secondary_rank"].append(j + 1)
            columns["expression_score"].append(metrics.expression_score())
            columns["combined_score"].append(metrics.combined_score())
    return pd.DataFrame(columns, columns=columns.keys())
