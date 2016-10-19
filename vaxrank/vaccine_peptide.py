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

from collections import namedtuple

import numpy as np

from .manufacturability import compute_manufacturability_scores

VaccinePeptideBase = namedtuple(
    "VaccinePeptide", [
        "mutant_protein_fragment",
        "epitope_predictions",
        "mutant_epitope_score",
        "wildtype_epitope_score",
        "manufacturability_scores"])

class VaccinePeptide(VaccinePeptideBase):
    """
    VaccinePeptide combines the sequence information of MutantProteinFragment
    with MHC binding predictions for subsequences of the protein fragment.
    """
    def __new__(
            cls,
            mutant_protein_fragment,
            epitope_predictions,
            min_epitope_score=0):
        wildtype_epitope_score = sum(
            p.logistic_score()
            for p in epitope_predictions
            if not p.overlaps_mutation)

        mutant_epitope_score = sum(
            p.logistic_score()
            for p in epitope_predictions
            if p.overlaps_mutation)

        return VaccinePeptideBase.__new__(
            cls,
            mutant_protein_fragment=mutant_protein_fragment,
            epitope_predictions=epitope_predictions,
            mutant_epitope_score=mutant_epitope_score,
            wildtype_epitope_score=wildtype_epitope_score,
            manufacturability_scores=compute_manufacturability_scores(
                mutant_protein_fragment.amino_acids))

    def lexicographic_sort_key(self):
        """
        Create tuple of scores so that candidates get sorted lexicographically
        by multiple criteria. Make sure to make the wildtype epitope
        score positive (since we want fewer wildtype epitopes) but the others
        negative (since we want more of them).
        """
        # since we're sorting in decreasing order, numbers which we want
        # to be larger must have their signs flipped
        essential_score_tuple = (
            # Sum of normalized MHC binding affinities of subsequences
            -self.mutant_epitope_score,
            # Number of reads supporting the variant
            -self.mutant_protein_fragment.n_alt_reads
        )
        manufacturability_score_tuple = tuple(self.manufacturability_scores)
        extra_score_tuple = (
            # Number of reads supporting the particular protein sequence
            # sequence we're using for this vaccine peptide. Currently
            # all vaccine peptides are drawn from the same larger sequence
            # so this score shouldn't change.
            -self.mutant_protein_fragment.n_alt_reads_supporting_protein_sequence,
            # Minimize the sum of non-mutant MHC binding scores
            self.wildtype_epitope_score,
            # All else being equal, we prefer to maximize the number of
            # mutant amino acids
            -self.mutant_protein_fragment.n_mutant_amino_acids,
            # If nothing else can serve as a tie break then try to center
            # the mutation in the vaccine peptide.
            -self.mutant_protein_fragment.mutation_distance_from_edge
        )
        return (
            essential_score_tuple +
            manufacturability_score_tuple +
            extra_score_tuple
        )

    @property
    def expression_score(self):
        return np.sqrt(self.mutant_protein_fragment.n_alt_reads)

    @property
    def combined_score(self):
        return self.expression_score * self.mutant_epitope_score
