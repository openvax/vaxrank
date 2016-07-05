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

class VaccinePeptideMetrics(namedtuple("VaccinePeptideMetrics", [
        "mutant_epitope_score",
        "wildtype_epitope_score",
        "n_mutant_residues",
        "mutation_distance_from_edge",
        "n_alt_reads",
        "n_alt_reads_supporting_protein_sequence"])):

    @classmethod
    def from_epitope_predictions(
            cls,
            epitope_predictions,
            mutant_protein_fragment):
        wildtype_epitope_score = sum(
            p.logistic_score()
            for p in epitope_predictions
            if not p.overlaps_mutation)
        mutant_epitope_score = sum(
            p.logistic_score()
            for p in epitope_predictions
            if p.overlaps_mutation)
        n_supporting = mutant_protein_fragment.n_alt_reads_supporting_protein_sequence
        return VaccinePeptideMetrics(
            mutant_epitope_score=mutant_epitope_score,
            wildtype_epitope_score=wildtype_epitope_score,
            n_mutant_residues=mutant_protein_fragment.n_mutant_amino_acids,
            mutation_distance_from_edge=mutant_protein_fragment.mutation_distance_from_edge,
            n_alt_reads=mutant_protein_fragment.n_alt_reads,
            n_alt_reads_supporting_protein_sequence=n_supporting)

    def lexicographic_sort_key(self):
        """
        Create tuple of scores so that candidates get sorted lexicographically
        by multiple criteria. Make sure to make the wildtype epitope
        score positive (since we want fewer wildtype epitopes) but the others
        negative (since we want more of them).
        """
        return (
            -self.mutant_epitope_score,
            -self.n_alt_reads,
            -self.n_alt_reads_supporting_protein_sequence,
            self.wildtype_epitope_score,
            -self.n_mutant_residues,
            -self.mutation_distance_from_edge
        )

    def expression_score(self):
        return np.sqrt(self.n_alt_reads)

    def combined_score(self):
        return self.expression_score() * self.mutant_epitope_score
