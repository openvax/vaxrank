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

from .manufacturability import ManufacturabilityScores

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
            if not p.overlaps_mutation or p.occurs_in_reference)
        mutant_epitope_score = sum(
            p.logistic_score()
            for p in epitope_predictions
            if p.overlaps_mutation and not p.occurs_in_reference)
        return VaccinePeptideBase.__new__(
            cls,
            mutant_protein_fragment=mutant_protein_fragment,
            epitope_predictions=epitope_predictions,
            mutant_epitope_score=mutant_epitope_score,
            wildtype_epitope_score=wildtype_epitope_score,
            manufacturability_scores=ManufacturabilityScores.from_amino_acids(
                mutant_protein_fragment.amino_acids))

    def peptide_synthesis_difficulty_score_tuple(
            self,
            hydropathy_threshold_high_priority=1.5,
            hydropathy_threshold_low_priority=0):
        """
        Generates a tuple of scores used for lexicographic sorting of vaccine
        peptides.

        The most important criterion for choosing a vaccine peptide is to
        get the mean hydropathy of the C-terminal residues below 1.5
        and also to ensure that no window of amino acids within
        the sequence has a mean hydropathy score > 1.5 (using
        AA values from Table 2 of Kyte & Doolittle 1982).

        If there are multiple vaccine peptides all of whose subsequence
        windows satisfy the GRAVY (mean hydropathy) < 1.5 constraint then
        let's optimize the terminal amino acids to exclude ones known to
        make solid phase synthesis difficult. Additionally we try to
        minimize the number of cysteine residues to prevent the formation
        of disulfide bonds.

        If there are multiple vaccine peptides without difficult terminal
        residues then we return to trying to minimize the GRAVY score, this
        time trying to push it below 0.

        If we are succsessful at getting the max GRAVY score of any
        subsequence below 0 then also try to eliminate N-terminal
        asparagine residues (not as harmful) and asparagine-proline bonds
        (known to dissociate easily).

        (Sort criteria determined through conversations with manufacturer)
        """
        cterm_7mer_gravy = self.manufacturability_scores.cterm_7mer_gravy_score
        max_7mer_gravy = self.manufacturability_scores.max_7mer_gravy_score

        # numbers we want to minimize, so a bigger number is worse

        return (
            # Priority I: C-terminal 7mer GRAVY score < 1.5
            max(0, cterm_7mer_gravy - hydropathy_threshold_high_priority),

            # Priority II: max 7mer GRAVY score < 1.5
            max(0, max_7mer_gravy - hydropathy_threshold_high_priority),

            # Priority III: avoid N-terminal Gln, Glu, Cys
            self.manufacturability_scores.difficult_n_terminal_residue,

            # Priority IV: avoid C-terminal Cys
            self.manufacturability_scores.c_terminal_cysteine,

            # Priority V: avoid C-terminal Pro
            self.manufacturability_scores.c_terminal_proline,

            # Priority VI: total number of Cys residues
            self.manufacturability_scores.cysteine_count,

            # Priority VII: minimize C-terminal 7mer GRAVY score (if > 0)
            max(0, cterm_7mer_gravy - hydropathy_threshold_low_priority),

            # Priority VIII: minimize max 7mer GRAVY score (if > 0)
            max(0, max_7mer_gravy - hydropathy_threshold_low_priority),

            # Priority IX: avoid N-terminal Asn
            self.manufacturability_scores.n_terminal_asparagine,

            # Priority X: avoid Asp-Pro bonds
            self.manufacturability_scores.asparagine_proline_bond_count,
        )

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
        manufacturability_score_tuple = self.peptide_synthesis_difficulty_score_tuple()
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

    def to_dict(self):
        return {
            "mutant_protein_fragment": self.mutant_protein_fragment,
            "epitope_predictions": self.epitope_predictions,
        }
