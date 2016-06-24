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


def clamp(i, last_pos):
    """
    Given a number between [-inf, inf], clamp it between [0,last_pos]

    Parameters
    ----------
    i : int

    last_pos : int
    """
    return min(max(i, 0), last_pos)

def generate_candidate_vaccine_peptides(
        seq,
        epitopes,
        mutation_start,
        mutation_end,
        epitope_scorer,
        result_length,
        padding):
    """
    Generate all possible vaccine peptides tiling over an amino acid sequence.

    Parameters
    ----------

    seq : str

    epitopes : list
        List of epitope records, each containing a nested list of per-allele
        binding predictions

    mutation_start : int
        Where in the given sequence is the first mutated residue?

    mutation_end : int
        Where in the given sequence is the last mutated residue?

    epitope_scorer : EpitopeScorer

    result_length : int
        How big of a substring are we looking to pull out as a vaccine peptide?

    padding : int
    """
    n = len(seq)

    if n <= result_length:
        # if source sequence is too short, just return whatever we have
        first_pos = 0
        n_candidates = 1
    elif n <= result_length + 2 * padding:
        # if the source sequence is too short for the full amount of requested
        # padding, then center as best as we can
        actual_combined_padding = n - result_length
        first_pos = actual_combined_padding / 2
        # if there are two equally good ways to center the insufficiently
        # padded sequence, try them both
        n_candidates = 1 if actual_combined_padding % 2 == 1 else 2
    else:
        first_pos = padding
        n_candidates = n - result_length - 2 * padding + 1

    # in case the mutation is at the beginning or end of the peptide,
    # make sure we cover it
    if mutation_start < first_pos:
        difference = first_pos - mutation_start
        first_pos = mutation_start
        n_candidates += difference

    if mutation_start > first_pos + result_length + n_candidates:
        difference = mutation_start - (first_pos + result_length + n_candidates)
        n_candidates += difference

    # we're going to lexically sort each peptide by four criteria:
    #   - average score of its mutated epitopes
    #   - negative average score of it's wildtype epitopes
    #   - number of mutated residues covered
    #   - distance from the edge of the spurce sequence
    candidate_peptides = []

    for peptide_start in xrange(first_pos, first_pos + n_candidates):

        peptide_seq = seq[peptide_start : peptide_start+result_length]

        # use this instead of 'result_length' just in case
        # we're dealing with a source sequence shorter than
        # the desired full length
        peptide_length = len(peptide_seq)
        peptide_end = peptide_start + peptide_length

        peptide_mutation_start = clamp(
            mutation_start - peptide_start, peptide_length)
        peptide_mutation_end = clamp(
            mutation_end - peptide_start, peptide_length)
        number_mutant_residues = peptide_mutation_end - peptide_mutation_start

        mutation_distance_from_edge = min(
            peptide_mutation_start,
            peptide_length - peptide_mutation_start)

        mutant_score = 0.0
        wildtype_score = 0.0
        for epitope in epitopes:
            epitope_start = epitope['EpitopeStart']
            epitope_end = epitope['EpitopeEnd']
            overlaps = (
                epitope_start >= peptide_start
                and epitope_end <= peptide_end
            )
            mutant = is_mutant_epitope(epitope, mutation_start, mutation_end)

            if overlaps:
                score = epitope_scorer.epitope_score(epitope)
                if is_mutant_epitope(epitope, mutation_start, mutation_end):
                    mutant_score += score
                else:
                    wildtype_score += score

        vaccine_peptide_record = {
            'VaccinePeptide' : peptide_seq,
            'VaccinePeptideMutationStart' : peptide_mutation_start,
            'VaccinePeptideMutationEnd' : peptide_mutation_end,
            'MutantEpitopeScore' : mutant_score,
            'WildtypeEpitopeScore' : wildtype_score,
            'VaccinePeptideStart' : peptide_start,
            'VaccinePeptideEnd' : peptide_end,
            'VaccinePeptideLength' : peptide_end - peptide_start,
            'NumMutantResidues' : number_mutant_residues,
            'MutationDistanceFromEdge' : mutation_distance_from_edge,
        }
        candidate_peptides.append(vaccine_peptide_record)

    return candidate_peptides
