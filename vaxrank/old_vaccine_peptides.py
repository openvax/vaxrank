from __future__ import absolute_import, print_function, division

import pandas as pd
import numpy as np
from isovar.cli_helpers import (
    variants_to_protein_sequences_dataframe_from_args,
)

from epitope_scoring import simple_ic50_epitope_scorer
from peptide_binding_measure import IC50_FIELD_NAME, PERCENTILE_RANK_FIELD_NAME
from immunogenicity import THYMIC_DELETION_FIELD_NAME


# To clarify the nomenclature (source seq vs. peptide vs. epitope)
# look at this example.
#
# Let's say we have a small protein with this 48 residue sequence:
#
#   TAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQR
#
# and a somatic variant results in the a single amino acid change
# of the first 'A' for a 'Q', then our mutated protein will be:
#
#   T[Q]ADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQR
#
# now either this full protein or some subset of it will become
# the 'SourceSequence' depending on how large of a vaccine peptide
# we're aiming to generate. If, for example, we only want 15-mer
# vaccine peptides, then the SourceSequence will be the union
# of all 15-mer windows containing the modified residue. In this
# case, since the change was near the start of the protein, there
# are only two such windows, yielding a SourceSequence of length 16:
#
#  T[Q]ADMAAQTTKHKWEA
#
# From this source sequence, we generate shorter substrings we
# expect to bind to MHC molecules (typically 8-11 residues).
# Assuming for simplicity that epitopes will all be 9-mers,
# the set of epitopes will be:
#
#  T[Q]ADMAAQT
#  [Q]ADMAAQTT
#  ADMAAQTTK
#  DMAAQTTKH
#  MAAQTTKHK
#  AAQTTKHKW
#  AQTTKHKWE
#
# We look at the scores for each of these epitopes to generate
# scores for the 15-mer peptides which contains these epitopes.
# So, the first peptide "TQADMAAQTTKHKWE" contains all but the
# the last epitope, and the second peptide "QADMAAQTTKHKWEA"
# contains all but the first epitope.
#

def is_mutant_epitope(epitope, mutation_start, mutation_end):
    """
    An epitope is considered mutant if it overlaps the mutated region
    of the source sequence and isn't similar to the thymically
    presented self epitopes.

    Parameters
    ----------

    epitope : dict

    mutation_start : int

    mutation_end : int
    """
    start = epitope['EpitopeStart']
    end = epitope['EpitopeEnd']
    overlaps = (start < mutation_end) and (end > mutation_start)
    if THYMIC_DELETION_FIELD_NAME in epitope:
        return not epitope[THYMIC_DELETION_FIELD_NAME] and overlaps
    else:
        return overlaps

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

def select_vaccine_peptide(
        seq,
        epitopes,
        mutation_start,
        mutation_end,
        epitope_scorer,
        result_length,
        padding):
    """
    Choose the best vaccine peptide from a longer amino acid sequence.
    The score of each candidate vaccine peptide is determined by
        - sum of mutant epitope scores
        - negative sum of wildtype epitope scores
        - total number of mutant residues in the vaccine peptide
        - mutation start placed closest to the center

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

    candidate_peptides = generate_candidate_vaccine_peptides(
        seq,
        epitopes,
        mutation_start,
        mutation_end,
        epitope_scorer,
        result_length,
        padding)

    def score_tuple(record):
        """
        Create tuple of scores so that candidates get sorted lexicographically
        by multiple criteria. Make sure to make the wildtype epitope
        score negative (since we want fewer wildtype epitopes)
        """
        return (
            record['MutantEpitopeScore'],
            -record['WildtypeEpitopeScore'],
            record['NumMutantResidues'],
            record['MutationDistanceFromEdge'],
        )
    candidate_peptides.sort(key=score_tuple, reverse=True)
    best = candidate_peptides[0]
    return best

def select_vaccine_peptides(
        source_peptides,
        epitope_scorer=simple_ic50_epitope_scorer,
        vaccine_peptide_length=31,
        padding=10):
    """
    Given a set of longer peptides and their associated predicted epitopes,
    find the best set of vaccine peptides overlapping each mutation.

    Parameters
    ----------

    source_peptides : list
        List of peptide source sequence records

    epitope_scorer : EpitopeScorer instance

    vaccine_peptide_length : int

    padding : int
        Maximum distance from edges of vaccine peptide where mutation can start.
    """

    results = []
    for peptide_record in source_peptides:
        seq = peptide_record['SourceSequence']
        assert len(seq) > 0, "Invalid empty peptide"
        mutation_start = peptide_record['MutationStart']
        mutation_end = peptide_record["MutationEnd"]
        epitopes = peptide_record['Epitopes']

        vaccine_peptide_record = select_vaccine_peptide(
            seq,
            epitopes,
            mutation_start,
            mutation_end,
            epitope_scorer,
            result_length=vaccine_peptide_length,
            padding=padding,
        )
        start_idx = vaccine_peptide_record['VaccinePeptideStart']

        assert start_idx >= 0
        # must overlap the mutation to some degree
        assert start_idx <= mutation_end

        # augment the vaccine peptide with all info that was attached to its
        # source sequence
        for k, v in peptide_record.iteritems():
            if k not in vaccine_peptide_record:
                vaccine_peptide_record[k] = v

        results.append(vaccine_peptide_record)

    # Make sure that sort is in descending order by vaccine peptide score.
    # When comparing across genes/mutations we only care about the mutant
    # epitope score.
    results.sort(key=lambda record: record['MutantEpitopeScore'], reverse=True)
    return results
