
"""
Each candidate vaccine peptide is a subsequence of a MutantProteinFragment
combined with several scores such as expression and MHC binding.
"""

VaccinePeptide = namedtuple("VaccinePeptide", (
    "source_mutant_protein_fragment",
    "offset",
    "amino_acids",
    "mutation_start_offset"
            'VaccinePeptideMutationStart': peptide_mutation_start,
            'VaccinePeptideMutationEnd': peptide_mutation_end,
            'MutantEpitopeScore': mutant_score,
            'WildtypeEpitopeScore': wildtype_score,
            'VaccinePeptideStart': peptide_start,
            'VaccinePeptideEnd': peptide_end,
            'VaccinePeptideLength': peptide_end - peptide_start,
            'NumMutantResidues': number_mutant_residues,
            'MutationDistanceFromEdge': mutation_distance_from_edge,
))
    for peptide_start in xrange(first_pos, first_pos + n_candidates):

        peptide_seq = seq[peptide_start:peptide_start + result_length]

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
                epitope_start >= peptide_start and epitope_end <= peptide_end
            )
            # mutant = is_mutant_epitope(epitope, mutation_start, mutation_end)

            if overlaps:
                score = epitope_scorer.epitope_score(epitope)
                if is_mutant_epitope(epitope, mutation_start, mutation_end):
                    mutant_score += score
                else:
                    wildtype_score += score

        vaccine_peptide_record = {

        }
        candidate_peptides.append(vaccine_peptide_record)

    return candidate_peptides
