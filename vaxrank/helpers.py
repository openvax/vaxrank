
def create_vaccine_fasta_lines(vaccine_peptide_records):
    fasta_lines = []
    for i, record in enumerate(vaccine_peptide_records):
        line = ">%d Gene=%s, Transcript=%s, Mut=%s (%d:%d), Score=%0.6f" % (
            i,
            record['Gene'],
            record['TranscriptId'],
            record['PeptideMutationInfo'],
            record['VaccinePeptideMutationStart'],
            record['VaccinePeptideMutationEnd'],
            record['MutantEpitopeScore']
        )
        fasta_lines.append(line)
        fasta_lines.append(record['VaccinePeptide'])
    return fasta_lines

def print_epitopes(source_sequences):
    print("\nEpitopes")
    print("--------\n")
    for record in sorted(source_sequences, key=lambda r: r['TranscriptId']):
        mut_start = record['MutationStart']
        mut_end = record['MutationEnd']
        print(">Gene=%s, Transcript=%s, Mut=%s (%d:%d)" % (
            record['Gene'],
            record['TranscriptId'],
            record['PeptideMutationInfo'],
            mut_start,
            mut_end,))
        print(record['SourceSequence'])
        for epitope in sorted(
                record['Epitopes'], key=lambda e: e['EpitopeStart']):
            overlap_start = epitope['EpitopeStart'] < mut_end
            overlap_end = epitope['EpitopeEnd'] > mut_start
            mutant = overlap_start and overlap_end
            print("\t", epitope['Epitope'], ("<-- MUTANT" if mutant else ""))
            for prediction in sorted(
                    epitope["MHC_Allele_Scores"],
                    key=lambda p: p['Allele']):
                print("\t\t", "allele = %s, IC50=%0.4f" % (
                    prediction['Allele'],
                    prediction[IC50_FIELD_NAME]))
