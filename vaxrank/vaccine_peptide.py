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

from collections import namedtuple

VaccinePeptide = namedtuple("VaccinePeptide", (
    "variant",  # varcode.Variant object
    "isovar_protein_sequence"

    ###
    # Translated protein sequence, aggregated from possibly multiple
    # synonymous coding sequences
    ###

    "amino_acids",
    # offsets of amino acids which differ due to the mutation
    "mutant_amino_acid_start_offset",
    "mutant_amino_acid_end_offset",
    "n_mutant_amino_acids",
    # offsets of codons containing mutant nucleotides, even if
    # synonymous with original reference sequence
    "mutant_codon_start_offset",
    "mutant_codon_end_offset",
    "n_mutant_codons",

    ###
    # RNA evidence
    ###

    # number of RNA reads fully spanning the cDNA sequence(s) from which we
    # translated this amino acid sequence.
    "n_rna_reads",

    ###
    # Properties which affect ranking or filtering of vaccine peptide
    ###

    # how many amino acids from the center is the first mutant amino acid
    "distance_from_center",
    # list of epitope predictions
    "epitope_predictions",
    # sum of normalized IC50s
    "combined_epitope_score",
))
