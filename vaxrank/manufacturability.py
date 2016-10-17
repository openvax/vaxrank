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

"""
Scoring functions for determing which sequences are easy to manufacture using
solid-phase synthesis.

For more information see: https://github.com/hammerlab/vaxrank/issues/2
"""

from __future__ import absolute_import, print_function, division

# Amino Acid Hydropathy Score
# Table 2 from Kyte and Doolittle"s
# "A Simple Method for Displaying the Hydropathic Character of a Protein"

hydropathy_dict = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3
}


def gravy_score(amino_acids):
    """
    Mean amino acid hydropathy averaged across residues of a peptide
    or protein sequence.
    """
    total = 0.0
    count = 0
    for amino_acid in amino_acids:
        if amino_acid == "X":
            # skip unknown amino acids
            continue
        elif amino_acid == "*":
            # just in case stop codons sneak into the amino acid sequence
            break
        total += hydropathy_dict[amino_acid]
        count += 1
    if count == 0:
        raise ValueError("Can't compute GRAVY score for '%s'" % amino_acids)
    return total / count

def max_kmer_gravy_score(amino_acids, k):
    """
    Returns max GRAVY score of any kmer in the amino acid sequence,
    used to determine if there are any extremely hydrophobic regions within a
    longer amino acid sequence.
    """
    return max(
        gravy_score(amino_acids[i:i + k])
        for i in range(len(amino_acids) - k))

def difficult_n_terminal_residue(amino_acids):
    """
    Is the N-terminus one of {Gln, Glu, Asn, Cys}
    """
    nterm_amino_acid = amino_acids[0]
    return nterm_amino_acid in {"Q", "E", "N", "C"}

def difficult_c_terminal_residue(amino_acids):
    """
    Is the C-terminus either Pro and Cys?
    """
    cterm_amino_acid = amino_acids[-1]
    return cterm_amino_acid in {"P", "C"}

def count_asparagine_proline_bonds(amino_acids):
    """
    Count the number of Asparagine/Asn/N-Proline/Pro/P bonds
    Problem with Asn-Pro bonds: can spontaneously cleave the peptide
    """
    return sum(
        amino_acids[i:i + 2] == "NP"
        for i in range(len(amino_acids) - 1))

def count_cysteine_residues(amino_acids):
    """
    How many cysteines are in the amino acid sequence?
    Problem with cysteine residues: They can form disulfide bonds across
    distant parts of the peptide,
    """
    return sum(amino_acid == "C" for amino_acid in amino_acids)
