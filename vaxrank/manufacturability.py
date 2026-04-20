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
Scoring functions for determining which sequences are easy to manufacture using
solid-phase synthesis.

For more information see: https://github.com/hammerlab/vaxrank/issues/2
"""

from collections import namedtuple

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
    if not amino_acids:
        return 0.0
    total = sum(
        hydropathy_dict[amino_acid] for amino_acid in amino_acids)
    return total / len(amino_acids)


def max_kmer_gravy_score(amino_acids, k):
    """
    Returns max GRAVY score of any kmer in the amino acid sequence,
    used to determine if there are any extremely hydrophobic regions within a
    longer amino acid sequence.
    """
    if len(amino_acids) < k:
        return gravy_score(amino_acids)
    return max(
        gravy_score(amino_acids[i:i + k])
        for i in range(len(amino_acids) - k + 1))


def max_7mer_gravy_score(amino_acids):
    return max_kmer_gravy_score(amino_acids, 7)


def cterm_kmer_gravy_score(amino_acids, k):
    """
    Mean hydropathy of last k residues on the C-terminus of the peptide.
    """
    n = len(amino_acids)
    return gravy_score(amino_acids[n - k:n])


def cterm_7mer_gravy_score(amino_acids):
    return cterm_kmer_gravy_score(amino_acids, 7)


def difficult_n_terminal_residue(amino_acids):
    """
    Is the N-terminus one of {Gln, Glu, Cys}?
    ---
    Priority I: avoid N-terminal Gln, Glu, Cys
    """
    return amino_acids[0] in {"Q", "E", "C"}


def c_terminal_proline(amino_acids):
    """
    Is the right-most (C-terminal) amino acid a proline?
    """
    return amino_acids[-1] == "P"


def c_terminal_cysteine(amino_acids):
    """
    Is the right-most (C-terminal) amino acid a cysteine?
    """
    return amino_acids[-1] == "C"


def n_terminal_asparagine(amino_acids):
    """
    Asparagine at the N-terminus of a peptide is also hard
    to synthesize, though not as bad as {Gln, Glu, Cys}
    """
    return amino_acids[0] == "N"


def aspartate_proline_bond_count(amino_acids):
    """
    Count the number of Aspartate/Asp/D-Proline/Pro/P bonds.
    Asp-Pro bonds are particularly susceptible to non-enzymatic hydrolysis
    under acidic conditions, which can spontaneously cleave the peptide.
    """
    return sum(
        amino_acids[i:i + 2] == "DP"
        for i in range(len(amino_acids) - 1))


def cysteine_count(amino_acids):
    """
    How many cysteines are in the amino acid sequence?
    Problem with cysteine residues: They can form disulfide bonds across
    distant parts of the peptide
    """
    return sum(amino_acid == "C" for amino_acid in amino_acids)


def combine_scoring_functions(*scoring_functions):
    """
    Given a list of scoring functions, make a namedtuple with
    fields of the same names. Returns the ManufacturabilityScores class.
    """
    names = [fn.__name__ for fn in scoring_functions]

    class ManufacturabilityScores(namedtuple('ManufacturabilityScores', names)):
        @classmethod
        def from_amino_acids(cls, amino_acids):
            return cls(*[fn(amino_acids) for fn in scoring_functions])

    return ManufacturabilityScores

ManufacturabilityScores = combine_scoring_functions(

    # GRAVY score of 7 residues closest to the C terminus
    cterm_7mer_gravy_score,

    # GRAVY score of any 7mer window in the peptide sequence
    max_7mer_gravy_score,

    # avoid N-terminal Gln, Glu, Cys
    difficult_n_terminal_residue,

    # avoid C-terminal Cys
    c_terminal_cysteine,

    # avoid C-terminal Pro
    c_terminal_proline,

    # total number of Cys residues
    cysteine_count,

    # avoid N-terminal Asn
    n_terminal_asparagine,

    # avoid Asp-Pro bonds
    aspartate_proline_bond_count,
)


# Registry of manufacturability rules used to build the lexicographic
# sort tuple in VaccinePeptide. Each rule takes (scores, thresholds) and
# returns one numeric value to minimize.
#
# Thresholds is a dict with keys:
#   max_c_terminal_hydropathy, min_kmer_hydropathy,
#   max_kmer_hydropathy_low_priority, max_kmer_hydropathy_high_priority


def _cysteine_count_rule(scores, _thresholds):
    return scores.cysteine_count


def _cterm_hydropathy_rule(scores, thresholds):
    return max(
        0,
        scores.cterm_7mer_gravy_score - thresholds["max_c_terminal_hydropathy"],
    )


def _max_hydropathy_high_rule(scores, thresholds):
    return max(
        0,
        scores.max_7mer_gravy_score - thresholds["max_kmer_hydropathy_high_priority"],
    )


def _difficult_n_terminal_rule(scores, _thresholds):
    return scores.difficult_n_terminal_residue


def _c_terminal_cysteine_rule(scores, _thresholds):
    return scores.c_terminal_cysteine


def _c_terminal_proline_rule(scores, _thresholds):
    return scores.c_terminal_proline


def _n_terminal_asparagine_rule(scores, _thresholds):
    return scores.n_terminal_asparagine


def _aspartate_proline_rule(scores, _thresholds):
    return scores.aspartate_proline_bond_count


def _max_hydropathy_low_rule(scores, thresholds):
    return max(
        0,
        scores.max_7mer_gravy_score - thresholds["max_kmer_hydropathy_low_priority"],
    )


def _min_hydropathy_rule(scores, thresholds):
    return max(
        0,
        thresholds["min_kmer_hydropathy"] - scores.max_7mer_gravy_score,
    )


MANUFACTURABILITY_RULE_REGISTRY = {
    "cysteine_count": _cysteine_count_rule,
    "cterm_hydropathy": _cterm_hydropathy_rule,
    "max_hydropathy_high": _max_hydropathy_high_rule,
    "difficult_n_terminal": _difficult_n_terminal_rule,
    "c_terminal_cysteine": _c_terminal_cysteine_rule,
    "c_terminal_proline": _c_terminal_proline_rule,
    "n_terminal_asparagine": _n_terminal_asparagine_rule,
    "aspartate_proline": _aspartate_proline_rule,
    "max_hydropathy_low": _max_hydropathy_low_rule,
    "min_hydropathy": _min_hydropathy_rule,
}

# Legacy ordering preserved byte-for-byte: this is the exact tuple order
# `peptide_synthesis_difficulty_score_tuple` produced before rules were
# made configurable. Do not reorder without updating parity tests.
DEFAULT_MANUFACTURABILITY_RULES = (
    "cysteine_count",
    "cterm_hydropathy",
    "max_hydropathy_high",
    "difficult_n_terminal",
    "c_terminal_cysteine",
    "c_terminal_proline",
    "n_terminal_asparagine",
    "aspartate_proline",
    "max_hydropathy_low",
    "min_hydropathy",
)


def compute_manufacturability_tuple(scores, thresholds, rules=None):
    """Apply the ordered rule list against a ManufacturabilityScores instance.

    ``rules`` defaults to ``DEFAULT_MANUFACTURABILITY_RULES`` (the legacy
    order). Each entry must be a key in ``MANUFACTURABILITY_RULE_REGISTRY``.
    """
    if rules is None:
        rules = DEFAULT_MANUFACTURABILITY_RULES
    tuple_values = []
    for rule_name in rules:
        try:
            fn = MANUFACTURABILITY_RULE_REGISTRY[rule_name]
        except KeyError:
            raise ValueError(
                f"Unknown manufacturability rule '{rule_name}'. "
                f"Available: {sorted(MANUFACTURABILITY_RULE_REGISTRY)}"
            ) from None
        tuple_values.append(fn(scores, thresholds))
    return tuple(tuple_values)
