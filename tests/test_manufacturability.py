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

import pytest

from vaxrank.manufacturability import (
    DEFAULT_MANUFACTURABILITY_RULES,
    MANUFACTURABILITY_RULE_REGISTRY,
    ManufacturabilityScores,
    compute_manufacturability_tuple,
)


DEFAULT_THRESHOLDS = {
    "max_c_terminal_hydropathy": 1.5,
    "min_kmer_hydropathy": 0.0,
    "max_kmer_hydropathy_low_priority": 1.5,
    "max_kmer_hydropathy_high_priority": 2.5,
}

def test_c_terminal_proline():
    scores = ManufacturabilityScores.from_amino_acids("A" * 6 + "P")
    assert scores.c_terminal_proline

    scores = ManufacturabilityScores.from_amino_acids("A" * 7)
    assert not scores.c_terminal_proline


def test_n_terminal_cysteine():
    scores = ManufacturabilityScores.from_amino_acids("C" + 6 * "A")
    assert scores.difficult_n_terminal_residue

    scores = ManufacturabilityScores.from_amino_acids(7 * "A")
    assert not scores.difficult_n_terminal_residue


def test_n_terminal_glutamic_acid():
    scores = ManufacturabilityScores.from_amino_acids("E" + 6 * "A")
    assert scores.difficult_n_terminal_residue

    scores = ManufacturabilityScores.from_amino_acids(7 * "A")
    assert not scores.difficult_n_terminal_residue


def test_n_terminal_glutamine():
    scores = ManufacturabilityScores.from_amino_acids("Q" + 6 * "A")
    assert scores.difficult_n_terminal_residue

    scores = ManufacturabilityScores.from_amino_acids(7 * "A")
    assert not scores.difficult_n_terminal_residue


def test_asp_pro_bond_count():
    scores = ManufacturabilityScores.from_amino_acids("A" * 7)
    assert scores.aspartate_proline_bond_count == 0

    scores = ManufacturabilityScores.from_amino_acids("DP" + "A" * 7 + "DP")
    assert scores.aspartate_proline_bond_count == 2


def test_cysteine_count():
    scores = ManufacturabilityScores.from_amino_acids("C" * 7)
    assert scores.cysteine_count == 7


def cterm_7mer_gravy_score():
    scores = ManufacturabilityScores.from_amino_acids("QLFY" + "A" * 7)
    # hydropathy of alanine is 1.8 from Kyte & Doolittle 1982
    assert scores.cterm_7mer_gravy_score == 1.8


def max_7mer_gravy_score():
    scores = ManufacturabilityScores.from_amino_acids("H" * 3 + "A" * 7)
    # hydropathy of alanine is 1.8, histidine is -3.2
    # from Kyte & Doolittle 1982
    assert scores.max_7mer_gravy_score == 1.8


def _legacy_manufacturability_tuple(scores, thresholds):
    """Inline replica of the legacy 10-entry tuple — parity baseline."""
    cterm_7mer_gravy = scores.cterm_7mer_gravy_score
    max_7mer_gravy = scores.max_7mer_gravy_score
    return (
        scores.cysteine_count,
        max(0, cterm_7mer_gravy - thresholds["max_c_terminal_hydropathy"]),
        max(0, max_7mer_gravy - thresholds["max_kmer_hydropathy_high_priority"]),
        scores.difficult_n_terminal_residue,
        scores.c_terminal_cysteine,
        scores.c_terminal_proline,
        scores.n_terminal_asparagine,
        scores.aspartate_proline_bond_count,
        max(0, max_7mer_gravy - thresholds["max_kmer_hydropathy_low_priority"]),
        max(0, thresholds["min_kmer_hydropathy"] - max_7mer_gravy),
    )


@pytest.mark.parametrize("seq", [
    "A" * 25,
    "C" * 7 + "A" * 18,
    "Q" + "A" * 23 + "P",
    "DP" + "A" * 21 + "NP",
    "H" * 3 + "A" * 22,
])
def test_default_rules_match_legacy_tuple(seq):
    scores = ManufacturabilityScores.from_amino_acids(seq)
    legacy = _legacy_manufacturability_tuple(scores, DEFAULT_THRESHOLDS)
    registry_tuple = compute_manufacturability_tuple(scores, DEFAULT_THRESHOLDS)
    assert legacy == registry_tuple


def test_rule_registry_covers_default_order():
    # Every rule in the default order must be in the registry.
    for rule_name in DEFAULT_MANUFACTURABILITY_RULES:
        assert rule_name in MANUFACTURABILITY_RULE_REGISTRY


def test_rule_registry_unknown_rule_raises():
    scores = ManufacturabilityScores.from_amino_acids("A" * 25)
    with pytest.raises(ValueError, match="Unknown manufacturability rule"):
        compute_manufacturability_tuple(scores, DEFAULT_THRESHOLDS, rules=["bogus"])


def test_rule_registry_reorders_and_filters():
    scores = ManufacturabilityScores.from_amino_acids("C" * 3 + "A" * 22)
    # Only two rules, reversed from default order.
    result = compute_manufacturability_tuple(
        scores,
        DEFAULT_THRESHOLDS,
        rules=["cterm_hydropathy", "cysteine_count"],
    )
    assert len(result) == 2
    # cysteine_count should now be at index 1
    assert result[1] == scores.cysteine_count


def test_n_terminal_methionine_flag():
    """N-terminal Met flagged as True; other residues flagged False.
    (Manufacturer concern — oxidation risk. See vaxrank#188.)"""
    assert ManufacturabilityScores.from_amino_acids("M" + "A" * 6).n_terminal_methionine
    assert not ManufacturabilityScores.from_amino_acids("A" * 7).n_terminal_methionine
    # Internal M doesn't count
    assert not ManufacturabilityScores.from_amino_acids("AMAAAAA").n_terminal_methionine


def test_n_terminal_methionine_is_opt_in_not_in_default_rules():
    """Legacy parity: n_terminal_methionine must NOT appear in the default
    rule tuple. Users opt in through manufacturability.rules config."""
    assert "n_terminal_methionine" not in DEFAULT_MANUFACTURABILITY_RULES
    assert "n_terminal_methionine" in MANUFACTURABILITY_RULE_REGISTRY


def test_n_terminal_methionine_rule_in_registry_drives_sort():
    """Using the rule via the registry: N-term M peptide gets a worse
    score than an otherwise-identical N-term A peptide."""
    met_scores = ManufacturabilityScores.from_amino_acids("M" + "A" * 24)
    ala_scores = ManufacturabilityScores.from_amino_acids("A" * 25)
    met_tuple = compute_manufacturability_tuple(
        met_scores, DEFAULT_THRESHOLDS, rules=["n_terminal_methionine"],
    )
    ala_tuple = compute_manufacturability_tuple(
        ala_scores, DEFAULT_THRESHOLDS, rules=["n_terminal_methionine"],
    )
    # Lexicographic: smaller tuple = better, so methionine peptide sorts worse
    assert met_tuple > ala_tuple
