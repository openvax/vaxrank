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

"""Tests for ``vaxrank.ranking`` — the rule registry that backs
``VaccinePeptide.lexicographic_sort_key``. Pinned to the legacy sort
tuple produced before ranking rules were extracted (Phase 1 of #199)."""

from types import SimpleNamespace

import pytest

from vaxrank.manufacturability import ManufacturabilityScores
from vaxrank.ranking import (
    DEFAULT_RANKING_RULES,
    MANUFACTURABILITY_SENTINEL,
    RANKING_RULE_REGISTRY,
    compute_ranking_tuple,
)


def _fake_peptide(
    amino_acids="A" * 25,
    n_alt_reads=7,
    n_alt_reads_supporting=5,
    n_mutant_amino_acids=3,
    mutation_distance_from_edge=9,
    target_epitope_score=0.73,
    self_epitope_score=0.42,
    manufacturability_rules=None,
    manufacturability_thresholds=None,
):
    """Minimal duck-typed VaccinePeptide — gives compute_ranking_tuple and
    the manufacturability expansion everything they read, nothing more."""
    fragment = SimpleNamespace(
        amino_acids=amino_acids,
        n_alt_reads=n_alt_reads,
        n_alt_reads_supporting_protein_sequence=n_alt_reads_supporting,
        n_mutant_amino_acids=n_mutant_amino_acids,
        mutation_distance_from_edge=mutation_distance_from_edge,
    )
    manufacturability_scores = ManufacturabilityScores.from_amino_acids(amino_acids)
    peptide = SimpleNamespace(
        mutant_protein_fragment=fragment,
        manufacturability_scores=manufacturability_scores,
        manufacturability_rules=manufacturability_rules,
        manufacturability_thresholds=manufacturability_thresholds or {},
        target_epitope_score=target_epitope_score,
        self_epitope_score=self_epitope_score,
    )

    # Mirror VaccinePeptide.peptide_synthesis_difficulty_score_tuple so the
    # sentinel expansion in compute_ranking_tuple gets a real tuple back.
    def _difficulty_tuple(rules=None, **thresholds):
        from vaxrank.manufacturability import compute_manufacturability_tuple

        filled = {
            "max_c_terminal_hydropathy": 1.5,
            "min_kmer_hydropathy": 0.0,
            "max_kmer_hydropathy_low_priority": 1.5,
            "max_kmer_hydropathy_high_priority": 2.5,
        }
        filled.update(thresholds)
        return compute_manufacturability_tuple(manufacturability_scores, filled, rules=rules)

    peptide.peptide_synthesis_difficulty_score_tuple = _difficulty_tuple
    return peptide


def _legacy_sort_tuple(peptide):
    """Inline replica of the pre-refactor VaccinePeptide.lexicographic_sort_key —
    kept verbatim so the parity test can't accidentally drift with the
    registry.  If the default order ever changes intentionally, update
    BOTH this and DEFAULT_RANKING_RULES together."""
    fragment = peptide.mutant_protein_fragment
    essential = (
        -round(peptide.target_epitope_score, 6),
        -fragment.n_alt_reads,
    )
    manufacturability = peptide.peptide_synthesis_difficulty_score_tuple(
        rules=peptide.manufacturability_rules,
        **peptide.manufacturability_thresholds,
    )
    extra = (
        -fragment.n_alt_reads_supporting_protein_sequence,
        round(peptide.self_epitope_score, 6),
        -fragment.n_mutant_amino_acids,
        -fragment.mutation_distance_from_edge,
    )
    return essential + manufacturability + extra


@pytest.mark.parametrize("seq", [
    "A" * 25,
    "C" * 7 + "A" * 18,
    "Q" + "A" * 23 + "P",
    "DP" + "A" * 21 + "NP",
    "H" * 3 + "A" * 22,
    "M" + "A" * 24,
])
@pytest.mark.parametrize("epitope_score", [0.0, 0.5, 1.2])
@pytest.mark.parametrize("n_alt_reads", [0, 3, 42])
def test_default_rules_match_legacy_tuple(seq, epitope_score, n_alt_reads):
    """Default registry output byte-matches the pre-refactor sort tuple
    across representative peptides, epitope scores, and read counts."""
    peptide = _fake_peptide(
        amino_acids=seq,
        n_alt_reads=n_alt_reads,
        target_epitope_score=epitope_score,
    )
    assert compute_ranking_tuple(peptide) == _legacy_sort_tuple(peptide)


def test_registry_covers_every_default_rule():
    for rule_name in DEFAULT_RANKING_RULES:
        assert rule_name in RANKING_RULE_REGISTRY


def test_manufacturability_sentinel_is_registered():
    """The sentinel string must sit in the registry so
    validators / introspection don't flag it as unknown, even though
    its value is None and expansion is handled specially."""
    assert MANUFACTURABILITY_SENTINEL in RANKING_RULE_REGISTRY
    assert RANKING_RULE_REGISTRY[MANUFACTURABILITY_SENTINEL] is None


def test_unknown_rule_raises():
    peptide = _fake_peptide()
    with pytest.raises(ValueError, match="Unknown ranking rule"):
        compute_ranking_tuple(peptide, rules=["bogus_rule"])


def test_custom_rules_reorder_and_subset():
    """Users can narrow the sort tuple — lower-priority rules drop out,
    and order in the list drives lexicographic precedence."""
    peptide = _fake_peptide()
    # Two-element tuple, sort first by n_alt_reads then by target_epitope_score.
    result = compute_ranking_tuple(
        peptide, rules=["n_alt_reads", "target_epitope_score"]
    )
    fragment = peptide.mutant_protein_fragment
    assert result == (
        -fragment.n_alt_reads,
        -round(peptide.target_epitope_score, 6),
    )


def test_custom_rules_can_drop_manufacturability():
    peptide = _fake_peptide()
    result = compute_ranking_tuple(
        peptide, rules=["target_epitope_score", "n_alt_reads"]
    )
    assert len(result) == 2
    # None of the manufacturability fields should leak through
    assert all(isinstance(v, (int, float)) for v in result)


def test_manufacturability_sentinel_respects_peptide_manufacturability_rules():
    """When a peptide uses a subset of manufacturability rules, the
    sentinel should expand to that subset, not the 10-rule default."""
    peptide = _fake_peptide(
        manufacturability_rules=("cysteine_count", "cterm_hydropathy"),
    )
    result = compute_ranking_tuple(peptide)
    # 2 (essential) + 2 (manufacturability subset) + 4 (extra) = 8
    assert len(result) == 8
