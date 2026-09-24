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

"""Tests for the variant-level sort in ``ranked_vaccine_peptides`` —
verifies the three-level tiebreak (combined_score, n_rna_alt,
target_epitope_score) introduced for vaxrank#151 and refined in #394."""

from types import SimpleNamespace

import pytest

from vaxrank.core_logic import ranked_vaccine_peptides


def _fake_peptide(
        combined_score,
        n_alt_reads,
        target_epitope_score,
        n_rna_alt=None):
    return SimpleNamespace(
        combined_score=combined_score,
        target_epitope_score=target_epitope_score,
        mutant_protein_fragment=SimpleNamespace(
            n_alt_reads=n_alt_reads,
            n_rna_alt=n_alt_reads if n_rna_alt is None else n_rna_alt,
            rna_evidence_method="rna_alignment",
            rna_evidence_subject="fragments",
        ),
    )


def test_primary_sort_by_combined_score():
    v_low = "v_low"
    v_high = "v_high"
    d = {
        v_low: [_fake_peptide(combined_score=1.0, n_alt_reads=100, target_epitope_score=0.5)],
        v_high: [_fake_peptide(combined_score=9.0, n_alt_reads=10, target_epitope_score=0.9)],
    }
    ranked = ranked_vaccine_peptides(d)
    assert [variant for variant, _ in ranked] == [v_high, v_low]


def test_tiebreak_uses_independent_rna_evidence_when_combined_score_ties():
    """A duplicate-heavy locus does not win an otherwise tied score."""
    v_duplicate_heavy = "v_duplicate_heavy"
    v_clean = "v_clean"
    d = {
        v_duplicate_heavy: [_fake_peptide(
            combined_score=0.0,
            n_alt_reads=40,
            n_rna_alt=9,
            target_epitope_score=0.5,
        )],
        v_clean: [_fake_peptide(
            combined_score=0.0,
            n_alt_reads=25,
            n_rna_alt=16,
            target_epitope_score=0.5,
        )],
    }
    ranked = ranked_vaccine_peptides(d)
    assert [variant for variant, _ in ranked] == [v_clean, v_duplicate_heavy]


def test_tiebreak_by_target_epitope_score_when_higher_tiers_tie():
    """When combined_score and n_rna_alt both tie, target_epitope_score
    breaks the tie — higher MHC binding score wins."""
    v_low_epitope = "v_low_epitope"
    v_high_epitope = "v_high_epitope"
    d = {
        v_low_epitope: [_fake_peptide(combined_score=0.0, n_alt_reads=10, target_epitope_score=0.1)],
        v_high_epitope: [_fake_peptide(combined_score=0.0, n_alt_reads=10, target_epitope_score=0.9)],
    }
    ranked = ranked_vaccine_peptides(d)
    assert [variant for variant, _ in ranked] == [v_high_epitope, v_low_epitope]


def test_empty_peptide_list_sorts_last():
    """An absent construct cannot outrank a valid negative-score construct."""
    v_empty = "v_empty"
    v_normal = "v_normal"
    d = {
        v_empty: [],
        v_normal: [_fake_peptide(combined_score=-1.0, n_alt_reads=5, target_epitope_score=0.2)],
    }
    ranked = ranked_vaccine_peptides(d)
    assert [variant for variant, _ in ranked] == [v_normal, v_empty]


def test_all_zeros_keep_insertion_order():
    """Python's sort is stable, so variants tied on every key preserve
    their original dict insertion order."""
    v1, v2, v3 = "v1", "v2", "v3"
    zero = _fake_peptide(combined_score=0.0, n_alt_reads=0, target_epitope_score=0.0)
    d = {v1: [zero], v2: [zero], v3: [zero]}
    ranked = ranked_vaccine_peptides(d)
    assert [variant for variant, _ in ranked] == [v1, v2, v3]


@pytest.mark.parametrize("unavailable", [
    {"rna_evidence_subject": "reads"},
    {"rna_evidence_method": "rna_depth_x_vaf"},
    {"rna_evidence_method": ""},
    {"rna_evidence_subject": ""},
    {"n_rna_alt": None},
    {"n_rna_alt": float("nan")},
    None,
])
def test_incomparable_rna_skips_the_tier_for_the_whole_score_tie(unavailable):
    high_rna = _fake_peptide(1.0, 100, 0.1)
    low_rna = _fake_peptide(1.0, 1, 0.5)
    other = _fake_peptide(1.0, 0, 0.9)
    if unavailable is None:
        other.mutant_protein_fragment = None
    else:
        for key, value in unavailable.items():
            setattr(other.mutant_protein_fragment, key, value)
    ranked = ranked_vaccine_peptides({
        "high_rna": [high_rna], "low_rna": [low_rna], "other": [other],
    })
    assert [source for source, _ in ranked] == ["other", "low_rna", "high_rna"]
    if unavailable is not None:
        for key, value in unavailable.items():
            assert getattr(other.mutant_protein_fragment, key) is value


def test_unrelated_score_group_does_not_disable_comparable_rna_tiebreak():
    high_rna = _fake_peptide(1.0, 10, 0.1)
    low_rna = _fake_peptide(1.0, 0, 0.9)
    other = _fake_peptide(2.0, 0, 0.1)
    other.mutant_protein_fragment = None
    ranked = ranked_vaccine_peptides({
        "low_rna": [low_rna], "high_rna": [high_rna], "other": [other],
    })
    assert [source for source, _ in ranked] == ["other", "high_rna", "low_rna"]


def test_final_ranking_preserves_selected_window_and_alternates():
    selected = _fake_peptide(1.0, 10, 0.1)
    alternate = _fake_peptide(2.0, 10, 0.9)
    other = _fake_peptide(1.5, 10, 0.5)
    windows = [selected, alternate]
    ranked = ranked_vaccine_peptides({"windows": windows, "other": [other]})
    assert [source for source, _ in ranked] == ["other", "windows"]
    assert ranked[1][1] is windows
    assert windows == [selected, alternate]
