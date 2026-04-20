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

"""Tests for configurable manufacturability rules and combined_score on
``VaccinePeptide``. Separate from ``test_epitope_prediction.py`` because
these don't require an MHC predictor."""

from types import SimpleNamespace

import numpy as np
import pytest

from vaxrank.manufacturability import DEFAULT_MANUFACTURABILITY_RULES
from vaxrank.vaccine_peptide import VaccinePeptide


def _make_fragment(seq="A" * 25, n_alt_reads=9):
    """Minimal duck-typed MutantProteinFragment for VaccinePeptide tests."""
    return SimpleNamespace(
        amino_acids=seq,
        n_alt_reads=n_alt_reads,
        n_alt_reads_supporting_protein_sequence=n_alt_reads,
        n_mutant_amino_acids=1,
        mutation_distance_from_edge=5,
    )


def _make_mutant_epitope(ic50=10.0):
    """Duck-typed EpitopePrediction that overlaps mutation, not in reference."""
    return SimpleNamespace(
        ic50=ic50,
        percentile_rank=0.5,
        overlaps_mutation=True,
        occurs_in_reference=False,
        logistic_epitope_score=lambda **kw: 0.9,
    )


def test_default_manufacturability_tuple_length():
    """Default sort tuple length = 10 (legacy behavior, unchanged)."""
    vp = VaccinePeptide(
        mutant_protein_fragment=_make_fragment(),
        epitope_predictions=[_make_mutant_epitope()],
    )
    tup = vp.peptide_synthesis_difficulty_score_tuple()
    assert len(tup) == 10
    assert len(DEFAULT_MANUFACTURABILITY_RULES) == 10


def test_custom_rules_shorten_tuple():
    vp = VaccinePeptide(
        mutant_protein_fragment=_make_fragment(),
        epitope_predictions=[_make_mutant_epitope()],
        manufacturability_rules=["cysteine_count", "cterm_hydropathy"],
    )
    tup = vp.peptide_synthesis_difficulty_score_tuple(
        rules=vp.manufacturability_rules,
    )
    assert len(tup) == 2


def test_combined_score_default_mode_is_legacy():
    fragment = _make_fragment(n_alt_reads=9)
    vp = VaccinePeptide(
        mutant_protein_fragment=fragment,
        epitope_predictions=[_make_mutant_epitope()],
    )
    # Legacy: sqrt(9) * mutant_epitope_score
    assert vp.expression_score == pytest.approx(np.sqrt(9))
    assert vp.combined_score == pytest.approx(np.sqrt(9) * vp.mutant_epitope_score)


def test_combined_score_reads_times_epitope_mode():
    fragment = _make_fragment(n_alt_reads=9)
    vp = VaccinePeptide(
        mutant_protein_fragment=fragment,
        epitope_predictions=[_make_mutant_epitope()],
        combined_score_mode="reads_times_epitope",
    )
    assert vp.expression_score == 9.0
    assert vp.combined_score == pytest.approx(9.0 * vp.mutant_epitope_score)


def test_combined_score_epitope_only_mode():
    fragment = _make_fragment(n_alt_reads=9)
    vp = VaccinePeptide(
        mutant_protein_fragment=fragment,
        epitope_predictions=[_make_mutant_epitope()],
        combined_score_mode="epitope_only",
    )
    assert vp.expression_score == 1.0
    assert vp.combined_score == pytest.approx(vp.mutant_epitope_score)


def test_vaccine_peptide_roundtrip_preserves_custom_config():
    """to_dict / from_dict preserves non-default rules + score mode."""
    fragment = _make_fragment()
    vp = VaccinePeptide(
        mutant_protein_fragment=fragment,
        epitope_predictions=[_make_mutant_epitope()],
        manufacturability_rules=["cysteine_count"],
        combined_score_mode="epitope_only",
    )
    d = vp.to_dict()
    assert d["manufacturability_rules"] == ["cysteine_count"]
    assert d["combined_score_mode"] == "epitope_only"


def test_vaccine_peptide_roundtrip_omits_defaults():
    """to_dict doesn't include None/default rule + mode fields."""
    fragment = _make_fragment()
    vp = VaccinePeptide(
        mutant_protein_fragment=fragment,
        epitope_predictions=[_make_mutant_epitope()],
    )
    d = vp.to_dict()
    assert "manufacturability_rules" not in d
    assert "combined_score_mode" not in d
