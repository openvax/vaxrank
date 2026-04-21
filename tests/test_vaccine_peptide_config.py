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

import argparse
import os
import tempfile
from types import SimpleNamespace

import numpy as np
import pytest

from vaxrank.cli.vaccine_config_args import vaccine_config_from_args
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
    # expression_score is an honest metric; mode only affects combined_score
    assert vp.expression_score == pytest.approx(np.sqrt(9))
    assert vp.combined_score == pytest.approx(9.0 * vp.mutant_epitope_score)


def test_combined_score_epitope_only_mode():
    fragment = _make_fragment(n_alt_reads=9)
    vp = VaccinePeptide(
        mutant_protein_fragment=fragment,
        epitope_predictions=[_make_mutant_epitope()],
        combined_score_mode="epitope_only",
    )
    assert vp.expression_score == pytest.approx(np.sqrt(9))
    assert vp.combined_score == pytest.approx(vp.mutant_epitope_score)


def test_vaccine_peptide_rejects_unknown_combined_score_mode():
    with pytest.raises(ValueError, match="combined_score_mode"):
        VaccinePeptide(
            mutant_protein_fragment=_make_fragment(),
            epitope_predictions=[_make_mutant_epitope()],
            combined_score_mode="garbage",
        )


def test_end_to_end_yaml_rules_through_vaccine_config_to_peptide():
    """YAML with manufacturability.rules loads through load_vaxrank_config,
    then vaccine_config_from_args, producing a tuple on VaccineConfig.
    A VaccinePeptide built with those rules gets a sort tuple of the
    expected length."""
    yaml_content = """
manufacturability:
  rules:
    - cysteine_count
    - cterm_hydropathy
    - aspartate_proline
vaccine_peptides:
  combined_score_mode: epitope_only
"""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".yaml", delete=False
    ) as f:
        f.write(yaml_content)
        config_path = f.name

    try:
        args = argparse.Namespace(
            config=config_path,
            vaccine_peptide_length=None,
            padding_around_mutation=None,
            max_vaccine_peptides_per_variant=None,
            num_epitopes_per_vaccine_peptide=None,
            config_set_overrides=None,
            config_expr_overrides=None,
        )
        vc = vaccine_config_from_args(args)
        # Type contract: VaccineConfig declares Optional[tuple[str, ...]]
        assert vc.manufacturability_rules == (
            "cysteine_count", "cterm_hydropathy", "aspartate_proline",
        )
        assert isinstance(vc.manufacturability_rules, tuple)
        assert vc.combined_score_mode == "epitope_only"

        vp = VaccinePeptide(
            mutant_protein_fragment=_make_fragment(),
            epitope_predictions=[_make_mutant_epitope()],
            manufacturability_rules=vc.manufacturability_rules,
            combined_score_mode=vc.combined_score_mode,
        )
        tup = vp.peptide_synthesis_difficulty_score_tuple(
            rules=vp.manufacturability_rules,
        )
        assert len(tup) == 3
        # combined_score under epitope_only is just the epitope score
        assert vp.combined_score == pytest.approx(vp.mutant_epitope_score)
    finally:
        os.unlink(config_path)


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
    assert "ranking_rules" not in d


def test_vaccine_peptide_ranking_rules_tuple_coercion():
    """List of ranking rule names coerced to tuple (same contract as
    manufacturability_rules)."""
    vp = VaccinePeptide(
        mutant_protein_fragment=_make_fragment(),
        epitope_predictions=[_make_mutant_epitope()],
        ranking_rules=["mutant_epitope_score", "n_alt_reads"],
    )
    assert vp.ranking_rules == ("mutant_epitope_score", "n_alt_reads")
    assert isinstance(vp.ranking_rules, tuple)
    # Sort tuple length reflects the custom 2-rule list.
    assert len(vp.lexicographic_sort_key()) == 2


def test_vaccine_peptide_custom_ranking_rules_in_to_dict():
    vp = VaccinePeptide(
        mutant_protein_fragment=_make_fragment(),
        epitope_predictions=[_make_mutant_epitope()],
        ranking_rules=["mutant_epitope_score"],
    )
    d = vp.to_dict()
    assert d["ranking_rules"] == ["mutant_epitope_score"]


def test_end_to_end_yaml_ranking_rules():
    """YAML with vaccine_peptides.ranking_rules loads through the config
    pipeline to a VaccineConfig tuple and threads into VaccinePeptide."""
    yaml_content = """
vaccine_peptides:
  ranking_rules:
    - mutant_epitope_score
    - n_alt_reads
    - manufacturability
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        f.write(yaml_content)
        config_path = f.name
    try:
        args = argparse.Namespace(
            config=config_path,
            vaccine_peptide_length=None,
            padding_around_mutation=None,
            max_vaccine_peptides_per_variant=None,
            num_epitopes_per_vaccine_peptide=None,
            config_set_overrides=None,
            config_expr_overrides=None,
        )
        vc = vaccine_config_from_args(args)
        assert vc.ranking_rules == (
            "mutant_epitope_score", "n_alt_reads", "manufacturability",
        )
        assert isinstance(vc.ranking_rules, tuple)
        # Custom rules drop the extra_score_tuple tiers (n_alt_reads_supporting,
        # wildtype_epitope_score, n_mutant_amino_acids, mutation_distance_from_edge).
        vp = VaccinePeptide(
            mutant_protein_fragment=_make_fragment(),
            epitope_predictions=[_make_mutant_epitope()],
            ranking_rules=vc.ranking_rules,
        )
        # 2 essential + 10 manufacturability defaults = 12 entries
        assert len(vp.lexicographic_sort_key()) == 12
    finally:
        os.unlink(config_path)


def test_vaccine_config_rejects_unknown_ranking_rule():
    from vaxrank.vaccine_config import VaccineConfig

    with pytest.raises(ValueError, match="Unknown ranking rule"):
        VaccineConfig(ranking_rules=("not_a_real_rule",))
