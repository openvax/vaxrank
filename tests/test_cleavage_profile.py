"""Exact processing evidence and coverage, without a biological verdict."""

from dataclasses import FrozenInstanceError, replace

from mhctools.cleavage import CleavageInput, CleavageModel, CleavageResult, CleavageSite
import pytest

from vaxrank.cleavage_profile import (
    CleavageProfile, profile_from_cleavage_result, profile_from_residue_scores,
)
from vaxrank.native_serialization import from_native_json, to_native_json


MODEL = CleavageModel(
    "fixture", "1", "proteasome", "", "", ("proteasomal",),
    "quantitative_model", ("synthetic-unit-test",), "synthetic",
    "Not biological evidence", "fixture_probability", "dimensionless")


def profile(sequence="ACDEFG", scores=(.1, .2, .3, .4, .5, 0.), **kw):
    return profile_from_residue_scores(CleavageInput(sequence), MODEL, scores, **kw)


def test_every_internal_bond_retained_and_endpoint_not_counted():
    p = profile()
    assert [s.bond for s in p.sites] == [1, 2, 3, 4, 5]
    assert p.terminal_output == 0.
    assert p.raw_residue_scores == (.1, .2, .3, .4, .5, 0.)
    assert p.unassessed_bonds == ()
    whole = p.interval_evidence(0, 6)
    assert whole.internal_sites == p.sites
    assert whole.n_boundary is whole.c_boundary is None
    assert whole.n_boundary_status == whole.c_boundary_status == "sequence_endpoint"


def test_internal_cuts_and_ligand_boundaries_remain_separate():
    p = profile()
    evidence = p.interval_evidence(1, 5)
    assert [s.bond for s in evidence.internal_sites] == [2, 3, 4]
    assert evidence.n_boundary.bond == 1 and evidence.c_boundary.bond == 5
    assert evidence.unassessed_internal_bonds == ()
    junction = p.interval_evidence(3, 3)
    assert junction.internal_sites == ()
    assert junction.n_boundary == junction.c_boundary == p.sites[2]


@pytest.mark.parametrize("scores", [[], [.1], [.1] * 7, [.1] * 5 + [True],
                                  [.1] * 5 + [float("nan")], [.1] * 5 + [float("inf")],
                                  [.1] * 5 + [-.1], [.1] * 5 + [1.1], [.1] * 5 + [".2"]])
def test_invalid_arrays_are_not_clean_profiles(scores):
    with pytest.raises(ValueError):
        profile(scores=scores)


def test_exact_p1_padding_coordinates():
    p = profile("A" * 20, [.1] * 20, upstream_context=8, downstream_context=8)
    assert p.padded_bonds == tuple(range(1, 9)) + tuple(range(13, 20))
    assert p.reason_codes == ("terminal_context_padded",)


def test_no_internal_bond_is_not_a_negative_cleavage_prediction():
    p = profile("A", [0.])
    assert p.status == "no_observations"
    assert p.reason_codes == ("sequence_has_no_internal_bonds",)
    assert not p.sites and p.terminal_output == 0.


def test_partial_motif_rule_preserves_native_semantics_and_unknown_bonds():
    motif = replace(MODEL, evidence="motif_rule", score_name=None, score_units=None)
    result = CleavageResult(CleavageInput("ACDEFG"), motif,
                            (CleavageSite(5, "not_matched", "limited terminal rule"),))
    p = profile_from_cleavage_result(result)
    assert p.status == "observations_returned"
    assert p.unassessed_bonds == (1, 2, 3, 4)
    evidence = p.interval_evidence(0, 6)
    assert evidence.unassessed_internal_bonds == (1, 2, 3, 4)
    assert evidence.internal_sites[0].score is None
    assert p.raw_residue_scores == () and p.terminal_output is None


def test_unsupported_chemistry_remains_explicit():
    result = CleavageResult(CleavageInput("ACDE", n_term="acetylated"), MODEL,
                            unsupported_reason="model excludes acetylation")
    p = profile_from_cleavage_result(result)
    assert p.status == "unassessed" and p.unassessed_bonds == (1, 2, 3)
    assert p.reason_codes == ("model excludes acetylation",)


def test_native_roundtrip_preserves_all_nested_types_and_provenance():
    p = profile(backend_version="3.41.0", settings=(("human_only", "False"),),
                model_asset_sha256="a" * 64, upstream_context=8, downstream_context=8)
    restored = from_native_json(to_native_json(p), CleavageProfile)
    assert restored == p
    assert isinstance(restored.peptide, CleavageInput)
    assert isinstance(restored.model, CleavageModel)
    assert all(isinstance(site, CleavageSite) for site in restored.sites)
    interval = p.interval_evidence(1, 5)
    assert from_native_json(to_native_json(interval), type(interval)) == interval
    with pytest.raises(FrozenInstanceError):
        restored.status = "unassessed"


@pytest.mark.parametrize("updates", [
    {"status": "unassessed"}, {"padded_bonds": (0,)}, {"padded_bonds": (6,)},
    {"padded_bonds": (1, 1)}, {"raw_residue_scores": (.9,) * 6},
    {"model_asset_sha256": "not-a-digest"}, {"settings": (("a", "1"), ("a", "2"))},
])
def test_inconsistent_serialized_profiles_fail_validation(updates):
    with pytest.raises(ValueError):
        replace(profile(), **updates)


def test_repeated_motifs_retain_distinct_context_bonds():
    p = profile("ACDACD", (.1, .2, .3, .7, .8, 0.))
    assert [s.score for s in p.interval_evidence(0, 3).internal_sites] == [.1, .2]
    assert [s.score for s in p.interval_evidence(3, 6).internal_sites] == [.7, .8]


def test_source_offsets_and_terminal_chemistry_survive_roundtrip():
    result = CleavageResult(
        CleavageInput("ACDE", c_term="amidated", source_id="parent", source_start=37), MODEL,
        unsupported_reason="unsupported chemistry")
    p = profile_from_cleavage_result(result)
    assert from_native_json(to_native_json(p), CleavageProfile) == p


@pytest.mark.parametrize("score", [-2., 400.])
def test_native_quantitative_scores_are_not_forced_into_probability_bounds(score):
    model = replace(MODEL, score_name="native_activity", score_units="native_units")
    result = CleavageResult(CleavageInput("ACDE"), model,
                            (CleavageSite(2, "scored", "fixture native output", score),))
    p = profile_from_cleavage_result(result)
    assert p.sites[0].score == score
    assert p.unassessed_bonds == (1, 3)
    assert p.raw_residue_scores == ()
