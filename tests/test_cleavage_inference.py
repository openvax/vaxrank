"""Full-context orchestration with explicit failures and real motif semantics."""

from unittest.mock import patch
from types import SimpleNamespace

import mhctools
import pytest
from mhctools.cleavage import CleavageInput
from mhctools.peptidases import get_cleavage_model

from vaxrank.cleavage_inference import audit_pepsickle_inputs, audit_peptidase_inputs

from .test_cleavage_profile import MODEL


def identity():
    return patch("vaxrank.cleavage_inference._pepsickle_model_identity",
                 return_value=(MODEL, "a" * 64, ""))


def test_contexts_deduplicated_without_losing_occurrences_or_source_offsets():
    inputs = [CleavageInput("ACDE", source_id="a"),
              CleavageInput("ACDE", source_id="b", source_start=40),
              CleavageInput("ACDEKK", source_id="final")]
    with identity(), patch.object(mhctools.Pepsickle, "cleavage_probs_many", return_value={
            "ACDE": [.1, .2, .3, 0.], "ACDEKK": [.8, .7, .6, .5, .4, 0.]}) as run:
        profiles = audit_pepsickle_inputs(inputs)
    run.assert_called_once_with(("ACDE", "ACDEKK"))
    assert [p.peptide for p in profiles] == inputs
    assert profiles[0].sites[0].score != profiles[2].sites[0].score
    assert profiles[1].peptide.source_start == 40


def test_chemistry_and_missing_model_never_enter_sequence_only_backend():
    inputs = [CleavageInput("ACDE", n_term="acetylated"), CleavageInput("ACDE")]
    with patch("vaxrank.cleavage_inference._pepsickle_model_identity", return_value=(MODEL, "", "")), \
            patch.object(mhctools.Pepsickle, "cleavage_probs_many") as run:
        profiles = audit_pepsickle_inputs(inputs)
    run.assert_not_called()
    assert [p.status for p in profiles] == ["unassessed", "unassessed"]
    assert profiles[0].reason_codes == ("unsupported_terminal_chemistry",)
    assert profiles[1].reason_codes == ("pepsickle_model_assets_unavailable",)


def test_bad_array_does_not_destroy_other_sequence_evidence():
    with identity(), patch.object(mhctools.Pepsickle, "cleavage_probs_many", return_value={
            "ACDE": [float("nan")] * 4, "ACDEK": [.1] * 5}):
        profiles = audit_pepsickle_inputs([CleavageInput("ACDE"), CleavageInput("ACDEK")])
    assert [p.status for p in profiles] == ["failed", "observations_returned"]
    assert profiles[0].unassessed_bonds == (1, 2, 3)


def test_batch_failure_is_explicit_not_no_cuts():
    with identity(), patch.object(mhctools.Pepsickle, "cleavage_probs_many", side_effect=RuntimeError("failed")):
        p, = audit_pepsickle_inputs([CleavageInput("ACDE")])
    assert p.status == "failed" and p.error_message == "failed"
    assert p.sites == ()


def test_real_cpn_rule_assesses_only_the_exposed_terminal_bond():
    # This tests the released API's stated motif behavior, not serum kinetics.
    predictor = get_cleavage_model("cpn-basic")
    peptide = CleavageInput("ACDEKK", source_id="JLF-shaped-fixture")
    p, repeated = audit_peptidase_inputs([peptide, peptide], predictor)
    assert p is repeated
    assert p.model.evidence == "motif_rule"
    assert [(s.bond, s.status, s.score) for s in p.sites] == [(5, "matched", None)]
    assert p.unassessed_bonds == (1, 2, 3, 4)
    assert p.raw_residue_scores == ()


def test_canonical_peptidase_chemistry_gate_is_preserved():
    p, = audit_peptidase_inputs([CleavageInput("ACDEK", c_term="amidated")],
                                get_cleavage_model("cpn-basic"))
    assert p.status == "unassessed" and p.reason_codes


def test_unreadable_weights_are_explicit_and_do_not_crash_the_audit(tmp_path):
    weights = tmp_path / "trained_model_dict.pickle"
    weights.write_bytes(b"fixture")
    with patch("vaxrank.cleavage_inference.util.find_spec", return_value=SimpleNamespace(
            origin=str(tmp_path / "__init__.py"))), \
            patch("vaxrank.cleavage_inference.metadata.version", return_value="fixture"), \
            patch("pathlib.Path.open", side_effect=PermissionError("weights not readable")), \
            patch.object(mhctools.Pepsickle, "cleavage_probs_many") as run:
        p, = audit_pepsickle_inputs([CleavageInput("ACDE")])
    run.assert_not_called()
    assert p.status == "unassessed" and p.error_message == "weights not readable"


@pytest.mark.parametrize("threshold", [None, True, float("nan"), float("inf"), -1., 2., ".5"])
def test_invalid_model_settings_fail_before_loading_assets(threshold):
    with pytest.raises(ValueError, match="settings"):
        audit_pepsickle_inputs([CleavageInput("ACDE")], threshold=threshold)


def test_empty_input_does_not_load_models():
    with patch("vaxrank.cleavage_inference._pepsickle_model_identity") as identify:
        assert audit_pepsickle_inputs([]) == ()
    identify.assert_not_called()
