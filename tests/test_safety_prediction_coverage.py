"""Candidate coverage is a request/response check, independent of risk scoring."""

from dataclasses import replace
import json

import pandas as pd
import pytest
from topiary import TopiaryPredictor

from tests.test_near_self import _index
from tests.test_safety_assessment import _antigen, _Reference, _row
from vaxrank import MHCRequest
from vaxrank.native_serialization import from_native_json, to_native_json
from vaxrank.safety_assessment import (
    SAFETY_REASON_INCOMPLETE_PREDICTIONS,
    SAFETY_REASON_MISSING_PREDICTION_SCOPE,
    SafetyAssessmentError,
    WindowSafetyAssessment,
    assess_antigen_safety,
    assess_vaccine_antigen_window,
    safety_assessment_from_prediction_frame,
)
from vaxrank.safety_policy import SafetyEnforcementPolicy, SafetyPolicyError, evaluate_antigen_safety_policy


REQUESTS = tuple(
    MHCRequest(kind, model, version, allele, length)
    for kind, model, version in (
        ("pMHC_affinity", "mhcflurry", "2.1.4"),
        ("pMHC_affinity", "mhcflurry", "2.2.0"),
        ("pMHC_affinity", "alternate", "1"),
        ("pMHC_presentation", "presentation", "2.1.4"),
    )
    for allele in ("HLA-A*02:01", "HLA-B*07:02")
    for length in (9, 10)
)


def _rows(requests=REQUESTS, sequence=None):
    sequence = sequence or _antigen().amino_acids
    return [
        dict(_row(sequence[offset:offset + request.peptide_length], offset, kind=request.kind,
                  predictor=request.predictor_name, allele=request.allele),
             predictor_version=request.predictor_version)
        for request in requests
        for offset in range(len(sequence) - request.peptide_length + 1)
    ]


def _window(rows, requests=REQUESTS):
    return safety_assessment_from_prediction_frame(
        pd.DataFrame(rows), antigen=_antigen(), reference_proteome=_Reference(),
        prediction_requests=requests,
    )


def _decision(window, mode="enforce"):
    return evaluate_antigen_safety_policy(
        assess_antigen_safety(window, _index()), SafetyEnforcementPolicy(mode=mode),
    )


@pytest.mark.parametrize("omission", ["allele", "length", "kind", "model", "model_version", "position", "all"])
def test_missing_requested_evidence_fails_even_when_returned_rows_are_consistent(omission):
    rows = _rows()
    omit = {
        "allele": lambda row: row["allele"] == "HLA-B*07:02",
        "length": lambda row: row["peptide_length"] == 10,
        "kind": lambda row: row["kind"] == "pMHC_presentation",
        "model": lambda row: row["prediction_method_name"] == "alternate",
        "model_version": lambda row: row["predictor_version"] == "2.2.0",
        "position": lambda row: row["peptide_offset"] == 3,
        "all": lambda row: True,
    }[omission]
    window = _window([row for row in rows if not omit(row)])
    assert not window.has_complete_prediction_coverage
    assert window.prediction_coverage_gaps
    audit = _decision(window, "audit")
    assert audit.finding == "incomplete"
    assert audit.action == "allow"
    assert not audit.has_complete_coverage
    assert SAFETY_REASON_INCOMPLETE_PREDICTIONS in audit.reason_codes
    with pytest.raises(SafetyPolicyError) as error:
        _decision(window)
    assert error.value.decision.action == "error"
    assert error.value.decision.assessment.window_assessment == window
    if omission == "position":
        assert all(gap["missing_window_offsets"] == (3,) for gap in window.prediction_coverage_gaps)


@pytest.mark.parametrize("percentile", [None, float("nan"), float("inf")])
def test_nonfinite_percentile_is_a_gap_even_when_other_metrics_and_models_exist(percentile):
    rows = _rows()
    rows[3]["percentile_rank"] = percentile
    window = _window(rows)
    gap, = window.prediction_coverage_gaps
    assert gap["missing_window_offsets"] == ()
    assert gap["nonfinite_window_offsets"] == (3,)
    with pytest.raises(SafetyPolicyError):
        _decision(window)


def test_original_sparse_reproduction_and_historical_scope_are_not_complete():
    sparse = [_row(_antigen().amino_acids[:9], 0, kind=kind)
              for kind in ("pMHC_affinity", "pMHC_presentation")]
    window = _window(sparse, ())
    assert window.prediction_coverage_reason_codes == (SAFETY_REASON_MISSING_PREDICTION_SCOPE,)
    with pytest.raises(SafetyPolicyError):
        _decision(window)
    # A complete-looking historical output still cannot declare its own scope.
    with pytest.raises(SafetyPolicyError):
        _decision(_window(_rows(), ()))
    disabled = _decision(window, "off")
    assert disabled.finding == "disabled" and disabled.action == "allow"
    assert not disabled.has_complete_coverage


def test_complete_request_passes_and_remains_complete_after_lossless_replay():
    window = _window(_rows())
    assert window.has_complete_prediction_coverage
    decision = _decision(window)
    assert decision.finding == "safe" and decision.has_complete_coverage
    assert WindowSafetyAssessment.from_json(window.to_json()) == window
    replay = from_native_json(to_native_json(window), WindowSafetyAssessment)
    assert _decision(replay) == decision
    with pytest.raises(ValueError, match="requires an incomplete finding"):
        replace(decision, assessment=replace(decision.assessment,
                window_assessment=replace(window, prediction_requests=())))


def test_report_and_replay_retain_requests_and_exact_gaps():
    window = _window(_rows()[1:])
    replay = from_native_json(to_native_json(window), WindowSafetyAssessment)
    assert replay.prediction_coverage_gaps == window.prediction_coverage_gaps
    payload = json.loads(json.dumps(_decision(replay, "audit").to_report_dict()))
    saved = payload["assessment"]["window_assessment"]
    gap, = saved["prediction_coverage_gaps"]
    assert gap["missing_window_offsets"] == [0]
    assert gap["request"]["allele"] == "HLA-A*02:01"
    assert len(saved["prediction_requests"]) == len(REQUESTS)
    # Dropping a new field from an old-style native archive cannot imply success.
    legacy = json.loads(to_native_json(window))
    del legacy["prediction_requests"]
    historical = from_native_json(json.dumps(legacy), WindowSafetyAssessment)
    assert not historical.has_complete_prediction_coverage


def test_complete_evidence_keeps_audit_and_enforcement_actions_distinct():
    rows = _rows()
    for row in rows:
        if row["peptide_offset"] == 11 and row["peptide_length"] == 9:
            row["percentile_rank"] = 0.5  # A complete, strong non-target occurrence.
    window = _window(rows)
    enforce = _decision(window)
    audit = _decision(window, "audit")
    assert enforce.finding == audit.finding == "unsafe"
    assert enforce.action == "repair" and audit.action == "allow"
    assert enforce.has_complete_coverage and audit.has_complete_coverage
    assert "threshold_non_target_ligand" in enforce.reason_codes


def test_mixed_class_lengths_are_checked_per_allele_without_cross_product():
    requests = (
        MHCRequest("pMHC_affinity", "class-I", "1", "HLA-A*02:01", 9),
        MHCRequest("pMHC_affinity", "class-I", "1", "HLA-A*02:01", 10),
        MHCRequest("pMHC_affinity", "class-II", "1", "HLA-DRB1*01:01", 15),
    )
    window = _window(_rows(requests), requests)
    assert window.has_complete_prediction_coverage
    assert {len(ligand.peptide) for ligand in window.ligands} == {9, 10, 15}


@pytest.mark.parametrize("field,value", [
    ("allele", "HLA-C*04:01"), ("kind", "pMHC_stability"),
    ("prediction_method_name", "unrequested"), ("predictor_version", "99"),
])
def test_output_outside_request_is_rejected(field, value):
    rows = _rows()
    rows[0][field] = value
    with pytest.raises(SafetyAssessmentError, match="outside the declared request"):
        _window(rows)


def test_high_level_api_freezes_scope_before_model_runs():
    requests = list(REQUESTS)

    class Predictor(TopiaryPredictor):
        def __init__(self):
            pass

        def predict_from_named_sequences(self, sequences):
            requests.clear()
            return pd.DataFrame(_rows()[:1])

    window = assess_vaccine_antigen_window(
        Predictor(), _antigen(), reference_proteome=_Reference(), prediction_requests=requests,
    )
    assert len(window.prediction_requests) == len(REQUESTS)
    assert not window.has_complete_prediction_coverage
    with pytest.raises(SafetyPolicyError):
        _decision(window)


@pytest.mark.parametrize("length", [0, -1, True, 9.5])
def test_request_requires_positive_integer_lengths(length):
    with pytest.raises(ValueError, match="positive integer"):
        replace(REQUESTS[0], peptide_length=length)


def test_duplicate_scope_and_unrequested_length_fail_explicitly():
    with pytest.raises(ValueError, match="Duplicate"):
        _window([], (REQUESTS[0], REQUESTS[0]))
    with pytest.raises(SafetyAssessmentError, match="outside the declared request"):
        _window(_rows(), tuple(r for r in REQUESTS if r.peptide_length == 9))
    request = replace(REQUESTS[0], allele="HLA-A0201")
    assert request.allele == "HLA-A*02:01"


def test_too_long_lengths_do_not_invent_windows_or_make_empty_output_complete():
    requests = (REQUESTS[0], replace(REQUESTS[0], peptide_length=25))
    assert _window(_rows(requests), requests).has_complete_prediction_coverage
    window = _window([], (replace(REQUESTS[0], peptide_length=25),))
    assert window.prediction_coverage_gaps == ()
    assert not window.has_complete_prediction_coverage


def test_repeated_peptide_sequence_does_not_fill_a_missing_occurrence():
    antigen = replace(_antigen(), amino_acids="A" * 20)
    rows = _rows(sequence=antigen.amino_acids)
    rows = [row for row in rows if row["peptide_offset"] != 4]
    window = safety_assessment_from_prediction_frame(
        pd.DataFrame(rows), antigen=antigen, reference_proteome=_Reference(),
        prediction_requests=REQUESTS,
    )
    assert all(gap["missing_window_offsets"] == (4,) for gap in window.prediction_coverage_gaps)
    assert not window.has_complete_prediction_coverage


def test_requested_positions_are_relative_to_the_scanned_window():
    antigen = _antigen()
    window = safety_assessment_from_prediction_frame(
        pd.DataFrame(_rows(sequence=antigen.amino_acids[2:17])),
        antigen=antigen, reference_proteome=_Reference(),
        prediction_requests=REQUESTS, window_start=2, window_end=17,
    )
    assert window.has_complete_prediction_coverage
    assert window.ligands[0].antigen_start_offset == 2


def test_existing_context_audit_request_import_and_archives_remain_supported():
    from vaxrank.context_audit import MHCRequest as ContextRequest

    assert ContextRequest is MHCRequest
    historical = json.loads(to_native_json(REQUESTS[0]))
    historical["__class__"]["__module__"] = "vaxrank.context_audit"
    assert from_native_json(json.dumps(historical), MHCRequest) == REQUESTS[0]
