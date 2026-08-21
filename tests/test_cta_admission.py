import json

import pandas as pd
import pytest

from vaxrank.cta_admission import (
    CTA_REASON_CANONICAL_EXPRESSION_PASS,
    CTA_REASON_EXPLICIT_OVERRIDE,
    CTA_REASON_EXPRESSION_BELOW_THRESHOLD,
    CTA_REASON_NONCANONICAL_HELD_OUT,
    CTAAdmissionAssessment,
    CTAAdmissionError,
    CTAAdmissionPolicy,
    CTAOverrideEvidence,
    PatientTumorExpressionEvidence,
    assess_cta_antigen,
)
from vaxrank.vaccine_antigen import (
    ATTESTATION_ADMITTED,
    ATTESTATION_HELD_OUT,
    ATTESTATION_OVERRIDDEN,
    TumorSpecificityEvidence,
)


PRAME = "ENSG00000185686"
HELD_OUT = "ENSG00000279804"
OTHER_CTA = "ENSG00000199999"


def _frame():
    rows = []
    for gene_id, symbol, canonical in (
        (PRAME, "PRAME", True),
        (HELD_OUT, "PRAMEF18", False),
        (OTHER_CTA, "OTHERCTA", False),
    ):
        rows.append({
            "Ensembl_Gene_ID": gene_id,
            "Symbol": symbol,
            "passes_filters": True,
            "never_expressed": not canonical,
            "specificity_status": (
                "canonical_default" if canonical else "canonical_low_expression"
            ),
            "specificity_action": (
                "include_default" if canonical else "exclude_default"
            ),
            "specificity_source_anchor": "derived:oncoref.cta.passes_filters",
            "specificity_rationale": "fixture",
            "restriction": "TESTIS",
            "restriction_confidence": "HIGH",
            "safety_flags": "",
            "source_databases": "CTpedia",
            "Canonical_Transcript_ID": f"ENST_{symbol}",
        })
    return pd.DataFrame(rows)


def _patch_oncoref(monkeypatch, *, frame=None, canonical=None, unfiltered=None):
    monkeypatch.setattr(
        "oncoref.cta.cta_evidence", lambda: _frame() if frame is None else frame
    )
    monkeypatch.setattr(
        "oncoref.cta.cta_gene_ids", lambda: {PRAME} if canonical is None else canonical
    )
    monkeypatch.setattr(
        "oncoref.cta.cta_unfiltered_gene_ids",
        lambda: {PRAME, HELD_OUT, OTHER_CTA}
        if unfiltered is None else unfiltered,
    )


def _expression(gene_id=PRAME, value=5.0, unit="TPM"):
    return PatientTumorExpressionEvidence(
        gene_id=gene_id,
        sample_id="patient-tumor-1",
        value=value,
        unit=unit,
        evidence_source="RNA-seq workflow",
        evidence_version="1.2.3",
        assay="bulk RNA-seq",
    )


def _policy(minimum=5.0, unit="TPM"):
    return CTAAdmissionPolicy(
        min_tumor_expression=minimum,
        expression_unit=unit,
    )


def _assess(gene_id=PRAME, expression=None, policy=None, override=None):
    return assess_cta_antigen(
        amino_acids="ACDEFGHIKLMN",
        gene_id=gene_id,
        tumor_expression=expression or _expression(gene_id),
        policy=policy or _policy(),
        override_evidence=override,
    )


def test_canonical_cta_at_expression_threshold_is_admitted(monkeypatch):
    _patch_oncoref(monkeypatch)
    result = _assess()

    assert result.antigen.tumor_specificity.status == ATTESTATION_ADMITTED
    assert result.antigen.tumor_specificity.rationale_code == (
        CTA_REASON_CANONICAL_EXPRESSION_PASS
    )
    assert result.antigen.tumor_specificity.admits_construct
    assert result.antigen.targetable_mask.overlaps(0, 1)
    assert result.antigen.gene_name == "PRAME"
    assert result.antigen.transcript_ids == ("ENST_PRAME",)
    assert result.antigen.self_reference_excluded_gene_ids == tuple(sorted({
        PRAME, HELD_OUT, OTHER_CTA
    }))
    assert result.reference_evidence.canonical_default
    assert result.reference_evidence.oncoref_version == "1.8.174"
    records = result.antigen.tumor_specificity.evidence_records
    assert [record.evidence_kind for record in records] == [
        "oncoref_cta_membership", "patient_tumor_expression"
    ]
    assert records[1].numeric_value == records[1].threshold == 5.0

    payload = result.to_report_dict()
    assert json.loads(json.dumps(payload))["antigen"]["gene_name"] == "PRAME"
    assert CTAAdmissionAssessment.from_json(result.to_json()) == result


def test_expression_below_threshold_is_held_out_even_with_override(monkeypatch):
    _patch_oncoref(monkeypatch)
    result = _assess(
        expression=_expression(value=4.99),
        override=CTAOverrideEvidence("review", "2026-08", "clinical rationale"),
    )
    attestation = result.antigen.tumor_specificity
    assert attestation.status == ATTESTATION_HELD_OUT
    assert attestation.rationale_code == CTA_REASON_EXPRESSION_BELOW_THRESHOLD
    assert not attestation.admits_construct
    assert all(
        record.evidence_kind != "cta_admission_override"
        for record in attestation.evidence_records
    )


def test_noncanonical_cta_requires_explicit_versioned_override(monkeypatch):
    _patch_oncoref(monkeypatch)
    held = _assess(gene_id=HELD_OUT)
    assert held.antigen.tumor_specificity.status == ATTESTATION_HELD_OUT
    assert held.antigen.tumor_specificity.rationale_code == (
        CTA_REASON_NONCANONICAL_HELD_OUT
    )

    override = CTAOverrideEvidence(
        evidence_source="tumor board",
        evidence_version="protocol-7",
        rationale="orthogonal tumor-specificity evidence",
        evidence_id="review-123",
    )
    admitted = _assess(gene_id=HELD_OUT, override=override)
    attestation = admitted.antigen.tumor_specificity
    assert attestation.status == ATTESTATION_OVERRIDDEN
    assert attestation.rationale_code == CTA_REASON_EXPLICIT_OVERRIDE
    assert attestation.override_reason == override.rationale
    assert attestation.requires_review
    assert attestation.evidence_records[-1].evidence_kind == (
        "cta_admission_override"
    )


def test_gene_outside_unfiltered_oncoref_universe_is_not_a_cta(monkeypatch):
    _patch_oncoref(monkeypatch)
    with pytest.raises(CTAAdmissionError, match="not in oncoref"):
        _assess(gene_id="ENSG_NOT_CTA")


@pytest.mark.parametrize(
    "expression, policy, message",
    [
        (_expression(gene_id=HELD_OUT), _policy(), "gene IDs differ"),
        (_expression(unit="nTPM"), _policy(unit="TPM"), "units do not match"),
    ],
)
def test_expression_must_match_cta_gene_and_policy_units(
    monkeypatch, expression, policy, message
):
    _patch_oncoref(monkeypatch)
    with pytest.raises(ValueError, match=message):
        _assess(expression=expression, policy=policy)


@pytest.mark.parametrize("value", [-1, float("nan"), float("inf")])
def test_invalid_quantitative_expression_fails(value):
    with pytest.raises(ValueError, match="finite and non-negative"):
        _expression(value=value)


def test_expression_policy_rejects_zero_threshold():
    with pytest.raises(ValueError, match="greater than zero"):
        _policy(minimum=0)


def test_oncoref_missing_ambiguous_or_untyped_evidence_fails_closed(monkeypatch):
    _patch_oncoref(monkeypatch, canonical=set(), unfiltered=set())
    with pytest.raises(CTAAdmissionError, match="empty CTA reference"):
        _assess()

    duplicate = pd.concat([_frame(), _frame().iloc[[0]]], ignore_index=True)
    _patch_oncoref(monkeypatch, frame=duplicate)
    with pytest.raises(CTAAdmissionError, match="Expected one"):
        _assess()

    malformed = _frame().astype({"passes_filters": object})
    malformed.loc[0, "passes_filters"] = "unknown"
    _patch_oncoref(monkeypatch, frame=malformed)
    with pytest.raises(CTAAdmissionError, match="not boolean"):
        _assess()


def test_generic_attestation_evidence_rejects_unversioned_quantitation():
    with pytest.raises(ValueError, match="requires units"):
        TumorSpecificityEvidence(
            evidence_kind="expression",
            evidence_source="fixture",
            numeric_value=1.0,
        )


def test_live_oncoref_prame_regression_uses_two_distinct_cta_sets():
    from oncoref.cta import cta_gene_ids, cta_unfiltered_gene_ids

    result = _assess()
    assert PRAME in cta_gene_ids()
    assert PRAME in cta_unfiltered_gene_ids()
    assert len(cta_unfiltered_gene_ids()) > len(cta_gene_ids())
    assert result.antigen.tumor_specificity.status == ATTESTATION_ADMITTED
    assert set(result.antigen.self_reference_excluded_gene_ids) == {
        gene_id.split(".")[0] for gene_id in cta_unfiltered_gene_ids()
    }
