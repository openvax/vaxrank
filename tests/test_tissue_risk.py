import hashlib
import json

import pandas as pd
import pytest

from vaxrank.tissue_risk import (
    DEFAULT_HPA_LEVELS,
    DEFAULT_HPA_RELIABILITIES,
    PRAME_GENE_ID,
    HPAReferenceProvenance,
    TissueRiskDerivationError,
    build_hpa_reference_provenance,
    derive_tissue_risk_proteins,
    resolve_tissue_risk_policy,
)


def _provenance(verification_status="verified"):
    return HPAReferenceProvenance(
        oncoref_version="1.8.174",
        oncoref_data_version="5.23.17",
        oncoref_source_matrix_version="5.22.10",
        source_name="hpa_normal_tissue",
        hpa_version="v23",
        source_path="/data/normal_tissue.tsv",
        source_bytes=100,
        source_sha256="a" * 64,
        verification_status=verification_status,
    )


def _frame():
    return pd.DataFrame([
        # Same non-CTA gene qualifies through Medium and High evidence.
        ["ENSG_SAFE.4", "SAFE", "heart muscle", "cardiomyocytes",
         "Medium", "Approved"],
        ["ENSG_SAFE", "SAFE", "lung", "pneumocytes",
         "High", "Enhanced"],
        # Low is deliberately below the vaxrank risk threshold even though
        # oncoref's broader CTA detection rule treats Low as protein detected.
        ["ENSG_LOW", "LOW", "liver", "hepatocytes", "Low", "Enhanced"],
        # Medium with uncertain antibody evidence does not qualify.
        ["ENSG_UNCERTAIN", "UNCERTAIN", "pancreas", "exocrine cells",
         "Medium", "Uncertain"],
        # Missing evidence remains separately visible.
        ["ENSG_MISSING_LEVEL", "MISSING_LEVEL", "lung", "macrophages",
         None, "Supported"],
        ["ENSG_MISSING_REL", "MISSING_REL", "liver", "hepatocytes",
         "High", None],
        [None, "MISSING_GENE", "heart muscle", "cardiomyocytes",
         "High", "Supported"],
        # PRAME would otherwise qualify but must be removed by source gene.
        [PRAME_GENE_ID, "PRAME", "lung", "pneumocytes",
         "High", "Enhanced"],
        # Same tissue evidence from a non-CTA source remains eligible. This is
        # the source-specific foundation for later shared-peptide retention.
        ["ENSG_SHARED", "SHARED", "lung", "pneumocytes",
         "High", "Enhanced"],
        # A strong row outside selected safety tissues is not considered.
        ["ENSG_TESTIS", "TESTIS", "testis", "germ cells",
         "High", "Enhanced"],
    ], columns=[
        "Gene", "Gene name", "Tissue", "Cell type", "Level", "Reliability"
    ])


def test_default_policy_is_resolved_from_oncoref_not_copied():
    from oncoref.cta_tissues import SAFETY_TISSUE_GROUPS

    policy = resolve_tissue_risk_policy()
    assert {group.name for group in policy.tissue_groups} == set(
        SAFETY_TISSUE_GROUPS
    )
    assert {
        group.name: set(group.tissues) for group in policy.tissue_groups
    } == SAFETY_TISSUE_GROUPS
    assert set(policy.allowed_levels) == set(DEFAULT_HPA_LEVELS)
    assert set(policy.allowed_reliabilities) == set(DEFAULT_HPA_RELIABILITIES)


def test_medium_high_reliability_rule_and_uncertainty_counts(monkeypatch):
    monkeypatch.setattr(
        "oncoref.cta.cta_unfiltered_gene_ids",
        lambda: {PRAME_GENE_ID, "ENSG_OTHER_CTA"},
    )
    result = derive_tissue_risk_proteins(
        hpa_frame=_frame(),
        provenance=_provenance(),
    )

    assert [protein.gene_id for protein in result.proteins] == [
        "ENSG_SAFE", "ENSG_SHARED"
    ]
    safe = result.proteins[0]
    assert {evidence.tissue_group for evidence in safe.evidence} == {
        "heart", "lung"
    }
    assert result.excluded_cta_source_gene_ids == (PRAME_GENE_ID,)
    assert result.cta_unfiltered_gene_count == 2
    assert result.coverage.total_rows == 10
    assert result.coverage.selected_tissue_rows == 9
    assert result.coverage.missing_gene_id_rows == 1
    assert result.coverage.missing_level_rows == 1
    assert result.coverage.missing_reliability_rows == 1
    assert result.coverage.nonqualifying_level_rows == 1
    assert result.coverage.nonqualifying_reliability_rows == 1
    assert result.coverage.qualifying_rows == 4
    assert result.coverage.cta_source_rows_excluded == 1

    payload = result.to_report_dict()
    assert json.loads(json.dumps(payload)) == payload
    assert {row["gene_id"] for row in result.evidence_rows()} == {
        "ENSG_SAFE", "ENSG_SHARED"
    }


def test_prame_is_in_pinned_oncoref_unfiltered_cta_universe():
    from oncoref.cta import cta_unfiltered_gene_ids

    assert PRAME_GENE_ID in cta_unfiltered_gene_ids()


def test_custom_group_names_still_resolve_through_oncoref(monkeypatch):
    monkeypatch.setattr(
        "oncoref.cta.cta_unfiltered_gene_ids", lambda: {"ENSG_CTA"}
    )
    result = derive_tissue_risk_proteins(
        hpa_frame=_frame(),
        provenance=_provenance(),
        tissue_group_names=("heart",),
    )
    assert [protein.gene_id for protein in result.proteins] == ["ENSG_SAFE"]
    assert [group.name for group in result.policy.tissue_groups] == ["heart"]
    with pytest.raises(ValueError, match="Unknown oncoref"):
        resolve_tissue_risk_policy(("not-a-real-group",))


def test_hpa_schema_and_conflicting_gene_symbols_fail_closed(monkeypatch):
    monkeypatch.setattr(
        "oncoref.cta.cta_unfiltered_gene_ids", lambda: {"ENSG_CTA"}
    )
    with pytest.raises(TissueRiskDerivationError, match="missing columns"):
        derive_tissue_risk_proteins(
            hpa_frame=pd.DataFrame({"Gene": ["ENSG1"]}),
            provenance=_provenance(),
        )

    conflicting = pd.DataFrame([
        ["ENSG1", "OLD", "lung", "cells", "High", "Enhanced"],
        ["ENSG1", "NEW", "lung", "cells", "High", "Enhanced"],
    ], columns=[
        "Gene", "Gene name", "Tissue", "Cell type", "Level", "Reliability"
    ])
    with pytest.raises(TissueRiskDerivationError, match="conflicting symbols"):
        derive_tissue_risk_proteins(
            hpa_frame=conflicting,
            provenance=_provenance(),
        )


def test_injected_hpa_frame_cannot_claim_implicit_provenance():
    with pytest.raises(ValueError, match="explicit provenance"):
        derive_tissue_risk_proteins(hpa_frame=_frame())


@pytest.mark.parametrize(
    "status", ["verified", "checksum_mismatch", "manifest_unavailable"]
)
def test_hpa_provenance_preserves_integrity_state(status):
    assert _provenance(status).verification_status == status


def test_build_hpa_provenance_hashes_and_verifies_exact_file(monkeypatch, tmp_path):
    source = tmp_path / "normal_tissue.tsv"
    source.write_bytes(b"exact hpa bytes")
    monkeypatch.setattr(
        "oncoref.reference_data.ensure", lambda name, version: source
    )
    monkeypatch.setattr(
        "oncoref.reference_data.verify", lambda name, version: True
    )

    result = build_hpa_reference_provenance()
    assert result.source_path == str(source)
    assert result.source_bytes == len(b"exact hpa bytes")
    assert result.source_sha256 == hashlib.sha256(b"exact hpa bytes").hexdigest()
    assert result.verification_status == "verified"


def test_build_hpa_provenance_never_calls_unverified_data_verified(
    monkeypatch, tmp_path
):
    from oncoref.reference_data import ReferenceDataError

    source = tmp_path / "normal_tissue.tsv"
    source.write_bytes(b"hpa")
    monkeypatch.setattr(
        "oncoref.reference_data.ensure", lambda name, version: source
    )
    monkeypatch.setattr(
        "oncoref.reference_data.verify", lambda name, version: False
    )
    assert (
        build_hpa_reference_provenance().verification_status
        == "checksum_mismatch"
    )

    def no_manifest(name, version):
        raise ReferenceDataError("no manifest")

    monkeypatch.setattr("oncoref.reference_data.verify", no_manifest)
    assert (
        build_hpa_reference_provenance().verification_status
        == "manifest_unavailable"
    )


def test_empty_cta_universe_is_rejected_as_unsafe_dependency_state(monkeypatch):
    monkeypatch.setattr("oncoref.cta.cta_unfiltered_gene_ids", lambda: set())
    with pytest.raises(ValueError, match="cannot be empty"):
        derive_tissue_risk_proteins(
            hpa_frame=_frame(),
            provenance=_provenance(),
        )
