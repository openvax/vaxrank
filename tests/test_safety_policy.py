import json
from dataclasses import replace

import pytest

from vaxrank import (
    DEFAULT_NEAR_SELF_COMPARATOR,
    SAFETY_MODE_AUDIT,
    SAFETY_MODE_ENFORCE,
    SAFETY_MODE_OFF,
    AminoAcidInterval,
    AntigenSafetyAssessment,
    AntigenSafetyPolicyDecision,
    EmittedSafetyLigand,
    RiskLigand,
    RiskLigandIndexProvenance,
    RiskLigandPrediction,
    RiskLigandSource,
    SafetyEnforcementPolicy,
    SafetyPrediction,
    SafetyPredictionCoverage,
    SafetyPolicyError,
    TargetableMask,
    TissueRiskEvidence,
    TumorSpecificityAttestation,
    VaccineAntigen,
    WindowSafetyAssessment,
    evaluate_antigen_safety_policy,
    near_self_queries_for_window,
)
from vaxrank.near_self import (
    NEAR_SELF_REASON_PROTEIN_RESOLUTION,
    NearSelfAssessment,
    NearSelfQuery,
    NearestRiskLigandHit,
)
from vaxrank.safety_assessment import SAFETY_REASON_NO_PREDICTIONS
from vaxrank.safety_policy import (
    SAFETY_ACTION_ALLOW,
    SAFETY_ACTION_ERROR,
    SAFETY_ACTION_REPAIR,
    SAFETY_FINDING_DISABLED,
    SAFETY_FINDING_INCOMPLETE,
    SAFETY_FINDING_UNSAFE,
    SAFETY_MODEL_COMBINATION_ALL,
    SAFETY_REASON_DISABLED,
    SAFETY_REASON_EXACT_SELF,
    SAFETY_REASON_INCOMPLETE_NEAR_SELF,
    SAFETY_REASON_MISSING_PREDICTION_KIND,
    SAFETY_REASON_NEAR_SELF_CONSERVATIVE,
    SAFETY_REASON_NEAR_SELF_DISTANCE,
    SAFETY_REASON_NON_TARGET,
)
from vaxrank.vaccine_antigen import (
    ANTIGEN_KIND_MUTATION,
    ATTESTATION_ADMITTED,
    SelfReferenceMatch,
)


@pytest.mark.parametrize(
    "target_peptide,risk_peptide,expected_reason",
    (
        ("AAAAAAAAI", "AAAAAAAAI", SAFETY_REASON_NEAR_SELF_DISTANCE),
        ("AAAAAAAAI", "AAAAAAAAL", SAFETY_REASON_NEAR_SELF_DISTANCE),
        ("AAAAAAAII", "AAAAAAALL", SAFETY_REASON_NEAR_SELF_CONSERVATIVE),
        ("AAAAAAADW", "AAAAAAAWC", None),
        ("AAAAAAVWY", "AAAAAACCA", None),
    ),
)
def test_default_near_self_gate_covers_distance_and_conservative_cases(
    target_peptide,
    risk_peptide,
    expected_reason,
):
    tissue_evidence = TissueRiskEvidence(
        tissue_group="heart",
        hpa_tissue="heart muscle",
        cell_type="cardiomyocytes",
        level="High",
        reliability="Enhanced",
    )
    risk_source = RiskLigandSource(
        source_sequence_name="risk-protein",
        protein_sequence_sha256="a" * 64,
        gene_id="ENSG000001",
        gene_name="RISK1",
        transcript_ids=("ENST000001",),
        protein_ids=("ENSP000001",),
        protein_offsets=(4,),
        tissue_evidence=(tissue_evidence,),
    )
    risk_ligand = RiskLigand(
        peptide=risk_peptide,
        allele="HLA-A*02:01",
        predictions=(RiskLigandPrediction(
            kind="pMHC_affinity",
            predictor_name="netMHCpan",
            predictor_version="4.2",
            allele="HLA-A*02:01",
            score=0.99,
            value=20.0,
            percentile_rank=0.1,
        ),),
        inclusion_kinds=("pMHC_affinity",),
        sources=(risk_source,),
    )
    substitutions = DEFAULT_NEAR_SELF_COMPARATOR.substitutions(
        target_peptide,
        risk_peptide,
    )
    near_self = NearSelfAssessment(
        query=NearSelfQuery(
            peptide=target_peptide,
            allele="HLA-A*02:01",
        ),
        normalized_allele="HLA-A*02:01",
        nearest_distance=len(substitutions),
        nearest_hits=(NearestRiskLigandHit(
            risk_ligand=risk_ligand,
            hamming_distance=len(substitutions),
            substitutions=substitutions,
            substitution_score_sum=sum(
                substitution.matrix_score for substitution in substitutions
            ),
        ),),
        prediction_coverage_gaps=(),
        reason_codes=(),
        comparator=DEFAULT_NEAR_SELF_COMPARATOR.provenance,
        risk_index_cache_identity_sha256="e" * 64,
    )

    reasons = SafetyEnforcementPolicy().near_self_risk_reason_codes(near_self)

    assert reasons == ((expected_reason,) if expected_reason else ())
    if target_peptide == "AAAAAAADW":
        override = SafetyEnforcementPolicy(
            additional_conservative_pairs=(("W", "C"),),
        )
        assert override.near_self_risk_reason_codes(near_self) == (
            SAFETY_REASON_NEAR_SELF_CONSERVATIVE,
        )


def test_prediction_gate_supports_union_intersection_and_missing_coverage():
    affinity = SafetyPrediction(
        kind="pMHC_affinity",
        predictor_name="netMHCpan",
        predictor_version="4.2",
        allele="HLA-A*02:01",
        percentile_rank=0.5,
    )
    presentation = SafetyPrediction(
        kind="pMHC_presentation",
        predictor_name="MHCflurry-presentation",
        predictor_version="2.1.4",
        allele="HLA-A*02:01",
        percentile_rank=2.0,
    )
    union_policy = SafetyEnforcementPolicy()
    intersection_policy = SafetyEnforcementPolicy(
        model_combination=SAFETY_MODEL_COMBINATION_ALL,
    )

    assert union_policy.prediction_gate_is_met((affinity, presentation))
    assert not intersection_policy.prediction_gate_is_met((affinity, presentation))
    assert union_policy.threshold_predictions((affinity, presentation)) == (affinity,)
    assert union_policy.missing_prediction_kinds((affinity, presentation)) == ()
    assert union_policy.missing_prediction_kinds((affinity,)) == (
        "pMHC_presentation",
    )


def test_off_audit_and_enforce_separate_findings_from_actions():
    antigen = VaccineAntigen(
        kind=ANTIGEN_KIND_MUTATION,
        amino_acids="CCCCCCCCCAAAAAAAAIGGGGGGGGG",
        targetable_mask=TargetableMask((
            AminoAcidInterval(17, 18),
            AminoAcidInterval(26, 27),
        )),
        tumor_specificity=TumorSpecificityAttestation(
            status=ATTESTATION_ADMITTED,
            evidence_kind="somatic_variant",
            evidence_source="test",
            patient_specific=True,
        ),
        gene_name="TP53",
        source_identifier="TP53:p.R175H",
    )
    affinity = SafetyPrediction(
        kind="pMHC_affinity",
        predictor_name="netMHCpan",
        predictor_version="4.2",
        allele="HLA-A*02:01",
        percentile_rank=0.5,
    )
    presentation = SafetyPrediction(
        kind="pMHC_presentation",
        predictor_name="MHCflurry-presentation",
        predictor_version="2.1.4",
        allele="HLA-A*02:01",
        percentile_rank=0.5,
    )
    window = WindowSafetyAssessment(
        antigen=antigen,
        window_start_offset=0,
        window_end_offset=27,
        ligands=(
            EmittedSafetyLigand(
                peptide="CCCCCCCCC",
                window_start_offset=0,
                window_end_offset=9,
                antigen_start_offset=0,
                antigen_end_offset=9,
                overlaps_targetable=False,
                self_reference_match=SelfReferenceMatch(
                    peptide="CCCCCCCCC",
                    occurs=False,
                    antigen_kind=ANTIGEN_KIND_MUTATION,
                ),
                predictions=(affinity, presentation),
            ),
            EmittedSafetyLigand(
                peptide="AAAAAAAAI",
                window_start_offset=9,
                window_end_offset=18,
                antigen_start_offset=9,
                antigen_end_offset=18,
                overlaps_targetable=True,
                self_reference_match=SelfReferenceMatch(
                    peptide="AAAAAAAAI",
                    occurs=False,
                    antigen_kind=ANTIGEN_KIND_MUTATION,
                ),
                predictions=(affinity, presentation),
            ),
            EmittedSafetyLigand(
                peptide="GGGGGGGGG",
                window_start_offset=18,
                window_end_offset=27,
                antigen_start_offset=18,
                antigen_end_offset=27,
                overlaps_targetable=True,
                self_reference_match=SelfReferenceMatch(
                    peptide="GGGGGGGGG",
                    occurs=True,
                    antigen_kind=ANTIGEN_KIND_MUTATION,
                ),
                predictions=(affinity, presentation),
            ),
        ),
        coverage=(
            SafetyPredictionCoverage(
                kind="pMHC_affinity",
                predictor_name="netMHCpan",
                predictor_version="4.2",
                alleles=("HLA-A*02:01",),
                peptide_lengths=(9,),
                prediction_count=3,
            ),
            SafetyPredictionCoverage(
                kind="pMHC_presentation",
                predictor_name="MHCflurry-presentation",
                predictor_version="2.1.4",
                alleles=("HLA-A*02:01",),
                peptide_lengths=(9,),
                prediction_count=3,
            ),
        ),
    )
    risk_source = RiskLigandSource(
        source_sequence_name="risk-protein",
        protein_sequence_sha256="a" * 64,
        gene_id="ENSG000001",
        gene_name="RISK1",
        transcript_ids=("ENST000001",),
        protein_ids=("ENSP000001",),
        protein_offsets=(4,),
        tissue_evidence=(TissueRiskEvidence(
            tissue_group="heart",
            hpa_tissue="heart muscle",
            cell_type="cardiomyocytes",
            level="High",
            reliability="Enhanced",
        ),),
    )
    risk_ligand = RiskLigand(
        peptide="AAAAAAAAL",
        allele="HLA-A*02:01",
        predictions=(RiskLigandPrediction(
            kind="pMHC_affinity",
            predictor_name="netMHCpan",
            predictor_version="4.2",
            allele="HLA-A*02:01",
            score=0.99,
            value=20.0,
            percentile_rank=0.1,
        ),),
        inclusion_kinds=("pMHC_affinity",),
        sources=(risk_source,),
    )
    substitutions = DEFAULT_NEAR_SELF_COMPARATOR.substitutions(
        "AAAAAAAAI",
        "AAAAAAAAL",
    )
    near_self = NearSelfAssessment(
        query=near_self_queries_for_window(window)[0],
        normalized_allele="HLA-A*02:01",
        nearest_distance=1,
        nearest_hits=(NearestRiskLigandHit(
            risk_ligand=risk_ligand,
            hamming_distance=1,
            substitutions=substitutions,
            substitution_score_sum=2,
        ),),
        prediction_coverage_gaps=(),
        reason_codes=(),
        comparator=DEFAULT_NEAR_SELF_COMPARATOR.provenance,
        risk_index_cache_identity_sha256="e" * 64,
    )
    provenance = RiskLigandIndexProvenance(
        genome_release="110",
        species="Homo sapiens",
        protein_resolution_policy="all_protein_coding_isoforms",
        tissue_risk_policy_sha256="b" * 64,
        hpa_source_sha256="c" * 64,
        cta_unfiltered_gene_ids_sha256="d" * 64,
        predictor_identities=("pMHC_affinity|netMHCpan|4.2",),
        cache_identity_sha256="e" * 64,
    )
    assessment = AntigenSafetyAssessment(
        window_assessment=window,
        risk_index_provenance=provenance,
        near_self_assessments=(near_self,),
    )

    off = evaluate_antigen_safety_policy(
        assessment,
        SafetyEnforcementPolicy(mode=SAFETY_MODE_OFF),
    )
    audit = evaluate_antigen_safety_policy(
        assessment,
        SafetyEnforcementPolicy(mode=SAFETY_MODE_AUDIT),
    )
    enforce = evaluate_antigen_safety_policy(
        assessment,
        SafetyEnforcementPolicy(mode=SAFETY_MODE_ENFORCE),
    )

    assert off.finding == SAFETY_FINDING_DISABLED
    assert off.action == SAFETY_ACTION_ALLOW
    assert off.reason_codes == (SAFETY_REASON_DISABLED,)
    assert not off.has_complete_coverage
    assert audit.finding == SAFETY_FINDING_UNSAFE
    assert audit.action == SAFETY_ACTION_ALLOW
    assert audit.would_require_repair
    assert audit.has_complete_coverage
    assert enforce.finding == SAFETY_FINDING_UNSAFE
    assert enforce.action == SAFETY_ACTION_REPAIR
    assert enforce.would_require_repair
    assert set(enforce.reason_codes) == {
        SAFETY_REASON_NON_TARGET,
        SAFETY_REASON_EXACT_SELF,
        SAFETY_REASON_NEAR_SELF_DISTANCE,
    }
    assert json.loads(json.dumps(enforce.to_report_dict())) == enforce.to_report_dict()
    assert AntigenSafetyPolicyDecision.from_json(enforce.to_json()) == enforce
    assert SafetyEnforcementPolicy().mode == SAFETY_MODE_ENFORCE

    target_ligand = window.ligands[1]
    incomplete_ligand = replace(target_ligand, predictions=(affinity,))
    incomplete_window = WindowSafetyAssessment(
        antigen=antigen,
        window_start_offset=0,
        window_end_offset=27,
        ligands=(incomplete_ligand,),
        coverage=(SafetyPredictionCoverage(
            kind="pMHC_affinity",
            predictor_name="netMHCpan",
            predictor_version="4.2",
            alleles=("HLA-A*02:01",),
            peptide_lengths=(9,),
            prediction_count=1,
        ),),
    )
    incomplete_assessment = AntigenSafetyAssessment(
        window_assessment=incomplete_window,
        risk_index_provenance=provenance,
        near_self_assessments=(replace(
            near_self,
            query=near_self_queries_for_window(incomplete_window)[0],
        ),),
    )
    incomplete_audit = evaluate_antigen_safety_policy(
        incomplete_assessment,
        SafetyEnforcementPolicy(mode=SAFETY_MODE_AUDIT),
    )
    assert incomplete_audit.finding == SAFETY_FINDING_INCOMPLETE
    assert incomplete_audit.action == SAFETY_ACTION_ALLOW
    assert incomplete_audit.would_require_repair
    assert SAFETY_REASON_MISSING_PREDICTION_KIND in incomplete_audit.reason_codes

    with pytest.raises(SafetyPolicyError) as error:
        evaluate_antigen_safety_policy(incomplete_assessment)
    assert error.value.decision.finding == SAFETY_FINDING_INCOMPLETE
    assert error.value.decision.action == SAFETY_ACTION_ERROR
    assert error.value.decision.would_require_repair

    incomplete_near_self = replace(
        near_self,
        reason_codes=(NEAR_SELF_REASON_PROTEIN_RESOLUTION,),
    )
    incomplete_near_self_assessment = AntigenSafetyAssessment(
        window_assessment=window,
        risk_index_provenance=provenance,
        near_self_assessments=(incomplete_near_self,),
    )
    incomplete_near_self_audit = evaluate_antigen_safety_policy(
        incomplete_near_self_assessment,
        SafetyEnforcementPolicy(mode=SAFETY_MODE_AUDIT),
    )
    assert incomplete_near_self_audit.finding == SAFETY_FINDING_INCOMPLETE
    assert (
        SAFETY_REASON_INCOMPLETE_NEAR_SELF
        in incomplete_near_self_audit.reason_codes
    )
    with pytest.raises(SafetyPolicyError):
        evaluate_antigen_safety_policy(incomplete_near_self_assessment)

    empty_window = WindowSafetyAssessment(
        antigen=antigen,
        window_start_offset=0,
        window_end_offset=27,
        ligands=(),
        coverage=(),
        reason_codes=(SAFETY_REASON_NO_PREDICTIONS,),
    )
    empty_assessment = AntigenSafetyAssessment(
        window_assessment=empty_window,
        risk_index_provenance=provenance,
        near_self_assessments=(),
    )
    with pytest.raises(SafetyPolicyError) as empty_error:
        evaluate_antigen_safety_policy(empty_assessment)
    assert empty_error.value.decision.action == SAFETY_ACTION_ERROR
def test_safety_prediction_rejects_out_of_range_percentiles():
    with pytest.raises(ValueError, match="0..100"):
        SafetyPrediction(
            kind="pMHC_affinity",
            predictor_name="netMHCpan",
            predictor_version="4.2",
            allele="HLA-A*02:01",
            percentile_rank=101.0,
        )
