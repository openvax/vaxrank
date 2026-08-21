import json

import pytest

from vaxrank import (
    AminoAcidInterval,
    AntigenSafetyAssessment,
    EmittedSafetyLigand,
    PatientHLARiskLigandIndex,
    ProteinResolutionCoverage,
    RiskLigand,
    RiskLigandCoverage,
    RiskLigandIndexProvenance,
    RiskLigandPrediction,
    RiskLigandSelectionPolicy,
    RiskLigandSource,
    SafetyPrediction,
    SafetyPredictionCoverage,
    TargetableMask,
    TissueRiskEvidence,
    TumorSpecificityAttestation,
    VaccineAntigen,
    WindowSafetyAssessment,
    assess_antigen_safety,
    near_self_queries_for_window,
)
from vaxrank.near_self import NEAR_SELF_REASON_QUERY_OUTSIDE_INDEX
from vaxrank.vaccine_antigen import (
    ANTIGEN_KIND_MUTATION,
    ATTESTATION_ADMITTED,
    SelfReferenceMatch,
)


def test_public_antigen_safety_composes_complete_window_and_near_self_evidence():
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
    affinity_a = SafetyPrediction(
        kind="pMHC_affinity",
        predictor_name="netMHCpan",
        predictor_version="4.2",
        allele="HLA-A*02:01",
        percentile_rank=0.1,
    )
    affinity_b = SafetyPrediction(
        kind="pMHC_affinity",
        predictor_name="netMHCpan",
        predictor_version="4.2",
        allele="HLA-B*07:02",
        percentile_rank=0.2,
    )
    presentation_a = SafetyPrediction(
        kind="pMHC_presentation",
        predictor_name="MHCflurry-presentation",
        predictor_version="2.1.4",
        allele="HLA-A*02:01",
        percentile_rank=0.3,
    )
    window = WindowSafetyAssessment(
        antigen=antigen,
        window_start_offset=0,
        window_end_offset=len(antigen.amino_acids),
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
                predictions=(affinity_a,),
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
                predictions=(affinity_a, affinity_b, presentation_a),
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
                predictions=(affinity_a,),
            ),
        ),
        coverage=(
            SafetyPredictionCoverage(
                kind="pMHC_affinity",
                predictor_name="netMHCpan",
                predictor_version="4.2",
                alleles=("HLA-A*02:01", "HLA-B*07:02"),
                peptide_lengths=(9,),
                prediction_count=4,
            ),
            SafetyPredictionCoverage(
                kind="pMHC_presentation",
                predictor_name="MHCflurry-presentation",
                predictor_version="2.1.4",
                alleles=("HLA-A*02:01",),
                peptide_lengths=(9,),
                prediction_count=1,
            ),
        ),
        genome_release="110",
    )

    tissue_evidence = TissueRiskEvidence(
        tissue_group="heart",
        hpa_tissue="heart muscle",
        cell_type="cardiomyocytes",
        level="High",
        reliability="Enhanced",
    )
    risk_source = RiskLigandSource(
        source_sequence_name="risk_000001_ENSG000001",
        protein_sequence_sha256="a" * 64,
        gene_id="ENSG000001",
        gene_name="RISK1",
        transcript_ids=("ENST000001",),
        protein_ids=("ENSP000001",),
        protein_offsets=(4,),
        tissue_evidence=(tissue_evidence,),
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
    risk_index = PatientHLARiskLigandIndex(
        alleles=("HLA-A*02:01",),
        peptide_lengths=(9,),
        selection_policy=RiskLigandSelectionPolicy(),
        ligands=(risk_ligand,),
        protein_resolution_coverage=ProteinResolutionCoverage(
            requested_gene_count=1,
            resolved_gene_count=1,
            unresolved_gene_ids=(),
            protein_coding_transcript_count=1,
            invalid_sequence_transcript_count=0,
            distinct_sequence_count=1,
        ),
        coverage=RiskLigandCoverage(
            prediction_row_count=1,
            configured_prediction_row_count=1,
            passing_prediction_row_count=1,
            retained_ligand_count=1,
        ),
        provenance=RiskLigandIndexProvenance(
            genome_release="110",
            species="Homo sapiens",
            protein_resolution_policy="all_protein_coding_isoforms",
            tissue_risk_policy_sha256="b" * 64,
            hpa_source_sha256="c" * 64,
            cta_unfiltered_gene_ids_sha256="d" * 64,
            predictor_identities=("pMHC_affinity|netMHCpan|4.2",),
            cache_identity_sha256="e" * 64,
        ),
    )

    queries = near_self_queries_for_window(window)
    assert len(window.target_ligands) == 1
    assert [query.allele for query in queries] == [
        "HLA-A*02:01",
        "HLA-B*07:02",
    ]
    assert {query.source_offset for query in queries} == {9}
    assert len({query.target_id for query in queries}) == 1

    result = assess_antigen_safety(window, risk_index)
    by_allele = {
        assessment.normalized_allele: assessment
        for assessment in result.near_self_assessments
    }
    assert by_allele["HLA-A*02:01"].nearest_distance == 1
    assert by_allele["HLA-A*02:01"].nearest_hits[0].substitutions[0].matrix_score == 2
    assert (
        NEAR_SELF_REASON_QUERY_OUTSIDE_INDEX
        in by_allele["HLA-B*07:02"].reason_codes
    )
    assert not result.has_complete_near_self_coverage
    assert result.risk_index_provenance == risk_index.provenance

    rows = result.near_self_rows()
    assert len(rows) == 2
    risk_row = next(row for row in rows if row["risk_gene_id"])
    assert risk_row["risk_gene_id"] == "ENSG000001"
    assert risk_row["tissue_group"] == "heart"
    assert risk_row["substitutions"] == "8:I>L:2"
    missing_row = next(row for row in rows if not row["risk_gene_id"])
    assert missing_row["risk_peptide"] == ""
    assert NEAR_SELF_REASON_QUERY_OUTSIDE_INDEX in missing_row["reason_codes"]
    assert all(
        not isinstance(value, (dict, list, tuple))
        for row in rows
        for value in row.values()
    )

    payload = result.to_report_dict()
    assert json.loads(json.dumps(payload)) == payload
    assert AntigenSafetyAssessment.from_json(result.to_json()) == result
    with pytest.raises(ValueError, match="every target ligand/allele"):
        AntigenSafetyAssessment(
            window_assessment=window,
            risk_index_provenance=risk_index.provenance,
            near_self_assessments=result.near_self_assessments[:1],
        )


def test_repeated_target_peptide_positions_remain_distinct_public_queries():
    antigen = VaccineAntigen(
        kind=ANTIGEN_KIND_MUTATION,
        amino_acids="AAAAAAAAI" * 2,
        targetable_mask=TargetableMask((
            AminoAcidInterval(8, 9),
            AminoAcidInterval(17, 18),
        )),
        tumor_specificity=TumorSpecificityAttestation(
            status=ATTESTATION_ADMITTED,
            evidence_kind="somatic_variant",
            evidence_source="test",
            patient_specific=True,
        ),
        source_identifier="repeated-target",
    )
    prediction = SafetyPrediction(
        kind="pMHC_affinity",
        predictor_name="netMHCpan",
        predictor_version="4.2",
        allele="HLA-A*02:01",
        percentile_rank=0.1,
    )
    self_match = SelfReferenceMatch(
        peptide="AAAAAAAAI",
        occurs=False,
        antigen_kind=ANTIGEN_KIND_MUTATION,
    )
    window = WindowSafetyAssessment(
        antigen=antigen,
        window_start_offset=0,
        window_end_offset=18,
        ligands=(
            EmittedSafetyLigand(
                peptide="AAAAAAAAI",
                window_start_offset=0,
                window_end_offset=9,
                antigen_start_offset=0,
                antigen_end_offset=9,
                overlaps_targetable=True,
                self_reference_match=self_match,
                predictions=(prediction,),
            ),
            EmittedSafetyLigand(
                peptide="AAAAAAAAI",
                window_start_offset=9,
                window_end_offset=18,
                antigen_start_offset=9,
                antigen_end_offset=18,
                overlaps_targetable=True,
                self_reference_match=self_match,
                predictions=(prediction,),
            ),
        ),
        coverage=(SafetyPredictionCoverage(
            kind="pMHC_affinity",
            predictor_name="netMHCpan",
            predictor_version="4.2",
            alleles=("HLA-A*02:01",),
            peptide_lengths=(9,),
            prediction_count=2,
        ),),
    )

    queries = near_self_queries_for_window(window)

    assert [query.source_offset for query in queries] == [0, 9]
    assert len({query.target_id for query in queries}) == 2


@pytest.mark.parametrize(
    "missing_field",
    ("kind", "predictor_name", "predictor_version", "allele"),
)
def test_safety_prediction_requires_complete_model_and_allele_provenance(
    missing_field,
):
    values = {
        "kind": "pMHC_affinity",
        "predictor_name": "netMHCpan",
        "predictor_version": "4.2",
        "allele": "HLA-A*02:01",
    }
    values[missing_field] = ""

    with pytest.raises(ValueError, match="predictor, version, and allele"):
        SafetyPrediction(**values)


def test_safety_prediction_normalizes_nonhuman_alleles_through_mhcgnomes():
    prediction = SafetyPrediction(
        kind="pMHC_affinity",
        predictor_name="netMHCpan",
        predictor_version="4.2",
        allele="H-2-Kb",
    )
    coverage = SafetyPredictionCoverage(
        kind="pMHC_affinity",
        predictor_name="netMHCpan",
        predictor_version="4.2",
        alleles=("H-2-Kb",),
        peptide_lengths=(8,),
        prediction_count=1,
    )

    assert prediction.allele == "H2-K*b"
    assert coverage.alleles == ("H2-K*b",)
