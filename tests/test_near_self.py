import json

import pytest

from vaxrank.near_self import (
    BLOSUM62,
    DEFAULT_NEAR_SELF_COMPARATOR,
    NEAR_SELF_REASON_NO_RISK_LIGANDS,
    NEAR_SELF_REASON_PREDICTION_COVERAGE,
    NEAR_SELF_REASON_PROTEIN_RESOLUTION,
    NEAR_SELF_REASON_QUERY_OUTSIDE_INDEX,
    AminoAcidSubstitutionMatrix,
    NearSelfAssessment,
    NearSelfError,
    NearSelfQuery,
    assess_near_self,
    assess_near_self_queries,
)
from vaxrank.risk_ligand import (
    PatientHLARiskLigandIndex,
    ProteinResolutionCoverage,
    RiskCoverageGap,
    RiskLigand,
    RiskLigandCoverage,
    RiskLigandIndexProvenance,
    RiskLigandPrediction,
    RiskLigandSelectionPolicy,
    RiskLigandSource,
)
from vaxrank.tissue_risk import TissueRiskEvidence


EVIDENCE = TissueRiskEvidence(
    tissue_group="heart",
    hpa_tissue="heart muscle",
    cell_type="cardiomyocytes",
    level="High",
    reliability="Enhanced",
)


def _source(gene_id, source_name):
    return RiskLigandSource(
        source_sequence_name=source_name,
        protein_sequence_sha256=(gene_id[-1].lower() * 64),
        gene_id=gene_id,
        gene_name=f"GENE_{gene_id[-1]}",
        transcript_ids=(f"T_{gene_id[-1]}",),
        protein_ids=(f"P_{gene_id[-1]}",),
        protein_offsets=(4,),
        tissue_evidence=(EVIDENCE,),
    )


def _ligand(peptide, allele, gene_id, source_name):
    return RiskLigand(
        peptide=peptide,
        allele=allele,
        predictions=(RiskLigandPrediction(
            kind="pMHC_affinity",
            predictor_name="netMHCpan",
            predictor_version="4.2",
            allele=allele,
            score=0.99,
            value=20.0,
            percentile_rank=0.1,
        ),),
        inclusion_kinds=("pMHC_affinity",),
        sources=(_source(gene_id, source_name),),
    )


def _index(*, gaps=(), unresolved=(), invalid_sequences=0):
    ligands = (
        _ligand("AAAAAAAAL", "HLA-A*02:01", "ENSG1", "risk_1"),
        _ligand("AAAAAAAAV", "HLA-A*02:01", "ENSG2", "risk_2"),
        _ligand("CCCCCCCCC", "HLA-A*02:01", "ENSG3", "risk_3"),
        _ligand("AAAAAAAAI", "HLA-B*07:02", "ENSG4", "risk_4"),
        _ligand("AAAAAAAAII", "HLA-A*02:01", "ENSG5", "risk_5"),
    )
    resolved_count = 5
    return PatientHLARiskLigandIndex(
        alleles=("HLA-A*02:01", "HLA-B*07:02"),
        peptide_lengths=(9, 10),
        selection_policy=RiskLigandSelectionPolicy(),
        ligands=ligands,
        protein_resolution_coverage=ProteinResolutionCoverage(
            requested_gene_count=resolved_count + len(unresolved),
            resolved_gene_count=resolved_count,
            unresolved_gene_ids=tuple(unresolved),
            protein_coding_transcript_count=resolved_count + invalid_sequences,
            invalid_sequence_transcript_count=invalid_sequences,
            distinct_sequence_count=resolved_count,
        ),
        coverage=RiskLigandCoverage(
            prediction_row_count=100,
            configured_prediction_row_count=100,
            passing_prediction_row_count=5,
            retained_ligand_count=5,
            missing_combinations=tuple(gaps),
        ),
        provenance=RiskLigandIndexProvenance(
            genome_release="110",
            species="Homo sapiens",
            protein_resolution_policy="all_protein_coding_isoforms",
            tissue_risk_policy_sha256="a" * 64,
            hpa_source_sha256="b" * 64,
            cta_unfiltered_gene_ids_sha256="c" * 64,
            predictor_identities=("pMHC_affinity|netMHCpan|4.2",),
            cache_identity_sha256="d" * 64,
        ),
    )


def test_blosum62_known_scores_and_provenance_are_stable():
    assert BLOSUM62.score("L", "I") == 2
    assert BLOSUM62.score("I", "V") == 3
    assert BLOSUM62.score("D", "E") == 2
    assert BLOSUM62.score("W", "W") == 11
    assert BLOSUM62.score("W", "C") == -2
    assert BLOSUM62.score("I", "L") == BLOSUM62.score("L", "I")
    assert BLOSUM62.source_url == "https://doi.org/10.1073/pnas.89.22.10915"
    assert BLOSUM62.sha256 == (
        "3a565180d86e4df38549766d42f2c4e20153083c7c1a638609ecefffcd801641"
    )


def test_nearest_search_is_same_allele_equal_length_and_retains_all_ties():
    result = assess_near_self("AAAAAAAAI", "A*02:01", _index())

    assert result.normalized_allele == "HLA-A*02:01"
    assert result.nearest_distance == 1
    assert [hit.risk_ligand.peptide for hit in result.nearest_hits] == [
        "AAAAAAAAL", "AAAAAAAAV"
    ]
    assert [hit.substitutions[0].matrix_score for hit in result.nearest_hits] == [
        2, 3
    ]
    assert all(hit.substitutions[0].position == 8 for hit in result.nearest_hits)
    assert {
        hit.risk_ligand.sources[0].gene_id for hit in result.nearest_hits
    } == {"ENSG1", "ENSG2"}
    assert result.reason_codes == ()
    assert result.has_complete_coverage

    payload = result.to_report_dict()
    assert json.loads(json.dumps(payload))["nearest_distance"] == 1
    assert NearSelfAssessment.from_json(result.to_json()) == result


def test_exact_same_allele_risk_ligand_has_distance_zero():
    result = assess_near_self("AAAAAAAAI", "HLA-B*07:02", _index())
    assert result.nearest_distance == 0
    assert len(result.nearest_hits) == 1
    assert result.nearest_hits[0].substitutions == ()
    assert result.nearest_hits[0].substitution_score_sum == 0


def test_batch_results_preserve_query_identity_and_source_coordinates():
    queries = (
        NearSelfQuery(
            target_id="candidate-1",
            peptide="AAAAAAAAI",
            allele="A*02:01",
            source_name="TP53",
            source_offset=12,
        ),
        NearSelfQuery(
            target_id="candidate-2",
            peptide="CCCCCCCCC",
            allele="A*02:01",
            source_name="KRAS",
            source_offset=3,
        ),
    )
    results = assess_near_self_queries(queries, _index())
    assert tuple(result.query for result in results) == queries
    assert results[0].query.source_offset == 12
    assert results[1].nearest_distance == 0


def test_prediction_and_protein_resolution_gaps_remain_visible_with_hits():
    gap = RiskCoverageGap(
        allele="HLA-A*02:01",
        peptide_length=9,
        kind="pMHC_presentation",
        expected_window_count=100,
        observed_window_count=99,
    )
    result = assess_near_self(
        "AAAAAAAAI",
        "HLA-A*02:01",
        _index(gaps=(gap,), unresolved=("ENSG_MISSING",), invalid_sequences=1),
    )
    assert result.nearest_hits
    assert result.prediction_coverage_gaps == (gap,)
    assert set(result.reason_codes) == {
        NEAR_SELF_REASON_PREDICTION_COVERAGE,
        NEAR_SELF_REASON_PROTEIN_RESOLUTION,
    }
    assert not result.has_complete_coverage


def test_query_outside_index_is_not_reported_as_safe_empty_result():
    result = assess_near_self("AAAAAAAA", "HLA-A*02:01", _index())
    assert result.nearest_distance is None
    assert result.nearest_hits == ()
    assert set(result.reason_codes) == {
        NEAR_SELF_REASON_NO_RISK_LIGANDS,
        NEAR_SELF_REASON_QUERY_OUTSIDE_INDEX,
    }
    assert not result.has_complete_coverage


def test_configured_group_with_no_risk_ligands_is_explicit_but_complete():
    # There are no B*07:02 10-mers in the synthetic retained risk set, but the
    # index declares that group configured and has no coverage gaps.
    result = assess_near_self("AAAAAAAAAA", "HLA-B*07:02", _index())
    assert result.reason_codes == (NEAR_SELF_REASON_NO_RISK_LIGANDS,)
    assert result.has_complete_coverage


@pytest.mark.parametrize("peptide", ["", "AAAAAAAAX"])
def test_invalid_target_peptides_fail(peptide):
    with pytest.raises((NearSelfError, ValueError), match="canonical"):
        assess_near_self(peptide, "HLA-A*02:01", _index())


def test_invalid_allele_fails_through_mhcgnomes():
    with pytest.raises(NearSelfError, match="Invalid MHC allele"):
        assess_near_self("AAAAAAAAI", "not-an-allele", _index())


def test_custom_matrix_must_be_complete_and_symmetric():
    rows = list(BLOSUM62.rows)
    changed = list(rows[0])
    changed[1] += 1
    rows[0] = tuple(changed)
    with pytest.raises(ValueError, match="symmetric"):
        AminoAcidSubstitutionMatrix(
            name="broken",
            version="1",
            alphabet=BLOSUM62.alphabet,
            rows=tuple(rows),
            citation="test",
            source_url="https://example.test",
        )


def test_comparator_failure_is_not_an_empty_success():
    class BrokenComparator:
        provenance = DEFAULT_NEAR_SELF_COMPARATOR.provenance

        def distances(self, target, risk_peptides):
            raise RuntimeError("broken")

    with pytest.raises(NearSelfError, match="comparison failed"):
        assess_near_self(
            "AAAAAAAAI", "HLA-A*02:01", _index(), comparator=BrokenComparator()
        )
