import json
from types import SimpleNamespace

import pandas as pd
import pytest
from topiary import TopiaryPredictor

from vaxrank.risk_ligand import (
    PROTEIN_RESOLUTION_ALL_ISOFORMS,
    PROTEIN_RESOLUTION_LONGEST_ISOFORM,
    RISK_COMBINATION_INTERSECTION,
    PatientHLARiskLigandIndex,
    RiskLigandError,
    RiskLigandSelectionPolicy,
    build_patient_hla_risk_ligand_index,
    resolve_tissue_risk_protein_sequences,
    risk_ligand_index_from_prediction_frame,
)
from vaxrank.tissue_risk import (
    HPAEvidenceCoverage,
    HPAReferenceProvenance,
    SafetyTissueGroup,
    TissueRiskAssessment,
    TissueRiskEvidence,
    TissueRiskPolicy,
    TissueRiskProtein,
)


EVIDENCE = TissueRiskEvidence(
    tissue_group="lung",
    hpa_tissue="lung",
    cell_type="pneumocytes",
    level="High",
    reliability="Enhanced",
)


def _tissue_risk():
    return TissueRiskAssessment(
        policy=TissueRiskPolicy(
            tissue_groups=(SafetyTissueGroup("lung", ("lung",)),),
        ),
        proteins=(
            TissueRiskProtein("ENSG1", "G1", (EVIDENCE,)),
            TissueRiskProtein("ENSG2", "G2", (EVIDENCE,)),
            TissueRiskProtein("ENSG3", "G3", (EVIDENCE,)),
        ),
        coverage=HPAEvidenceCoverage(
            total_rows=3,
            selected_tissue_rows=3,
            missing_gene_id_rows=0,
            missing_level_rows=0,
            missing_reliability_rows=0,
            nonqualifying_level_rows=0,
            nonqualifying_reliability_rows=0,
            qualifying_rows=3,
            cta_source_rows_excluded=0,
        ),
        hpa_provenance=HPAReferenceProvenance(
            oncoref_version="1.8.174",
            oncoref_data_version="5.23.17",
            oncoref_source_matrix_version="5.22.10",
            source_name="hpa_normal_tissue",
            hpa_version="v23",
            source_path="/data/normal_tissue.tsv",
            source_bytes=100,
            source_sha256="a" * 64,
            verification_status="verified",
        ),
        cta_unfiltered_gene_count=2,
        cta_unfiltered_gene_ids_sha256="b" * 64,
    )


class _Transcript:
    def __init__(self, transcript_id, sequence, protein_id=None, coding=True):
        self.transcript_id = transcript_id
        self.id = transcript_id
        self.protein_sequence = sequence
        self.protein_id = protein_id or f"P_{transcript_id}"
        self.is_protein_coding = coding


class _Genome:
    release = 110
    species = SimpleNamespace(latin_name="Homo sapiens")

    def __init__(self):
        self.genes = {
            "ENSG1": SimpleNamespace(transcripts=(
                _Transcript("T1", "ACDEFGHIKLMN"),
                _Transcript("T1_DUP", "ACDEFGHIKLMN"),
                _Transcript("T2", "ACDEFGHIKLAA"),
                _Transcript("T_BAD", "ACDEFGHIKLMX"),
                _Transcript("T_NONCODING", "ACDEFGHIKLMN", coding=False),
            )),
            "ENSG2": SimpleNamespace(transcripts=(
                _Transcript("T3", "TTTACDEFGHIKVV"),
            )),
        }

    def gene_by_id(self, gene_id):
        return self.genes[gene_id]


def _protein_sequences(policy=PROTEIN_RESOLUTION_ALL_ISOFORMS):
    return resolve_tissue_risk_protein_sequences(
        _tissue_risk(), _Genome(), resolution_policy=policy
    )


def _row(source, peptide, offset, *, kind="pMHC_affinity",
         allele="HLA-A*02:01", percentile=0.5, predictor="mhcflurry"):
    return {
        "source_sequence_name": source,
        "peptide": peptide,
        "peptide_length": len(peptide),
        "peptide_offset": offset,
        "kind": kind,
        "allele": allele,
        "percentile_rank": percentile,
        "score": 1.0 - percentile / 100.0,
        "affinity": 50.0 if kind == "pMHC_affinity" else float("nan"),
        "value": 50.0 if kind == "pMHC_affinity" else float("nan"),
        "prediction_method_name": predictor,
        "predictor_version": "2.1.4",
    }


def test_resolve_all_isoforms_deduplicates_sequences_and_exposes_gaps():
    result = _protein_sequences()
    assert result.resolution_policy == PROTEIN_RESOLUTION_ALL_ISOFORMS
    assert result.genome_release == "110"
    assert result.species == "Homo sapiens"
    assert result.coverage.requested_gene_count == 3
    assert result.coverage.resolved_gene_count == 2
    assert result.coverage.unresolved_gene_ids == ("ENSG3",)
    assert result.coverage.protein_coding_transcript_count == 5
    assert result.coverage.invalid_sequence_transcript_count == 1
    assert result.coverage.distinct_sequence_count == 3
    g1_shared = next(
        protein for protein in result.proteins
        if protein.gene_id == "ENSG1" and protein.amino_acids == "ACDEFGHIKLMN"
    )
    assert g1_shared.transcript_ids == ("T1", "T1_DUP")


def test_longest_isoform_policy_is_explicit_performance_tradeoff():
    result = _protein_sequences(PROTEIN_RESOLUTION_LONGEST_ISOFORM)
    assert len(result.proteins) == 2
    assert {protein.gene_id for protein in result.proteins} == {"ENSG1", "ENSG2"}


def _prediction_frame(proteins):
    names = proteins.named_sequences()
    rows = []
    for source, sequence in names.items():
        for offset in range(len(sequence) - 9 + 1):
            peptide = sequence[offset:offset + 9]
            affinity_rank = 0.5 if peptide == "ACDEFGHIK" else 50.0
            presentation_rank = 50.0
            if peptide == "ACDEFGHIK":
                presentation_rank = 2.0
            elif peptide == "CDEFGHIKL":
                affinity_rank = 5.0
                presentation_rank = 1.0
            rows.extend([
                _row(source, peptide, offset, percentile=affinity_rank),
                _row(
                    source,
                    peptide,
                    offset,
                    kind="pMHC_presentation",
                    percentile=presentation_rank,
                    predictor="mhcflurry-presentation",
                ),
            ])
    return pd.DataFrame(rows)


def test_default_union_selects_top_one_percent_and_preserves_shared_sources():
    proteins = _protein_sequences()
    frame = _prediction_frame(proteins)
    result = risk_ligand_index_from_prediction_frame(
        frame,
        protein_sequences=proteins,
        alleles=("A*02:01",),
        peptide_lengths=(9,),
    )

    assert result.alleles == ("HLA-A*02:01",)
    assert [ligand.peptide for ligand in result.ligands] == [
        "ACDEFGHIK", "CDEFGHIKL"
    ]
    shared = result.ligands[0]
    assert shared.inclusion_kinds == ("pMHC_affinity",)
    assert {source.gene_id for source in shared.sources} == {"ENSG1", "ENSG2"}
    assert len(shared.predictions) == 2
    presentation_only = result.ligands[1]
    assert presentation_only.inclusion_kinds == ("pMHC_presentation",)
    assert result.coverage.missing_combinations == ()
    assert result.coverage.retained_ligand_count == 2
    assert result.for_allele_and_length("A*02:01", 9) == result.ligands

    payload = result.to_report_dict()
    assert json.loads(json.dumps(payload)) == payload
    assert PatientHLARiskLigandIndex.from_json(result.to_json()) == result


def test_intersection_requires_every_configured_prediction_kind():
    proteins = _protein_sequences()
    policy = RiskLigandSelectionPolicy(combination=RISK_COMBINATION_INTERSECTION)
    result = risk_ligand_index_from_prediction_frame(
        _prediction_frame(proteins),
        protein_sequences=proteins,
        alleles=("HLA-A*02:01",),
        peptide_lengths=(9,),
        selection_policy=policy,
    )
    assert result.ligands == ()


def test_missing_model_coverage_is_explicit_and_cache_is_deterministic():
    proteins = _protein_sequences()
    frame = _prediction_frame(proteins)
    affinity_only = frame[frame["kind"] == "pMHC_affinity"]
    first = risk_ligand_index_from_prediction_frame(
        affinity_only,
        protein_sequences=proteins,
        alleles=("HLA-A*02:01", "HLA-B*07:02"),
        peptide_lengths=(9, 10),
    )
    assert len(first.coverage.missing_combinations) == 7
    assert {
        (gap.allele, gap.peptide_length, gap.kind)
        for gap in first.coverage.missing_combinations
    } >= {
        ("HLA-A*02:01", 9, "pMHC_presentation"),
        ("HLA-B*07:02", 9, "pMHC_affinity"),
        ("HLA-A*02:01", 10, "pMHC_affinity"),
    }
    second = risk_ligand_index_from_prediction_frame(
        affinity_only.iloc[::-1],
        protein_sequences=proteins,
        alleles=("HLA-B*07:02", "A*02:01"),
        peptide_lengths=(10, 9),
    )
    assert (
        first.provenance.cache_identity_sha256
        == second.provenance.cache_identity_sha256
    )


def test_partial_model_output_is_reported_as_a_coverage_gap():
    proteins = _protein_sequences()
    frame = _prediction_frame(proteins)
    partial = frame.drop(frame.index[0])
    result = risk_ligand_index_from_prediction_frame(
        partial,
        protein_sequences=proteins,
        alleles=("HLA-A*02:01",),
        peptide_lengths=(9,),
    )
    assert len(result.coverage.missing_combinations) == 1
    gap = result.coverage.missing_combinations[0]
    assert gap.kind == "pMHC_affinity"
    assert gap.expected_window_count == gap.observed_window_count + 1


def test_missing_percentile_rank_is_unusable_model_coverage():
    proteins = _protein_sequences()
    frame = _prediction_frame(proteins)
    unique_affinity_row = frame.index[
        (frame["peptide"] == "DEFGHIKVV")
        & (frame["kind"] == "pMHC_affinity")
    ][0]
    frame.loc[unique_affinity_row, "percentile_rank"] = float("nan")
    result = risk_ligand_index_from_prediction_frame(
        frame,
        protein_sequences=proteins,
        alleles=("HLA-A*02:01",),
        peptide_lengths=(9,),
    )
    assert len(result.coverage.missing_combinations) == 1
    assert result.coverage.missing_combinations[0].missing_window_count == 1


@pytest.mark.parametrize(
    "mutator, message",
    [
        (lambda row: row.update(source_sequence_name="WRONG"), "source"),
        (lambda row: row.update(peptide_offset=99), "source protein"),
        (lambda row: row.update(peptide_length=8), "length"),
        (lambda row: row.update(allele="HLA-B*07:02"), "unconfigured allele"),
        (lambda row: row.update(percentile_rank=-1), "percentile rank"),
        (lambda row: row.update(predictor_version=""), "predictor name, and version"),
    ],
)
def test_malformed_or_unexpected_prediction_rows_fail_closed(mutator, message):
    proteins = _protein_sequences()
    row = _prediction_frame(proteins).iloc[0].to_dict()
    mutator(row)
    with pytest.raises(RiskLigandError, match=message):
        risk_ligand_index_from_prediction_frame(
            pd.DataFrame([row]),
            protein_sequences=proteins,
            alleles=("HLA-A*02:01",),
            peptide_lengths=(9,),
        )


def test_conflicting_duplicate_model_prediction_fails_closed():
    proteins = _protein_sequences()
    row = _prediction_frame(proteins).iloc[0].to_dict()
    conflicting = dict(row, percentile_rank=0.9)
    with pytest.raises(RiskLigandError, match="Conflicting"):
        risk_ligand_index_from_prediction_frame(
            pd.DataFrame([row, conflicting]),
            protein_sequences=proteins,
            alleles=("HLA-A*02:01",),
            peptide_lengths=(9,),
        )


def test_high_level_predictor_error_is_not_empty_success():
    class BrokenPredictor(TopiaryPredictor):
        def __init__(self):
            pass

        def predict_from_named_sequences(self, sequences):
            raise RuntimeError("unavailable")

    with pytest.raises(RiskLigandError, match="prediction failed"):
        build_patient_hla_risk_ligand_index(
            BrokenPredictor(),
            _protein_sequences(),
            alleles=("HLA-A*02:01",),
            peptide_lengths=(9,),
        )


def test_invalid_allele_uses_mhcgnomes_and_fails_closed():
    with pytest.raises(RiskLigandError, match="Invalid MHC allele"):
        risk_ligand_index_from_prediction_frame(
            pd.DataFrame(),
            protein_sequences=_protein_sequences(),
            alleles=("not-an-allele",),
            peptide_lengths=(9,),
        )
