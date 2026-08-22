import json

import pandas as pd
import pytest
from topiary import TopiaryPredictor

from vaxrank.safety_assessment import (
    SAFETY_REASON_NO_PREDICTIONS,
    ConstructBoundary,
    EmittedSafetyLigand,
    SafetyAssessmentError,
    SafetyPrediction,
    SafetyPredictionCoverage,
    WindowSafetyAssessment,
    assess_vaccine_antigen_window,
    safety_assessment_from_prediction_frame,
)
from vaxrank.vaccine_antigen import (
    ANTIGEN_KIND_MUTATION,
    ATTESTATION_ADMITTED,
    AminoAcidInterval,
    TargetableMask,
    TumorSpecificityAttestation,
    VaccineAntigen,
)


def _antigen(excluded_gene_ids=()):
    return VaccineAntigen(
        kind=ANTIGEN_KIND_MUTATION,
        amino_acids="ACDEFGHIKLMNPQRSTVWY",
        targetable_mask=TargetableMask((AminoAcidInterval(8, 10),)),
        tumor_specificity=TumorSpecificityAttestation(
            status=ATTESTATION_ADMITTED,
            evidence_kind="somatic_variant",
            evidence_source="test",
            patient_specific=True,
            rationale_code="test_admission",
        ),
        self_reference_excluded_gene_ids=excluded_gene_ids,
        gene_name="GENE",
    )


def _row(peptide, offset, *, predictor="mhcflurry", kind="pMHC_affinity",
         allele="HLA-A*02:01", score=0.01, value=40000.0,
         percentile_rank=99.0):
    return {
        "peptide": peptide,
        "peptide_length": len(peptide),
        "peptide_offset": offset,
        "allele": allele,
        "affinity": value if kind == "pMHC_affinity" else float("nan"),
        "value": value if kind == "pMHC_affinity" else float("nan"),
        "score": score,
        "percentile_rank": percentile_rank,
        "prediction_method_name": predictor,
        "predictor_version": "2.1.4",
        "source_sequence_name": "GENE",
        "kind": kind,
    }


def test_safety_prediction_from_prediction_row_is_public_and_normalizes_values():
    row = _row(
        "ACDEFGHIK",
        0,
        predictor="mhcflurry-presentation",
        kind="pMHC_presentation",
        allele="H-2-Kb",
        score="0.42",
        value=float("nan"),
        percentile_rank="1.7",
    )

    prediction = SafetyPrediction.from_prediction_row(row)

    assert prediction.kind == "pMHC_presentation"
    assert prediction.allele == "H2-K*b"
    assert prediction.score == 0.42
    assert prediction.value is None
    assert prediction.percentile_rank == 1.7


def test_safety_prediction_from_prediction_row_rejects_non_numeric_values():
    row = _row("ACDEFGHIK", 0, score="strong")

    with pytest.raises(ValueError, match="not numeric"):
        SafetyPrediction.from_prediction_row(row)


class _Reference:
    def __init__(self, peptides=(), excluded_gene_ids=()):
        self.peptides = set(peptides)
        self.excluded_gene_ids = frozenset(excluded_gene_ids)

    def contains(self, peptide):
        return peptide in self.peptides


def test_complete_inventory_retains_weak_non_target_and_all_prediction_kinds():
    antigen = _antigen()
    peptide = antigen.amino_acids[:8]
    rows = [
        _row(peptide, 0, predictor="mhcflurry", score=0.00001),
        _row(
            peptide,
            0,
            predictor="mhcflurry-presentation",
            kind="pMHC_presentation",
            score=0.00002,
            value=None,
        ),
        _row(peptide, 0, predictor="netmhcpan", score=0.00003),
    ]
    result = safety_assessment_from_prediction_frame(
        pd.DataFrame(rows),
        antigen=antigen,
        reference_proteome=_Reference({peptide}),
    )

    assert len(result.ligands) == 1
    ligand = result.ligands[0]
    assert ligand.is_non_target
    assert ligand.occurs_in_self_reference
    assert len(ligand.predictions) == 3
    assert {prediction.kind for prediction in ligand.predictions} == {
        "pMHC_affinity", "pMHC_presentation"
    }
    assert ligand.predictions[-1].value is None
    assert {item.predictor_name for item in result.coverage} == {
        "mhcflurry", "mhcflurry-presentation", "netmhcpan"
    }


def test_window_coordinates_target_mask_and_boundary_crossing_are_exact():
    antigen = _antigen()
    # Scan antigen [4, 16); peptide local [3, 9) => antigen [7, 13),
    # overlapping target [8, 10) and crossing only the boundary at 8.
    window = antigen.amino_acids[4:16]
    peptide = window[3:9]
    result = safety_assessment_from_prediction_frame(
        pd.DataFrame([_row(peptide, 3)]),
        antigen=antigen,
        reference_proteome=_Reference(),
        window_start=4,
        window_end=16,
        construct_boundaries=(
            ConstructBoundary(3, "LEFT", "PEPTIDE"),
            ConstructBoundary(8, "PEPTIDE", "RIGHT"),
            ConstructBoundary(9, "RIGHT", "TAIL"),
        ),
        genome_release="110",
    )

    ligand = result.ligands[0]
    assert ligand.window_start_offset == 3
    assert ligand.window_end_offset == 9
    assert ligand.antigen_start_offset == 7
    assert ligand.antigen_end_offset == 13
    assert ligand.overlaps_targetable
    assert [boundary.offset for boundary in ligand.crossed_construct_boundaries] == [8]
    assert ligand.crosses_construct_boundary
    assert not ligand.crosses_window_boundary
    assert not ligand.crosses_antigen_boundary
    assert ligand.self_reference_match.genome_release == "110"


def test_empty_predictor_output_is_explicit_not_reassuring_success():
    result = safety_assessment_from_prediction_frame(
        pd.DataFrame(),
        antigen=_antigen(),
        reference_proteome=_Reference(),
    )
    assert result.ligands == ()
    assert result.coverage == ()
    assert result.reason_codes == (SAFETY_REASON_NO_PREDICTIONS,)


@pytest.mark.parametrize(
    "mutator, message",
    [
        (lambda row: row.update(peptide_length=8), "length"),
        (lambda row: row.update(peptide_offset=99), "outside"),
        (lambda row: row.update(peptide="AAAAAAAAA"), "does not match"),
        (lambda row: row.update(peptide_offset=0.5), "integer"),
    ],
)
def test_malformed_prediction_coordinates_fail_closed(mutator, message):
    antigen = _antigen()
    row = _row(antigen.amino_acids[:9], 0)
    mutator(row)
    with pytest.raises(SafetyAssessmentError, match=message):
        safety_assessment_from_prediction_frame(
            pd.DataFrame([row]),
            antigen=antigen,
            reference_proteome=_Reference(),
        )


def test_duplicate_model_allele_leaf_fails_instead_of_silent_collapse():
    antigen = _antigen()
    row = _row(antigen.amino_acids[:9], 0)
    duplicate = dict(row, score=0.9)
    with pytest.raises(SafetyAssessmentError, match="Duplicate"):
        safety_assessment_from_prediction_frame(
            pd.DataFrame([row, duplicate]),
            antigen=antigen,
            reference_proteome=_Reference(),
        )


def test_nonfinite_prediction_numbers_are_explicit_missing_values():
    antigen = _antigen()
    row = _row(antigen.amino_acids[:9], 0)
    row.update(score=float("nan"), value=float("inf"), affinity=float("inf"))
    result = safety_assessment_from_prediction_frame(
        pd.DataFrame([row]),
        antigen=antigen,
        reference_proteome=_Reference(),
    )
    prediction = result.ligands[0].predictions[0]
    assert prediction.score is None
    assert prediction.value is None


def test_high_level_api_never_swallows_predictor_errors():
    class BrokenPredictor(TopiaryPredictor):
        def __init__(self):
            pass

        def predict_from_named_sequences(self, _sequences):
            raise RuntimeError("model unavailable")

    with pytest.raises(SafetyAssessmentError, match="MHC safety prediction failed"):
        assess_vaccine_antigen_window(
            BrokenPredictor(),
            _antigen(),
            reference_proteome=_Reference(),
        )


def test_high_level_api_checks_reference_exclusions_and_serializes():
    antigen = _antigen(("ENSG000001.2",))

    class StubPredictor(TopiaryPredictor):
        def __init__(self):
            pass

        def predict_from_named_sequences(self, sequences):
            sequence = next(iter(sequences.values()))
            return pd.DataFrame([_row(sequence[:9], 0)])

    with pytest.raises(ValueError, match="exclusions do not match"):
        assess_vaccine_antigen_window(
            StubPredictor(),
            antigen,
            reference_proteome=_Reference(),
        )

    result = assess_vaccine_antigen_window(
        StubPredictor(),
        antigen,
        reference_proteome=_Reference(excluded_gene_ids={"ENSG000001"}),
    )
    payload = result.to_report_dict()
    assert json.loads(json.dumps(payload)) == payload
    assert payload["ligands"][0]["peptide"] == antigen.amino_acids[:9]
    assert payload["coverage"][0]["prediction_count"] == 1
    row = result.prediction_rows()[0]
    assert row["peptide"] == antigen.amino_acids[:9]
    assert row["predictor_name"] == "mhcflurry"
    assert all(not isinstance(value, (dict, list, tuple)) for value in row.values())
    assert WindowSafetyAssessment.from_json(result.to_json()) == result


def test_construct_boundary_must_be_internal_to_window():
    antigen = _antigen()
    with pytest.raises(ValueError, match="outside"):
        safety_assessment_from_prediction_frame(
            pd.DataFrame(),
            antigen=antigen,
            reference_proteome=_Reference(),
            construct_boundaries=(ConstructBoundary(len(antigen.amino_acids)),),
        )


def test_high_level_api_rejects_prediction_rows_from_another_source():
    class WrongSourcePredictor(TopiaryPredictor):
        def __init__(self):
            pass

        def predict_from_named_sequences(self, sequences):
            sequence = next(iter(sequences.values()))
            return pd.DataFrame([
                dict(_row(sequence[:9], 0), source_sequence_name="OTHER")
            ])

    with pytest.raises(SafetyAssessmentError, match="source"):
        assess_vaccine_antigen_window(
            WrongSourcePredictor(),
            _antigen(),
            reference_proteome=_Reference(),
        )


def test_result_rejects_coverage_that_omits_prediction_evidence():
    antigen = _antigen()
    peptide = antigen.amino_acids[:8]
    prediction = SafetyPrediction(
        kind="pMHC_affinity",
        predictor_name="mhcflurry",
        predictor_version="2.1.4",
        allele="HLA-A*02:01",
        score=0.1,
    )
    ligand = EmittedSafetyLigand(
        peptide=peptide,
        window_start_offset=0,
        window_end_offset=8,
        antigen_start_offset=0,
        antigen_end_offset=8,
        overlaps_targetable=False,
        self_reference_match=antigen.self_reference_match(peptide, False),
        predictions=(prediction,),
    )
    with pytest.raises(ValueError, match="omits"):
        WindowSafetyAssessment(
            antigen=antigen,
            window_start_offset=0,
            window_end_offset=len(antigen.amino_acids),
            ligands=(ligand,),
            coverage=(),
        )

    with pytest.raises(ValueError, match="disagrees"):
        WindowSafetyAssessment(
            antigen=antigen,
            window_start_offset=0,
            window_end_offset=len(antigen.amino_acids),
            ligands=(ligand,),
            coverage=(SafetyPredictionCoverage(
                kind="pMHC_affinity",
                predictor_name="mhcflurry",
                predictor_version="2.1.4",
                alleles=("HLA-B*07:02",),
                peptide_lengths=(8,),
                prediction_count=1,
            ),),
        )
