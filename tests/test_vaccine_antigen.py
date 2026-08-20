# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

from types import SimpleNamespace

import pytest

from vaxrank.candidate_epitope import CandidateEpitope
from vaxrank.epitope_logic import predict_epitopes
from vaxrank.vaccine_antigen import (
    ANTIGEN_KIND_CTA,
    ANTIGEN_KIND_MUTATION,
    ATTESTATION_ADMITTED,
    ATTESTATION_HELD_OUT,
    ATTESTATION_OVERRIDDEN,
    AminoAcidInterval,
    SelfReferenceMatch,
    SelfReferenceSource,
    TargetableMask,
    TumorSpecificityAttestation,
    VaccineAntigen,
)
from vaxrank.vaccine_peptide import VaccinePeptide


def admitted_attestation():
    return TumorSpecificityAttestation(
        status=ATTESTATION_ADMITTED,
        evidence_kind="patient_tumor_expression",
        evidence_source="test fixture",
        patient_specific=True,
        rationale_code="test_admission",
    )


def cta_antigen(sequence="ACDEFGHIKL"):
    return VaccineAntigen(
        kind=ANTIGEN_KIND_CTA,
        amino_acids=sequence,
        targetable_mask=TargetableMask((AminoAcidInterval(0, len(sequence)),)),
        tumor_specificity=admitted_attestation(),
        self_reference_excluded_gene_ids=("ENSG00000185686.4",),
        gene_name="PRAME",
        gene_id="ENSG00000185686.4",
    )


def test_targetable_interval_validation_and_deletion_parity():
    with pytest.raises(ValueError, match="non-negative"):
        AminoAcidInterval(-1, 1)
    with pytest.raises(ValueError, match="must not precede"):
        AminoAcidInterval(2, 1)

    deletion = AminoAcidInterval(4, 4)
    assert deletion.overlaps(0, 8)
    assert not deletion.overlaps(4, 8)


def test_targetable_mask_rejects_ambiguous_or_out_of_bounds_intervals():
    with pytest.raises(ValueError, match="sorted"):
        TargetableMask((AminoAcidInterval(5, 6), AminoAcidInterval(1, 2)))
    with pytest.raises(ValueError, match="must not overlap"):
        TargetableMask((AminoAcidInterval(1, 4), AminoAcidInterval(3, 5)))
    with pytest.raises(ValueError, match="extends beyond"):
        VaccineAntigen(
            kind=ANTIGEN_KIND_CTA,
            amino_acids="ACDEF",
            targetable_mask=TargetableMask((AminoAcidInterval(0, 6),)),
            tumor_specificity=admitted_attestation(),
        )
    with pytest.raises(ValueError, match="non-canonical amino acids: J"):
        VaccineAntigen(
            kind=ANTIGEN_KIND_CTA,
            amino_acids="ACDEFGHIJL",
            targetable_mask=TargetableMask((AminoAcidInterval(0, 10),)),
            tumor_specificity=admitted_attestation(),
        )
    with pytest.raises(ValueError, match="requires targetable content"):
        VaccineAntigen(
            kind=ANTIGEN_KIND_CTA,
            amino_acids="ACDEF",
            targetable_mask=TargetableMask(),
            tumor_specificity=admitted_attestation(),
        )


def test_attestation_requires_known_status_evidence_and_override_reason():
    with pytest.raises(ValueError, match="Unknown tumor-specificity status"):
        TumorSpecificityAttestation("unknown", "expression", "fixture")
    with pytest.raises(ValueError, match="kind and source are required"):
        TumorSpecificityAttestation(ATTESTATION_ADMITTED, "", "fixture")
    with pytest.raises(ValueError, match="requires an override reason"):
        TumorSpecificityAttestation(
            ATTESTATION_OVERRIDDEN, "clinical_review", "fixture"
        )


def test_mutation_adapter_preserves_target_span_and_source_ids():
    transcript = SimpleNamespace(
        gene_id="ENSG00000123456.9",
        transcript_id="ENST000001",
        protein_id="ENSP000001",
    )
    fragment = SimpleNamespace(
        amino_acids="ACDEFGHIKL",
        mutant_amino_acid_start_offset=3,
        mutant_amino_acid_end_offset=4,
        supporting_reference_transcripts=[transcript],
        gene_name="GENE1",
        variant="1:100 A>G",
    )

    antigen = VaccineAntigen.from_mutant_protein_fragment(fragment)

    assert antigen.kind == ANTIGEN_KIND_MUTATION
    assert antigen.interval_is_targetable(0, 4)
    assert not antigen.interval_is_targetable(4, 8)
    assert antigen.gene_id == "ENSG00000123456"
    assert antigen.transcript_ids == ("ENST000001",)
    assert antigen.protein_ids == ("ENSP000001",)
    assert antigen.self_reference_excluded_gene_ids == ()
    assert antigen.tumor_specificity.admits_construct


def test_self_reference_match_rejects_excluded_or_inconsistent_sources():
    prame_source = SelfReferenceSource(gene_id="ENSG00000185686")
    with pytest.raises(ValueError, match="cannot include excluded genes"):
        SelfReferenceMatch(
            peptide="SIINFEKL",
            occurs=True,
            antigen_kind=ANTIGEN_KIND_CTA,
            excluded_gene_ids=("ENSG00000185686",),
            sources=(prame_source,),
            source_provenance_complete=True,
        )
    with pytest.raises(ValueError, match="must agree with its source list"):
        SelfReferenceMatch(
            peptide="SIINFEKL",
            occurs=True,
            antigen_kind=ANTIGEN_KIND_CTA,
            source_provenance_complete=True,
        )
    with pytest.raises(ValueError, match="non-canonical amino acids: X"):
        SelfReferenceMatch(
            peptide="SIINFEKX",
            occurs=False,
            antigen_kind=ANTIGEN_KIND_CTA,
        )


def test_vaccine_peptide_uses_targetable_mask_and_antigen_self_result():
    antigen = cta_antigen()
    fragment = SimpleNamespace(amino_acids=antigen.amino_acids)
    target = CandidateEpitope(
        sequence="ACDEFGHI",
        source_sequence=antigen.amino_acids,
        overlaps_targetable=True,
        occurs_in_reference=True,
        self_reference_match=antigen.self_reference_match("ACDEFGHI", False),
    )
    non_target = CandidateEpitope(
        sequence="CDEFGHIK",
        source_sequence=antigen.amino_acids,
        offset=1,
        overlaps_targetable=False,
        occurs_in_reference=False,
        self_reference_match=antigen.self_reference_match("CDEFGHIK", False),
    )
    exact_self = CandidateEpitope(
        sequence="DEFGHIKL",
        source_sequence=antigen.amino_acids,
        offset=2,
        overlaps_targetable=True,
        occurs_in_reference=False,
        self_reference_match=antigen.self_reference_match("DEFGHIKL", True),
    )

    vaccine_peptide = VaccinePeptide(
        mutant_protein_fragment=fragment,
        antigen=antigen,
        epitopes=[target, non_target, exact_self],
    )

    assert vaccine_peptide.target_epitopes == [target]
    assert vaccine_peptide.non_target_epitopes == [non_target]
    assert vaccine_peptide.self_epitopes == [exact_self]


def test_vaccine_peptide_fails_closed_for_held_out_or_mismatched_policy():
    held_out = VaccineAntigen(
        kind=ANTIGEN_KIND_CTA,
        amino_acids="ACDEFGHIKL",
        targetable_mask=TargetableMask((AminoAcidInterval(0, 10),)),
        tumor_specificity=TumorSpecificityAttestation(
            status=ATTESTATION_HELD_OUT,
            evidence_kind="missing_patient_expression",
            evidence_source="test fixture",
        ),
    )
    fragment = SimpleNamespace(amino_acids=held_out.amino_acids)
    with pytest.raises(ValueError, match="held-out antigen"):
        VaccinePeptide(mutant_protein_fragment=fragment, antigen=held_out)

    antigen = cta_antigen()
    wrong_policy = SelfReferenceMatch(
        peptide="ACDEFGHI",
        occurs=False,
        antigen_kind=ANTIGEN_KIND_MUTATION,
    )
    epitope = CandidateEpitope(
        sequence="ACDEFGHI",
        overlaps_targetable=True,
        self_reference_match=wrong_policy,
    )
    with pytest.raises(ValueError, match="does not match antigen kind"):
        VaccinePeptide(
            mutant_protein_fragment=SimpleNamespace(
                amino_acids=antigen.amino_acids
            ),
            antigen=antigen,
            epitopes=[epitope],
        )


def test_prediction_refuses_non_mutation_antigen_until_admission_path_exists():
    antigen = cta_antigen()
    fragment = SimpleNamespace(amino_acids=antigen.amino_acids)

    with pytest.raises(NotImplementedError, match="issue #303"):
        predict_epitopes(
            mhc_predictor=None,
            protein_fragment=fragment,
            antigen=antigen,
        )


def test_vaccine_antigen_json_roundtrip_preserves_policy():
    antigen = cta_antigen()

    restored = VaccineAntigen.from_json(antigen.to_json())

    assert restored == antigen
    assert restored.self_reference_excluded_gene_ids == ("ENSG00000185686",)
