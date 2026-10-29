"""Exercise released reference data alongside the stubbed CTA unit tests."""

import oncoref
from packaging.version import Version
import pytest
from oncoref import cta
from oncoref.version import DATA_VERSION, SOURCE_MATRIX_VERSION
from pyensembl import EnsemblRelease

from vaxrank.cta_admission import (
    CTAAdmissionAssessment,
    CTAAdmissionPolicy,
    CTA_REASON_NONCANONICAL_HELD_OUT,
    PatientTumorExpressionEvidence,
    assess_cta_antigen,
    resolve_cta_reference_evidence,
)
from vaxrank.cta_expression import CTAExpressionResult
from .test_cta_expression import run
from .test_cta_identity import Annotation


def test_public_default_cta_admission_keeps_extended_targets_explicit():
    core = cta.cta_df().set_index("Symbol")
    assert "PRAME" in cta.cta_gene_names()
    assert "GPX5" in cta.cta_extended_gene_names()
    assert "GPX5" not in cta.cta_gene_names()
    genome = EnsemblRelease(93)
    excluded = tuple(sorted(oncoref.cta_annotation_gene_ids(genome, unfiltered=True)))

    for symbol, canonical in (("PRAME", True), ("GPX5", False), ("SPAG9", False)):
        resolution = resolve_cta_reference_evidence(core.loc[symbol, "Ensembl_Gene_ID"], genome=genome)
        evidence = resolution.evidence
        assert evidence.canonical_default is canonical
        assert evidence.oncoref_version == oncoref.__version__
        assert evidence.oncoref_data_version == DATA_VERSION
        assert evidence.oncoref_source_matrix_version == SOURCE_MATRIX_VERSION
        assert resolution.self_reference_excluded_gene_ids == excluded
        assert evidence.gene_id in excluded


@pytest.mark.parametrize("symbol", ["GPX5", "SPAG9"])
def test_public_noncanonical_cta_stays_held_out_and_replays_without_reference_queries(
    monkeypatch, symbol,
):
    # GPX5 is extended-only; SPAG9 is a newly nominated, filtered-out record.
    # Neither source nomination nor extended membership admits a construct.
    row = cta.cta_df().set_index("Symbol").loc[symbol]
    assessment = assess_cta_antigen(
        # Synthetic sequence: this test checks reference admission, not source
        # protein identity or peptide presentation.
        amino_acids="ACDEFGHIKLMNPQRSTVWY",
        gene_id=row.Ensembl_Gene_ID,
        genome=EnsemblRelease(93),
        tumor_expression=PatientTumorExpressionEvidence(
            gene_id=row.Ensembl_Gene_ID, sample_id="synthetic-patient",
            value=10, unit="TPM", evidence_source="synthetic fixture",
            evidence_version="1", assay="synthetic RNA-seq"),
        policy=CTAAdmissionPolicy(2, "TPM"),
    )
    assert assessment.antigen.tumor_specificity.rationale_code == CTA_REASON_NONCANONICAL_HELD_OUT
    assert not assessment.antigen.tumor_specificity.admits_construct
    saved = assessment.to_json()

    def forbidden(*args, **kwargs):
        raise AssertionError("Saved CTA admission must not re-resolve a reference")

    monkeypatch.setattr(cta, "cta_evidence", forbidden)
    monkeypatch.setattr(cta, "cta_gene_ids", forbidden)
    monkeypatch.setattr(cta, "cta_unfiltered_gene_ids", forbidden)
    monkeypatch.setattr(oncoref, "__version__", "future-reference")
    assert CTAAdmissionAssessment.from_json(saved) == assessment


def test_identical_rbmy_transcripts_retain_independent_expression_and_source_identity(
    tmp_path, monkeypatch,
):
    genes = ("ENSG00000169800", "ENSG00000226941")
    # The corrected public protein group is not a gene alias or evidence that
    # independently measured patient features may be pooled by the consumer.
    assert {oncoref.proteoform_key(gid, scope="cta") for gid in genes} == {"RBMY1F/J"}
    genome = Annotation(112, fixture_name="rbmy-ensembl-112.json")
    transcripts = ("ENST00000303766", "ENST00000250831")
    result = run(tmp_path, genome, list(zip(transcripts, (3, 7))), level="transcript")
    assert [decision.status for decision in result.decisions] == ["admitted", "admitted"]
    antigens = result.admitted_antigens
    assert antigens[0].amino_acids == antigens[1].amino_acids
    assert len(antigens[0].amino_acids) == 496
    assert len({antigen.source_identifier for antigen in antigens}) == 2
    assert result.annotation_reference.self_reference_excluded_gene_ids == tuple(sorted(genes))
    for decision, gene, transcript, value in zip(result.decisions, genes, transcripts, (3, 7)):
        assert decision.gene_identity["source_gene_id"] == gene
        assert decision.gene_identity["canonical_gene_id"] == gene
        assert decision.assessment.antigen.gene_id == gene
        assert decision.assessment.antigen.transcript_ids == (transcript,)
        assert decision.assessment.tumor_expression.transcript_id == transcript
        assert decision.assessment.tumor_expression.value == value
        evidence = decision.assessment.reference_evidence
        assert evidence.gene_id == gene
        assert evidence.oncoref_version == oncoref.__version__
        assert evidence.oncoref_data_version == DATA_VERSION
        assert evidence.oncoref_source_matrix_version == SOURCE_MATRIX_VERSION

    path = tmp_path / "admission.json"
    result.save(path)

    def forbidden(*args, **kwargs):
        raise AssertionError("Saved gene-specific CTA evidence queried a current reference")

    monkeypatch.setattr(oncoref, "resolve_gene_identity", forbidden)
    monkeypatch.setattr(oncoref, "cta_annotation_gene_identities", forbidden)
    monkeypatch.setattr(oncoref, "proteoform_key", forbidden)
    monkeypatch.setattr(cta, "cta_evidence", forbidden)
    monkeypatch.setattr(cta, "cta_unfiltered_gene_ids", forbidden)
    monkeypatch.setattr(oncoref, "__version__", "future-reference")
    monkeypatch.setattr("oncoref.version.DATA_VERSION", "future-data")
    assert CTAExpressionResult.load(path) == result


def test_released_species_cta_snapshots_keep_evidence_and_offline_identity(tmp_path):
    from oncoref.species_antigens import CTAEvidenceSnapshot

    for collection, taxon in ((oncoref.human, 9606), (oncoref.canine, 9615),
                               (oncoref.mouse, 10090)):
        snapshot = collection.cta_snapshot()
        assert snapshot.taxon_id == taxon
        path = tmp_path / f'{taxon}.json'
        snapshot.write(path)
        replay = CTAEvidenceSnapshot.read(path, expected_snapshot_id=snapshot.snapshot_id)
        assert replay.as_dict() == snapshot.as_dict()
        assert replay.candidates() == snapshot.candidates()
        assert replay.evidence() == snapshot.evidence()

    dog = oncoref.canine.cta_snapshot()
    assert dog.coverage['status'] == 'partial'
    assert dog.gene_ids(reference='refseq106')
    assert any(c['mapping']['status'] == 'ambiguous' for c in dog.candidates())
    assert any(e['tissue'] == 'kidney' for e in dog.evidence())
    mouse = oncoref.mouse.cta_snapshot()
    if Version(oncoref.__version__) < Version('1.22.1'):
        assert mouse.coverage['status'] == 'unavailable'
    else:
        assert mouse.coverage['status'] == 'partial'
        assert {candidate['id'] for candidate in mouse.candidates()} == {'P1A'}
        assert any(row['tissue_role'] == 'somatic' for row in mouse.evidence())
        assert all(row['sample_count'] is None and row['donor_count'] is None
                   for row in mouse.evidence())
