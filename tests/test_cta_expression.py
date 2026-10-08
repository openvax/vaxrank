"""Expression admission preserves measurements, occurrences and decisions."""

from dataclasses import replace
import hashlib
from types import SimpleNamespace

import pandas as pd
import pytest

from tests.test_cta_admission import _frame, _patch_oncoref, PRAME, HELD_OUT, OTHER_CTA
from vaxrank.cta_admission import CTAAdmissionPolicy, PatientTumorExpressionEvidence, assess_cta_antigen
from vaxrank.cta_expression import (
    CTAExpressionInput, CTAExpressionResult, CTATargetExclusionPolicy,
    admit_cta_expression,
)


MAGEA3 = OTHER_CTA
MAGEA4 = "ENSG00000147381"
NOT_CTA = "ENSG00000111111"
PRAME_TRANSCRIPT = "ENST00000398741"
ALTERNATE_TRANSCRIPT = "ENST00000398742"
MAGEA3_TRANSCRIPT = "ENST00000398743"
MAGEA4_TRANSCRIPT = "ENST00000398744"
HELD_OUT_TRANSCRIPT = "ENST00000398745"
UNKNOWN_TRANSCRIPT = "ENST00000399999"
SEQUENCE = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMN"


@pytest.mark.parametrize('suffix,separator', [('tsv', '\t'), ('csv', ',')])
def test_duplicate_expression_headers_fail_before_measurement_selection(tmp_path, genome, suffix, separator):
    path = tmp_path / ('duplicate.' + suffix)
    path.write_text(separator.join(['feature_id', 'patient_tpm', 'patient_tpm']) + '\n' +
                    separator.join([PRAME, '3', '30']) + '\n')
    with pytest.raises(ValueError, match='unique column names'):
        admit_cta_expression(path, input_contract=contract(),
                             admission_policy=CTAAdmissionPolicy(2, 'TPM'), genome=genome)


def contract(level="gene", **changes):
    return replace(CTAExpressionInput(
        measurement_level=level, id_column="feature_id", value_column="patient_tpm",
        expression_unit="TPM", sample_id="patient-001",
        evidence_source="Salmon", evidence_version="1.10.3", assay="bulk RNA-seq"),
        **changes)


@pytest.fixture
def genome(monkeypatch):
    frame = _frame()
    frame.loc[frame.Ensembl_Gene_ID == MAGEA3, "Symbol"] = "MAGEA3"
    frame.loc[frame.Ensembl_Gene_ID == MAGEA3, "never_expressed"] = False
    row = frame.loc[frame.Ensembl_Gene_ID == MAGEA3].iloc[0].copy()
    row["Ensembl_Gene_ID"], row["Symbol"] = MAGEA4, "MAGEA4"
    frame = pd.concat([frame, pd.DataFrame([row])], ignore_index=True)
    canonical = {PRAME, MAGEA3, MAGEA4}
    ids = {PRAME: PRAME_TRANSCRIPT, MAGEA3: MAGEA3_TRANSCRIPT,
           MAGEA4: MAGEA4_TRANSCRIPT, HELD_OUT: HELD_OUT_TRANSCRIPT}
    frame["Canonical_Transcript_ID"] = frame.Ensembl_Gene_ID.map(ids)
    _patch_oncoref(monkeypatch, frame=frame, canonical=canonical,
                   unfiltered=set(frame.Ensembl_Gene_ID))
    class Transcripts(dict):
        def __call__(self):
            return list(self.values())
    transcripts = Transcripts({tid: SimpleNamespace(gene_id=gene, transcript_id=tid, protein_id=f"protein-{tid}",
                   protein_sequence=SEQUENCE, transcript_version=3)
                   for gene, tid in ids.items()})
    # Same gene, different isoform and sequence: do not select the canonical
    # isoform merely because it is the OncoRef representative.
    transcripts[ALTERNATE_TRANSCRIPT] = SimpleNamespace(
        gene_id=PRAME, protein_id="alternate-protein", protein_sequence=SEQUENCE[::-1],
        transcript_version=2)
    transcripts[ALTERNATE_TRANSCRIPT].transcript_id = ALTERNATE_TRANSCRIPT
    genes = {r.Ensembl_Gene_ID: SimpleNamespace(gene_id=r.Ensembl_Gene_ID,
             gene_name=r.Symbol, contig="22", start=1, end=100, strand="+",
             transcripts=[t for t in transcripts.values() if t.gene_id == r.Ensembl_Gene_ID])
             for r in frame.itertuples()}
    def lookup(gid):
        if gid not in genes:
            raise ValueError(gid)
        return genes[gid]
    return SimpleNamespace(
        species=SimpleNamespace(latin_name="homo_sapiens"), reference_name="GRCh38",
        annotation_name="ensembl", annotation_version=93, release=93,
        gene_ids=lambda: list(genes), gene_by_id=lookup,
        transcript_ids=lambda: list(transcripts), transcript_by_id=transcripts.__getitem__,
        transcripts=transcripts)


def run(tmp_path, genome, rows, *, level="gene", suffix="tsv", policy=None,
        input_contract=None, exclusion_policy=None):
    path = tmp_path / f"expression.{suffix}"
    pd.DataFrame(rows, columns=["feature_id", "patient_tpm"]).to_csv(
        path, sep="\t" if suffix == "tsv" else ",", index=False)
    return admit_cta_expression(
        path, input_contract=input_contract or contract(level),
        admission_policy=policy or CTAAdmissionPolicy(2, "TPM"), genome=genome,
        exclusion_policy=exclusion_policy)


@pytest.mark.parametrize("suffix", ["tsv", "csv"])
def test_gene_measurements_select_reference_proteins_and_preserve_unknowns(tmp_path, genome, suffix):
    result = run(tmp_path, genome, [
        (PRAME + ".8", 2), (HELD_OUT, 9), (MAGEA3, 100),
        (MAGEA4, 3), (NOT_CTA, 99)], suffix=suffix)
    assert [d.status for d in result.decisions] == [
        "admitted", "held_out", "excluded_target", "admitted", "not_cta"]
    assert result.input_sha256 == hashlib.sha256(
        (tmp_path / f"expression.{suffix}").read_bytes()).hexdigest()
    assert result.reference_assembly == "GRCh38"
    assert result.annotation_name == "ensembl"
    assert result.annotation_version == "93"
    first = result.admitted_antigens[0]
    assert first.gene_id == PRAME
    assert first.transcript_ids == (PRAME_TRANSCRIPT,)
    assert first.protein_ids == ("protein-" + PRAME_TRANSCRIPT,)
    assert first.species == "Homo sapiens"
    evidence = first.tumor_specificity.evidence_records[1]
    assert evidence.numeric_value == 2
    assert dict(evidence.details) == {
        "assay": "bulk RNA-seq", "measurement_level": "gene",
        "input_identifier": PRAME + ".8", "input_sha256": result.input_sha256}
    assert dict(first.source_metadata)["sequence_selection"] == "oncoref_canonical"
    # Identical proteins keep independent measurements and source identities.
    assert len(result.admitted_antigens) == 2
    assert result.admitted_antigens[0].amino_acids == result.admitted_antigens[1].amino_acids
    assert result.admitted_antigens[0].source_identifier != result.admitted_antigens[1].source_identifier
    assert not hasattr(first, "n_rna_alt")


def test_transcript_measurement_selects_only_measured_isoform(tmp_path, genome):
    result = run(tmp_path, genome, [(ALTERNATE_TRANSCRIPT + ".2", 4),
                                  (PRAME_TRANSCRIPT + ".3", 0)], level="transcript")
    assert [d.status for d in result.decisions] == ["admitted", "held_out"]
    antigen = result.admitted_antigens[0]
    assert antigen.amino_acids == SEQUENCE[::-1]
    assert antigen.transcript_ids == (ALTERNATE_TRANSCRIPT,)
    evidence = result.decisions[0].assessment.tumor_expression
    assert evidence.measurement_level == "transcript"
    assert evidence.transcript_id == ALTERNATE_TRANSCRIPT + ".2"
    assert evidence.value == 4
    assert dict(antigen.source_metadata)["sequence_selection"] == "measured_transcript"


def test_missing_expression_is_unknown_and_zero_is_measured(tmp_path, genome):
    result = run(tmp_path, genome, [(PRAME, None), (MAGEA4, 0)])
    assert result.decisions[0].status == "unknown_expression"
    assert result.decisions[0].value is None
    assert result.decisions[0].assessment is None
    assert result.decisions[0].reference_resolution.evidence.gene_id == PRAME
    zero = result.decisions[1]
    assert zero.status == "held_out"
    assert zero.assessment.tumor_expression.value == 0
    assert not result.admitted_antigens


def test_exclusion_exceptions_do_not_change_negative_reference(tmp_path, genome):
    rows = [(MAGEA3, 100), (MAGEA4, 3)]
    default = run(tmp_path, genome, rows)
    allow_all = run(tmp_path, genome, rows, exclusion_policy=CTATargetExclusionPolicy((), ()))
    assert len(default.admitted_antigens) == 1
    assert len(allow_all.admitted_antigens) == 2
    excluded = default.decisions[0].reference_resolution
    assert excluded == allow_all.decisions[0].reference_resolution
    assert default.admitted_antigens[0].self_reference_excluded_gene_ids == (
        allow_all.admitted_antigens[0].self_reference_excluded_gene_ids)
    assert MAGEA3 in excluded.self_reference_excluded_gene_ids


@pytest.mark.parametrize("value", [-1, float("inf"), -float("inf"), "broken", "True"])
def test_invalid_abundance_fails_instead_of_becoming_unknown(tmp_path, genome, value):
    with pytest.raises(ValueError, match="measurement"):
        run(tmp_path, genome, [(PRAME, value)])


@pytest.mark.parametrize("identifier", [None, "PRAME", "ENST00000398741", "ENSG1", ""])
def test_invalid_gene_namespace_does_not_silently_drop_rows(tmp_path, genome, identifier):
    with pytest.raises(ValueError, match="Ensembl gene ID"):
        run(tmp_path, genome, [(identifier, 4)])


def test_duplicate_stable_feature_requires_explicit_aggregation(tmp_path, genome):
    with pytest.raises(ValueError, match="Duplicate normalized"):
        run(tmp_path, genome, [(PRAME + ".7", 2), (PRAME + ".8", 3)])


def test_multiple_transcripts_are_not_aggregated_into_gene_expression(tmp_path, genome):
    result = run(tmp_path, genome, [(PRAME_TRANSCRIPT, 1), (ALTERNATE_TRANSCRIPT, 1)],
                 level="transcript")
    assert [d.status for d in result.decisions] == ["held_out", "held_out"]
    assert [d.assessment.tumor_expression.value for d in result.decisions] == [1, 1]


@pytest.mark.parametrize("level,identifier", [("gene", PRAME), ("transcript", UNKNOWN_TRANSCRIPT)])
def test_unresolved_transcripts_remain_auditable(tmp_path, genome, level, identifier):
    del genome.transcripts[PRAME_TRANSCRIPT]
    result = run(tmp_path, genome, [(identifier, 9)], level=level)
    assert result.decisions[0].status == "unresolved_transcript"
    assert not result.admitted_antigens


def test_non_coding_transcripts_and_version_mismatches_are_held_out(tmp_path, genome):
    genome.transcripts[ALTERNATE_TRANSCRIPT].protein_sequence = None
    result = run(tmp_path, genome, [(ALTERNATE_TRANSCRIPT + ".2", 5),
                                  (PRAME_TRANSCRIPT + ".4", 5)], level="transcript")
    assert [d.status for d in result.decisions] == ["no_protein_sequence", "transcript_version_mismatch"]
    assert not result.admitted_antigens


def test_missing_annotation_errors_propagate(tmp_path, genome):
    def unavailable():
        raise ValueError("Missing annotation cache")
    genome.transcript_ids = unavailable
    with pytest.raises(ValueError, match="Missing annotation cache"):
        run(tmp_path, genome, [(PRAME, 3)])


def test_inconsistent_selected_gene_fails(tmp_path, genome):
    genome.transcripts[PRAME_TRANSCRIPT].gene_id = MAGEA4
    with pytest.raises(ValueError, match="different gene"):
        run(tmp_path, genome, [(PRAME, 3)])


def test_units_are_not_converted_or_inferred(tmp_path, genome):
    with pytest.raises(ValueError, match="units"):
        run(tmp_path, genome, [(PRAME, 3)], policy=CTAAdmissionPolicy(2, "FPKM"))


def test_replay_keeps_policies_decisions_and_sources_without_upstream_reads(tmp_path, genome, monkeypatch):
    result = run(tmp_path, genome, [(PRAME, 3), (MAGEA3, 10), (MAGEA4, None)])
    saved = tmp_path / "admission.json"
    result.save(saved)
    (tmp_path / "expression.tsv").unlink()
    def forbidden(*args, **kwargs):
        raise AssertionError("Saved admission must not query upstream")
    monkeypatch.setattr("oncoref.cta.cta_evidence", forbidden)
    monkeypatch.setattr("oncoref.cta.cta_unfiltered_gene_ids", forbidden)
    genome.transcript_by_id = forbidden
    reloaded = CTAExpressionResult.load(saved)
    assert reloaded == result
    assert reloaded.admitted_antigens == result.admitted_antigens
    assert reloaded.exclusion_policy == CTATargetExclusionPolicy()


@pytest.mark.parametrize("field", ["id_column", "value_column", "expression_unit", "sample_id",
                                   "evidence_source", "evidence_version", "assay"])
def test_contract_requires_declared_provenance(field):
    with pytest.raises(ValueError, match=field):
        contract(**{field: ""})


@pytest.mark.parametrize("kwargs", [
    {"gene_patterns": "MAGE*"}, {"gene_exceptions": "MAGEA4"},
    {"gene_patterns": ("",)}, {"gene_exceptions": ("MAGEA*",)},
])
def test_exclusion_policy_rejects_ambiguous_configuration(kwargs):
    with pytest.raises(ValueError):
        CTATargetExclusionPolicy(**kwargs)


def test_transcript_evidence_cannot_admit_another_isoform(genome):
    measurement = PatientTumorExpressionEvidence(
        gene_id=PRAME, sample_id="patient", value=3, unit="TPM", assay="RNA-seq",
        evidence_source="Salmon", evidence_version="1", measurement_level="transcript",
        transcript_id=ALTERNATE_TRANSCRIPT)
    with pytest.raises(ValueError, match="measured transcript"):
        assess_cta_antigen(amino_acids=SEQUENCE, gene_id=PRAME, genome=genome, tumor_expression=measurement,
                           policy=CTAAdmissionPolicy(2, "TPM"), transcript_ids=(PRAME_TRANSCRIPT,))


def test_expression_antigen_reaches_existing_peptide_and_mrna_builders(tmp_path, genome):
    from mhctools import RandomBindingPredictor
    from vaxrank.core_logic import vaccine_peptides_for_antigen
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.peptide import PeptideConstructConfig, assemble_peptide_constructs
    from vaxrank.mrna import RNAConstructConfig, assemble_mrna_constructs

    result = run(tmp_path, genome, [(PRAME, 5)])
    antigen = result.admitted_antigens[0]
    # Synthetic prediction smoke: a constant DSL score checks consumption of
    # the admitted source, without asserting a biological binding outcome.
    windows = vaccine_peptides_for_antigen(
        antigen=antigen, mhc_predictor=RandomBindingPredictor(["HLA-A*02:01"]),
        epitope_config=EpitopeConfig(score_expr="1", min_epitope_score=0),
        vaccine_peptide_length=25)
    assert windows
    ranked = [(antigen, windows)]
    peptides = assemble_peptide_constructs(ranked, options=PeptideConstructConfig(
        mode="minimal_epitope", min_antigen_length_aa=8))
    mrna = assemble_mrna_constructs(ranked, options=RNAConstructConfig(
        antigen_content="minimal_epitope", signal_peptide=None, include_mitd=False,
        optimize_linkers=False))
    assert peptides and mrna
    assert all(window.mutant_protein_fragment is None for window in windows)
    assert all(window.antigen.gene_id == PRAME for window in windows)
    assert all(window.antigen.tumor_specificity == antigen.tumor_specificity for window in windows)
    assert all(peptide.sequence in SEQUENCE for peptide in peptides)


@pytest.mark.parametrize("suffix", ["csv", "tsv"])
@pytest.mark.parametrize("level", ["gene", "transcript"])
@pytest.mark.parametrize("n_rows", [0, 1, 3])
def test_expression_tables_of_any_size_use_the_public_loader(tmp_path, genome, suffix, level, n_rows):
    features = ([PRAME, MAGEA4, HELD_OUT] if level == "gene" else
                [PRAME_TRANSCRIPT, ALTERNATE_TRANSCRIPT, MAGEA4_TRANSCRIPT])
    result = run(tmp_path, genome, [(identity, 3) for identity in features[:n_rows]],
                 suffix=suffix, level=level)
    assert len(result.decisions) == n_rows
    assert [d.input_identifier for d in result.decisions] == features[:n_rows]
    assert len(result.admitted_antigens) == (min(n_rows, 2) if level == "gene" else n_rows)
