"""Final-context assessment and actual emitted-file regressions."""

import csv
from dataclasses import replace
import json

import pandas as pd
import pytest
from topiary import TopiaryPredictor

from vaxrank.construct_audit import (
    ConstructAudit,
    audit_construct_sequence,
    load_construct_sequences,
    save_construct_sequences,
    write_construct_audits,
)
from vaxrank.construct_sequence import (
    ConstructChemicalModification,
    ConstructPlacement,
    ConstructSequenceEdit,
)
from vaxrank.mrna import RNAConstruct, write_mrna_outputs
from vaxrank.native_serialization import from_native_json, to_native_json
from vaxrank.peptide import (
    PeptideConstructConfig,
    peptide_construct_from_sequence,
    write_peptide_outputs,
)

from .test_construct_sequence import DOCUMENTED, JLF, NATIVE, RATIONALE, construct, native_antigen
from .test_reference_proteome import create_mock_genome, create_mock_transcript


class RecordingPredictor(TopiaryPredictor):
    """Software-contract fixture, not a biological prediction cache."""
    def __init__(self, failure=False, empty=False):
        self.calls = []
        self.failure = failure
        self.empty = empty

    def predict_from_named_sequences(self, sequences):
        self.calls.append(dict(sequences))
        if self.failure:
            raise RuntimeError("fixture model unavailable")
        if self.empty:
            return pd.DataFrame()
        rows = []
        for name, sequence in sequences.items():
            for offset in range(len(sequence) - 8):
                rows.append(dict(
                    peptide=sequence[offset:offset + 9], peptide_offset=offset,
                    peptide_length=9, kind="pMHC_affinity", allele="HLA-B*27:05",
                    prediction_method_name="software_fixture", predictor_version="1",
                    source_sequence_name=name, score=0.5, value=100.0, percentile_rank=1.0))
        return pd.DataFrame(rows)


def jlf_construct():
    return construct(JLF, edits=(ConstructSequenceEdit(22, 22, "KK", DOCUMENTED, RATIONALE),))


def test_final_context_is_scanned_including_new_tail_boundary_and_self_sources():
    record = jlf_construct()
    self_peptide = JLF[-9:]
    genome = create_mock_genome([
        create_mock_transcript("SELF_TX", "M" + self_peptide, gene_id="ENSG_SELF")],
        species_name="Homo sapiens", release=114)
    predictor = RecordingPredictor()
    audit = audit_construct_sequence(record, predictor, genome=genome)
    assert predictor.calls == [{record.cache_identity: JLF}]
    assert record.native_sequence == NATIVE
    assert audit.mhc_status == "predictions_returned"
    assert len(audit.mhc_assessment.ligands) == 16
    tail = audit.mhc_assessment.ligands[-1]
    assert tail.peptide == self_peptide
    assert tail.crosses_construct_boundary
    assert not tail.overlaps_targetable
    assert tail.self_reference_match.source_provenance_complete
    assert tail.self_reference_match.sources[0].gene_id == "ENSG_SELF"
    assert set(audit.unassessed_processing_compartments) == {
        "proteasomal", "extracellular_serum", "endolysosomal"}
    assert from_native_json(to_native_json(audit), ConstructAudit) == audit


def test_shared_cta_non_cta_sources_survive_final_construct_assessment():
    cta_gene = "ENSG00000185686"
    native = replace(native_antigen(kind="CTA"), self_reference_excluded_gene_ids=(cta_gene,))
    record = replace(jlf_construct(), native_antigen=native)
    shared = JLF[:9]
    genome = create_mock_genome([
        create_mock_transcript("CTA_TX", shared, gene_id=cta_gene),
        create_mock_transcript("SELF_TX", shared, gene_id="ENSG_SELF")])
    audit = audit_construct_sequence(record, RecordingPredictor(), genome=genome)
    match = audit.mhc_assessment.ligands[0].self_reference_match
    assert [source.gene_id for source in match.sources] == ["ENSG_SELF"]
    assert match.occurs


def test_unavailable_chemical_prediction_never_calls_backend():
    record = replace(jlf_construct(), chemical_modifications=(
        ConstructChemicalModification(24, 24, "amidation", DOCUMENTED),))
    predictor = RecordingPredictor()
    audit = audit_construct_sequence(record, predictor)
    assert audit.mhc_status == "unassessed"
    assert audit.reason_codes == ("unsupported_chemical_modifications",)
    assert predictor.calls == []
    assert audit.mhc_assessment is None


def test_no_model_missing_reference_empty_output_and_failure_stay_distinct():
    record = jlf_construct()
    assert audit_construct_sequence(record).reason_codes == ("mhc_predictor_not_requested",)
    assert audit_construct_sequence(record, RecordingPredictor()).reason_codes == (
        "self_reference_unavailable",)
    genome = create_mock_genome([])
    assert audit_construct_sequence(record, RecordingPredictor(empty=True), genome=genome).mhc_status == (
        "no_predictions")
    failed = audit_construct_sequence(record, RecordingPredictor(failure=True), genome=genome)
    assert failed.mhc_status == "failed"
    assert failed.error_message
    with pytest.raises(ValueError, match="Processing remains unassessed"):
        replace(failed, unassessed_processing_compartments=())


def test_actual_construct_input_and_report_outputs_preserve_provenance_and_escape_html(tmp_path):
    record = replace(jlf_construct(), name="<script>not markup</script>")
    input_path = tmp_path / "input.json"
    save_construct_sequences([record], input_path)
    assert load_construct_sequences(input_path) == [record]
    audit = audit_construct_sequence(record)
    json_path, html_path = tmp_path / "audit.json", tmp_path / "audit.html"
    write_construct_audits([audit], json_path=json_path, html_path=html_path)
    assert from_native_json(json_path.read_text(), list) == [audit]
    html = html_path.read_text()
    assert "<script>" not in html
    assert "&lt;script&gt;" in html
    for expected in (JLF, NATIVE, "reported", "solubility", "extracellular_serum", "unassessed"):
        assert expected in html
    input_path.write_text(to_native_json(["not a construct"]))
    with pytest.raises(ValueError, match="input.json.*Record 1"):
        load_construct_sequences(input_path)


def test_actual_manufacturing_fasta_manifest_and_order_form_roundtrip(tmp_path):
    record = jlf_construct()
    product = peptide_construct_from_sequence(record)
    fasta, manifest, order = (tmp_path / name for name in ("vaccine.fa", "manifest.json", "order.csv"))
    write_peptide_outputs([product], fasta, manifest, order)
    assert JLF in fasta.read_text()
    item = json.loads(manifest.read_text())[0]
    assert item["sequence"] == JLF
    restored = from_native_json(json.dumps(item["construct_provenance"]), tuple)
    assert restored == (ConstructPlacement(record, 0),)
    with order.open() as handle:
        row = next(csv.DictReader(handle))
    assert row["sequence"] == JLF
    assert from_native_json(row["construct_provenance_json"], tuple) == restored


def test_chemical_order_form_uses_recorded_modification_not_writer_defaults(tmp_path):
    record = replace(jlf_construct(), chemical_modifications=(
        ConstructChemicalModification(24, 24, "C-terminal amidation", DOCUMENTED),))
    product = peptide_construct_from_sequence(record)
    write_peptide_outputs([product], tmp_path / "seq.fa", tmp_path / "seq.json", tmp_path / "seq.csv")
    with (tmp_path / "seq.csv").open() as handle:
        row = next(csv.DictReader(handle))
    assert row["c_terminal_modification"] == "C-terminal amidation"
    assert row["sequence"] == JLF
    assert "amidation" in row["displayed_sequence"]
    with pytest.raises(ValueError, match="writer options"):
        write_peptide_outputs([product], tmp_path / "bad.fa", tmp_path / "bad.json",
                              options=PeptideConstructConfig(n_terminal_acetylation=True))
    assert not (tmp_path / "bad.fa").exists()


def test_modified_product_cannot_export_stale_provenance_or_lose_sidecar(tmp_path):
    product = peptide_construct_from_sequence(jlf_construct())
    with pytest.raises(ValueError, match="sidecar"):
        write_peptide_outputs([product], tmp_path / "seq.fa")
    product.sequence = NATIVE
    with pytest.raises(ValueError, match="does not match"):
        write_peptide_outputs([product], tmp_path / "seq.fa", tmp_path / "seq.json")
    assert not (tmp_path / "seq.fa").exists()


def test_mrna_manifest_records_encoded_changes_and_checks_actual_translation(tmp_path):
    native = native_antigen("MACDE", ((2, 3),))
    record = replace(construct("MACDEK", native=native, edits=(
        ConstructSequenceEdit(5, 5, "K", DOCUMENTED),)), modality="mrna")
    coding = "ATGGCTTGTGATGAAAAATAA"
    product = RNAConstruct(
        name="test", antigen_names=["DYNC1H1"], sequence=coding,
        cds_aa="MACDEK", cds_nt=coding, no_polya_nt=coding, full_nt=coding,
        construct_placements=(ConstructPlacement(record, 0),))
    manifest = tmp_path / "mrna.json"
    write_mrna_outputs([product], tmp_path / "mrna", manifest_path=manifest)
    item = json.loads(manifest.read_text())[0]
    assert item["cds"]["aa"] == "MACDEK"
    assert from_native_json(json.dumps(item["construct_provenance"]), tuple) == product.construct_placements
    product.cds_nt = "ATGAAATAA"
    with pytest.raises(ValueError, match="actual CDS translation"):
        write_mrna_outputs([product], tmp_path / "bad", manifest_path=tmp_path / "bad.json")
    assert not (tmp_path / "bad").exists()


@pytest.mark.parametrize("changed", ["cds_nt", "full_nt", "no_polya_nt", "sequence"])
def test_mrna_all_emitted_views_must_match_provenance(tmp_path, changed):
    native = native_antigen("MACDE", ((2, 3),))
    record = replace(construct(native=native, sequence="MACDE"), modality="mrna")
    coding = "ATGGCTTGTGATGAATAA"
    product = RNAConstruct(
        name="test", antigen_names=["x"], sequence=coding, cds_aa="MACDE",
        cds_nt=coding, no_polya_nt=coding, full_nt=coding,
        construct_placements=(ConstructPlacement(record, 0),))
    setattr(product, changed, getattr(product, changed)[:-3])
    with pytest.raises(ValueError):
        write_mrna_outputs([product], tmp_path / "out", manifest_path=tmp_path / "out.json")
    assert not (tmp_path / "out").exists()


def test_pathlike_fasta_suffix_rejection_is_preserved(tmp_path):
    with pytest.raises(ValueError, match="looks like a FASTA"):
        write_mrna_outputs([], tmp_path / "not-a-directory.fasta")


def test_chemical_termini_cannot_be_embedded_inside_a_longer_product():
    from vaxrank.construct_sequence import validate_construct_placements

    record = replace(jlf_construct(), chemical_modifications=(
        ConstructChemicalModification(0, 0, "N-acetylation", DOCUMENTED),))
    with pytest.raises(ValueError, match="terminus is internal"):
        validate_construct_placements((ConstructPlacement(record, 1),), "M" + JLF, "peptide")
