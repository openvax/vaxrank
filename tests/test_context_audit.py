"""Software contracts for full contexts, not synthetic biological validation."""

from dataclasses import replace
from unittest.mock import patch

from Bio.Data import CodonTable
import mhctools
import pandas as pd
import pytest
from mhctools.peptidases import get_cleavage_model

from vaxrank.context_audit import MHCRequest, SequenceContextAudit, audit_sequence_contexts
from vaxrank.context_report import write_sequence_context_audits
from vaxrank.construct_sequence import ConstructChemicalModification, ConstructPlacement, ConstructSequenceEdit
from vaxrank.mrna import RNAConstruct
from vaxrank.native_serialization import from_native_json, to_native_json
from vaxrank.sequence_context import (
    SequenceContext, final_sequence_context, native_sequence_context, translated_sequence_context,
)

from .test_cleavage_inference import identity
from .test_construct_audit import RecordingPredictor, jlf_construct
from .test_construct_sequence import DOCUMENTED, JLF, NATIVE, construct, native_antigen


REQUEST = MHCRequest("pMHC_affinity", "software_fixture", "1", "HLA-B*27:05", 9)
FREE = dict(n_term="free", c_term="free", chemistry_basis="Conditional unmodified free-terminal analysis, not manufacturing evidence")


def product_fixture():
    record = replace(construct(), modality="mrna")
    protein = "M" + NATIVE + "GG" + NATIVE + "K"
    codons = {aa: codon for codon, aa in CodonTable.unambiguous_dna_by_id[1].forward_table.items()}
    nt = "".join(codons[aa] for aa in protein) + "TAA"
    return RNAConstruct("multi-antigen", ["first", "second"], "CCC" + nt + "TTTAAA",
                        cds_aa=protein, cds_nt=nt, no_polya_nt="CCC" + nt + "TTT",
                        full_nt="CCC" + nt + "TTTAAA", poly_a_nt="AAA",
                        elements={"utr_5p": {"nt": "CCC"}, "utr_3p": {"nt": "TTT"},
                                  "linkers_per_junction": [{"aa": "GG", "name": "fixture-linker"}]},
                        construct_placements=(ConstructPlacement(record, 1), ConstructPlacement(record, 25)))


def scores(sequences):
    return {s: ([.2 if s == NATIVE else .8] * (len(s) - 1) + [0.]) for s in sequences}


def run(contexts, predictor=None, requests=(REQUEST,), **kwargs):
    with identity(), patch.object(mhctools.Pepsickle, "cleavage_probs_many", side_effect=scores):
        return audit_sequence_contexts(contexts, mhc_predictor=predictor,
                                       mhc_requests=requests, pepsickle=True, **kwargs)


def test_distinct_context_scores_and_repeated_source_provenance_are_not_merged():
    native = native_sequence_context(native_antigen())
    final = final_sequence_context(jlf_construct(), **FREE)
    repeated = replace(final, name="independent manufacturing occurrence")
    predictor = RecordingPredictor()
    with identity(), patch.object(mhctools.Pepsickle, "cleavage_probs_many", side_effect=scores) as infer:
        audits = audit_sequence_contexts(iter([native, final, repeated]), mhc_predictor=predictor,
                                        mhc_requests=[REQUEST], pepsickle=True)
    infer.assert_called_once_with((NATIVE, JLF))
    assert predictor.calls == [{"context_0": NATIVE, "context_1": JLF}]
    assert audits[0].profiles[0].sites[0].score == .2
    assert audits[1].profiles[0].sites[0].score == .8
    assert audits[1].profiles[0].peptide.source_id != audits[2].profiles[0].peptide.source_id
    assert [len(a.ligands) for a in audits] == [14, 16, 16]
    assert all(a.coverage[0][2] == () for a in audits)
    assert audits[1].context.construct.native_residue_offsets[-2:] == (None, None)
    assert from_native_json(to_native_json(audits), tuple) == audits


def test_complete_product_includes_linkers_both_occurrences_and_full_provenance():
    product = product_fixture()
    context = translated_sequence_context(product, DOCUMENTED)
    audit, = run([context], RecordingPredictor())
    assert len(context.sequence) == 48 and len(audit.ligands) == 40
    assert context.sequence[23:25] == "GG"
    assert [t.start for t in context.targets] == [8, 32]
    assert len(context.source_antigens) == 2
    assert [ligand.start for ligand in audit.ligands if ligand.peptide == NATIVE[:9]] == [1, 25]
    assert "molecular_termini_not_established" in audit.profiles[0].reason_codes
    assert len(audit.profiles[0].sites) == 47
    assert from_native_json(to_native_json(audit), SequenceContextAudit) == audit
    product.elements["linkers_per_junction"][0]["name"] = "changed later"
    assert context.product.elements["linkers_per_junction"][0]["name"] == "fixture-linker"


@pytest.mark.parametrize("field,value", [("cds_nt", "ATGAAA"), ("cds_aa", "MKK"),
    ("no_polya_nt", "AAA"), ("full_nt", "AAA"), ("sequence", "AAA")])
def test_invalid_complete_product_views_cannot_enter_audit(field, value):
    product = product_fixture()
    setattr(product, field, value)
    with pytest.raises(ValueError):
        translated_sequence_context(product, DOCUMENTED)


def test_mutable_assembly_record_tampering_is_rechecked_before_export(tmp_path):
    audit, = run([translated_sequence_context(product_fixture(), DOCUMENTED)], RecordingPredictor())
    audit.context.product.full_nt = "AAA"
    with pytest.raises(ValueError, match="emitted full"):
        write_sequence_context_audits([audit], json_path=tmp_path / "invalid.json")
    assert not (tmp_path / "invalid.json").exists()


def test_cpn_never_treats_native_window_end_as_an_exposed_terminal_lysine():
    native = native_sequence_context(native_antigen(JLF))
    final = final_sequence_context(jlf_construct(), **FREE)
    audits = run([native, final], peptidase_predictors=[get_cleavage_model("cpn-basic")])
    assert audits[0].profiles[1].status == "unassessed"
    assert audits[0].profiles[1].reason_codes == ("native_window_molecular_termini_unestablished",)
    assert audits[1].profiles[1].sites[0].bond == 23
    assert audits[1].profiles[1].sites[0].status == "matched"
    assert audits[1].profiles[1].sites[0].score is None


def test_encoded_antigen_segment_does_not_claim_exposed_termini():
    context = final_sequence_context(replace(construct(), modality="mrna"), **FREE)
    audit, = run([context], peptidase_predictors=[get_cleavage_model("cpn-basic")])
    assert audit.profiles[1].reason_codes == ("encoded_segment_molecular_termini_unestablished",)


@pytest.mark.parametrize("chemistry", [{}, {"c_term": "amidated", "chemistry_basis": "conditional"}])
def test_unknown_or_unsupported_final_chemistry_is_not_silently_stripped(chemistry):
    context = final_sequence_context(jlf_construct(), **chemistry)
    with identity(), patch.object(mhctools.Pepsickle, "cleavage_probs_many") as backend:
        audit, = audit_sequence_contexts([context], pepsickle=True,
                                        peptidase_predictors=[get_cleavage_model("cpn-basic")])
    backend.assert_not_called()
    assert all(p.status == "unassessed" for p in audit.profiles)


def test_documented_noncanonical_modification_cannot_be_overridden_by_free_assumption():
    record = replace(jlf_construct(), chemical_modifications=(
        ConstructChemicalModification(2, 3, "D-amino acid", DOCUMENTED),))
    predictor = RecordingPredictor()
    context = final_sequence_context(record, **FREE)
    with identity(), patch.object(mhctools.Pepsickle, "cleavage_probs_many") as backend:
        audit, = audit_sequence_contexts([context], mhc_predictor=predictor,
                                        mhc_requests=[REQUEST], pepsickle=True)
    backend.assert_not_called()
    assert predictor.calls == []
    assert audit.mhc_status == audit.profiles[0].status == "unassessed"
    assert "unsupported_chemical_modifications" in audit.reason_codes


@pytest.mark.parametrize("kind", ["mutation", "CTA", "viral"])
def test_masks_and_all_ligand_boundary_internal_overlays_are_source_agnostic(kind):
    context = native_sequence_context(native_antigen(intervals=((0, 2), (7, 7), (20, 22)), kind=kind))
    audit, = run([context], RecordingPredictor())
    overlays = audit.overlays(audit.profiles[0])
    assert len(overlays["targets"]) == 3 and len(overlays["ligands"]) == 14
    assert overlays["targets"][0][1].n_boundary_status == "sequence_endpoint"
    junction = overlays["targets"][1][1]
    assert junction.internal_sites == () and junction.n_boundary.bond == junction.c_boundary.bond == 7
    assert overlays["targets"][2][1].c_boundary_status == "sequence_endpoint"
    last_ligand, targets, evidence = overlays["ligands"][-1]
    assert last_ligand.end == len(NATIVE) and targets[0].start == 20
    assert evidence.c_boundary_status == "sequence_endpoint"
    assert len(evidence.internal_sites) == 8


def test_edit_deletion_junction_masks_map_to_actual_final_coordinates():
    native = native_antigen(intervals=((7, 7),))
    record = construct(NATIVE[2:], native=native, edits=(ConstructSequenceEdit(0, 2, "", DOCUMENTED),))
    audit, = run([final_sequence_context(record, **FREE)], RecordingPredictor())
    assert [(t.start, t.end) for t in audit.context.targets] == [(5, 5)]
    assert audit.profiles[0].peptide.sequence == NATIVE[2:]


def test_partial_mhc_outputs_retain_missing_allele_length_and_occurrence_axes():
    class Partial(RecordingPredictor):
        def predict_from_named_sequences(self, sequences):
            return super().predict_from_named_sequences(sequences).iloc[::2]
    extra_allele = replace(REQUEST, allele="HLA-A*01:01")
    extra_length = replace(REQUEST, peptide_length=10)
    extra_model = replace(REQUEST, predictor_name="other_model")
    audit, = run([native_sequence_context(native_antigen())], Partial(),
                 requests=(REQUEST, extra_allele, extra_length, extra_model))
    assert audit.mhc_status == "predictions_returned"
    assert audit.coverage[0][2] == tuple(range(1, 14, 2))
    assert len(audit.coverage[1][2]) == len(audit.coverage[3][2]) == 14
    assert len(audit.coverage[2][2]) == 13
    assert audit.unrequested_predictions == ()


def test_unrequested_output_is_retained_not_filtered():
    audit, = run([native_sequence_context(native_antigen())], RecordingPredictor(),
                 requests=(replace(REQUEST, allele="HLA-A*01:01"),))
    assert len(audit.unrequested_predictions) == 14
    assert audit.coverage[0][1] == ()


@pytest.mark.parametrize("failure,empty,expected", [(True, False, "failed"), (False, True, "no_predictions")])
def test_failure_and_empty_output_are_not_no_cuts_or_no_binders(failure, empty, expected):
    audit, = run([native_sequence_context(native_antigen())], RecordingPredictor(failure=failure, empty=empty))
    assert audit.mhc_status == expected and audit.reason_codes
    assert len(audit.coverage[0][2]) == 14
    assert audit.profiles[0].status == "observations_returned"


@pytest.mark.parametrize("error", ["source", "coordinates", "peptide", "duplicates", "columns"])
def test_malformed_outputs_fail_closed_without_unrelated_processing_loss(error):
    class Broken(RecordingPredictor):
        def predict_from_named_sequences(self, sequences):
            frame = super().predict_from_named_sequences(sequences)
            if error == "source":
                frame.loc[0, "source_sequence_name"] = "unrequested"
            elif error == "coordinates":
                frame.loc[0, "peptide_offset"] = 1000
            elif error == "peptide":
                frame.loc[0, "peptide"] = "AAAAAAAAA"
            elif error == "duplicates":
                frame = pd.concat([frame, frame.iloc[:1]])
            else:
                frame = pd.concat([frame, frame[["peptide_offset"]]], axis=1)
            return frame
    audit, = run([native_sequence_context(native_antigen())], Broken())
    assert audit.mhc_status == "failed" and audit.error_message
    assert audit.profiles[0].sites and not audit.ligands


def test_readable_report_and_native_json_include_complete_graph_and_derived_overlays(tmp_path):
    context = replace(final_sequence_context(jlf_construct(), **FREE), name="<script>bad</script>")
    audit, = run([context], RecordingPredictor(), peptidase_predictors=[get_cleavage_model("cpn-basic")])
    json_path, html_path = tmp_path / "audit.json", tmp_path / "audit.html"
    write_sequence_context_audits([audit], json_path=json_path, html_path=html_path)
    restored = from_native_json(json_path.read_text(), dict)
    assert restored["audits"] == (audit,)
    assert restored["coverage"][0] == audit.coverage
    assert restored["overlays"][0][0] == audit.overlays(audit.profiles[0])
    html = html_path.read_text()
    assert "<script>bad</script>" not in html and "&lt;script&gt;bad&lt;/script&gt;" in html
    for text in (JLF, "L16", "23", "cpn-basic", "fixture", "HLA-B*27:05", "ENSG00000197102",
                 "NM_001376", "patient", "sequence", "unassessed", "not manufacturing evidence"):
        assert text.lower() in html.lower()


@pytest.mark.parametrize("change", [{"n_term": "free"}, {"kind": "invented"}, {"sequence": "KKK"}])
def test_context_rejects_ambiguous_or_false_identity(change):
    with pytest.raises(ValueError):
        replace(native_sequence_context(native_antigen()), **change)


def test_request_and_profile_identity_are_checked_before_export():
    context = native_sequence_context(native_antigen())
    with pytest.raises(ValueError, match="explicit"):
        audit_sequence_contexts([context], mhc_predictor=RecordingPredictor())
    with pytest.raises(ValueError, match="distinct"):
        audit_sequence_contexts([context], mhc_requests=[REQUEST, REQUEST])
    audit, = run([context], RecordingPredictor())
    with pytest.raises(ValueError, match="different context"):
        replace(audit, context=replace(context, name="different"))


def test_empty_audit_does_not_load_or_run_models():
    with patch("vaxrank.context_audit.audit_pepsickle_inputs") as model:
        assert audit_sequence_contexts([], pepsickle=True) == ()
    model.assert_not_called()


def test_enzyme_inference_deduplicates_chemistry_without_losing_source_occurrences():
    enzyme = get_cleavage_model("cpn-basic")
    context = final_sequence_context(jlf_construct(), **FREE)
    with patch.object(type(enzyme), "predict", autospec=True, side_effect=type(enzyme).predict) as predict:
        first, second = audit_sequence_contexts([context, replace(context, name="second")],
                                                peptidase_predictors=[enzyme])
    assert predict.call_count == 1
    assert first.profiles[0].sites == second.profiles[0].sites
    assert first.profiles[0].peptide.source_id != second.profiles[0].peptide.source_id


def test_complete_product_without_placements_is_not_claimed_to_have_no_targets():
    product = product_fixture()
    product.construct_placements = ()
    context = translated_sequence_context(product, DOCUMENTED)
    assert context.target_mapping_status == "unassessed_no_placements"
    assert context.targets == ()
    assert from_native_json(to_native_json(context), SequenceContext) == context


def test_multi_gene_cta_viral_product_never_invents_one_antigen_policy():
    product = product_fixture()
    first, second = product.construct_placements
    cta = replace(first.construct.native_antigen, kind="CTA", gene_name="CTA-gene", gene_id="ENSG_CTA")
    viral = replace(second.construct.native_antigen, kind="viral", gene_name="viral-gene", gene_id="VIRAL_ID",
                    species="virus", transcript_ids=("VIRAL_TX",), source_identifier="viral_source")
    product.construct_placements = (replace(first, construct=replace(first.construct, native_antigen=cta)),
                                    replace(second, construct=replace(second.construct, native_antigen=viral)))
    audit, = run([translated_sequence_context(product, DOCUMENTED)], RecordingPredictor())
    assert [(s.kind, s.gene_id, s.species) for s in audit.context.source_antigens] == [
        ("CTA", "ENSG_CTA", "human"), ("viral", "VIRAL_ID", "virus")]
    assert audit.context.native_antigen is None and audit.context.construct is None
    assert from_native_json(to_native_json(audit), SequenceContextAudit) == audit


def test_unresolved_mapping_remains_explicit_without_silently_claiming_zero_targets():
    record = replace(jlf_construct(), mapping_status="unresolved", mapping_reason="No source alignment",
                     native_start=None, native_end=None, sequence_edits=())
    context = final_sequence_context(record, **FREE)
    audit, = run([context], RecordingPredictor())
    assert context.target_mapping_status == "unresolved"
    assert audit.mhc_status == audit.profiles[0].status == "unassessed"
    assert context.source_antigens[0].gene_name == "DYNC1H1"


@pytest.mark.parametrize("change", [{"allele": "not-an-allele"}, {"allele": "HLA-A"},
    {"peptide_length": -1}, {"peptide_length": True}, {"kind": ""}])
def test_bad_mhc_request_fails_at_input(change):
    with pytest.raises(ValueError):
        replace(REQUEST, **change)
