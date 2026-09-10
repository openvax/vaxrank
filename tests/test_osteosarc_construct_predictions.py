"""Versioned real models, full-source context, exact overlays and offline reports."""

from dataclasses import replace
from hashlib import sha256
import socket

import pytest

from vaxrank.native_serialization import from_native_json, to_native_json
from tests.data.osteosarc.construct_audit.generate_predictions import canonical_audits
from tests.data.osteosarc.construct_audit.render_report import (
    comparison_rows, load_evidence, render, strong_predictions,
)
from .osteosarc_construct_helpers import DOCUMENTED, RNA


@pytest.fixture(autouse=True)
def no_network(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Cached Sid processing/self tests must remain offline")
    monkeypatch.setattr(socket.socket, "connect", forbidden)
    monkeypatch.setattr(socket, "create_connection", forbidden)


@pytest.fixture(scope="module")
def cached():
    return load_evidence()


def test_real_cache_has_full_reference_and_actual_model_provenance(cached):
    evidence, metadata, inventory = cached
    assert metadata["package_versions"]["isovar"] == "1.8.1"
    assert metadata["package_versions"]["mhctools"] == "3.41.0"
    assert metadata["model_assets_sha256"]
    assert len(metadata["executable_sha256"]) == 64
    reference = evidence["self_reference"]
    assert reference["release"] == 114 and len(reference["dataset_identity"]) == 64
    assert reference["annotated_protein_ids"] == 112375
    assert reference["protein_coding_transcripts"] == 89872
    assert len(reference["files"]) == 2
    catalog = evidence["cta_catalog"]
    assert len(catalog["gene_ids"]) == 390
    assert sha256(("\n".join(catalog["gene_ids"]) + "\n").encode()).hexdigest() == catalog["gene_ids_sha256"]
    assert len(inventory["source_ids"]) == 164
    assert "HLA-A*01:11N" not in metadata["alleles"]


def test_model_maps_use_full_original_regions_not_selected_fixture_contexts(cached):
    rows = {r["source_id"]: r for r in cached[0]["rna_results"]}
    assert len(rows) == 9
    for case in RNA["cases"]:
        row = rows[case["source_id"]]
        assert row["input_scope"] == "full_original_region"
        assert any(r.get("bam_sha256") == row["regional_bam_sha256"]
                   for r in case["source"]["acquisition_receipts"])
    assert rows["53f498a544883d51"]["n_alt_reads"] == 1560
    assert rows["1b66c15da594a3ef"]["n_alt_reads"] == 627
    assert rows["1b66c15da594a3ef"]["n_alt_fragments"] == 454
    assert not rows["0066232879babe83"]["passes_all_filters"]
    for sid in ("53f498a544883d51", "1120a096937e29e1", "8b8a9a02cbcbcf4b"):
        long_window = next(w for w in rows[sid]["windows"] if w["construct_id"] == "dync1h1-mrna-long")
        assert long_window["status"] == "absent_from_selected_rna_context"
    for audit in cached[0]["audits"]:
        assert all("rna_scope=full_original_region" in a.source_identifier for a in audit.context.source_antigens)


def test_every_requested_ligand_and_internal_bond_is_preserved(cached):
    audits = cached[0]["audits"]
    assert len(audits) == len({a.context.source_id for a in audits}) == 13
    for audit in audits:
        assert audit.mhc_status == "predictions_returned"
        assert len(audit.mhc_requests) == 40
        assert not audit.unrequested_predictions
        assert all(not missing for _, _, missing in audit.coverage)
        assert len(audit.ligands) == sum(max(0, len(audit.context.sequence) - k + 1) for k in (8, 9, 10, 11))
        assert all(len(ligand.predictions) == 10 for ligand in audit.ligands)
        profile = audit.profiles[0]
        assert profile.model.version == "0.1.3"
        assert profile.model_asset_sha256 == "96da1619131413722c6f4e21e6ff3435985073bceb61aec5fb02b32ed597affa"
        assert len(profile.raw_residue_scores) == len(audit.context.sequence)
        assert {s.bond for s in profile.sites} == set(range(1, len(audit.context.sequence)))
        assert profile.terminal_output == 0


def test_intended_target_overlays_are_not_single_residue_mutation_masks(cached):
    rows = comparison_rows(cached[0])
    for row in rows:
        audit = row["audit"]
        assert audit.context.sequence[row["target_start"]:row["target_end"]] == "KRFHATISF"
        assert row["target_end"] - row["target_start"] == 9
        assert all(target.end - target.start == 1 for target in audit.context.targets)
        profile, interval = row["target_profiles"][0]
        assert len(interval.internal_sites) == 8
        assert {s.bond for s in interval.internal_sites} == set(range(row["target_start"] + 1, row["target_end"]))
        assert all(s in profile.sites for s in interval.internal_sites)
    minimal = next(r for r in rows if r["audit"].context.name == "dync1h1-mrna-minimal")
    _, interval = minimal["target_profiles"][0]
    assert interval.n_boundary_status == interval.c_boundary_status == "sequence_endpoint"
    assert interval.n_boundary is interval.c_boundary is None


def test_real_non_target_self_flag_retains_all_source_transcripts(cached):
    evidence = cached[0]
    assert not evidence["self_matches"][DOCUMENTED["intended_ligand"]].occurs
    assert evidence["self_matches"][DOCUMENTED["intended_ligand"]].source_provenance_complete
    rows = comparison_rows(evidence)
    for name, offset in (("dync1h1-jlf-v2-v3", 12), ("dync1h1-mrna-long", 22)):
        row = next(r for r in rows if r["audit"].context.name == name)
        finding, = row["non_target_self_findings"]
        assert finding["ligand"].peptide == "DTGLKQAL" and finding["ligand"].start == offset
        prediction, = finding["predictions"]
        assert prediction.kind == "pMHC_presentation" and prediction.allele == "HLA-B*08:01"
        assert prediction.percentile_rank == .47  # Real model output, not biological truth.
        assert not finding["overlaps_intended"]
        assert len(finding["non_cta_sources"]) == 8
        assert {s.gene_id for s in finding["non_cta_sources"]} == {"ENSG00000197102"}
        assert all(s.transcript_id and s.protein_id and s.species == "homo_sapiens"
                   for s in finding["non_cta_sources"])


def test_dsl_strict_threshold_does_not_merge_prediction_kinds(cached):
    # Software-boundary sensitivity check, not new biological prediction data.
    audit = next(a for a in cached[0]["audits"] if a.context.name == "dync1h1-jlf-v2-v3")
    ligand = next(ligand for ligand in audit.ligands if ligand.peptide == "DTGLKQAL")
    changed = replace(ligand, predictions=tuple(
        replace(p, percentile_rank=.5) if p.kind == "pMHC_presentation" else p
        for p in ligand.predictions))
    assert strong_predictions([replace(audit, ligands=(changed,))]) == set()
    assert strong_predictions([replace(audit, ligands=(ligand,))])


def test_enzyme_scope_and_manufactured_chemistry_stay_conditional(cached):
    for audit in cached[0]["audits"]:
        enzyme = audit.profiles[1]
        if audit.context.kind == "native_context" or audit.context.construct.modality == "mrna":
            assert enzyme.status == "unassessed"
        else:
            assert "Conditional" in audit.context.chemistry_basis
            assert audit.context.construct.chemical_modifications == ()
            if audit.context.name == "dync1h1-jlf-v2-v3":
                assert [s.bond for s in enzyme.sites] == [23]
                assert enzyme.sites[0].score is None


def test_complete_evidence_graph_roundtrips(cached):
    evidence = cached[0]
    assert from_native_json(to_native_json(evidence), dict) == evidence


def test_real_observations_have_canonical_order_without_changing_values(cached):
    audit = cached[0]["audits"][0]
    reordered = replace(audit, ligands=tuple(
        replace(ligand, predictions=tuple(reversed(ligand.predictions)))
        for ligand in reversed(audit.ligands)))
    assert canonical_audits([reordered]) == (audit,)


def test_report_renders_all_contexts_and_unassessed_scopes(cached, tmp_path):
    render(tmp_path)
    html = (tmp_path / "sid-construct-audit.html").read_text()
    for audit in cached[0]["audits"]:
        assert audit.context.name in html
    for text in ("DTGLKQAL", "0.47", "ENSG00000197102", "Class II", "serum kinetics",
                 "JLF V1", "three-alternate-read", "CSBio", "164", "not a safety clearance"):
        assert text in html
    assert (tmp_path / "full-contexts.html").is_file()
    restored = from_native_json((tmp_path / "full-contexts.json").read_text(), dict)
    assert restored["audits"] == cached[0]["audits"]
