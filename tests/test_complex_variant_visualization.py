import json
from pathlib import Path
import struct
from xml.etree import ElementTree

from vaxrank.complex_variant_visualization import (
    _balanced_chunks,
    generate_complex_variant_results,
    render_complex_result_svg,
)


def test_summary_chunks_are_balanced_without_exceeding_page_capacity():
    chunks = _balanced_chunks(list(range(22)))
    assert [len(chunk) for chunk in chunks] == [6, 6, 5, 5]


def example_record(outcome="selected"):
    selected = outcome == "selected"
    return {
        "combined_score": 4.25,
        "decision_reason": "The result follows the evidence and selection gates.",
        "dna_only": "DNA does not establish the expressed mutant protein.",
        "id": "GENE1-complex-event",
        "interpretation": "Assembly establishes what Vaxrank may predict and rank.",
        "limitation": "Research use only.",
        "outcome": outcome,
        "protein_sequence": "MABCDEFGHIJK",
        "rna_evidence": {
            "alt_fragments": 7,
            "other_fragments": 0,
            "passes_vaxrank_filters": selected,
            "platform": "Long-read RNA",
            "ref_fragments": 12,
            "source": "T1 fixture",
            "top_protein_fragments": 6,
        },
        "sample": "T1",
        "selected_long_peptide": "ABCDEFGHIJK" if selected else None,
        "source_labels": ["fixture", "predictor"],
        "source_urls": ["https://example.org/evidence"],
        "summary_rna": "Assembly establishes a mutant protein.",
        "target_epitope_score": 2.5,
        "title": "GENE1 complex event",
        "top_epitopes": ([{
            "allele": "HLA-A*02:01",
            "epitope_score": 0.9,
            "ic50_nm": 32.1,
            "percentile_rank": 0.2,
            "sequence": "ABCDEFGHI",
        }] if selected else []),
        "variant": "chr1:g.100_120del",
        "variant_class": "Large deletion",
    }


def example_payload(records):
    return {
        "metadata": {
            "analysis_date": "2026-09-16",
            "hla_alleles": ["HLA-A*02:01"],
            "predictor_version": "4.2c",
            "software_versions": {
                "isovar": "1.17.0",
                "vaxrank": "3.18.2",
            },
        },
        "records": records,
    }


def test_result_svg_shows_selected_peptide_and_prediction():
    svg = render_complex_result_svg(example_record())
    ElementTree.fromstring(svg)
    assert '<rect width="1200" height="760" fill="#ffffff"/>' in svg
    assert "SELECTED" in svg
    assert "ABCDEFGHIJK" in svg
    assert "HLA-A*02:01" in svg
    assert "32.10" in svg


def test_result_svg_keeps_negative_result_explicit():
    svg = render_complex_result_svg(example_record("no_target_binder"))
    assert "NO TARGET BINDER" in svg
    assert "No peptide selected for a vaccine construct." in svg
    assert "No target epitope passed the configured filters." in svg


def test_result_svg_limits_table_to_five_rows_but_source_keeps_all():
    record = example_record()
    record["top_epitopes"] = [
        {
            "allele": "HLA-A*02:01",
            "epitope_score": 0.9 - index / 10,
            "ic50_nm": 32.1 + index,
            "percentile_rank": 0.2 + index,
            "sequence": "PEPTIDE%02d" % index,
        }
        for index in range(6)
    ]

    svg = render_complex_result_svg(record)

    assert "PEPTIDE04" in svg
    assert "PEPTIDE05" not in svg
    assert len(record["top_epitopes"]) == 6


def test_generate_compact_pdf_and_high_resolution_pages(tmp_path):
    records = [
        example_record("selected"),
        example_record("no_target_binder"),
        example_record("held_out"),
        example_record("no_translation"),
    ]
    for index, record in enumerate(records):
        record["id"] += "-%s" % index
        if record["outcome"] == "no_translation":
            record["protein_sequence"] = None
    source = tmp_path / "results.json"
    source.write_text(json.dumps(example_payload(records)))
    combined = tmp_path / "combined.pdf"
    run = generate_complex_variant_results(
        source,
        tmp_path / "runs",
        timestamp="2026-09-16T120000Z",
        combined_output=combined,
    )
    assert combined.read_bytes().startswith(b"%PDF")
    manifest = json.loads((run / "manifest.json").read_text())
    assert manifest["page_count"] == 6
    assert manifest["png_dimensions"] == {"width": 3600, "height": 2280}
    pngs = sorted(run.glob("*/result.png"))
    assert len(pngs) == 6
    assert struct.unpack(">II", pngs[0].read_bytes()[16:24]) == (3600, 2280)
    records = sorted(run.glob("*/record.json"))
    assert len(records) == 4


def test_generate_paginates_large_decision_matrix(tmp_path):
    records = [example_record() for _ in range(8)]
    for index, record in enumerate(records):
        record["id"] += "-%s" % index
    source = tmp_path / "results.json"
    source.write_text(json.dumps(example_payload(records)))

    run = generate_complex_variant_results(
        source,
        tmp_path / "runs",
        timestamp="2026-09-16T120001Z",
        png_scale=1,
    )

    manifest = json.loads((run / "manifest.json").read_text())
    assert manifest["page_count"] == 11
    assert [page["name"] for page in manifest["pages"][:2]] == [
        "summary", "summary-2"]
    assert len(list(run.glob("*/record.json"))) == 8


def test_committed_complex_inputs_cover_requested_cases_and_provenance():
    source = (
        Path(__file__).parents[1]
        / "examples" / "osteosarc_complex_results" / "source")
    antigens = json.loads((source / "assembled_antigens.json").read_text())
    rankings = json.loads((source / "assembled_rankings.json").read_text())
    sv_audits = json.loads((source / "sv_audits.json").read_text())
    audit = json.loads((source / "additional_candidate_audit.json").read_text())
    by_id = {record["id"]: record for record in antigens["antigens"]}
    ranked_by_id = {record["id"]: record for record in rankings["records"]}

    expected = {
        "GTF3C5-p-Glu503_Glu506del",
        "TECPR1-p-Thr259fs",
        "GLIS3-p-Ser775fs",
        "RNF213-p-Ile1070delLeu",
        "DYNC1H1-p-Val314Ile",
        "DYNC1H1-p-Gln3267His",
        "NAV2-p-Ala1809Val",
        "MAP2-compound-haplotype",
        "NTF3-compound-haplotype",
        "CD109-long-read-phased",
        "ZNF436-long-read-context",
    }
    assert set(by_id) == expected
    assert set(ranked_by_id) == expected
    assert all(record["evidence_gate"] == "pass" for record in by_id.values())
    assert all(record["source_metadata"] for record in by_id.values())
    assert by_id["CD109-long-read-phased"]["targetable_intervals"] == [
        [2, 3], [71, 72]]
    assert ranked_by_id["CD109-long-read-phased"][
        "selected_long_peptide"] is None
    assert "DRAH" in by_id["MAP2-compound-haplotype"]["amino_acids"]
    assert "DVSENY" in by_id["NTF3-compound-haplotype"]["amino_acids"]
    assert {
        record["id"] for record in sv_audits["records"]
    } == {
        "FOXO3-STRADA-CCDC47-unresolved",
        "PARD3B-CDKN2B-unresolved",
        "AMPH-internal-deletion-DNA-only",
        "chr21-unnamed-long-read-rich-unresolved",
        "KTN1-p-Leu340fs-unresolved",
        "GABBR1-SLC29A1-unresolved",
        "OTUD7A-FMN1-unresolved",
    }
    inventory = {
        record["id"]: record for record in sv_audits["inventory_only"]}
    assert inventory["MYO15B-chr17-75589953"][
        "published_vaccine_peptide"] == "AGRRAQAPTRVLGLAPP"
    assert inventory["DLG5--DLG5"]["class"] == "79.5-kb same-gene deletion"
    assert by_id["GTF3C5-p-Glu503_Glu506del"]["source_metadata"][
        "direct_alt_fragments_t2_ont"] == "128"
    audit_indels = {record["gene"]: record for record in audit["indels"]}
    assert audit_indels["KTN1"]["fragment_counts"]["2bc1fc291308debb"][
        "alt"] == 10
    audit_svs = {
        record["id"]: record for record in audit["structural_variants"]}
    assert audit_svs["OTUD7A--FMN1"]["tagged_ont_paths"][
        "2bc1fc291308debb"]["cell_umi_keys"] == 6
