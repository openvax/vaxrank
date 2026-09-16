import json
import struct
from xml.etree import ElementTree

from vaxrank.complex_variant_visualization import (
    generate_complex_variant_results,
    render_complex_result_svg,
)


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
