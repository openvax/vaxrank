import json
import struct
from xml.etree import ElementTree

from vaxrank.evidence_visualization import (
    generate_evidence_figures,
    render_evidence_svg,
)


def example_record():
    return {
        "id": "GENE1-GENE2-t2",
        "title": "GENE1 :: GENE2 · T2",
        "kind": "Matched-sample long-read rescue",
        "locus": "chr1:100 → chr2:200",
        "conclusion": "LONG-READ RESCUE",
        "evidence": {
            "dna": {"status": "detected", "detail": "Ten DNA fragments."},
            "short_read": {
                "status": "not_detected", "detail": "No matched call."},
            "short_read_assembly": {
                "status": "not_assessed", "detail": "No call to assemble."},
            "long_read": {
                "status": "detected", "detail": "Twenty-two spanning reads."},
        },
        "nucleotide": {
            "left_label": "GENE1",
            "right_label": "GENE2",
            "left_sequence": "AACCGGTTAACCGGTT",
            "right_sequence": "TTGGCCAATTGGCCAA",
        },
        "protein": {
            "sequence": None,
            "reason": "A coding frame has not been established.",
        },
        "interpretation": "Long reads recover an expressed DNA-supported junction.",
        "source_urls": ["https://example.org/evidence"],
    }


def test_evidence_svg_separates_modalities_and_withholds_protein():
    svg = render_evidence_svg(example_record())
    ElementTree.fromstring(svg)
    assert '<rect width="1200" height="760" fill="#ffffff"/>' in svg
    assert "SHORT-READ + ASSEMBLY" in svg
    assert "LONG-READ RNA" in svg
    assert "PROTEIN WITHHELD" in svg
    assert "GENE1" in svg and "GENE2" in svg


def test_generate_evidence_run_has_high_resolution_png(tmp_path):
    source = tmp_path / "evidence.json"
    source.write_text(json.dumps([example_record()]))
    run = generate_evidence_figures(
        source,
        tmp_path / "figures",
        formats=("svg", "png"),
        timestamp="2026-09-16T041500Z",
    )
    png, = run.glob("*/evidence-context.png")
    assert struct.unpack(">II", png.read_bytes()[16:24]) == (3600, 2280)
    manifest = json.loads((run / "manifest.json").read_text())
    assert manifest["png_dimensions"] == {"width": 3600, "height": 2280}
    record, = run.glob("*/record.json")
    assert json.loads(record.read_text())["conclusion"] == "LONG-READ RESCUE"
