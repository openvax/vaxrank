import json
import struct
from xml.etree import ElementTree

import pandas as pd
import pytest

import vaxrank.mutation_visualization as mutation_visualization
from vaxrank.cli.mutation_visualization import main
from vaxrank.mutation_visualization import (
    generate_mutation_figures,
    global_alignment,
    render_mutation_svg,
    render_transcript_svg,
)


def example_records():
    shared = {
        "predicted_effect": "p.AAK2del",
        "predicted_effect_gene_name": "REPEAT1",
        "predicted_effect_transcript_id": "ENST00000000001",
        "predicted_effect_original_protein_sequence": "QQAAKPKVVKPKQQ",
        "predicted_effect_mutant_protein_sequence": "QQVVKPKQQ",
        "predicted_effect_aa_mutation_start_offset": 2,
        "trimmed_reference_protein_sequence": "AAKPKVVKPK",
        "trimmed_predicted_mutant_protein_sequence": "VVKPK",
        "protein_sequence_mutation_start_idx": 0,
        "protein_sequence_mutation_end_idx": 1,
        "protein_sequence_gene_names": "RNA_GENE_A;RNA_GENE_B",
        "protein_sequence_gene_ids": "ENSGRNA1;ENSGRNA2",
        "protein_sequence_transcript_names": "RNA-201;RNA-202",
        "protein_sequence_transcript_ids": "ENSTRNA1;ENSTRNA2",
        "num_ref_fragments": 8,
        "num_other_fragments": 1,
        "figure_source_url": "https://example.org/variant",
        "reference_cdna_sequence": "AACCGGTTAAACCCGGG",
        "reference_cdna_variant_start": 8,
        "reference_cdna_variant_end": 9,
        "annotation_cdna_sequence": "AACCGGTTAATCCCGGG",
        "annotation_cdna_variant_start": 8,
        "annotation_cdna_variant_end": 9,
        "rna_assembled_cdna_sequence": "AACCGGTTAATCCCGGG",
        "rna_assembled_cdna_variant_start": 8,
        "rna_assembled_cdna_variant_end": 9,
    }
    return [
        {
            **shared,
            "variant": "chr1 g.10_14delAAKPK",
            "protein_sequence": "VVKPK",
            "num_alt_fragments": 11,
            "num_fragments_supporting_top_protein_sequence": 10,
            "figure_label": "CONFIRMED1",
        },
        {
            **shared,
            "variant": "chr1 g.20_24delAAKPK",
            "protein_sequence": "VIKPK",
            "num_alt_fragments": 7,
            "num_fragments_supporting_top_protein_sequence": 6,
            "figure_label": "REFINED1",
        },
        {
            **shared,
            "variant": "chr1 g.30_34delAAKPK",
            "protein_sequence": None,
            "protein_sequence_gene_names": None,
            "protein_sequence_gene_ids": None,
            "protein_sequence_transcript_names": None,
            "protein_sequence_transcript_ids": None,
            "num_alt_fragments": 0,
            "num_fragments_supporting_top_protein_sequence": 0,
            "figure_label": "WITHHELD1",
        },
    ]


def test_global_alignment_preserves_repeated_residues():
    left, right = global_alignment("AAKPKVVKPK", "VVKPK")
    assert left.replace("–", "") == "AAKPKVVKPK"
    assert right.replace("–", "") == "VVKPK"
    assert len(left) == len(right)


@pytest.mark.parametrize("index,state", [
    (0, "CONFIRMED"), (1, "REFINED"), (2, "WITHHELD")])
def test_svg_classifies_rna_assembly_outcome(index, state):
    svg = render_mutation_svg(example_records()[index])
    ElementTree.fromstring(svg)
    assert '<rect width="1200" height="700" fill="#ffffff"/>' in svg
    assert state in svg
    assert "sequence evidence, not a clinical recommendation" in svg
    assert "https://example.org/variant" in svg
    assert "ANNOTATION SOURCE" in svg
    assert "RNA-ASSEMBLY SOURCE" in svg


def test_missing_protein_with_alt_fragments_does_not_claim_zero_alt():
    record = {**example_records()[2], "num_alt_fragments": 4}
    svg = render_mutation_svg(record)
    assert "4 alternate RNA fragments observed; no protein assembled" in svg
    assert "No alternate RNA fragments" not in svg


def test_transcript_svg_compares_annotation_and_assembly():
    svg = render_transcript_svg(example_records()[0])
    ElementTree.fromstring(svg)
    assert "transcript nucleotide context" in svg
    assert "Reference transcript" in svg
    assert "RNA assembled" in svg
    assert "5-prime to 3-prime" in svg


def test_unassessed_rna_is_not_presented_as_zero_alt_evidence():
    record = {
        **example_records()[2],
        "figure_rna_status": "not_assessed",
        "rna_assembled_cdna_sequence": None,
    }
    svg = render_transcript_svg(record)
    assert "ANNOTATION ONLY" in svg
    assert "RNA evidence not assessed" in svg
    assert "No alternate RNA fragments" not in svg
    assert "ALT fragments  0" not in svg


def test_figure_keeps_annotation_and_multi_source_assembly_provenance_separate():
    svg = render_mutation_svg(example_records()[0])
    assert "genes  REPEAT1" in svg
    assert "genes  RNA_GENE_A;RNA_GENE_B" in svg
    assert "transcripts  ENSTRNA1;ENSTRNA2" in svg


def test_generate_timestamped_figure_run(tmp_path):
    input_csv = tmp_path / "isovar.csv"
    pd.DataFrame(example_records()).to_csv(input_csv, index=False)
    run = generate_mutation_figures(
        input_csv,
        tmp_path / "figures",
        variants=["refined1", "withheld1"],
        formats=("svg",),
        timestamp="2026-09-15T220000Z",
    )
    assert run.name == "2026-09-15T220000Z"
    manifest = json.loads((run / "manifest.json").read_text())
    assert manifest["created_utc"] == "2026-09-15T220000Z"
    assert [figure["label"] for figure in manifest["figures"]] == [
        "REFINED1", "WITHHELD1"]
    for figure in manifest["figures"]:
        directory = run / figure["directory"]
        assert (directory / "protein-context.svg").is_file()
        assert (directory / "transcript-context.svg").is_file()
        record = json.loads((directory / "record.json").read_text())
        assert record["figure_label"] == figure["label"]
        assert record["displayed_sequences"]["annotation_only"] == "VVKPK"
        assert record["assembly_outcome"] in {"refined", "withheld"}
        if record["assembly_outcome"] == "refined":
            assert record["protein_sequence_gene_names"] == "RNA_GENE_A;RNA_GENE_B"
            assert record["protein_sequence_transcript_ids"] == "ENSTRNA1;ENSTRNA2"
            assert record["annotation_source"]["gene_names"] == "REPEAT1"
            assert record["rna_assembly_source"] == {
                "gene_ids": "ENSGRNA1;ENSGRNA2",
                "gene_names": "RNA_GENE_A;RNA_GENE_B",
                "transcript_ids": "ENSTRNA1;ENSTRNA2",
                "transcript_names": "RNA-201;RNA-202",
            }

    with pytest.raises(FileExistsError, match="already exists"):
        generate_mutation_figures(
            input_csv, tmp_path / "figures", formats=("svg",),
            timestamp="2026-09-15T220000Z")


def test_command_prints_generated_directory(tmp_path, capsys):
    input_csv = tmp_path / "isovar.csv"
    pd.DataFrame(example_records()[:1]).to_csv(input_csv, index=False)
    output_root = tmp_path / "out"
    main([
        str(input_csv), "--output-root", str(output_root), "--format", "svg",
        "--timestamp", "2026-09-15T221500Z",
    ])
    assert capsys.readouterr().out.strip() == str(output_root / "2026-09-15T221500Z")


def test_high_resolution_png_is_three_times_svg_dimensions(tmp_path):
    input_csv = tmp_path / "isovar.csv"
    pd.DataFrame(example_records()[:1]).to_csv(input_csv, index=False)
    run = generate_mutation_figures(
        input_csv,
        tmp_path / "figures",
        formats=("png",),
        timestamp="2026-09-15T223000Z",
    )
    png, = run.glob("*/protein-context.png")
    data = png.read_bytes()
    assert data[:8] == b"\x89PNG\r\n\x1a\n"
    assert struct.unpack(">II", data[16:24]) == (3600, 2100)
    manifest = json.loads((run / "manifest.json").read_text())
    assert manifest["png_scale"] == 3
    assert manifest["png_dimensions"] == {"width": 3600, "height": 2100}


def test_missing_required_columns_are_reported(tmp_path):
    input_csv = tmp_path / "invalid.csv"
    pd.DataFrame([{"variant": "chr1 g.1A>T"}]).to_csv(input_csv, index=False)
    with pytest.raises(ValueError, match="missing required columns"):
        generate_mutation_figures(input_csv, tmp_path / "out", formats=("svg",))


def test_failed_render_does_not_leave_partial_run(tmp_path, monkeypatch):
    input_csv = tmp_path / "isovar.csv"
    pd.DataFrame(example_records()[:1]).to_csv(input_csv, index=False)

    def fail_pdf(_svg, _path):
        raise RuntimeError("PDF backend failed")

    monkeypatch.setattr(mutation_visualization, "_write_pdf", fail_pdf)
    output_root = tmp_path / "figures"
    with pytest.raises(RuntimeError, match="PDF backend failed"):
        generate_mutation_figures(
            input_csv, output_root, formats=("svg", "pdf"),
            timestamp="2026-09-15T230000Z")
    assert list(output_root.iterdir()) == []
