import json
from xml.etree import ElementTree

import pandas as pd
import pytest

import vaxrank.mutation_visualization as mutation_visualization
from vaxrank.cli.mutation_visualization import main
from vaxrank.mutation_visualization import (
    generate_mutation_figures,
    global_alignment,
    render_mutation_svg,
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
        "num_ref_fragments": 8,
        "num_other_fragments": 1,
        "figure_source_url": "https://example.org/variant",
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
        record = json.loads((directory / "record.json").read_text())
        assert record["figure_label"] == figure["label"]
        assert record["displayed_sequences"]["annotation_only"] == "VVKPK"
        assert record["assembly_outcome"] in {"refined", "withheld"}

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
