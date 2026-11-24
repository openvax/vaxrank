# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests for ``vaxrank compare`` run-to-run diff CLI (openvax/vaxrank#920)."""

import csv
import json
from pathlib import Path

import pytest

from vaxrank.cli.entry_point import main as vaxrank_main
from vaxrank.compare import (
    ComparePolicy,
    _kendall_tau_b,
    _spearman_rho,
    compare_runs,
    load_run_bundle,
    render_comparison_ascii,
    render_comparison_html,
    write_comparison_outputs,
)
from vaxrank.config import load_vaxrank_config


def _write_ranked_csv(path: Path, rows: list[dict]) -> None:
    fieldnames = [
        "source_rank", "window_rank", "source_kind", "variant_types",
        "rna_product", "source_identifier", "gene_name", "amino_acids",
        "chr", "pos", "ref", "alt", "combined_score", "expression_score",
        "target_epitope_score", "biological_evidence", "cta_targeting",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            full = {k: "" for k in fieldnames}
            full.update(r)
            writer.writerow(full)


def _write_run_bundle(
    root: Path,
    rows: list[dict],
    *,
    peptide_manifest: list[dict] | None = None,
    mrna_manifest: list[dict] | None = None,
    selection_policy: dict | None = None,
) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    _write_ranked_csv(root / "ranked_candidates.csv", rows)
    if peptide_manifest is not None:
        (root / "peptide").mkdir(parents=True, exist_ok=True)
        (root / "peptide" / "manifest.json").write_text(
            json.dumps(peptide_manifest, indent=2), encoding="utf-8"
        )
    if mrna_manifest is not None:
        (root / "mrna").mkdir(parents=True, exist_ok=True)
        (root / "mrna" / "manifest.json").write_text(
            json.dumps(mrna_manifest, indent=2), encoding="utf-8"
        )
    if selection_policy is not None:
        (root / "selection_policy.json").write_text(
            json.dumps(selection_policy, indent=2), encoding="utf-8"
        )
    return root


def test_compare_policy_validation_and_yaml(tmp_path):
    pol = ComparePolicy()
    assert pol.top_k_values == (5, 10, 20)
    assert pol.primary_top_k == 10
    assert pol.include_pdf is False

    # primary_top_k automatically included in top_k_values
    pol2 = ComparePolicy(top_k_values=(3, 5), primary_top_k=7)
    assert pol2.top_k_values == (3, 5, 7)

    with pytest.raises(ValueError, match="top_k_values"):
        ComparePolicy(top_k_values=())
    with pytest.raises(ValueError, match="expression_fold_change_threshold"):
        ComparePolicy(expression_fold_change_threshold=1.0)
    with pytest.raises(ValueError, match="Unknown compare policy"):
        ComparePolicy.resolve({"bogus": 1})

    cfg_path = tmp_path / "compare_policy.yaml"
    cfg_path.write_text(
        "compare:\n"
        "  top_k_values: [3, 5, 15]\n"
        "  primary_top_k: 5\n"
        "  score_relative_change_threshold: 0.15\n"
        "  expression_fold_change_threshold: 2.0\n"
        "  rank_shift_highlight_threshold: 2\n"
    )
    loaded = load_vaxrank_config(config_path=str(cfg_path))
    resolved = ComparePolicy.resolve(loaded["compare"])
    assert resolved.top_k_values == (3, 5, 15)
    assert resolved.primary_top_k == 5
    assert resolved.score_relative_change_threshold == 0.15
    assert resolved.expression_fold_change_threshold == 2.0
    assert resolved.rank_shift_highlight_threshold == 2


def test_rank_correlation_helpers():
    assert _spearman_rho([1, 2, 3, 4], [1, 2, 3, 4]) == pytest.approx(1.0)
    assert _spearman_rho([1, 2, 3, 4], [4, 3, 2, 1]) == pytest.approx(-1.0)
    assert _kendall_tau_b([1, 2, 3, 4], [1, 2, 3, 4]) == pytest.approx(1.0)
    assert _kendall_tau_b([1, 2, 3, 4], [4, 3, 2, 1]) == pytest.approx(-1.0)
    assert _spearman_rho([1], [1]) is None
    assert _kendall_tau_b([1], [1]) is None


def test_compare_runs_attribution_and_construct_diffs(tmp_path):
    run_a = _write_run_bundle(
        tmp_path / "dna_only",
        rows=[
            {
                "source_rank": 1, "window_rank": 1, "source_kind": "variant",
                "source_identifier": "TP53_p.R175H", "gene_name": "TP53",
                "amino_acids": "EVVRHCPHHERCSDSDGLAPPQHLI",
                "chr": "17", "pos": 7675088, "ref": "C", "alt": "T",
                "combined_score": 12.0, "expression_score": 4.0,
                "target_epitope_score": 3.0,
                "biological_evidence": "DNA-only proxy expression",
            },
            {
                "source_rank": 2, "window_rank": 1, "source_kind": "variant",
                "source_identifier": "FSGENE_p.L10fs", "gene_name": "FSGENE",
                "amino_acids": "PSWSLLDLKTGSCTSRWESAWSGLS",
                "chr": "1", "pos": 1000, "ref": "CA", "alt": "C",
                "combined_score": 9.0, "expression_score": 3.0,
                "target_epitope_score": 3.0,
                "biological_evidence": "DNA-only proxy expression",
            },
            {
                "source_rank": 3, "window_rank": 1, "source_kind": "variant",
                "source_identifier": "UNEXPRESSED_p.A1V", "gene_name": "UNEXPRESSED",
                "amino_acids": "MKLPVQRSTVWYACDEFGHIKLMNP",
                "chr": "2", "pos": 2000, "ref": "G", "alt": "A",
                "combined_score": 6.0, "expression_score": 2.0,
                "target_epitope_score": 3.0,
                "biological_evidence": "DNA-only proxy expression",
            },
            # Secondary window row (window_rank=2) must be ignored for primary target rank
            {
                "source_rank": 1, "window_rank": 2, "source_kind": "variant",
                "source_identifier": "TP53_p.R175H", "gene_name": "TP53",
                "amino_acids": "VVRHCPHHERCSDSDGLAPPQHLIR",
                "chr": "17", "pos": 7675088, "ref": "C", "alt": "T",
                "combined_score": 11.5, "expression_score": 4.0,
                "target_epitope_score": 2.8,
            },
        ],
        peptide_manifest=[
            {"modality": "peptide", "name": "peptide_001", "sequence": "EVVRHCPHHERCSDSDGLAPPQHLI",
             "antigen_names": ["TP53_p.R175H"], "pool": 1},
            {"modality": "peptide", "name": "peptide_002", "sequence": "PSWSLLDLKTGSCTSRWESAWSGLS",
             "antigen_names": ["FSGENE_p.L10fs"], "pool": 2},
            {"modality": "peptide", "name": "peptide_003", "sequence": "MKLPVQRSTVWYACDEFGHIKLMNP",
             "antigen_names": ["UNEXPRESSED_p.A1V"], "pool": 1},
        ],
        mrna_manifest=[
            {
                "modality": "mrna", "name": "mrna_001",
                "antigen_names": ["TP53_p.R175H", "FSGENE_p.L10fs"],
                "antigens": [
                    {"name": "TP53_p.R175H", "aa": "EVVRHCPHHERCSDSDGLAPPQHLI", "gene_names": ["TP53"]},
                    {"name": "FSGENE_p.L10fs", "aa": "PSWSLLDLKTGSCTSRWESAWSGLS", "gene_names": ["FSGENE"]},
                ],
            }
        ],
        selection_policy={"name": "openvax-v1", "sha256": "aaa111"},
    )

    run_b = _write_run_bundle(
        tmp_path / "rna_backed",
        rows=[
            {
                "source_rank": 1, "window_rank": 1, "source_kind": "variant",
                "source_identifier": "TP53_p.R175H", "gene_name": "TP53",
                "amino_acids": "EVVRHCPHHERCSDSDGLAPPQHLI",
                "chr": "chr17", "pos": 7675088, "ref": "C", "alt": "T",
                "combined_score": 18.0, "expression_score": 6.0,
                "target_epitope_score": 3.0,
                "biological_evidence": "RNA-supported 36 alt reads",
                "rna_product": "translated",
            },
            {
                "source_rank": 2, "window_rank": 1, "source_kind": "variant",
                "source_identifier": "KRAS_p.G12D", "gene_name": "KRAS",
                "amino_acids": "MTEYKLVVVGADGVGKSALTIQLIQ",
                "chr": "chr12", "pos": 25245350, "ref": "C", "alt": "T",
                "combined_score": 14.0, "expression_score": 5.0,
                "target_epitope_score": 2.8,
                "biological_evidence": "RNA-supported 25 alt reads",
                "rna_product": "translated",
            },
            {
                "source_rank": 3, "window_rank": 1, "source_kind": "variant",
                "source_identifier": "FSGENE_p.L10fs", "gene_name": "FSGENE",
                "amino_acids": "GPHDLLWAHSGRDRPWSGKESTGAI",
                "chr": "chr1", "pos": 1000, "ref": "CA", "alt": "C",
                "combined_score": 1.2, "expression_score": 4.8,
                "target_epitope_score": 0.25,
                "biological_evidence": "RNA-supported 23 alt reads",
                "rna_product": "translated",
            },
        ],
        peptide_manifest=[
            {"modality": "peptide", "name": "peptide_001", "sequence": "EVVRHCPHHERCSDSDGLAPPQHLI",
             "antigen_names": ["TP53_p.R175H"], "pool": 1},
            {"modality": "peptide", "name": "peptide_002", "sequence": "MTEYKLVVVGADGVGKSALTIQLIQ",
             "antigen_names": ["KRAS_p.G12D"], "pool": 2},
            {"modality": "peptide", "name": "peptide_003", "sequence": "GPHDLLWAHSGRDRPWSGKESTGAI",
             "antigen_names": ["FSGENE_p.L10fs"], "pool": 1},
        ],
        mrna_manifest=[
            {
                "modality": "mrna", "name": "mrna_001",
                "antigen_names": ["TP53_p.R175H", "KRAS_p.G12D", "FSGENE_p.L10fs"],
                "antigens": [
                    {"name": "TP53_p.R175H", "aa": "EVVRHCPHHERCSDSDGLAPPQHLI", "gene_names": ["TP53"]},
                    {"name": "KRAS_p.G12D", "aa": "MTEYKLVVVGADGVGKSALTIQLIQ", "gene_names": ["KRAS"]},
                    {"name": "FSGENE_p.L10fs", "aa": "GPHDLLWAHSGRDRPWSGKESTGAI", "gene_names": ["FSGENE"]},
                ],
            }
        ],
        selection_policy={"name": "tesla-recognition-v1", "sha256": "bbb222"},
    )

    comp = compare_runs(
        run_a, run_b,
        policy=ComparePolicy(top_k_values=(2, 3), primary_top_k=3),
        label_a="DNA-only",
        label_b="RNA-backed",
    )
    assert comp["schema"] == "vaxrank.comparison.v1"
    assert comp["policies_differ"] is True
    rp = comp["rank_preservation"]
    assert rp["shared_target_count"] == 2
    assert rp["only_in_a_count"] == 1
    assert rp["only_in_b_count"] == 1
    assert rp["spearman_rho"] == pytest.approx(1.0)
    assert rp["kendall_tau_b"] == pytest.approx(1.0)
    assert rp["top_k"]["3"]["overlap_count"] == 2
    assert rp["top_k"]["3"]["exact_sequence_matches"] == 1
    assert rp["top_k"]["3"]["shifted_window_matches"] == 1

    # Check attributions
    by_key = {row["target_key"]: row for row in comp["targets"]}
    fs_row = by_key["1:1000:CA>C"]
    assert fs_row["sequence_match"] == "shifted"
    assert "peptide_window_shift" in fs_row["attribution"]
    assert "expression_source_change_rna_vs_proxy" in fs_row["attribution"]
    assert "epitope_score_shift" in fs_row["attribution"]

    assert [r["gene_name"] for r in comp["gained_in_top_k"]] == ["KRAS"]
    assert "gained_with_rna_or_filter_inclusion" in comp["gained_in_top_k"][0]["attribution"]
    assert [r["gene_name"] for r in comp["lost_from_top_k"]] == ["UNEXPRESSED"]
    assert "dropped_without_rna_support" in comp["lost_from_top_k"][0]["attribution"]

    # Construct diffs
    pep_diff = comp["construct_diffs"]["peptide"]
    assert pep_diff["exact_match_count"] == 1
    assert pep_diff["shifted_window_count"] == 1
    assert pep_diff["added_count"] == 1
    assert pep_diff["dropped_count"] == 1

    mrna_diff = comp["construct_diffs"]["mrna"]
    assert mrna_diff["exact_match_count"] == 1
    assert mrna_diff["shifted_window_count"] == 1
    assert mrna_diff["added_count"] == 1
    assert mrna_diff["dropped_count"] == 0

    # Renderers + writers
    ascii_txt = render_comparison_ascii(comp)
    assert "Vaxrank Run Comparison" in ascii_txt
    assert "KRAS" in ascii_txt
    assert "UNEXPRESSED" in ascii_txt

    html_txt = render_comparison_html(comp)
    assert "<!DOCTYPE html>" in html_txt
    assert "peptide_window_shift" in html_txt

    out_dir = tmp_path / "diff_out"
    written = write_comparison_outputs(comp, out_dir)
    for kind in ("json", "csv", "ascii", "html"):
        assert Path(written[kind]).is_file()
    loaded_json = json.loads(Path(written["json"]).read_text())
    assert loaded_json["schema"] == "vaxrank.comparison.v1"


def test_compare_cli_via_vaxrank_main(tmp_path, capsys):
    csv_a = tmp_path / "run_a.csv"
    csv_b = tmp_path / "run_b.csv"
    _write_ranked_csv(csv_a, [
        {"source_rank": 1, "window_rank": 1, "source_identifier": "G1_p.A1V",
         "gene_name": "G1", "amino_acids": "AAAAAAAAAAAAAAAAAAAAAAAAA",
         "combined_score": 10.0, "expression_score": 2.0, "target_epitope_score": 5.0},
        {"source_rank": 2, "window_rank": 1, "source_identifier": "G2_p.B2L",
         "gene_name": "G2", "amino_acids": "BBBBBBBBBBBBBBBBBBBBBBBBB",
         "combined_score": 8.0, "expression_score": 2.0, "target_epitope_score": 4.0},
    ])
    _write_ranked_csv(csv_b, [
        {"source_rank": 1, "window_rank": 1, "source_identifier": "G2_p.B2L",
         "gene_name": "G2", "amino_acids": "BBBBBBBBBBBBBBBBBBBBBBBBB",
         "combined_score": 15.0, "expression_score": 5.0, "target_epitope_score": 3.0},
        {"source_rank": 2, "window_rank": 1, "source_identifier": "G1_p.A1V",
         "gene_name": "G1", "amino_acids": "AAAAAAAAAAAAAAAAAAAAAAAAA",
         "combined_score": 10.0, "expression_score": 2.0, "target_epitope_score": 5.0},
    ])

    out_dir = tmp_path / "cli_out"
    vaxrank_main([
        "compare", str(csv_a), str(csv_b),
        "--label-a", "Baseline", "--label-b", "Candidate",
        "--top-k", "1,2", "--primary-top-k", "1",
        "--output-dir", str(out_dir),
    ])
    stdout = capsys.readouterr().out
    assert "Baseline" in stdout
    assert "Candidate" in stdout
    assert (out_dir / "comparison_summary.json").is_file()
    assert (out_dir / "comparable_targets.csv").is_file()
    assert (out_dir / "comparison_report.html").is_file()

    # JSON stdout mode without --output-dir
    vaxrank_main(["compare", str(csv_a), str(csv_b), "--json"])
    payload = json.loads(capsys.readouterr().out)
    assert payload["rank_preservation"]["spearman_rho"] == pytest.approx(-1.0)


def test_load_run_bundle_errors_cleanly(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_run_bundle(tmp_path / "missing")
    empty = tmp_path / "empty_dir"
    empty.mkdir()
    with pytest.raises(ValueError, match="No ranked candidates CSV"):
        load_run_bundle(empty)
