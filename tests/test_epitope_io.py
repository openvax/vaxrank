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

"""
Tests for epitope_io: serialization, deserialization, and external imports.
"""

import os

import pytest

from vaxrank.epitope_prediction import EpitopePrediction
from vaxrank.epitope_io import (
    predictions_to_dataframe,
    save_predictions,
    load_predictions,
    load_pvacseq,
    load_lens,
    write_neoepitope_report,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "data", "epitope_fixtures")


def _make_prediction(**kwargs):
    defaults = dict(
        allele="HLA-A*02:01",
        peptide_sequence="SIINFEKL",
        wt_peptide_sequence="SIINFEKL",
        ic50=50.0,
        wt_ic50=5000.0,
        percentile_rank=0.5,
        prediction_method_name="test",
        overlaps_mutation=True,
        source_sequence="XXSIINFEKLXX",
        offset=2,
        occurs_in_reference=False,
    )
    defaults.update(kwargs)
    return EpitopePrediction(**defaults)


# ── Roundtrip ser/de ─────────────────────────────────────────────────────────

def test_predictions_to_dataframe():
    preds = [_make_prediction(), _make_prediction(allele="HLA-B*07:02", ic50=200.0)]
    df = predictions_to_dataframe(preds)
    assert len(df) == 2
    assert "allele" in df.columns
    assert "peptide_sequence" in df.columns


def test_save_and_load_csv_roundtrip(tmp_path):
    preds = [
        _make_prediction(),
        _make_prediction(allele="HLA-B*07:02", ic50=200.0, wt_ic50=None),
    ]
    path = tmp_path / "test.csv"
    save_predictions(preds, path)
    loaded = load_predictions(path)
    assert len(loaded) == 2
    assert loaded[0].allele == "HLA-A*02:01"
    assert loaded[0].ic50 == 50.0
    assert loaded[1].wt_ic50 is None


def test_save_and_load_tsv_roundtrip(tmp_path):
    preds = [_make_prediction()]
    path = tmp_path / "test.tsv"
    save_predictions(preds, path)
    loaded = load_predictions(path)
    assert len(loaded) == 1
    assert loaded[0].peptide_sequence == "SIINFEKL"


def test_roundtrip_preserves_values(tmp_path):
    original = _make_prediction(
        percentile_rank=1.23, occurs_in_reference=True, offset=5)
    path = tmp_path / "test.csv"
    save_predictions([original], path)
    loaded = load_predictions(path)[0]
    assert loaded.percentile_rank == pytest.approx(1.23)
    assert loaded.occurs_in_reference is True
    assert loaded.offset == 5


# ── pVACseq import ───────────────────────────────────────────────────────────

def test_load_pvacseq():
    path = os.path.join(DATA_DIR, "pvacseq_example.tsv")
    report_df, preds = load_pvacseq(path)
    assert len(preds) == 3
    assert len(report_df) == 3

    # Check first prediction
    p = preds[0]
    assert p.allele == "HLA-A*02:01"
    assert p.peptide_sequence == "SVVGSSSSS"
    assert p.ic50 == pytest.approx(76.11)
    assert p.wt_ic50 == pytest.approx(4520.30)
    assert p.percentile_rank == pytest.approx(0.5)
    assert p.overlaps_mutation is True

    # Check report DataFrame has expected columns
    assert "Gene name" in report_df.columns
    assert "Genomic variant" in report_df.columns
    assert "Tier" in report_df.columns

    # Check third prediction has occurs_in_reference=True (Ref Match=True)
    assert preds[2].occurs_in_reference is True


def test_load_pvacseq_missing_columns(tmp_path):
    bad_file = tmp_path / "bad.tsv"
    bad_file.write_text("Gene\tAA Change\nTP53\tG245S\n")
    with pytest.raises(ValueError, match="missing required columns"):
        load_pvacseq(bad_file)


# ── LENS import ──────────────────────────────────────────────────────────────

def test_load_lens_auto_prefers_mhcflurry():
    """Auto predictor should pick mhcflurry (canonical) when both are present."""
    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, preds = load_lens(path)
    assert len(preds) == 3
    assert len(report_df) == 3

    p = preds[0]
    # Allele should be normalized: LENS 'HLA-A02:01' -> 'HLA-A*02:01'
    assert p.allele == "HLA-A*02:01"
    assert p.peptide_sequence == "SVVGSSSSS"
    # mhcflurry affinity (95.4), not netmhcpan's 76.11
    assert p.ic50 == pytest.approx(95.4)
    # mhcflurry presentation percentile
    assert p.percentile_rank == pytest.approx(0.28)
    assert p.prediction_method_name == "mhcflurry"
    assert p.wt_peptide_sequence == ""
    assert p.wt_ic50 is None
    assert p.overlaps_mutation is True
    # pep_context is carried in source_sequence
    assert p.source_sequence == "AASVVGSSSSSGTR"

    # Report columns
    assert "Antigen source" in report_df.columns
    assert "Agretopicity" in report_df.columns
    assert "Predictors used" in report_df.columns
    # Per-predictor raw-value columns for every detected predictor
    assert "mhcflurry value (nM)" in report_df.columns
    assert "netmhcpan value (nM)" in report_df.columns
    assert report_df["Predictors used"].iloc[0] == "mhcflurry"
    assert report_df["mhcflurry value (nM)"].iloc[0] == pytest.approx(95.4)
    assert report_df["netmhcpan value (nM)"].iloc[0] == pytest.approx(76.11)


def test_load_lens_all_emits_per_predictor():
    """predictor='all' emits N predictions per (peptide, allele) pair."""
    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, preds = load_lens(path, predictor="all")
    # Fixture has 3 rows × 2 predictors = 6 predictions
    assert len(preds) == 6
    # Report DataFrame still has one row per (peptide, allele)
    assert len(report_df) == 3
    # Predictions alternate mhcflurry/netmhcpan (or in tool-sorted order)
    methods = {p.prediction_method_name for p in preds}
    assert methods == {"mhcflurry", "netmhcpan"}
    # Each (peptide, allele) pair should appear twice — once per predictor
    from collections import Counter
    pair_counts = Counter((p.peptide_sequence, p.allele) for p in preds)
    assert all(c == 2 for c in pair_counts.values())
    assert report_df["Predictors used"].iloc[0] == "mhcflurry,netmhcpan"


def test_load_lens_explicit_netmhcpan():
    """Explicit netmhcpan choice should use netmhcpan values."""
    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _, preds = load_lens(path, predictor="netmhcpan")
    p = preds[0]
    assert p.ic50 == pytest.approx(76.11)
    assert p.percentile_rank == pytest.approx(0.3)  # perc_rank_el
    assert p.prediction_method_name == "netmhcpan"


def test_load_lens_mhcflurry_only_v19():
    """A v1.9-dev style file (no netmhcpan cols) should load under auto."""
    path = os.path.join(DATA_DIR, "lens_v1.9_mhcflurry_only.tsv")
    report_df, preds = load_lens(path)
    assert len(preds) == 3

    # First row: ERV antigen, allele normalized, NA variant_coords falls through
    p = preds[0]
    assert p.allele == "HLA-A*01:01"
    assert p.peptide_sequence == "TAEFYQRY"
    assert p.ic50 == pytest.approx(254.15)
    assert p.prediction_method_name == "mhcflurry"
    # Antigen source captured
    assert report_df["Antigen source"].iloc[0] == "ERV"
    # When variant_coords / mut_aa_pos are NA, fall back to origin_descriptor
    assert "Hsap38" in report_df["Genomic variant"].iloc[0]
    # Only mhcflurry column present; no netmhcpan column at all
    assert "mhcflurry value (nM)" in report_df.columns
    assert "netmhcpan value (nM)" not in report_df.columns


def test_load_lens_explicit_netmhcpan_missing_raises():
    """Asking for netmhcpan on a mhcflurry-only file must raise."""
    path = os.path.join(DATA_DIR, "lens_v1.9_mhcflurry_only.tsv")
    with pytest.raises(ValueError, match="does not contain predictor 'netmhcpan'"):
        load_lens(path, predictor="netmhcpan")


def test_load_lens_missing_required_columns(tmp_path):
    bad_file = tmp_path / "bad.tsv"
    bad_file.write_text("gene_name\tpos\nTP53\t5\n")
    with pytest.raises(ValueError, match="missing required columns"):
        load_lens(bad_file)


def test_load_lens_no_predictor_columns(tmp_path):
    """File has peptide/allele but no recognized affinity column."""
    bad_file = tmp_path / "bad.tsv"
    bad_file.write_text(
        "peptide\tallele\tgene_name\nSIINFEKL\tHLA-A*02:01\tOVA\n")
    with pytest.raises(ValueError, match="no recognized predictor"):
        load_lens(bad_file)


def test_load_lens_detects_netmhcstabpan():
    """netmhcstabpan columns are detected as a pMHC_stability predictor."""
    path = os.path.join(DATA_DIR, "lens_v1.4_with_stability.tsv")
    report_df, preds = load_lens(path, predictor="all")
    methods = {p.prediction_method_name for p in preds}
    assert methods == {"mhcflurry", "netmhcpan", "netmhcstabpan"}
    # Stability value column is labeled with its unit
    assert "netmhcstabpan value (hours)" in report_df.columns


# ── DSL integration for external inputs ──────────────────────────────────────

def test_lens_dsl_combines_both_predictors(tmp_path):
    """A score_expr combining two predictors' affinities should average them."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, preds = load_lens(path, predictor="all")

    # 50/50 blend of the two predictors' logistic-normalized affinity scores
    cfg = EpitopeConfig(
        score_expr=(
            "affinity['mhcflurry'].logistic(350, 150) * 0.5 + "
            "affinity['netmhcpan'].logistic(350, 150) * 0.5"
        )
    )
    csv_path = tmp_path / "report.csv"
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path), epitope_config=cfg)

    import pandas as pd
    result = pd.read_csv(csv_path)
    assert len(result) == 3
    # Non-zero score means the DSL actually picked up both predictors
    assert (result["vaxrank_score"] > 0).any()


def test_lens_dsl_validation_catches_missing_predictor():
    """Formula referencing a predictor not in the file should error clearly."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_v1.9_mhcflurry_only.tsv")
    report_df, preds = load_lens(path)

    cfg = EpitopeConfig(
        score_expr="affinity['netmhcpan'].logistic(350, 150)"
    )
    with pytest.raises(ValueError, match="predictor 'netmhcpan'"):
        write_neoepitope_report(
            report_df, preds, csv_report_path="/tmp/_unused.csv",
            epitope_config=cfg)


def test_lens_default_scoring_matches_legacy(tmp_path):
    """With no filter_expr/score_expr, DSL-routed scoring should match the
    legacy EpitopePrediction.logistic_epitope_score output byte-for-byte."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, preds = load_lens(path)

    # Compute legacy scores by hand
    cfg = EpitopeConfig()
    legacy_scores = [
        round(p.logistic_epitope_score(
            midpoint=cfg.logistic_epitope_score_midpoint,
            width=cfg.logistic_epitope_score_width,
            ic50_cutoff=cfg.binding_affinity_cutoff,
            scoring_mode=cfg.scoring_mode,
            percentile_rank_cutoff=cfg.percentile_rank_cutoff,
        ), 6)
        for p in preds
    ]

    csv_path = tmp_path / "report.csv"
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path), epitope_config=cfg)

    import pandas as pd
    result = pd.read_csv(csv_path)
    # Re-order result rows back to prediction order so we can compare
    # element-wise (write_neoepitope_report sorts by score desc).
    original_order = [
        (p.peptide_sequence, p.allele) for p in preds
    ]
    result_scores_by_key = dict(zip(
        zip(result["Mutant peptide sequence"], result["Allele"]),
        result["vaxrank_score"],
    ))
    new_scores = [result_scores_by_key[key] for key in original_order]
    for legacy, new in zip(legacy_scores, new_scores):
        assert new == pytest.approx(legacy)


# ── Report generation from imports ───────────────────────────────────────────

def test_write_neoepitope_report_csv(tmp_path):
    path = os.path.join(DATA_DIR, "pvacseq_example.tsv")
    report_df, preds = load_pvacseq(path)
    csv_path = tmp_path / "report.csv"
    write_neoepitope_report(report_df, preds, csv_report_path=str(csv_path))
    assert csv_path.exists()

    import pandas as pd
    result = pd.read_csv(csv_path)
    assert "rank" in result.columns
    assert "vaxrank_score" in result.columns
    assert len(result) == 3
    # Should be sorted by score descending
    assert result["vaxrank_score"].iloc[0] >= result["vaxrank_score"].iloc[1]


def test_write_neoepitope_report_xlsx(tmp_path):
    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, preds = load_lens(path)
    xlsx_path = tmp_path / "report.xlsx"
    write_neoepitope_report(report_df, preds, excel_report_path=str(xlsx_path))
    assert xlsx_path.exists()


# ── Scoring modes ────────────────────────────────────────────────────────────

def test_scoring_mode_affinity():
    p = _make_prediction(ic50=100.0)
    score = p.logistic_epitope_score(scoring_mode="affinity")
    assert score > 0.5  # Strong binder


def test_scoring_mode_percentile_rank():
    p = _make_prediction(percentile_rank=0.5)
    score = p.logistic_epitope_score(scoring_mode="percentile_rank")
    assert score == pytest.approx(0.95)


def test_scoring_mode_percentile_rank_weak():
    p = _make_prediction(percentile_rank=10.0)
    score = p.logistic_epitope_score(scoring_mode="percentile_rank")
    assert score == 0.0


def test_scoring_mode_percentile_rank_none():
    p = _make_prediction(percentile_rank=None)
    score = p.logistic_epitope_score(scoring_mode="percentile_rank")
    assert score == 0.0


# ── CLI integration test ─────────────────────────────────────────────────────

def test_pvacseq_cli_csv_output(tmp_path):
    """Test the --input-pvacseq flag produces CSV output."""
    from vaxrank.cli.entry_point import main
    pvacseq_path = os.path.join(DATA_DIR, "pvacseq_example.tsv")
    csv_path = str(tmp_path / "output.csv")
    main([
        "--input-pvacseq", pvacseq_path,
        "--output-csv", csv_path,
    ])
    assert os.path.exists(csv_path)
    import pandas as pd
    df = pd.read_csv(csv_path)
    assert len(df) == 3
    assert "vaxrank_score" in df.columns


def test_lens_cli_csv_output(tmp_path):
    """Test the --input-lens flag produces CSV output."""
    from vaxrank.cli.entry_point import main
    lens_path = os.path.join(DATA_DIR, "lens_example.tsv")
    csv_path = str(tmp_path / "output.csv")
    main([
        "--input-lens", lens_path,
        "--output-csv", csv_path,
    ])
    assert os.path.exists(csv_path)
    import pandas as pd
    df = pd.read_csv(csv_path)
    assert len(df) == 3
