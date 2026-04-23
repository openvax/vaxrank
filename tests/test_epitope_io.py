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


def test_lens_dsl_validation_catches_missing_predictor(tmp_path):
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
            report_df, preds, csv_report_path=str(tmp_path / "unused.csv"),
            epitope_config=cfg)


def test_lens_dsl_validation_catches_missing_kind(tmp_path):
    """Formula using a Kind not present in the data (e.g. stability on a
    file without a stability predictor) should error clearly."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_v1.9_mhcflurry_only.tsv")
    report_df, preds = load_lens(path)
    cfg = EpitopeConfig(score_expr="stability['netmhcstabpan'].value")
    with pytest.raises(ValueError, match="references kind 'pMHC_stability'"):
        write_neoepitope_report(
            report_df, preds, csv_report_path=str(tmp_path / "unused.csv"),
            epitope_config=cfg)


def test_lens_dsl_validation_catches_missing_version(tmp_path):
    """Formula pinned to a specific predictor version not in the data should
    error clearly."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_v1.9_mhcflurry_only.tsv")
    report_df, preds = load_lens(path)
    cfg = EpitopeConfig(
        score_expr="affinity['mhcflurry', '9.9.9'].logistic(350, 150)"
    )
    with pytest.raises(ValueError, match="predictor version '9.9.9'"):
        write_neoepitope_report(
            report_df, preds, csv_report_path=str(tmp_path / "unused.csv"),
            epitope_config=cfg)


def test_lens_dsl_version_qualified_formula(tmp_path):
    """A formula pinned to a specific predictor version should match
    predictions carrying that version from LENS column sniffing."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, preds = load_lens(path, predictor="all")
    # Every LENS prediction should carry predictor_version from columns.
    assert all(p.predictor_version for p in preds)
    mhcflurry_version = next(
        p.predictor_version for p in preds if p.prediction_method_name == "mhcflurry")
    assert mhcflurry_version == "2.1.1"

    cfg = EpitopeConfig(
        score_expr=(
            "affinity['mhcflurry', '2.1.1'].logistic(350, 150)"
        )
    )
    csv_path = tmp_path / "report.csv"
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    assert (result["vaxrank_score"] > 0).any()


# ── Integration tests against real-LENS subsamples ───────────────────────────
# Fixtures in real_lens_subsets/ are 30-row stratified samples from actual
# LENS reports (v1.4-dev, v1.5.1, v1.9-dev) covering the full column surface
# (94 / 106 / 108 cols). They exercise the real schema including ERV / CTA /
# FUSION / SPLICE / INDEL antigen sources that synthetic fixtures miss.

REAL_LENS_DIR = os.path.join(DATA_DIR, "real_lens_subsets")


def test_real_lens_v14_detects_all_three_predictors():
    """v1.4-dev files emit mhcflurry + netmhcpan + netmhcstabpan. All three
    should be detected and emit predictions under predictor='all'."""
    path = os.path.join(REAL_LENS_DIR, "lens_v1.4_real_subset.tsv")
    report_df, preds = load_lens(path, predictor="all")
    methods = {p.prediction_method_name for p in preds}
    assert methods == {"mhcflurry", "netmhcpan", "netmhcstabpan"}
    # Every real LENS prediction carries a nonempty version.
    assert all(p.predictor_version for p in preds)
    # Report exposes per-tool columns with correct units.
    assert "mhcflurry value (nM)" in report_df.columns
    assert "netmhcpan value (nM)" in report_df.columns
    assert "netmhcstabpan value (hours)" in report_df.columns


def test_real_lens_v15_auto_picks_mhcflurry():
    """v1.5.1 has both mhcflurry and netmhcpan. Auto should pick mhcflurry
    (canonical); all mode should emit 2× predictions per report row."""
    path = os.path.join(REAL_LENS_DIR, "lens_v1.5_real_subset.tsv")
    report_df_auto, preds_auto = load_lens(path)
    assert all(p.prediction_method_name == "mhcflurry" for p in preds_auto)
    assert report_df_auto["Predictors used"].iloc[0] == "mhcflurry"

    report_df_all, preds_all = load_lens(path, predictor="all")
    # Under 'all' each report row yields up to N predictions (N = number of
    # detected predictors); some rows have NA for a predictor and get
    # fewer. Under 'auto' a row is dropped entirely if the canonical
    # predictor's value is NA, so report_df_all may have MORE rows than
    # report_df_auto — rows where mhcflurry is NA but netmhcpan has a value.
    assert len(preds_all) >= len(preds_auto)
    assert len(preds_all) <= 2 * len(report_df_all)
    assert len(report_df_all) >= len(report_df_auto)


def test_real_lens_v19_is_mhcflurry_only():
    """v1.9-dev emits only mhcflurry (netMHCpan was dropped)."""
    path = os.path.join(REAL_LENS_DIR, "lens_v1.9_real_subset.tsv")
    report_df, preds = load_lens(path, predictor="all")
    methods = {p.prediction_method_name for p in preds}
    assert methods == {"mhcflurry"}
    # Real data has HLA alleles in the un-asterisked form; loader normalizes.
    assert all(p.allele.startswith("HLA-") and "*" in p.allele for p in preds)


def test_real_lens_dsl_scoring_all_versions(tmp_path):
    """DSL scoring should produce nonzero scores for at least some rows in
    each real LENS fixture (with all predictors emitted). A 50/50 blend of
    mhcflurry + netmhcpan where both present; plain mhcflurry otherwise."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    blend = EpitopeConfig(score_expr=(
        "affinity['mhcflurry'].logistic(350, 150) * 0.5 + "
        "affinity['netmhcpan'].logistic(350, 150) * 0.5"))
    solo = EpitopeConfig(score_expr="affinity['mhcflurry'].logistic(350, 150)")

    cases = [
        ("lens_v1.4_real_subset.tsv", blend),
        ("lens_v1.5_real_subset.tsv", blend),
        ("lens_v1.9_real_subset.tsv", solo),  # no netmhcpan
    ]

    import pandas as pd
    for name, cfg in cases:
        path = os.path.join(REAL_LENS_DIR, name)
        report_df, preds = load_lens(path, predictor="all")
        csv_path = tmp_path / f"{name}.csv"
        write_neoepitope_report(
            report_df, preds, csv_report_path=str(csv_path),
            epitope_config=cfg)
        result = pd.read_csv(csv_path)
        assert len(result) == len(report_df), f"{name}: row count changed"
        assert (result["vaxrank_score"] > 0).any(), (
            f"{name}: no rows scored > 0 under {cfg.score_expr!r}")


def test_real_lens_v14_stability_formula(tmp_path):
    """v1.4 real fixture has netmhcstabpan — verify combined
    affinity + stability formula evaluates successfully."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(REAL_LENS_DIR, "lens_v1.4_real_subset.tsv")
    report_df, preds = load_lens(path, predictor="all")
    cfg = EpitopeConfig(score_expr=(
        "affinity['mhcflurry'].logistic(350, 150) * 0.7 + "
        "(stability['netmhcstabpan'].value / 24).clip(0, 1) * 0.3"))
    csv_path = tmp_path / "report.csv"
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    assert (result["vaxrank_score"] > 0).any()


def test_real_lens_version_qualified_auto_detects(tmp_path):
    """predictor_version sniffed from real LENS column names matches what
    a version-qualified DSL formula expects."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(REAL_LENS_DIR, "lens_v1.5_real_subset.tsv")
    report_df, preds = load_lens(path, predictor="all")
    versions = {(p.prediction_method_name, p.predictor_version) for p in preds}
    # Real v1.5.1 file: mhcflurry 2.1.1 + netmhcpan 4.1b
    assert ("mhcflurry", "2.1.1") in versions
    assert ("netmhcpan", "4.1b") in versions

    cfg = EpitopeConfig(score_expr=(
        "affinity['mhcflurry', '2.1.1'].logistic(350, 150)"))
    csv_path = tmp_path / "report.csv"
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    assert (result["vaxrank_score"] > 0).any()


def test_lens_dsl_stability_and_affinity_combined(tmp_path):
    """Formula combining pMHC_affinity and pMHC_stability Kinds should
    evaluate when both predictors are emitted."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_v1.4_with_stability.tsv")
    report_df, preds = load_lens(path, predictor="all")
    # Both kinds should be present
    kinds_seen = {
        "pMHC_affinity" if p.prediction_method_name in ("mhcflurry", "netmhcpan")
        else "pMHC_stability" for p in preds
    }
    assert kinds_seen == {"pMHC_affinity", "pMHC_stability"}

    cfg = EpitopeConfig(
        score_expr=(
            "affinity['mhcflurry'].logistic(350, 150) * 0.5 + "
            # stability: bound hours / 24h window, clipped to [0, 1]
            "(stability['netmhcstabpan'].value / 24).clip(0, 1) * 0.5"
        )
    )
    csv_path = tmp_path / "report.csv"
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    assert len(result) == 2
    assert (result["vaxrank_score"] > 0).all()


def test_lens_percentile_rank_scoring_matches_legacy(tmp_path):
    """With scoring_mode='percentile_rank' the DSL-routed path should match
    the legacy EpitopePrediction.logistic_epitope_score output element-wise."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, preds = load_lens(path)

    cfg = EpitopeConfig(scoring_mode="percentile_rank")
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
    result_scores_by_key = dict(zip(
        zip(result["Mutant peptide sequence"], result["Allele"]),
        result["vaxrank_score"],
    ))
    new_scores = [
        result_scores_by_key[(p.peptide_sequence, p.allele)] for p in preds
    ]
    for legacy, new in zip(legacy_scores, new_scores):
        assert new == pytest.approx(legacy)


def test_default_scoring_matches_legacy(tmp_path):
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
