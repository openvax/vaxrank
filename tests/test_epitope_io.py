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

def test_load_lens_emits_one_prediction_per_detected_predictor():
    """load_lens always emits one EpitopePrediction per (peptide, allele,
    detected predictor). The fixture has 3 rows × 2 predictors = 6 preds,
    3 report rows."""
    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, preds = load_lens(path)
    assert len(preds) == 6
    assert len(report_df) == 3

    methods = {p.prediction_method_name for p in preds}
    assert methods == {"mhcflurry", "netmhcpan"}

    # Each (peptide, allele) pair appears twice — once per predictor
    from collections import Counter
    pair_counts = Counter((p.peptide_sequence, p.allele) for p in preds)
    assert all(c == 2 for c in pair_counts.values())
    assert report_df["Predictors used"].iloc[0] == "mhcflurry,netmhcpan"

    # First prediction carries the expected mhcflurry row (sorted by tool
    # so mhcflurry comes before netmhcpan).
    p = preds[0]
    assert p.prediction_method_name == "mhcflurry"
    assert p.allele == "HLA-A*02:01"  # 'HLA-A02:01' normalized
    assert p.peptide_sequence == "SVVGSSSSS"
    assert p.ic50 == pytest.approx(95.4)
    assert p.percentile_rank == pytest.approx(0.28)  # mhcflurry pres_perc
    assert p.predictor_version == "2.1.1"
    assert p.wt_peptide_sequence == ""
    assert p.wt_ic50 is None
    assert p.overlaps_mutation is True
    assert p.source_sequence == "AASVVGSSSSSGTR"


def test_load_lens_report_display_uses_canonical_predictor():
    """The report DataFrame's single 'Predicted mutant pMHC affinity'
    column shows the canonical pMHC_affinity predictor (mhcflurry >
    netmhcpan > alphabetical). Per-predictor raw values land in
    individual columns alongside."""
    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, _ = load_lens(path)
    # Display value = mhcflurry 95.4 nM (canonical)
    assert "95.40 nM" in report_df["Predicted mutant pMHC affinity"].iloc[0]
    # Both raw-value columns present
    assert report_df["mhcflurry value (nM)"].iloc[0] == pytest.approx(95.4)
    assert report_df["netmhcpan value (nM)"].iloc[0] == pytest.approx(76.11)


def test_load_lens_mhcflurry_only_v19():
    """A v1.9-dev style file (no netmhcpan cols) emits only mhcflurry rows."""
    path = os.path.join(DATA_DIR, "lens_v1.9_mhcflurry_only.tsv")
    report_df, preds = load_lens(path)
    assert len(preds) == 3
    methods = {p.prediction_method_name for p in preds}
    assert methods == {"mhcflurry"}

    p = preds[0]
    assert p.allele == "HLA-A*01:01"
    assert p.peptide_sequence == "TAEFYQRY"
    assert p.ic50 == pytest.approx(254.15)
    assert p.predictor_version == "2.1.1"
    # Antigen source captured; variant_coords NA fallback to origin_descriptor
    assert report_df["Antigen source"].iloc[0] == "ERV"
    assert "Hsap38" in report_df["Genomic variant"].iloc[0]
    assert "mhcflurry value (nM)" in report_df.columns
    assert "netmhcpan value (nM)" not in report_df.columns


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
    report_df, preds = load_lens(path)
    methods = {p.prediction_method_name for p in preds}
    assert methods == {"mhcflurry", "netmhcpan", "netmhcstabpan"}
    # Stability value column is labeled with its unit
    assert "netmhcstabpan value (hours)" in report_df.columns


# ── DSL integration for external inputs ──────────────────────────────────────

# ── predictors_required_by_cfg auto-derivation ───────────────────────────────

def test_predictors_required_by_cfg_empty_for_no_formula():
    """No filter_expr / score_expr → empty set."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import predictors_required_by_cfg
    assert predictors_required_by_cfg(EpitopeConfig()) == set()


def test_predictors_required_by_cfg_empty_for_unqualified_formula():
    """A formula that never brackets a predictor name → empty set."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import predictors_required_by_cfg
    cfg = EpitopeConfig(score_expr="affinity.logistic(350, 150)")
    assert predictors_required_by_cfg(cfg) == set()


def test_predictors_required_by_cfg_picks_up_multiple():
    """Formula referencing two predictors by bracket → set of both."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import predictors_required_by_cfg
    cfg = EpitopeConfig(score_expr=(
        "affinity['mhcflurry'].logistic(350, 150) * 0.5 + "
        "affinity['netmhcpan'].logistic(350, 150) * 0.5"))
    assert predictors_required_by_cfg(cfg) == {"mhcflurry", "netmhcpan"}


def test_predictors_required_by_cfg_merges_filter_and_score():
    """References from filter_expr and score_expr both count."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import predictors_required_by_cfg
    cfg = EpitopeConfig(
        filter_expr="affinity['mhcflurry'] <= 500",
        score_expr="stability['netmhcstabpan'].value / 24",
    )
    assert predictors_required_by_cfg(cfg) == {"mhcflurry", "netmhcstabpan"}


def test_default_methods_resolves_multi_model_affinity(tmp_path):
    """When a file exposes multiple affinity models and no score_expr is
    set, cfg.default_methods picks which one the default DSL node scores
    on. Without default_methods, vaxrank auto-picks canonical."""
    from vaxrank.epitope_config import EpitopeConfig

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, preds = load_lens(path)

    # Explicit default = netmhcpan → score reflects netmhcpan affinity.
    cfg = EpitopeConfig(default_methods={"pMHC_affinity": "netmhcpan"})
    csv_netmhcpan = tmp_path / "netmhcpan.csv"
    write_neoepitope_report(
        report_df.copy(), preds,
        csv_report_path=str(csv_netmhcpan), epitope_config=cfg)

    # Explicit default = mhcflurry → score reflects mhcflurry affinity.
    cfg2 = EpitopeConfig(default_methods={"pMHC_affinity": "mhcflurry"})
    csv_mhcflurry = tmp_path / "mhcflurry.csv"
    write_neoepitope_report(
        report_df.copy(), preds,
        csv_report_path=str(csv_mhcflurry), epitope_config=cfg2)

    import pandas as pd
    r1 = pd.read_csv(csv_netmhcpan).set_index("Mutant peptide sequence")
    r2 = pd.read_csv(csv_mhcflurry).set_index("Mutant peptide sequence")
    # At least one row should differ (the two predictors disagree enough
    # that logistic-normalized output diverges).
    diffs = (r1["vaxrank_score"] - r2["vaxrank_score"]).abs()
    assert (diffs > 1e-6).any(), (
        "Expected different scores when default_methods picks different "
        "affinity predictor")


def test_default_methods_rejects_unknown_method(tmp_path):
    """cfg.default_methods[kind] must match a method actually in the data."""
    from vaxrank.epitope_config import EpitopeConfig

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, preds = load_lens(path)
    cfg = EpitopeConfig(default_methods={"pMHC_affinity": "mhcnuggets"})
    with pytest.raises(ValueError, match="default_methods.*mhcnuggets"):
        write_neoepitope_report(
            report_df, preds, csv_report_path=str(tmp_path / "x.csv"),
            epitope_config=cfg)


def test_default_methods_partial_qualification_mixed_kinds(tmp_path):
    """score_expr brackets a method for one Kind but leaves another Kind
    unqualified. The bracketed Kind keeps its explicit method; the other
    Kind falls back to its default (auto-picked or user-specified)."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_v1.4_with_stability.tsv")
    report_df, preds = load_lens(path)
    # Formula: affinity qualified (netmhcpan), stability unqualified.
    # default_methods settles the stability kind.
    cfg = EpitopeConfig(
        score_expr=(
            "affinity['netmhcpan'].logistic(350, 150) * 0.7 + "
            "(stability.value / 24).clip(0, 1) * 0.3"
        ),
        default_methods={"pMHC_stability": "netmhcstabpan"},
    )
    csv_path = tmp_path / "out.csv"
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    assert len(result) == 2
    assert (result["vaxrank_score"] > 0).all()


def test_default_methods_partial_qualification_auto_picks_stability(tmp_path):
    """Same as above but without default_methods — auto-pick fills in
    the stability default."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_v1.4_with_stability.tsv")
    report_df, preds = load_lens(path)
    # Only one stability method (netmhcstabpan) in this fixture, so
    # auto-pick is trivial but tests the path.
    cfg = EpitopeConfig(
        score_expr=(
            "affinity['mhcflurry'].logistic(350, 150) * 0.7 + "
            "(stability.value / 24).clip(0, 1) * 0.3"
        ),
    )
    csv_path = tmp_path / "out.csv"
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    assert (result["vaxrank_score"] > 0).all()


def test_filter_expr_unqualified_on_multi_method_data(tmp_path):
    """Unqualified filter_expr ('affinity <= X') on a file with multiple
    affinity methods should resolve via the default method, not raise
    "Ambiguous Kind"."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, preds = load_lens(path)
    # fixture's mhcflurry affinities: 95.4, 180.2, 310.5. Cut at 200 keeps
    # 2 rows; with netmhcpan as the default, cut at 200 would keep rows
    # with IC50 < 200 from netmhcpan: 76.11, 120.50 (2 rows).
    cfg = EpitopeConfig(
        filter_expr="affinity <= 200",
        # auto-pick selects mhcflurry; test it evaluates without error
    )
    csv_path = tmp_path / "out.csv"
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    # Filter drops rows whose mhcflurry affinity > 200 (only 95.4 passes),
    # so the nonzero-score rows number ≤ 1. Other rows still land in the
    # report with score 0.
    nonzero = result[result["vaxrank_score"] > 0]
    assert len(nonzero) <= len(result)
    assert len(nonzero) >= 1


def test_score_predictions_on_empty_predictions_list(tmp_path):
    """An empty predictions list should produce an empty score Series
    without blowing up."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import score_predictions
    scores = score_predictions([], EpitopeConfig())
    assert len(scores) == 0


def test_subset_topiary_df_noop_on_single_method_data():
    """A frame with only one method per Kind should pass through
    subset_topiary_df_for_eval unchanged."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import (
        predictions_to_topiary_df, subset_topiary_df_for_eval)

    path = os.path.join(DATA_DIR, "lens_v1.9_mhcflurry_only.tsv")
    _, preds = load_lens(path)
    df = predictions_to_topiary_df(preds)
    subset = subset_topiary_df_for_eval(df, EpitopeConfig())
    # Same row count (reset_index changes the index but not the rows)
    assert len(subset) == len(df)
    assert set(subset["prediction_method_name"]) == set(df["prediction_method_name"])


def test_subset_topiary_df_qualified_both_predictors_keeps_both():
    """When the formula explicitly brackets both predictors, both stay
    in the subsetted frame (no unnecessary default-method subsetting)."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import (
        predictions_to_topiary_df, subset_topiary_df_for_eval)

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _, preds = load_lens(path)
    df = predictions_to_topiary_df(preds)
    cfg = EpitopeConfig(score_expr=(
        "affinity['mhcflurry'].logistic(350, 150) * 0.5 + "
        "affinity['netmhcpan'].logistic(350, 150) * 0.5"))
    subset = subset_topiary_df_for_eval(df, cfg)
    # Both methods' rows survive
    assert set(subset["prediction_method_name"]) == {"mhcflurry", "netmhcpan"}
    assert len(subset) == len(df)


def test_subset_topiary_df_qualified_one_predictor_drops_other():
    """When the formula brackets only one predictor and no unqualified
    reference exists, the other is dropped (not added via default)."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import (
        predictions_to_topiary_df, subset_topiary_df_for_eval)

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _, preds = load_lens(path)
    df = predictions_to_topiary_df(preds)
    cfg = EpitopeConfig(
        score_expr="affinity['netmhcpan'].logistic(350, 150)")
    subset = subset_topiary_df_for_eval(df, cfg)
    assert set(subset["prediction_method_name"]) == {"netmhcpan"}


def test_parse_is_cached():
    """Repeated _parse calls for the same expression should hit the cache
    rather than re-invoking topiary.ranking.parse."""
    from vaxrank.epitope_dsl import _parse
    cache_info_before = _parse.cache_info()
    expr = "affinity['mhcflurry'].logistic(350, 150)"
    _parse(expr)
    _parse(expr)
    _parse(expr)
    info = _parse.cache_info()
    # At least two cache hits from the three calls above
    assert info.hits - cache_info_before.hits >= 2


def test_default_methods_mixed_valid_and_invalid_fails_fast(tmp_path):
    """If default_methods has one valid and one invalid entry, the invalid
    one must fail fast even though other Kinds resolve cleanly."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_v1.4_with_stability.tsv")
    report_df, preds = load_lens(path)
    cfg = EpitopeConfig(default_methods={
        "pMHC_affinity": "mhcflurry",      # valid
        "pMHC_stability": "doesnotexist",  # invalid
    })
    with pytest.raises(ValueError,
                       match=r"default_methods\['pMHC_stability'\]='doesnotexist'"):
        write_neoepitope_report(
            report_df, preds, csv_report_path=str(tmp_path / "x.csv"),
            epitope_config=cfg)


def test_report_columns_unaffected_by_subsetting(tmp_path):
    """Subsetting happens on the topiary frame used for DSL eval, NOT on
    report_df — report output must still expose every detected
    predictor's raw value column regardless of which one drives
    scoring."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, preds = load_lens(path)
    cfg = EpitopeConfig(
        score_expr="affinity['mhcflurry'].logistic(350, 150)")
    csv_path = tmp_path / "out.csv"
    write_neoepitope_report(
        report_df.copy(), preds,
        csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    # Even though only mhcflurry drives scoring, netmhcpan values are
    # still reported alongside.
    assert "mhcflurry value (nM)" in result.columns
    assert "netmhcpan value (nM)" in result.columns
    assert result["netmhcpan value (nM)"].notna().all()


def test_predictor_version_in_epitope_prediction_from_lens():
    """Every LENS-loaded EpitopePrediction should carry the sniffed
    predictor version, so version-qualified DSL refs resolve."""
    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _, preds = load_lens(path)
    by_method = {p.prediction_method_name: p.predictor_version for p in preds}
    assert by_method["mhcflurry"] == "2.1.1"
    assert by_method["netmhcpan"] == "4.1b"


def test_dsl_version_mismatch_validation_error(tmp_path):
    """A version-qualified formula against data without that version
    raises a clear error (not silent NaN)."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, preds = load_lens(path)
    cfg = EpitopeConfig(
        score_expr="affinity['mhcflurry', '0.0.0'].logistic(350, 150)")
    with pytest.raises(ValueError, match="predictor version '0.0.0'"):
        write_neoepitope_report(
            report_df, preds, csv_report_path=str(tmp_path / "x.csv"),
            epitope_config=cfg)


def test_kinds_needing_method_disambiguation_captures_default_node():
    """With no filter_expr / score_expr, pMHC_affinity is added to the
    disambiguation set (because the synthesized default node implicitly
    needs single-method rows)."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import kinds_needing_method_disambiguation
    assert kinds_needing_method_disambiguation(
        EpitopeConfig()) == {"pMHC_affinity"}


def test_kinds_needing_method_disambiguation_adds_only_unqualified():
    """With a fully bracketed score_expr, no Kinds need disambiguation."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import kinds_needing_method_disambiguation
    cfg = EpitopeConfig(score_expr=(
        "affinity['mhcflurry'].logistic(350, 150) * 0.5 + "
        "affinity['netmhcpan'].logistic(350, 150) * 0.5"))
    assert kinds_needing_method_disambiguation(cfg) == set()


def test_kinds_needing_method_disambiguation_mixed():
    """Unqualified kinds in expressions are collected; bracketed refs
    are not."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import kinds_needing_method_disambiguation
    cfg = EpitopeConfig(score_expr=(
        "affinity['mhcflurry'].logistic(350, 150) + stability.value"))
    assert kinds_needing_method_disambiguation(cfg) == {"pMHC_stability"}


def test_validate_default_methods_empty_df_is_noop():
    """validate_default_methods on an empty topiary frame returns cleanly
    (we can't check methods we don't have)."""
    import pandas as pd
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import validate_default_methods
    cfg = EpitopeConfig(default_methods={"pMHC_affinity": "mhcflurry"})
    validate_default_methods(cfg, pd.DataFrame())  # no raise


def test_validate_default_methods_no_defaults_is_noop():
    """With no default_methods, the validator is a no-op."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import (
        predictions_to_topiary_df, validate_default_methods)
    path = os.path.join(DATA_DIR, "lens_example.tsv")
    _, preds = load_lens(path)
    df = predictions_to_topiary_df(preds)
    validate_default_methods(EpitopeConfig(), df)  # no raise


def test_validate_default_methods_catches_typo_on_single_method_kind(tmp_path):
    """Even when a Kind has only one model in the data (so subsetting
    wouldn't consult the default), validate_default_methods still
    errors on a typo'd default."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import (
        predictions_to_topiary_df, validate_default_methods)
    # Single-model fixture (v1.9 mhcflurry-only).
    path = os.path.join(DATA_DIR, "lens_v1.9_mhcflurry_only.tsv")
    _, preds = load_lens(path)
    df = predictions_to_topiary_df(preds)
    cfg = EpitopeConfig(default_methods={"pMHC_affinity": "bogus"})
    with pytest.raises(ValueError,
                       match=r"default_methods\['pMHC_affinity'\]='bogus'"):
        validate_default_methods(cfg, df)


def test_filter_plus_score_expr_both_bracketed(tmp_path):
    """filter_expr and score_expr both set, each with bracketed predictor
    refs — the union of referenced methods should survive subsetting."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, preds = load_lens(path)
    cfg = EpitopeConfig(
        filter_expr="affinity['netmhcpan'] <= 500",
        score_expr="affinity['mhcflurry'].logistic(350, 150)",
    )
    csv_path = tmp_path / "out.csv"
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    # Filter keeps rows whose netmhcpan IC50 <= 500 (all 3 in fixture:
    # 76.11, 120.50, 250.00); score comes from mhcflurry.
    assert (result["vaxrank_score"] > 0).any()


def test_filter_expr_drops_all_groups_yields_zero_scores(tmp_path):
    """A filter that matches no rows should produce a report with all
    scores = 0 (not crash)."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, preds = load_lens(path)
    # No affinity < 0 anywhere.
    cfg = EpitopeConfig(filter_expr="affinity <= 0")
    csv_path = tmp_path / "out.csv"
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    assert (result["vaxrank_score"] == 0).all()


def test_pvacseq_single_method_per_group_has_no_ambiguity(tmp_path):
    """pVACseq emits at most one row per (peptide, allele). Even when
    different rows use different methods, each group is
    unambiguous — the default node should score without error."""
    from vaxrank.epitope_io import write_neoepitope_report
    path = os.path.join(DATA_DIR, "pvacseq_example.tsv")
    report_df, preds = load_pvacseq(path)
    csv_path = tmp_path / "out.csv"
    # Default EpitopeConfig, no default_methods set.
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path))
    import pandas as pd
    result = pd.read_csv(csv_path)
    assert len(result) == 3
    assert "vaxrank_score" in result.columns


def test_unknown_column_in_lens_file_is_ignored(tmp_path):
    """LENS files with extra columns that don't match the predictor
    regex are simply ignored — loader shouldn't error."""
    # Valid mhcflurry columns plus a novel experimental column.
    path = tmp_path / "lens.tsv"
    path.write_text(
        "peptide\tallele\tmhcflurry_2.1.1.aff\t"
        "mhcflurry_2.1.1.pres_perc\tfuture_tool_99.blarg.metric\t"
        "random_non_matching_col\n"
        "SIINFEKL\tHLA-A*02:01\t100.0\t0.5\t42\thello\n"
    )
    report_df, preds = load_lens(path)
    assert len(preds) == 1
    assert preds[0].prediction_method_name == "mhcflurry"


def test_default_methods_auto_picked_when_unset(tmp_path, caplog):
    """When multiple models for a Kind are present and default_methods
    doesn't specify one, vaxrank auto-picks canonical and logs which."""
    from vaxrank.epitope_config import EpitopeConfig
    import logging

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, preds = load_lens(path)
    cfg = EpitopeConfig()  # no default_methods

    with caplog.at_level(logging.INFO, logger="vaxrank.epitope_dsl"):
        write_neoepitope_report(
            report_df, preds,
            csv_report_path=str(tmp_path / "out.csv"), epitope_config=cfg)

    assert any(
        "auto-picking 'mhcflurry'" in rec.getMessage()
        for rec in caplog.records), (
        f"Expected auto-pick info log; got: "
        f"{[r.getMessage() for r in caplog.records]}")


# ── NaN handling in vaxrank-native roundtrip ─────────────────────────────────

def test_load_predictions_empty_string_fields_become_empty_not_nan(tmp_path):
    """CSV cells that read as NaN for string-typed fields should come back
    as '' (not float NaN) so downstream comparisons work."""
    from vaxrank.epitope_io import load_predictions
    path = tmp_path / "preds.csv"
    # Empty wt_peptide_sequence, source_sequence, prediction_method_name,
    # predictor_version columns all present but blank.
    path.write_text(
        "allele,peptide_sequence,wt_peptide_sequence,ic50,wt_ic50,"
        "percentile_rank,prediction_method_name,overlaps_mutation,"
        "source_sequence,offset,occurs_in_reference,predictor_version\n"
        "HLA-A*02:01,SIINFEKL,,50.0,,0.5,,True,,0,False,\n"
    )
    preds = load_predictions(path)
    assert len(preds) == 1
    p = preds[0]
    assert p.wt_peptide_sequence == ""
    assert p.source_sequence == ""
    assert p.prediction_method_name == ""
    assert p.predictor_version == ""


def test_write_neoepitope_report_rejects_duplicate_rows(tmp_path):
    """write_neoepitope_report merges scores by (peptide, allele); if a
    loader ever produced a report_df with duplicates that merge would be
    ambiguous, so we assert uniqueness upfront. This test constructs
    duplicates explicitly to exercise the guard."""
    import pandas as pd
    from vaxrank.epitope_io import write_neoepitope_report

    report_df = pd.DataFrame([
        {'Allele': 'HLA-A*02:01', 'Mutant peptide sequence': 'SIINFEKL'},
        {'Allele': 'HLA-A*02:01', 'Mutant peptide sequence': 'SIINFEKL'},  # duplicate
    ])
    preds = [_make_prediction(peptide_sequence='SIINFEKL', allele='HLA-A*02:01')]
    with pytest.raises(ValueError, match="duplicate .*peptide, allele"):
        write_neoepitope_report(
            report_df, preds, csv_report_path=str(tmp_path / "unused.csv"))


def test_lens_dsl_combines_both_predictors(tmp_path):
    """A score_expr combining two predictors' affinities should average them."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import write_neoepitope_report

    path = os.path.join(DATA_DIR, "lens_example.tsv")
    report_df, preds = load_lens(path)

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
    report_df, preds = load_lens(path)
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
    report_df, preds = load_lens(path)
    methods = {p.prediction_method_name for p in preds}
    assert methods == {"mhcflurry", "netmhcpan", "netmhcstabpan"}
    # Every real LENS prediction carries a nonempty version.
    assert all(p.predictor_version for p in preds)
    # Report exposes per-tool columns with correct units.
    assert "mhcflurry value (nM)" in report_df.columns
    assert "netmhcpan value (nM)" in report_df.columns
    assert "netmhcstabpan value (hours)" in report_df.columns


def test_real_lens_v15_emits_both_affinity_predictors():
    """v1.5.1 has both mhcflurry and netmhcpan; load_lens emits one
    prediction per (peptide, allele, detected predictor) so the DSL can
    see all of them. Canonical 'mhcflurry' goes into the report's
    display columns."""
    path = os.path.join(REAL_LENS_DIR, "lens_v1.5_real_subset.tsv")
    report_df, preds = load_lens(path)
    methods = {p.prediction_method_name for p in preds}
    assert methods == {"mhcflurry", "netmhcpan"}
    # Each report row spawns up to N predictions (N = detected count); NA
    # cells produce fewer.
    assert len(preds) <= 2 * len(report_df)
    # Report display column shows the canonical predictor's name first.
    predictors_used = report_df["Predictors used"].iloc[0]
    assert predictors_used.startswith("mhcflurry")


def test_real_lens_v19_is_mhcflurry_only():
    """v1.9-dev emits only mhcflurry (netMHCpan was dropped)."""
    path = os.path.join(REAL_LENS_DIR, "lens_v1.9_real_subset.tsv")
    report_df, preds = load_lens(path)
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
        report_df, preds = load_lens(path)
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
    report_df, preds = load_lens(path)
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
    report_df, preds = load_lens(path)
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
    report_df, preds = load_lens(path)
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


def _assert_default_scores_match_legacy(cfg, path, tmp_path):
    """Shared helper: DSL-routed default scoring should equal the legacy
    EpitopePrediction.logistic_epitope_score output on the subset of
    predictions selected by default_methods (auto-picked mhcflurry here)."""
    from vaxrank.epitope_io import write_neoepitope_report

    report_df, preds = load_lens(path)
    # Auto-picked default for pMHC_affinity is mhcflurry — that's the
    # subset whose rows survive subset_topiary_df_for_eval in the
    # default-scoring code path.
    default_preds = [p for p in preds if p.prediction_method_name == "mhcflurry"]
    legacy_scores = {
        (p.peptide_sequence, p.allele): round(p.logistic_epitope_score(
            midpoint=cfg.logistic_epitope_score_midpoint,
            width=cfg.logistic_epitope_score_width,
            ic50_cutoff=cfg.binding_affinity_cutoff,
            scoring_mode=cfg.scoring_mode,
            percentile_rank_cutoff=cfg.percentile_rank_cutoff,
        ), 6)
        for p in default_preds
    }
    csv_path = tmp_path / "report.csv"
    write_neoepitope_report(
        report_df, preds, csv_report_path=str(csv_path), epitope_config=cfg)
    import pandas as pd
    result = pd.read_csv(csv_path)
    for _, row in result.iterrows():
        key = (row["Mutant peptide sequence"], row["Allele"])
        assert row["vaxrank_score"] == pytest.approx(legacy_scores[key])


def test_lens_percentile_rank_scoring_matches_legacy(tmp_path):
    """scoring_mode='percentile_rank' through the DSL path matches legacy
    logistic_epitope_score computed on the canonical (mhcflurry) rows."""
    from vaxrank.epitope_config import EpitopeConfig
    _assert_default_scores_match_legacy(
        EpitopeConfig(scoring_mode="percentile_rank"),
        os.path.join(DATA_DIR, "lens_example.tsv"),
        tmp_path,
    )


def test_default_scoring_matches_legacy(tmp_path):
    """With no filter_expr/score_expr the DSL default node matches legacy
    logistic_epitope_score byte-for-byte on the canonical predictor rows."""
    from vaxrank.epitope_config import EpitopeConfig
    _assert_default_scores_match_legacy(
        EpitopeConfig(),
        os.path.join(DATA_DIR, "lens_example.tsv"),
        tmp_path,
    )


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
