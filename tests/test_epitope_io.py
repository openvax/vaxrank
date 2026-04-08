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
    load_epitopes,
    detect_format,
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
    preds = load_pvacseq(path)
    assert len(preds) == 3

    # Check first prediction
    p = preds[0]
    assert p.allele == "HLA-A*02:01"
    assert p.peptide_sequence == "SVVGSSSSS"
    assert p.ic50 == pytest.approx(76.11)
    assert p.wt_ic50 == pytest.approx(4520.30)
    assert p.percentile_rank == pytest.approx(0.5)
    assert p.overlaps_mutation is True

    # Check third prediction has occurs_in_reference=True (Ref Match=True)
    assert preds[2].occurs_in_reference is True


def test_load_pvacseq_missing_columns(tmp_path):
    bad_file = tmp_path / "bad.tsv"
    bad_file.write_text("Gene\tAA Change\nTP53\tG245S\n")
    with pytest.raises(ValueError, match="missing required columns"):
        load_pvacseq(bad_file)


# ── LENS import ──────────────────────────────────────────────────────────────

def test_load_lens():
    path = os.path.join(DATA_DIR, "lens_example.tsv")
    preds = load_lens(path)
    assert len(preds) == 3

    p = preds[0]
    assert p.allele == "HLA-A*02:01"
    assert p.peptide_sequence == "SVVGSSSSS"
    assert p.ic50 == pytest.approx(76.11)
    assert p.wt_ic50 == pytest.approx(4520.30)
    assert p.wt_peptide_sequence == "SVVGGSSSS"
    # Should prefer EL percentile rank
    assert p.percentile_rank == pytest.approx(0.3)


def test_load_lens_missing_columns(tmp_path):
    bad_file = tmp_path / "bad.tsv"
    bad_file.write_text("gene_name\tpos\nTP53\t5\n")
    with pytest.raises(ValueError, match="missing required columns"):
        load_lens(bad_file)


# ── Format auto-detection ────────────────────────────────────────────────────

def test_detect_format_pvacseq():
    path = os.path.join(DATA_DIR, "pvacseq_example.tsv")
    assert detect_format(path) == "pvacseq"


def test_detect_format_lens():
    path = os.path.join(DATA_DIR, "lens_example.tsv")
    assert detect_format(path) == "lens"


def test_detect_format_vaxrank(tmp_path):
    path = tmp_path / "vaxrank.csv"
    save_predictions([_make_prediction()], path)
    assert detect_format(path) == "vaxrank"


def test_detect_format_unknown(tmp_path):
    bad_file = tmp_path / "unknown.csv"
    bad_file.write_text("foo,bar\n1,2\n")
    with pytest.raises(ValueError, match="Cannot detect"):
        detect_format(bad_file)


def test_load_epitopes_auto_pvacseq():
    path = os.path.join(DATA_DIR, "pvacseq_example.tsv")
    preds = load_epitopes(path)
    assert len(preds) == 3


def test_load_epitopes_auto_lens():
    path = os.path.join(DATA_DIR, "lens_example.tsv")
    preds = load_epitopes(path)
    assert len(preds) == 3


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
