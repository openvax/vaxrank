# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests for the topiary-DSL construction in vaxrank.epitope_dsl.

The parity tests (``test_default_score_matches_legacy_*``) verify that the
default synthesized score node produces byte-identical output to the legacy
``EpitopePrediction.logistic_epitope_score``. This is the back-compat
contract that makes the DSL migration safe for existing users and fixtures.
"""

import math

import numpy as np
import pandas as pd
import pytest
from topiary.ranking import EvalContext, apply_filter

from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dsl import (
    _logistic_normalizer,
    build_filter_node,
    build_score_node,
)
from vaxrank.epitope_prediction import EpitopePrediction

from .common import eq_


def _predictions_df(rows):
    """Build a minimal topiary-shaped predictions DataFrame from row dicts."""
    base = {
        "source_sequence_name": "seq",
        "peptide_offset": 0,
        "kind": "pMHC_affinity",
        "prediction_method_name": "test",
        "predictor_version": "1.0",
    }
    records = []
    for i, r in enumerate(rows):
        rec = dict(base)
        rec.update(r)
        rec.setdefault("peptide_offset", i)
        rec["peptide_length"] = len(rec["peptide"])
        rec["value"] = rec.get("value", rec.get("affinity", np.nan))
        records.append(rec)
    return pd.DataFrame(records)


def _legacy_score(ic50, percentile_rank, cfg):
    pred = EpitopePrediction(
        allele="HLA-A*02:01",
        peptide_sequence="AAAAAAAAA",
        wt_peptide_sequence="AAAAAAAAA",
        ic50=ic50,
        wt_ic50=ic50,
        percentile_rank=percentile_rank,
        prediction_method_name="test",
        overlaps_mutation=True,
        source_sequence="AAAAAAAAA",
        offset=0,
        occurs_in_reference=False,
    )
    return pred.logistic_epitope_score(
        midpoint=cfg.logistic_epitope_score_midpoint,
        width=cfg.logistic_epitope_score_width,
        ic50_cutoff=cfg.binding_affinity_cutoff,
        scoring_mode=cfg.scoring_mode,
        percentile_rank_cutoff=cfg.percentile_rank_cutoff,
    )


def _eval_score(cfg, ic50_values, percentile_rank_values=None):
    rows = []
    for i, ic50 in enumerate(ic50_values):
        row = {
            "peptide": f"PEP{i:05d}",
            "allele": "HLA-A*02:01",
            "affinity": ic50,
            "value": ic50,
        }
        if percentile_rank_values is not None:
            row["percentile_rank"] = percentile_rank_values[i]
        rows.append(row)
    df = _predictions_df(rows)
    node = build_score_node(cfg)
    return node.eval(EvalContext(df)).reindex(EvalContext(df).group_index)


def test_logistic_normalizer_matches_legacy_formula():
    # Pre-5.0 vaxrank used 1/(1+exp(-midpoint/width)) as the normalizer so
    # the score equals 1.0 at ic50=0. Verify the helper reproduces that.
    eq_(_logistic_normalizer(350, 150), 1.0 / (1.0 + math.exp(-350 / 150)))


def test_default_score_matches_legacy_affinity_mode():
    cfg = EpitopeConfig()  # defaults: midpoint=350, width=150, cutoff=5000
    ic50s = [0.0, 50.0, 150.0, 350.0, 500.0, 1000.0, 4999.0, 5000.0, 10000.0]
    scores = _eval_score(cfg, ic50s)
    for ic50, dsl_score in zip(ic50s, scores):
        legacy = _legacy_score(ic50, None, cfg)
        assert float(dsl_score) == pytest.approx(legacy, abs=1e-12), (
            f"DSL score {float(dsl_score)} != legacy {legacy} at ic50={ic50}"
        )


def test_default_score_matches_legacy_percentile_rank_mode():
    cfg = EpitopeConfig(scoring_mode="percentile_rank", percentile_rank_cutoff=2.0)
    # Vaxrank stores percentile_rank as 0..100 with the "< cutoff" check.
    ranks = [0.0, 0.5, 1.0, 1.999, 2.0, 2.5, 10.0]
    ic50s = [100.0] * len(ranks)
    scores = _eval_score(cfg, ic50s, percentile_rank_values=ranks)
    for rank, dsl_score in zip(ranks, scores):
        legacy = _legacy_score(100.0, rank, cfg)
        assert float(dsl_score) == pytest.approx(legacy, abs=1e-12), (
            f"DSL score {float(dsl_score)} != legacy {legacy} at rank={rank}"
        )


def test_default_score_at_cutoff_is_zero_affinity_mode():
    cfg = EpitopeConfig(binding_affinity_cutoff=500.0)
    scores = _eval_score(cfg, [499.0, 500.0, 501.0])
    # Legacy semantics: ic50 >= cutoff → 0. So 500.0 itself is zeroed.
    assert float(scores.iloc[0]) > 0
    eq_(float(scores.iloc[1]), 0.0)
    eq_(float(scores.iloc[2]), 0.0)


def test_no_default_filter_node():
    """Default path has no filter node; the score node's cutoff mask +
    min_epitope_score gate reproduces legacy behavior. This preserves
    tests that set ``min_epitope_score=0`` expecting "keep every peptide"."""
    cfg = EpitopeConfig(binding_affinity_cutoff=500.0)
    assert build_filter_node(cfg) is None


def test_custom_filter_expr_is_parsed():
    cfg = EpitopeConfig(filter_expr="affinity <= 250")
    df = _predictions_df(
        [
            {"peptide": "AAAAAAAAA", "allele": "HLA-A*02:01", "affinity": 100.0},
            {"peptide": "BBBBBBBBB", "allele": "HLA-A*02:01", "affinity": 300.0},
        ]
    )
    filtered = apply_filter(df, build_filter_node(cfg))
    eq_(set(filtered["peptide"]), {"AAAAAAAAA"})


def test_custom_score_expr_is_parsed():
    cfg = EpitopeConfig(score_expr="affinity.logistic(350, 150)")
    # User-provided expression: no normalizer → raw sigmoid values.
    scores = _eval_score(cfg, [0.0])
    expected = 1.0 / (1.0 + math.exp(-350 / 150))
    assert float(scores.iloc[0]) == pytest.approx(expected, abs=1e-12)


def test_invalid_filter_expr_raises_at_build_time():
    cfg = EpitopeConfig(filter_expr="not a valid dsl expression")
    with pytest.raises(Exception):
        build_filter_node(cfg)


def test_unknown_column_error_is_actionable():
    cfg = EpitopeConfig(filter_expr="column('no_such_column') <= 1")
    df = _predictions_df(
        [{"peptide": "AAAAAAAAA", "allele": "HLA-A*02:01", "affinity": 100.0}]
    )
    with pytest.raises(ValueError, match="no_such_column"):
        apply_filter(df, build_filter_node(cfg))
