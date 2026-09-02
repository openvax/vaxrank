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
``pre-3.0 logistic_epitope_score``. This is the back-compat
contract that makes the DSL migration safe for existing users and fixtures.
"""

import math

import numpy as np
import pandas as pd
import pytest
import topiary
from topiary.ranking import EvalContext, apply_filter

from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dsl import build_filter_node, build_score_node

from ._legacy_score_reference import legacy_score_one
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
    return legacy_score_one(
        ic50, percentile_rank,
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


def test_invalid_filter_expr_raises_at_config_construction():
    # Eager validation: malformed expressions should fail at EpitopeConfig
    # construction, not later when build_filter_node() is called. This
    # means a bad YAML config fails at load, not mid-pipeline.
    with pytest.raises(ValueError, match="Invalid filter_expr"):
        EpitopeConfig(filter_expr="not a valid dsl expression")


def test_invalid_score_expr_raises_at_config_construction():
    with pytest.raises(ValueError, match="Invalid score_expr"):
        EpitopeConfig(score_expr="this is not valid")


def test_unknown_column_error_is_actionable():
    # column(IDENT) parses fine; the "column missing from df" error fires
    # when apply_filter actually runs.
    cfg = EpitopeConfig(filter_expr="column(no_such_column) <= 1")
    df = _predictions_df(
        [{"peptide": "AAAAAAAAA", "allele": "HLA-A*02:01", "affinity": 100.0}]
    )
    with pytest.raises(ValueError, match="no_such_column"):
        apply_filter(df, build_filter_node(cfg))


def test_default_score_parity_across_alleles():
    """Multi-allele score parity: one score per (peptide, allele) group tuple,
    each matching legacy per-row scorer."""
    cfg = EpitopeConfig()
    rows = [
        {"peptide": "AAAAAAAAA", "allele": "HLA-A*02:01", "affinity": 100.0},
        {"peptide": "AAAAAAAAA", "allele": "HLA-B*07:02", "affinity": 750.0},
        {"peptide": "BBBBBBBBB", "allele": "HLA-A*02:01", "affinity": 2000.0},
        {"peptide": "BBBBBBBBB", "allele": "HLA-B*07:02", "affinity": 6000.0},
    ]
    df = _predictions_df(rows)
    ctx = EvalContext(df)
    scores = (
        build_score_node(cfg).eval(ctx).reindex(ctx.group_index).fillna(0.0)
    )
    for row in df.to_dict("records"):
        group_key = (
            row["source_sequence_name"], row["peptide"],
            row["peptide_offset"], row["allele"],
        )
        legacy = _legacy_score(row["affinity"], None, cfg)
        assert float(scores[group_key]) == pytest.approx(legacy, abs=1e-12)


def test_multi_method_default_score_raises_ambiguity_error():
    """Topiary 5.0 explicitly rejects unqualified ``affinity`` when the df
    contains multiple prediction_method_name values. The error points at
    ``affinity['modelname']`` as the remedy. Vaxrank's default wraps a
    single mhctools predictor so this case doesn't fire in practice, but
    the behavior is part of the back-compat contract if a user passes a
    pre-built multi-model TopiaryPredictor and doesn't set score_expr."""
    cfg = EpitopeConfig()
    df = _predictions_df(
        [
            {
                "peptide": "AAAAAAAAA", "allele": "HLA-A*02:01",
                "peptide_offset": 0, "affinity": 100.0,
                "prediction_method_name": "netmhcpan",
            },
            {
                "peptide": "AAAAAAAAA", "allele": "HLA-A*02:01",
                "peptide_offset": 0, "affinity": 100.0,
                "prediction_method_name": "mhcflurry",
            },
        ]
    )
    with pytest.raises(ValueError, match="Ambiguous"):
        build_score_node(cfg).eval(EvalContext(df))


def test_multi_method_resolves_with_qualified_affinity():
    """Users with multiple methods can disambiguate via
    ``score_expr="affinity['netmhcpan'].logistic(350, 150)"`` — verify the
    qualified form picks out the right rows."""
    cfg = EpitopeConfig(score_expr="affinity['netmhcpan'].logistic(350, 150)")
    df = _predictions_df(
        [
            {
                "peptide": "AAAAAAAAA", "allele": "HLA-A*02:01",
                "peptide_offset": 0, "affinity": 100.0,
                "prediction_method_name": "netmhcpan",
            },
            {
                "peptide": "AAAAAAAAA", "allele": "HLA-A*02:01",
                "peptide_offset": 0, "affinity": 50000.0,
                "prediction_method_name": "mhcflurry",
            },
        ]
    )
    ctx = EvalContext(df)
    scores = build_score_node(cfg).eval(ctx).reindex(ctx.group_index)
    # Only the netmhcpan ic50=100 contributes; expected = raw sigmoid at ic50=100
    expected = 1.0 / (1.0 + math.exp((100.0 - 350.0) / 150.0))
    assert float(scores.iloc[0]) == pytest.approx(expected, abs=1e-12)


# ---- candidate_epitopes_from_rows + epitopes_to_topiary_df ---------


def _row(peptide='SIINFEKL', allele='HLA-A*02:01', ic50=50.0,
         wt_peptide=None, wt_ic50=None, predictor='mhcflurry',
         version='2.1.1', source='AAAASIINFEKLCCCC', offset=4,
         percentile_rank=0.5, overlaps_mutation=True,
         occurs_in_reference=False):
    """Concise test helper: build one row dict suitable for
    ``candidate_epitopes_from_rows``. Mirrors the per-row shape that
    ``predict_epitopes`` and the LENS / pVACseq loaders emit.
    ``wt_ic50=None`` skips the WT pairing entirely."""
    from mhctools.pred import Prediction
    mutant = Prediction(
        kind='pMHC_affinity', predictor_name=predictor,
        predictor_version=version, allele=allele, peptide=peptide,
        value=ic50, score=0.0, percentile_rank=percentile_rank,
    )
    wt = None
    if wt_ic50 is not None:
        wt = Prediction(
            kind='pMHC_affinity', predictor_name=predictor,
            predictor_version=version, allele=allele,
            peptide=wt_peptide if wt_peptide is not None else peptide,
            value=wt_ic50, score=0.0, percentile_rank=None,
        )
    return {
        'peptide': peptide, 'source': source, 'offset': offset,
        'mutant': mutant, 'wt': wt,
        'overlaps_mutation': overlaps_mutation,
        'occurs_in_reference': occurs_in_reference,
    }


def test_grouping_collapses_multi_allele_rows():
    """Per-(peptide, source, offset) grouping → one CandidateEpitope.
    Multi-allele rows for the same position collapse into one
    CandidateEpitope whose mutant context carries N predictions."""
    from vaxrank.candidate_epitope import candidate_epitopes_from_rows

    epitopes = candidate_epitopes_from_rows([
        _row(allele='HLA-A*02:01', ic50=50.0),
        _row(allele='HLA-B*07:02', ic50=200.0),
        _row(allele='HLA-C*03:04', ic50=800.0),
    ])
    assert len(epitopes) == 1
    e = epitopes[0]
    assert e.sequence == 'SIINFEKL'
    assert e.offset == 4
    assert sorted(e.alleles_for('pMHC_affinity')) == [
        'HLA-A*02:01', 'HLA-B*07:02', 'HLA-C*03:04']
    assert e.overlaps_mutation is True


def test_grouping_separates_distinct_peptides():
    """Two rows with different peptides → two Epitopes."""
    from vaxrank.candidate_epitope import candidate_epitopes_from_rows

    epitopes = candidate_epitopes_from_rows([
        _row(peptide='SIINFEKL', offset=4, wt_ic50=None),
        _row(peptide='SIINFEKM', offset=12, wt_ic50=None),
    ])
    assert len(epitopes) == 2
    assert {e.sequence for e in epitopes} == {
        'SIINFEKL', 'SIINFEKM'}


def test_wt_comparator_built_when_peptides_differ():
    """A parallel WT comparator context is built when each row carries
    a WT prediction whose peptide differs from the mutant."""
    from vaxrank.candidate_epitope import candidate_epitopes_from_rows

    epitopes = candidate_epitopes_from_rows([
        _row(peptide='SIINFEKL', wt_peptide='SIINFEKM',
             ic50=50.0, wt_ic50=500.0, allele='HLA-A*02:01'),
        _row(peptide='SIINFEKL', wt_peptide='SIINFEKM',
             ic50=100.0, wt_ic50=400.0, allele='HLA-B*07:02'),
    ])
    assert len(epitopes) == 1
    e = epitopes[0]
    assert e.wt is not None
    assert e.wt.sequence == 'SIINFEKM'
    wt_preds = e.wt.predictions_for('pMHC_affinity')
    assert {p.allele for p in wt_preds} == {'HLA-A*02:01', 'HLA-B*07:02'}
    assert {p.value for p in wt_preds} == {500.0, 400.0}


def test_self_wt_is_dropped():
    """A WT Prediction whose peptide equals the mutant peptide is
    dropped as a meaningless self-comparator. The CandidateEpitope has no WT
    context."""
    from mhctools.pred import Prediction
    from vaxrank.candidate_epitope import candidate_epitopes_from_rows

    mutant = Prediction(
        kind='pMHC_affinity', predictor_name='mhcflurry',
        predictor_version='2.1.1', allele='HLA-A*02:01',
        peptide='SIINFEKL', value=50.0, score=0.0, percentile_rank=0.5)
    self_wt = Prediction(
        kind='pMHC_affinity', predictor_name='mhcflurry',
        predictor_version='2.1.1', allele='HLA-A*02:01',
        peptide='SIINFEKL', value=50.0, score=0.0, percentile_rank=None)
    epitopes = candidate_epitopes_from_rows([{
        'peptide': 'SIINFEKL', 'source': 'AAAASIINFEKLCCCC', 'offset': 4,
        'mutant': mutant, 'wt': self_wt, 'overlaps_mutation': False,
    }])
    assert epitopes[0].wt is None
    assert epitopes[0].overlaps_mutation is False


def test_anonymous_wt_kept_when_peptide_empty():
    """LENS / pVACseq inputs sometimes carry a WT IC50 without a WT
    peptide sequence. An empty WT peptide is NOT treated as a
    self-match — the WT context is retained so the IC50 signal
    isn't lost."""
    from mhctools.pred import Prediction
    from vaxrank.candidate_epitope import candidate_epitopes_from_rows

    mutant = Prediction(
        kind='pMHC_affinity', predictor_name='mhcflurry',
        predictor_version='2.1.1', allele='HLA-A*02:01',
        peptide='SIINFEKL', value=50.0, score=0.0, percentile_rank=0.5)
    anon_wt = Prediction(
        kind='pMHC_affinity', predictor_name='mhcflurry',
        predictor_version='2.1.1', allele='HLA-A*02:01',
        peptide='', value=2500.0, score=0.0, percentile_rank=None)
    epitopes = candidate_epitopes_from_rows([{
        'peptide': 'SIINFEKL', 'source': 'AAAASIINFEKLCCCC', 'offset': 4,
        'mutant': mutant, 'wt': anon_wt, 'overlaps_mutation': True,
    }])
    assert epitopes[0].wt is not None
    wt_leaf = epitopes[0].wt.predictions_for('pMHC_affinity')[0]
    assert wt_leaf.value == 2500.0


def test_flags_or_across_rows():
    """``overlaps_mutation`` / ``occurs_in_reference`` OR across rows
    in the same (peptide, source, offset) group. If any row says
    True, the group is True — flags are position-level, not
    leaf-level."""
    from vaxrank.candidate_epitope import candidate_epitopes_from_rows

    epitopes = candidate_epitopes_from_rows([
        _row(allele='HLA-A*02:01', overlaps_mutation=False,
             occurs_in_reference=False),
        _row(allele='HLA-B*07:02', overlaps_mutation=True,
             occurs_in_reference=True),
    ])
    assert len(epitopes) == 1
    assert epitopes[0].overlaps_mutation is True
    assert epitopes[0].occurs_in_reference is True


def test_epitopes_to_topiary_df_emits_one_row_per_prediction():
    """Each leaf ``mhctools.Prediction`` in a CandidateEpitope's mutant
    context becomes one frame row."""
    from vaxrank.epitope_dsl import epitopes_to_topiary_df
    from vaxrank.candidate_epitope import candidate_epitopes_from_rows

    epitopes = candidate_epitopes_from_rows([
        _row(allele='HLA-A*02:01', ic50=50.0),
        _row(allele='HLA-B*07:02', ic50=200.0),
    ])
    df = epitopes_to_topiary_df(epitopes)
    assert len(df) == 2
    assert set(df.columns) >= {
        'peptide', 'allele', 'value', 'affinity', 'percentile_rank',
        'kind', 'prediction_method_name', 'predictor_version',
        'source_sequence_name', 'peptide_offset', 'peptide_length',
        'score', 'n_flank', 'c_flank',
    }
    assert set(df['allele']) == {'HLA-A*02:01', 'HLA-B*07:02'}
    assert (df['kind'] == 'pMHC_affinity').all()


def test_epitopes_to_topiary_df_schema_pinned():
    """Round-trip pin: rows → CandidateEpitope → frame must produce the
    canonical topiary schema (columns + dtype-stable values) that
    ``apply_filter`` / ``score_predictions`` consume."""
    from vaxrank.epitope_dsl import epitopes_to_topiary_df
    from vaxrank.candidate_epitope import candidate_epitopes_from_rows

    epitopes = candidate_epitopes_from_rows([
        _row(allele='HLA-A*02:01', ic50=50.0, percentile_rank=0.3),
        _row(allele='HLA-B*07:02', ic50=200.0, percentile_rank=1.5),
        _row(peptide='SIINFEKM', allele='HLA-A*02:01',
             ic50=800.0, offset=12),
    ])
    df = epitopes_to_topiary_df(epitopes)

    # One row per leaf prediction; alleles are preserved.
    assert len(df) == 3
    assert set(df['allele']) == {'HLA-A*02:01', 'HLA-B*07:02'}
    # ``value`` carries the IC50 verbatim; ``affinity`` is the alias.
    assert sorted(df['value']) == [50.0, 200.0, 800.0]
    assert (df['value'] == df['affinity']).all()
    # Multi-peptide inputs land in distinct rows.
    assert set(df['peptide']) == {'SIINFEKL', 'SIINFEKM'}
    # Every leaf inherits its parent CandidateEpitope's offset.
    by_peptide = df.set_index('peptide')['peptide_offset'].to_dict()
    assert by_peptide['SIINFEKL'] == 4
    assert by_peptide['SIINFEKM'] == 12


def test_allele_free_evidence_still_produces_per_allele_scores():
    """Processing-only evidence must still produce per-allele scores.

    A peptide whose only evidence is allele-free has one group, keyed on an
    empty allele, and the patient's genotype is not in the frame at all — so
    it would score 0 and be dropped, which is the failure mode of
    vaxrank#295. topiary takes the genotype and adds the groups
    (openvax/topiary#182, #219); vaxrank used to duplicate the rows itself.

    The assertion is on the scores rather than on the frame's shape: which
    rows exist is topiary's business and changed once already, while "a
    processing-only candidate scores against both alleles" is the property
    that has to hold either way.
    """
    from mhctools.pred import Prediction

    from vaxrank.candidate_epitope import CandidateEpitope
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import (
        attach_per_allele_scores, epitopes_to_topiary_df,
    )

    alleles = ("HLA-A*02:01", "HLA-B*07:02")
    epitope = CandidateEpitope(
        sequence="SIINFEKL", source_sequence="XXSIINFEKLXX", offset=2,
        prediction_id="processing-only",
        patient_alleles=alleles,
        predictions=(Prediction(
            kind="antigen_processing", predictor_name="mhcflurry",
            predictor_version="2.1.1", allele="", peptide="SIINFEKL",
            value=None, score=0.77),))

    # The canonical leaf stays allele-free, on the object and in the frame:
    # vaxrank states what was predicted and lets topiary add the groups.
    [leaf] = epitope.predictions_for("antigen_processing", predictor="mhcflurry")
    assert leaf.allele == ""
    assert set(epitopes_to_topiary_df([epitope])["allele"]) == {""}

    [scored] = attach_per_allele_scores(
        [epitope],
        EpitopeConfig(score_expr="processing[mhcflurry].score"))
    assert set(scored.per_allele_scores) == set(alleles)
    assert scored.epitope_score > 0


def test_pipeline_and_external_frames_each_group_on_their_own_identity():
    """Both frame shapes name their group keys; neither is left to inference.

    What topiary infers has changed twice (sample_name being prepended, then
    the 5.17.1 group-index fixes). A frame whose grouping shifts underneath it
    merges candidates that are not the same candidate, so vaxrank states the
    identity for each shape instead of relying on the default.
    """
    import pandas as pd

    from vaxrank.epitope_dsl import (
        PIPELINE_GROUP_COLUMNS, PREDICTION_GROUP_COLUMNS,
        prediction_group_columns,
    )

    common = {"peptide": ["SIINFEKL"], "peptide_offset": [0],
              "allele": ["HLA-A*02:01"]}
    # A TopiaryPredictor frame: identity is the source sequence name.
    pipeline = pd.DataFrame({**common, "source_sequence_name": ["ovalbumin"]})
    assert prediction_group_columns(pipeline) == list(PIPELINE_GROUP_COLUMNS)

    # A vaxrank / external-reader frame: identity is the provenance key, and
    # it wins even when a source name is also present.
    external = pd.DataFrame({
        **common, "source_sequence_name": ["ovalbumin"],
        "prediction_id": ["lens:abc"]})
    assert prediction_group_columns(external) == list(
        PREDICTION_GROUP_COLUMNS)

    # Neither identity present: refuse rather than group on sequence alone.
    with pytest.raises(ValueError, match="no usable score-group identity"):
        prediction_group_columns(pd.DataFrame(common))


def test_kind_support_is_forwarded_to_the_dsl():
    """A caller holding a predictor must not make topiary guess.

    ``mhc_dependence`` guards resolve from kind_support when supplied and by
    scanning rows otherwise. vaxrank's pipeline path holds the predictor that
    knows the answer, so leaving it unset made every such guard a guess —
    the same defect topiary fixed on its own --filter-by path.
    """
    from unittest import mock

    import pandas as pd
    import topiary.ranking

    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import score_predictions

    df = pd.DataFrame([{
        "prediction_id": "p1", "source_sequence_name": "ctx",
        "peptide": "SIINFEKL", "peptide_offset": 0, "peptide_length": 8,
        "allele": "HLA-A*02:01", "n_flank": "", "c_flank": "",
        "prediction_method_name": "mhcflurry", "predictor_version": "2.1.1",
        "kind": "pMHC_affinity", "value": 50.0, "affinity": 50.0,
        "percentile_rank": 0.4, "score": 0.0,
    }])
    support = {"mhcflurry": {
        "pMHC_affinity": {"mhc_dependence": "single_allele",
                          "mhc_class": "I"}}}

    real_context = topiary.ranking.EvalContext
    seen = {}

    def spy(frame, **kwargs):
        seen.update(kwargs)
        return real_context(frame, **kwargs)

    with mock.patch.object(topiary.ranking, "EvalContext", spy):
        scores = score_predictions(
            [], EpitopeConfig(), topiary_df=df, kind_support=support)

    assert seen["kind_support"] is support
    # The group identity is named too, not inferred.
    assert seen["group_keys"][0] == "prediction_id"
    assert len(scores) == 1


def test_external_readers_leave_kind_support_unset():
    """No predictor means no kind_support — absent, not fabricated."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import read_lens_report

    import tempfile
    import os

    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "lens.tsv")
        with open(path, "w") as handle:
            handle.write(
                "peptide\tallele\tpep_context\tmhcflurry_2.1.1.aff\n"
                "SIINFEKL\tHLA-A02:01\tXXSIINFEKLXX\t50\n")
        report = read_lens_report(path, epitope_config=EpitopeConfig())

    # Scores still land; dependence resolves from the rows, which is correct
    # for a frame no predictor produced.
    [epitope] = report.epitopes
    assert epitope.per_allele_scores


def test_canonical_method_preference_is_topiary_s_not_a_local_copy():
    """vaxrank must not carry its own idea of which model is canonical.

    It used to, and the two had already drifted: vaxrank ordered
    ``netmhcstabpan`` ahead of ``netmhcpan_el`` / ``netmhcpan_ba`` where
    topiary orders it last. The same frame could therefore resolve to
    different models depending on which table was consulted — exactly the
    disagreement topiary's own docstring warns consumers about.
    """
    import vaxrank.epitope_dsl as dsl

    assert not hasattr(dsl, "_CANONICAL_METHOD_PREFERENCE")
    assert not hasattr(dsl, "_auto_pick_canonical_method")


def test_default_methods_resolution_delegates_but_keeps_the_user_s_choice():
    """topiary picks the canonical model; an explicit config entry overrides.

    The auto-pick is topiary's to define. Which entry the *user* wrote is
    vaxrank's config, and has to win — that is the half topiary has no
    opinion about.
    """
    import pandas as pd
    from topiary.ranking import resolve_default_methods as topiary_resolve

    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import resolve_default_methods

    def row(method, kind="pMHC_affinity", value=50.0):
        return {
            "prediction_id": "p", "source_sequence_name": "ctx",
            "peptide": "SIINFEKL", "peptide_offset": 0, "peptide_length": 8,
            "allele": "HLA-A*02:01", "n_flank": "", "c_flank": "",
            "prediction_method_name": method, "predictor_version": "1",
            "kind": kind, "value": value, "affinity": value,
            "percentile_rank": None, "score": 0.0,
        }

    frame = pd.DataFrame([row("mhcflurry"), row("netmhcpan", value=60.0)])

    # With nothing configured, the answer is topiary's, verbatim.
    assert resolve_default_methods(EpitopeConfig(), frame) == dict(
        topiary_resolve(frame))

    # An explicit entry overrides the canonical pick.
    configured = resolve_default_methods(
        EpitopeConfig(default_methods={"pMHC_affinity": "netmhcpan"}), frame)
    assert configured["pMHC_affinity"] == "netmhcpan"


def test_invalid_default_methods_names_the_config_key():
    """The error must say which setting to edit, not only what is missing.

    topiary reports what the data lacks — the half it can know. Which config
    key produced it is the half it cannot, and is what the user needs.
    """
    import pandas as pd
    import pytest

    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import validate_default_methods

    frame = pd.DataFrame([{
        "prediction_id": "p", "source_sequence_name": "ctx",
        "peptide": "SIINFEKL", "peptide_offset": 0, "peptide_length": 8,
        "allele": "HLA-A*02:01", "n_flank": "", "c_flank": "",
        "prediction_method_name": "mhcflurry", "predictor_version": "1",
        "kind": "pMHC_affinity", "value": 50.0, "affinity": 50.0,
        "percentile_rank": None, "score": 0.0,
    }])
    cfg = EpitopeConfig(default_methods={"pMHC_affinity": "mhcnuggets"})
    with pytest.raises(
            ValueError,
            match=r"default_methods\['pMHC_affinity'\]='mhcnuggets'"):
        validate_default_methods(cfg, frame)


def test_canonical_pick_is_announced_only_when_the_user_did_not_choose():
    """Don't tell someone to set an entry they have already set.

    The log exists to make an implicit choice explicit. For a kind the user
    configured there is no implicit choice, and saying "no default_methods
    entry names one" is simply false.
    """
    import logging

    import pandas as pd

    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import resolve_default_methods

    def row(method, value):
        return {
            "prediction_id": "p", "source_sequence_name": "ctx",
            "peptide": "SIINFEKL", "peptide_offset": 0, "peptide_length": 8,
            "allele": "HLA-A*02:01", "n_flank": "", "c_flank": "",
            "prediction_method_name": method, "predictor_version": "1",
            "kind": "pMHC_affinity", "value": value, "affinity": value,
            "percentile_rank": None, "score": 0.0,
        }

    frame = pd.DataFrame([row("mhcflurry", 50.0), row("netmhcpan", 60.0)])
    messages = []

    class _Capture(logging.Handler):
        def emit(self, record):
            messages.append(record.getMessage())

    logger = logging.getLogger("vaxrank.epitope_dsl")
    handler = _Capture()
    logger.addHandler(handler)
    previous = logger.level
    logger.setLevel(logging.INFO)
    try:
        configured = resolve_default_methods(
            EpitopeConfig(default_methods={"pMHC_affinity": "netmhcpan"}),
            frame)
        assert configured == {"pMHC_affinity": "netmhcpan"}
        assert not messages

        messages.clear()
        auto = resolve_default_methods(EpitopeConfig(), frame)
        assert auto == {"pMHC_affinity": "mhcflurry"}
        assert any("canonical pick" in message for message in messages)
    finally:
        logger.removeHandler(handler)
        logger.setLevel(previous)


def test_typo_in_a_default_methods_kind_is_rejected():
    """A kind the data does not have is a config error, not a no-op."""
    import pandas as pd
    import pytest

    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import validate_default_methods

    frame = pd.DataFrame([{
        "prediction_id": "p", "source_sequence_name": "ctx",
        "peptide": "SIINFEKL", "peptide_offset": 0, "peptide_length": 8,
        "allele": "HLA-A*02:01", "n_flank": "", "c_flank": "",
        "prediction_method_name": "mhcflurry", "predictor_version": "1",
        "kind": "pMHC_affinity", "value": 50.0, "affinity": 50.0,
        "percentile_rank": None, "score": 0.0,
    }])
    with pytest.raises(ValueError, match=r"pMHC_affinty"):
        validate_default_methods(
            EpitopeConfig(default_methods={"pMHC_affinty": "mhcflurry"}),
            frame)


@pytest.mark.parametrize("kind", sorted(
    k for k, v in topiary.KIND_MHC_DEPENDENCE.items() if v == "none"))
def test_every_allele_free_kind_reaches_every_patient_allele(kind):
    """Allele-free evidence must score against the patient's alleles.

    A prediction that describes the peptide rather than a peptide-MHC pair
    carries no allele, so it forms a group keyed on the empty string and
    contributes to no per-allele score — the candidate scores 0 and is
    dropped (openvax/topiary#182). vaxrank sidesteps that by projecting such
    evidence across the patient's alleles when it builds the frame.

    That projection used to name ``antigen_processing`` literally, because
    the set of allele-free kinds lived in a private topiary table. There are
    five, so the other four were silently dropped from scoring: mhctools
    ships netcleave, deeptap, eramer and proteasome_predictor, all of which
    emit them. No LENS fixture carries one, which is why nothing caught it.
    """
    from mhctools import Prediction

    from vaxrank.candidate_epitope import CandidateEpitope
    from vaxrank.epitope_dsl import attach_per_allele_scores

    alleles = ("HLA-A*02:01", "HLA-B*07:02")
    epitope = CandidateEpitope(
        sequence="SIINFEKL", offset=0, patient_alleles=alleles,
        predictions=(Prediction(
            kind=kind, score=0.8, peptide="SIINFEKL", allele="",
            predictor_name="mhcflurry", predictor_version="2.1.1"),))

    [scored] = attach_per_allele_scores(
        [epitope],
        EpitopeConfig(score_expr="%s[mhcflurry].score" % kind))

    assert dict(scored.per_allele_scores) == {a: 0.8 for a in alleles}


def test_a_candidate_is_never_scored_against_another_candidates_alleles():
    """Each candidate's genotype is its own, not the frame's union.

    LENS emits one row per (peptide, allele) that passed its own threshold,
    so candidates in one file arrive having been reported against different
    alleles — every LENS fixture here holds between two and eight distinct
    sets. Declaring the union frame-wide would add a group pairing each
    peptide with alleles it was never scored against.

    The union is invisible under the default score expression, because it
    contains an affinity term and the invented groups carry no affinity
    rows, so they read NaN and score nothing. An expression reading only
    peptide-level evidence puts a real number in every one of them, which
    is what this pins.
    """
    from mhctools.pred import Prediction

    from vaxrank.candidate_epitope import CandidateEpitope
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_dsl import attach_per_allele_scores

    def candidate(name, peptide, allele):
        return CandidateEpitope(
            sequence=peptide, offset=0, prediction_id=name,
            patient_alleles=(allele,),
            predictions=(
                Prediction(
                    kind="pMHC_affinity", predictor_name="mhcflurry",
                    predictor_version="2.1.1", allele=allele, peptide=peptide,
                    value=50.0, score=None),
                Prediction(
                    kind="proteasome_cleavage", predictor_name="mhcflurry",
                    predictor_version="2.1.1", allele="", peptide=peptide,
                    value=None, score=0.8)))

    epitopes = [
        candidate("a", "SIINFEKLA", "HLA-A*02:01"),
        candidate("b", "AAAAAAAAA", "HLA-B*07:02"),
    ]

    scored = attach_per_allele_scores(
        epitopes,
        EpitopeConfig(score_expr="proteasome_cleavage[mhcflurry].score"))

    # Each candidate scores against its own allele only. Under a frame-wide
    # union both would carry both alleles — four scores, two of them for
    # pairings no model ever saw.
    assert [dict(e.per_allele_scores) for e in scored] == [
        {"HLA-A*02:01": pytest.approx(0.8)},
        {"HLA-B*07:02": pytest.approx(0.8)},
    ]


@pytest.mark.parametrize("stored", ["nan", "NaN", "", "   ", None, np.nan])
def test_a_pin_matching_no_named_version_is_rejected(stored):
    """Validation must not accept a pin that selects nothing.

    A cell that records no version is NaN, None, blank — or the literal
    string "nan", which survives dropna() and a truthiness test. That last
    shape once counted as a real version in vaxrank's resolver and dropped
    the rows carrying it from scoring.

    Here the consequence is the shape this validator exists to prevent:
    ``affinity['netmhcpan', 'nan']`` would pass validation and then select
    nothing, so the run reports a score of zero for a formula that looked
    checked.
    """
    from vaxrank.epitope_dsl import validate_dsl_against_predictions

    df = pd.DataFrame([{
        "prediction_id": "p", "source_sequence_name": "s",
        "peptide": "SIINFEKL", "peptide_offset": 0, "peptide_length": 8,
        "allele": "HLA-A*02:01", "prediction_method_name": "netmhcpan",
        "predictor_version": stored, "kind": "pMHC_affinity",
        "value": 50.0, "score": None}])

    with pytest.raises(ValueError, match="predictor version"):
        validate_dsl_against_predictions(
            EpitopeConfig(score_expr="affinity['netmhcpan', 'nan'].value"),
            [], topiary_df=df)

UNSTATED_SPELLINGS = [
    pytest.param(None, id="None"),          # mhctools
    pytest.param(np.nan, id="NaN"),         # a frame column
    pytest.param("", id="empty"),           # vaxrank's own builders
    pytest.param("nan", id="nan-string"),   # anything through astype(str)
]


def _affinity(version, value):
    from mhctools.pred import Prediction

    return Prediction(
        kind="pMHC_affinity", predictor_name="netmhcpan",
        predictor_version=version, allele="HLA-A*02:01",
        peptide="SIINFEKL", value=value, score=None)


def _candidate(*predictions):
    from vaxrank.candidate_epitope import CandidateEpitope

    return CandidateEpitope(
        sequence="SIINFEKL", offset=0, prediction_id="x",
        patient_alleles=("HLA-A*02:01",), predictions=predictions)


@pytest.mark.parametrize("missing", UNSTATED_SPELLINGS)
def test_every_spelling_of_no_version_reaches_the_frame_the_same_way(missing):
    """One meaning, one representation.

    Four producers spell "this prediction has no version" four ways, and
    they used to reach the frame as three different values:
    ``p.predictor_version or ""`` passes NaN straight through, because NaN
    is truthy, and passes the string "nan" through as if it named a version.
    """
    from vaxrank.epitope_dsl import epitopes_to_topiary_df

    frame = epitopes_to_topiary_df([_candidate(_affinity(missing, 50.0))])

    assert list(frame["predictor_version"]) == [""]


@pytest.mark.parametrize("second", UNSTATED_SPELLINGS)
@pytest.mark.parametrize("first", UNSTATED_SPELLINGS)
def test_two_spellings_of_no_version_are_one_bucket(first, second):
    """Unversioned predictions from two producers form one group.

    The nested store keys on ``predictor_version`` directly, so each
    spelling opened its own bucket — one predictor's evidence split in two,
    where unqualified access resolves to a single version and sees one half.
    """
    epitope = _candidate(_affinity(first, 50.0), _affinity(second, 60.0))

    buckets = epitope.predictions["pMHC_affinity"]["netmhcpan"]
    assert list(buckets) == [""]
    assert len(buckets[""]) == 2


@pytest.mark.parametrize("missing", UNSTATED_SPELLINGS)
def test_a_named_version_beside_an_unnamed_one_does_not_raise(missing):
    """Mixed version types must not break candidate construction.

    ``predictions_flat`` sorts on ``predictor_version``, so a predictor
    reporting a version on one record and nothing on another raised
    ``TypeError: '<' not supported between instances of 'NoneType' and
    'str'`` from the constructor, before any scoring happened.
    """
    epitope = _candidate(_affinity("4.1b", 50.0), _affinity(missing, 60.0))

    # The named version keeps its own bucket; only the unnamed collapse.
    assert sorted(epitope.predictions["pMHC_affinity"]["netmhcpan"]) == [
        "", "4.1b"]
    assert len(epitope.predictions_flat()) == 2


def _processing(allele, score):
    from mhctools.pred import Prediction

    return Prediction(
        kind="antigen_processing", predictor_name="mhcflurry",
        predictor_version="2.1.1", allele=allele, peptide="SIINFEKL",
        value=None, score=score)


@pytest.mark.parametrize("unstated", UNSTATED_SPELLINGS)
def test_an_unstated_allele_never_becomes_a_patient_allele(unstated):
    """A blank allele means peptide-level evidence, not a genotype entry.

    ``__post_init__`` derives ``patient_alleles`` from the alleles its
    predictions carry, filtered on truthiness — and two spellings of blank
    are truthy. NaN sorted against the real allele names and raised
    ``TypeError``; the string "nan" became a patient allele that no model
    ever scored, which is a per-allele score for a peptide-MHC pair that was
    never predicted.
    """
    from vaxrank.candidate_epitope import CandidateEpitope

    epitope = CandidateEpitope(
        sequence="SIINFEKL", offset=0, prediction_id="x",
        patient_alleles=("HLA-A*02:01",),
        predictions=(_processing(unstated, 0.8),))

    assert epitope.patient_alleles == ("HLA-A*02:01",)


@pytest.mark.parametrize("unstated", UNSTATED_SPELLINGS)
def test_an_unstated_allele_reaches_the_frame_as_blank(unstated):
    """One spelling in the frame, so one group.

    topiary keys its groups on the allele column, so a stringified null
    there is a group of its own — allele-free evidence split across buckets
    that peptide_view cannot broadcast from.
    """
    from vaxrank.candidate_epitope import CandidateEpitope
    from vaxrank.epitope_dsl import epitopes_to_topiary_df

    epitope = CandidateEpitope(
        sequence="SIINFEKL", offset=0, prediction_id="x",
        patient_alleles=("HLA-A*02:01",),
        predictions=(_processing(unstated, 0.8),))

    assert list(epitopes_to_topiary_df([epitope])["allele"]) == [""]


@pytest.mark.parametrize("unstated", UNSTATED_SPELLINGS)
def test_a_named_allele_beside_an_unstated_one_does_not_raise(unstated):
    """Mixed allele types must not break candidate construction.

    Same shape as the version axis: both are sorted, and a producer that
    spells blank differently from vaxrank's own builders made ``sorted``
    compare a float to a string.
    """
    from vaxrank.candidate_epitope import CandidateEpitope

    epitope = CandidateEpitope(
        sequence="SIINFEKL", offset=0, prediction_id="x",
        patient_alleles=(),
        predictions=(_processing("HLA-A*02:01", 0.9),
                     _processing(unstated, 0.8)))

    assert epitope.patient_alleles == ("HLA-A*02:01",)
    assert len(epitope.predictions_flat()) == 2


@pytest.mark.parametrize("unstated", UNSTATED_SPELLINGS)
def test_the_record_itself_is_normalized_not_just_the_keys(unstated):
    """Every reader of ``p.allele`` must see a stated value or blank.

    The first pass applied the rule only where vaxrank keyed or ordered on
    these fields, which fixed the buckets and the sorts and nothing else —
    around sixty call sites read ``p.allele`` directly, and every one that
    asks ``if p.allele`` to mean "is this allele-scoped" is wrong for NaN
    and for the string "nan", both truthy. ``alleles_for`` is one: it
    returned a phantom allele rather than the empty tuple.
    """
    from vaxrank.candidate_epitope import CandidateEpitope

    epitope = CandidateEpitope(
        sequence="SIINFEKL", offset=0, prediction_id="x",
        patient_alleles=("HLA-A*02:01",),
        predictions=(_processing(unstated, 0.8),))

    [leaf] = epitope.predictions_flat()
    assert leaf.allele == ""
    assert epitope.alleles_for("antigen_processing") == ()


@pytest.mark.parametrize("unstated", UNSTATED_SPELLINGS)
def test_a_caller_supplied_frame_cannot_add_a_phantom_allele(unstated):
    """per_allele_scores must reject an unstated allele from any frame.

    Normalizing predictions on the way into a CandidateEpitope does not
    cover this path: pVACseq passes its own topiary frame, which never went
    through that boundary. topiary collapses a stringified null in a group
    key to a real null, and NaN is truthy, so a bare ``if allele`` wrote it
    into a map whose contract is per patient allele.
    """
    from vaxrank.candidate_epitope import CandidateEpitope
    from vaxrank.epitope_dsl import attach_per_allele_scores

    frame = pd.DataFrame([{
        "prediction_id": "x", "source_sequence_name": "s",
        "peptide": "SIINFEKL", "peptide_offset": 0, "peptide_length": 8,
        "allele": unstated, "prediction_method_name": "mhcflurry",
        "predictor_version": "2.1.1", "kind": "antigen_processing",
        "value": None, "score": 0.8}])
    epitope = CandidateEpitope(
        sequence="SIINFEKL", offset=0, prediction_id="x",
        patient_alleles=("HLA-A*02:01",), predictions=())

    [scored] = attach_per_allele_scores(
        [epitope],
        EpitopeConfig(
            score_expr="peptide_view(processing[mhcflurry].score)"),
        topiary_df=frame)

    assert dict(scored.per_allele_scores) == {
        "HLA-A*02:01": pytest.approx(0.8)}


def test_an_unavailable_column_fails_when_the_config_is_read(tmp_path):
    """Naming a column the input lacks must fail early and name the setting.

    topiary raises on this too, with a better "did you mean" list than
    vaxrank would generate — but at eval time, mid-run, and without saying
    which config key produced the reference. That is the half topiary
    cannot know and the half the user needs (#388).
    """
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import read_lens_report

    path = tmp_path / "annotated.tsv"
    path.write_text(
        "peptide\tallele\tpep_context\ttpm\tmhcflurry_2.1.1.aff\n"
        "SIINFEKL\tHLA-A02:01\tXXSIINFEKLXX\t50\t20\n")

    config = EpitopeConfig(
        score_expr=(
            "affinity.value.logistic_normalized(350,150) * (nope > 1)"))

    with pytest.raises(ValueError) as excinfo:
        read_lens_report(path, epitope_config=config)

    message = str(excinfo.value)
    assert "epitopes.score_expr" in message      # which key to edit
    assert "nope" in message                     # what was wrong
    assert "gene_tpm" in message                 # what is available


def test_a_column_the_input_does_carry_is_accepted(tmp_path):
    """The check must not reject a working expression.

    ``gene_tpm`` is on a LENS frame, so a filter naming it is valid — the
    point is to reject references the input cannot answer, not to narrow
    what an expression may name.
    """
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import read_lens_report

    path = tmp_path / "annotated.tsv"
    path.write_text(
        "peptide\tallele\tpep_context\ttpm\tmhcflurry_2.1.1.aff\n"
        "SIINFEKL\tHLA-A02:01\tXXSIINFEKLXX\t50\t20\n")

    loaded = read_lens_report(
        path,
        epitope_config=EpitopeConfig(
            score_expr=(
                "affinity.value.logistic_normalized(350,150) * "
                "(gene_tpm > 1)")))

    assert list(loaded.epitopes)


def test_the_failing_setting_is_named_even_when_both_are_set(tmp_path):
    """filter_expr and score_expr are distinguished, not conflated."""
    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_io import read_lens_report

    path = tmp_path / "annotated.tsv"
    path.write_text(
        "peptide\tallele\tpep_context\ttpm\tmhcflurry_2.1.1.aff\n"
        "SIINFEKL\tHLA-A02:01\tXXSIINFEKLXX\t50\t20\n")

    config = EpitopeConfig(
        filter_expr="missing_here > 0",
        score_expr="affinity.value.logistic_normalized(350,150)")

    with pytest.raises(ValueError, match="epitopes.filter_expr"):
        read_lens_report(path, epitope_config=config)
