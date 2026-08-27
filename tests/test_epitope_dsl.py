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


def test_allele_free_evidence_is_projected_onto_patient_alleles():
    """Processing-only evidence must still produce per-allele scores.

    topiary 5.18's ``peptide_view()`` broadcasts an allele-free value across
    allele groups that already exist, but a peptide whose only evidence is
    allele-free has no such groups — its single group is keyed on an empty
    allele, and the patient's genotype is not in the frame at all
    (openvax/topiary#182). Dropping this projection in favor of
    ``peptide_view()`` would silently score such candidates 0 and drop them,
    which is the failure mode of vaxrank#295.
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

    # The canonical leaf stays allele-free on the object ...
    [leaf] = epitope.predictions_for("antigen_processing", predictor="mhcflurry")
    assert leaf.allele == ""

    # ... and is projected onto both patient alleles in the evaluation frame.
    frame = epitopes_to_topiary_df([epitope])
    assert set(frame["allele"]) == set(alleles)

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
