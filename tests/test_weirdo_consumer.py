"""Actual released category scores remain opt-in and retain their semantics."""

from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from mhctools import Kind, Prediction
from topiary import TopiaryResult
from vaxrank.candidate_epitope import CandidateEpitope
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.epitope_dsl import (
    attach_per_allele_scores, epitopes_to_topiary_df, score_predictions,
)


def test_released_weirdo_changes_only_opted_in_scores_and_replays(tmp_path, monkeypatch):
    pytest.importorskip('weirdo')
    from mhctools import WEIRDO

    model = WEIRDO()
    peptides = ['ACDEFGHIK', 'MTMDKSELV']
    observations = [result.foreignness for result in model.predict(peptides)]
    scores = {p.peptide: p.score for p in observations}
    high, low = max(scores, key=scores.get), min(scores, key=scores.get)
    assert high != low
    allele = 'HLA-A*02:01'
    candidates = [CandidateEpitope(
        sequence=p.peptide, prediction_id=p.peptide, patient_alleles=(allele,),
        predictions=(p, Prediction(
            kind=Kind.pMHC_affinity, peptide=p.peptide, allele=allele,
            value=50. if p.peptide == low else 500., score=0.5,
            predictor_name='binding_fixture', predictor_version='1')),
    ) for p in observations]
    default = EpitopeConfig(min_epitope_score=0)
    baseline = score_predictions(candidates, default)
    without = [replace(e, predictions=tuple(p for p in e.predictions_flat()
                                           if p.kind != Kind.foreignness)) for e in candidates]
    pd.testing.assert_series_equal(
        baseline.xs(allele, level='allele'),
        score_predictions(without, default).xs(allele, level='allele'))
    assert [e.per_allele_scores for e in attach_per_allele_scores(candidates, default)] == [
        e.per_allele_scores for e in attach_per_allele_scores(without, default)]
    assert baseline.idxmax()[1] == low

    opted_in = EpitopeConfig(score_expr='peptide_view(foreignness[weirdo].score)', min_epitope_score=0)
    selected = score_predictions(candidates, opted_in)
    assert selected.idxmax()[1] == high
    np.testing.assert_allclose(selected, [scores[key[1]] for key in selected.index])
    scored = attach_per_allele_scores(candidates, opted_in)
    dataset = EpitopeDataset.from_predictions(scored, config=opted_in)
    frame = dataset.result.df
    rows = frame.loc[frame.kind.eq(Kind.foreignness)]
    assert rows.allele.eq('').all() and len(rows) == 2
    assert rows.predictor_version.eq(model.predictor_version).all()
    assert all(c['score_semantics'] == 'human-self-complement-v1' for c in rows.measurement_context)

    path = tmp_path / 'native.tsv'
    dataset.save(path)
    monkeypatch.setattr(WEIRDO, 'predict', lambda *args: pytest.fail('replay must not run the model'))
    restored = EpitopeDataset.load(path)
    pd.testing.assert_frame_equal(epitopes_to_topiary_df(restored.epitopes), frame)
    pd.testing.assert_series_equal(score_predictions(restored.epitopes, restored.config), selected)
    assert [e.per_allele_scores for e in restored.epitopes] == [e.per_allele_scores for e in scored]


def test_released_weirdo_table_import_preserves_mutant_and_comparator_context():
    pytest.importorskip('weirdo')
    from mhctools import WEIRDO
    from topiary import TopiaryPredictor

    model = WEIRDO()
    frame = TopiaryPredictor(models=model).predict_from_named_peptides({'p': 'ACDEFGHIK'})
    comparator = model.predict(['ACDEFGHIL'])[0].foreignness
    frame['wt_peptide'] = comparator.peptide
    frame['wt_score'] = comparator.score
    frame['wt_predictor_version'] = comparator.predictor_version
    frame['wt_measurement_context'] = [comparator.measurement_context.to_dict()]
    dataset = EpitopeDataset.from_topiary(TopiaryResult(frame), sample_name='test-patient')
    candidate, = dataset.epitopes
    prediction, = candidate.predictions_flat()
    wt, = candidate.wt.predictions_flat()
    assert prediction.measurement_context.to_dict() == frame.measurement_context.iloc[0]
    assert wt.measurement_context == comparator.measurement_context
    assert prediction.predictor_version == wt.predictor_version == model.predictor_version
    assert prediction.allele == wt.allele == ''
