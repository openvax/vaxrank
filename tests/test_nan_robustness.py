# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Regression tests for #289.

Two layered guarantees:

1. ``finite_prediction_value`` coerces NaN / Inf to ``None`` so
   ``Prediction.value`` never carries NaN coming out of the topiary
   frame. The topiary frame legitimately emits NaN in the ``value``
   column for non-affinity kinds (``pMHC_presentation`` carries its
   probability in ``score``, not ``value``); the producer must
   normalize these on the way in or every downstream consumer pays.

2. ``cli.entry_point._to_json_nan_tolerant`` renders NaN / Inf as
   JSON ``null`` rather than crashing the writer. Defense in depth:
   even if a future producer leaks a NaN, the whole report doesn't
   get lost on the last step.
"""
import json
import math
from unittest.mock import patch

import pandas as pd
import pytest


# ---- finite_prediction_value -------------------------------------------


def test_finite_prediction_value_passes_through_finite():
    from vaxrank import finite_prediction_value

    assert finite_prediction_value(0.0) == 0.0
    assert finite_prediction_value(123.456) == 123.456
    assert finite_prediction_value(-1e9) == -1e9


def test_finite_prediction_value_coerces_nan_and_inf():
    from vaxrank import finite_prediction_value

    assert finite_prediction_value(float('nan')) is None
    assert finite_prediction_value(float('inf')) is None
    assert finite_prediction_value(float('-inf')) is None


def test_finite_prediction_value_passes_none():
    from vaxrank import finite_prediction_value

    assert finite_prediction_value(None) is None


def test_finite_prediction_value_accepts_ints_and_numeric_strings():
    from vaxrank import finite_prediction_value

    assert finite_prediction_value(0) == 0.0
    assert finite_prediction_value(42) == 42.0
    assert finite_prediction_value("42.5") == 42.5


def test_finite_prediction_value_rejects_non_numeric_input():
    from vaxrank import finite_prediction_value

    with pytest.raises(ValueError, match="not numeric"):
        finite_prediction_value("strong")


def test_prediction_integer_accepts_pandas_style_integral_float_only():
    from vaxrank import prediction_integer

    assert prediction_integer(9.0, "Peptide length") == 9
    with pytest.raises(ValueError, match="Peptide length must be an integer"):
        prediction_integer(9.5, "Peptide length")
    with pytest.raises(ValueError, match="Peptide length must be an integer"):
        prediction_integer(True, "Peptide length")


# ---- predict_epitopes: presentation rows don't leak NaN ----------------


def _stub_topiary_predictions_df():
    """Build a minimal topiary-shaped frame with both affinity (real
    IC50 in ``value``) and presentation (``value=NaN``) rows for the
    same peptide × allele. Mirrors what mhctools.MHCflurry emits via
    TopiaryPredictor on real input."""
    rows = []
    # Affinity row — has IC50 in 'value' and 'affinity'
    rows.append({
        'peptide': 'KLQGHSAPV', 'peptide_length': 9,
        'peptide_offset': 0, 'allele': 'HLA-A*02:01',
        'kind': 'pMHC_affinity',
        'value': 50.0, 'affinity': 50.0,
        'score': 0.5, 'percentile_rank': 0.5,
        'prediction_method_name': 'mhcflurry',
        'predictor_version': '2.1.1',
        'source_sequence_name': 'gene',
    })
    # Presentation row — value/affinity NaN; score carries the signal
    rows.append({
        'peptide': 'KLQGHSAPV', 'peptide_length': 9,
        'peptide_offset': 0, 'allele': 'HLA-A*02:01',
        'kind': 'pMHC_presentation',
        'value': float('nan'), 'affinity': float('nan'),
        'score': 0.42, 'percentile_rank': 1.7,
        'prediction_method_name': 'mhcflurry',
        'predictor_version': '2.1.1',
        'source_sequence_name': 'gene',
    })
    return pd.DataFrame(rows)


def test_predict_epitopes_does_not_emit_nan_value_for_presentation():
    """Regression for #289: a topiary frame with NaN in ``value``
    for ``pMHC_presentation`` rows must yield ``Prediction.value=None``,
    not ``Prediction.value=NaN``."""
    from topiary import TopiaryPredictor
    from varcode import Variant

    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.epitope_logic import predict_epitopes
    from vaxrank.mutant_protein_fragment import MutantProteinFragment

    fake_df = _stub_topiary_predictions_df()

    class _StubTopiary(TopiaryPredictor):
        def __init__(self):
            pass

        def predict_from_named_sequences(self, _named_sequences):
            return fake_df

    fragment = MutantProteinFragment(
        variant=Variant('X', '1000', 'C', 'A', ensembl=None),
        gene_name='gene',
        amino_acids='KLQGHSAPVLDVIVNCDESLLASSD',
        mutant_amino_acid_start_offset=12,
        mutant_amino_acid_end_offset=13,
        n_overlapping_reads=10, n_alt_reads=5, n_ref_reads=5,
        n_alt_reads_supporting_protein_sequence=5,
        supporting_reference_transcripts=[])

    cfg = EpitopeConfig(
        min_epitope_score=0,
        default_methods={'pMHC_affinity': 'mhcflurry'})
    with patch('vaxrank.epitope_logic.ReferenceProteome'):
        epitopes = predict_epitopes(
            mhc_predictor=_StubTopiary(),
            protein_fragment=fragment,
            epitope_config=cfg,
            genome=None)

    # Every leaf Prediction must have a finite or None ``value`` —
    # never NaN.
    nan_records = []
    for e in epitopes:
        for p in e.predictions_flat():
            if isinstance(p.value, float) and math.isnan(p.value):
                nan_records.append(p)
    assert nan_records == [], (
        f'Expected no NaN in Prediction.value; found {len(nan_records)}: '
        f'{nan_records[:3]}')

    # And the presentation prediction must come through as None (we
    # dropped the NaN, we didn't replace it with the affinity value).
    presentation_preds = [
        p for e in epitopes for p in e.predictions_flat()
        if p.kind == 'pMHC_presentation']
    assert presentation_preds, "expected pMHC_presentation leaf in output"
    assert all(p.value is None for p in presentation_preds), (
        "pMHC_presentation rows should have value=None, not NaN")


# ---- _to_json_nan_tolerant ---------------------------------------------


def test_json_writer_renders_nan_as_null():
    """Defense-in-depth: if a NaN slips past the producer guards, the
    JSON writer renders it as ``null`` rather than crashing."""
    from vaxrank.cli.entry_point import _to_json_nan_tolerant

    payload = {'finite': 1.5, 'nan': float('nan'),
               'inf': float('inf'), 'neg_inf': float('-inf'),
               'nested': {'also_nan': float('nan')}}
    out = _to_json_nan_tolerant(payload)
    parsed = json.loads(out)
    assert parsed['finite'] == 1.5
    assert parsed['nan'] is None
    assert parsed['inf'] is None
    assert parsed['neg_inf'] is None
    assert parsed['nested']['also_nan'] is None


def test_json_writer_handles_list_of_floats():
    """Lists of floats with embedded NaN also serialize cleanly."""
    from vaxrank.cli.entry_point import _to_json_nan_tolerant

    payload = {'values': [1.0, float('nan'), 2.0, float('inf')]}
    parsed = json.loads(_to_json_nan_tolerant(payload))
    assert parsed['values'] == [1.0, None, 2.0, None]


def test_json_writer_round_trips_real_payload():
    """A finite payload that doesn't trigger the NaN path should
    survive the writer with the same content."""
    from vaxrank.cli.entry_point import _to_json_nan_tolerant

    payload = {'a': 1, 'b': 'string', 'c': [1, 2, 3], 'd': {'e': 4.5}}
    parsed = json.loads(_to_json_nan_tolerant(payload))
    assert parsed == payload


# ---- stdlib json baseline (for documentation; not a real test) ---------


def test_simplejson_default_rejects_nan_as_documented():
    """Sanity: confirms the bug premise. The stdlib + simplejson
    *default* path rejects NaN — which is why we wrapped with
    ``ignore_nan=True``. This test fails the day simplejson changes
    its default; if it ever does, we can remove the wrapper."""
    import simplejson
    with pytest.raises(ValueError, match='Out of range'):
        simplejson.dumps({'x': float('nan')})
