import hashlib
import math

import pytest

from vaxrank.file_digest import file_sha256
from vaxrank.prediction_input import prediction_row_fields


def test_file_sha256_matches_hashlib_across_chunk_boundaries(tmp_path):
    payload = bytes(range(256)) * 17
    path = tmp_path / "payload.bin"
    path.write_bytes(payload)
    expected = hashlib.sha256(payload).hexdigest()
    assert file_sha256(path) == expected
    assert file_sha256(str(path), chunk_size=7) == expected


def test_file_sha256_of_empty_file(tmp_path):
    path = tmp_path / "empty"
    path.write_bytes(b"")
    assert file_sha256(path) == hashlib.sha256(b"").hexdigest()


def test_file_sha256_missing_file_raises_oserror(tmp_path):
    with pytest.raises(OSError):
        file_sha256(tmp_path / "missing")


def test_prediction_row_fields_prefers_affinity_and_defaults_kind():
    fields = prediction_row_fields(
        {"affinity": 50.0, "value": 7.0, "score": float("nan"),
         "prediction_method_name": "m", "predictor_version": "1",
         "allele": "HLA-A*02:01", "percentile_rank": 0.5},
        default_kind="pMHC_affinity",
    )
    assert fields == dict(
        kind="pMHC_affinity", predictor_name="m", predictor_version="1",
        allele="HLA-A*02:01", score=None, value=50.0, percentile_rank=0.5)


def test_prediction_row_fields_falls_back_to_value_and_blank_kind():
    fields = prediction_row_fields({"affinity": math.nan, "value": 7, "kind": ""})
    assert fields["value"] == 7.0
    assert fields["kind"] == ""
    assert fields["predictor_name"] == ""


def test_prediction_row_fields_rejects_non_numeric_score():
    with pytest.raises(ValueError):
        prediction_row_fields({"score": "high"})
