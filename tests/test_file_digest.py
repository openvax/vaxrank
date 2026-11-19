import hashlib
import math

import pytest

from vaxrank.file_digest import (
    canonical_json_sha256,
    file_sha256,
    is_sha256_digest,
    text_sha256,
)
from vaxrank.prediction_input import prediction_row_fields
from vaxrank.report import _format_fixed_width_table, _to_roman


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


@pytest.mark.parametrize('chunk_size', [0, -1, True, 1.5, '7'])
def test_file_sha256_rejects_invalid_chunk_size(tmp_path, chunk_size):
    path = tmp_path / 'nonempty'
    path.write_bytes(b'content must not receive an empty-file digest')
    with pytest.raises(ValueError, match='positive integer'):
        file_sha256(path, chunk_size=chunk_size)


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


@pytest.mark.parametrize("value", [
    "0" * 64,
    "0123456789abcdef" * 4,
    hashlib.sha256(b"vaxrank").hexdigest(),
])
def test_is_sha256_digest_accepts_lowercase_hex64(value):
    assert is_sha256_digest(value) is True


@pytest.mark.parametrize("value", [
    "",
    "a" * 63,
    "a" * 65,
    "A" * 64,
    "g" * 64,
    "0" * 63 + " ",
    None,
    12345,
    b"0" * 64,
])
def test_is_sha256_digest_rejects_invalid_values(value):
    assert is_sha256_digest(value) is False


def test_text_and_canonical_json_sha256():
    assert text_sha256("hello") == hashlib.sha256(b"hello").hexdigest()
    assert canonical_json_sha256({"b": 1, "a": 2}) == canonical_json_sha256({"a": 2, "b": 1})
    with pytest.raises(ValueError):
        canonical_json_sha256({"bad": float("nan")})


@pytest.mark.parametrize("number,expected", [
    (1, "I"),
    (4, "IV"),
    (9, "IX"),
    (14, "XIV"),
    (42, "XLII"),
    (99, "XCIX"),
    (400, "CD"),
    (944, "CMXLIV"),
    (1999, "MCMXCIX"),
    (3999, "MMMCMXCIX"),
])
def test_to_roman(number, expected):
    assert _to_roman(number) == expected


@pytest.mark.parametrize("bad", [0, -1, 4000, True, False, 1.5, "1"])
def test_to_roman_rejects_invalid(bad):
    with pytest.raises(ValueError):
        _to_roman(bad)


def test_format_fixed_width_table_matches_two_line_ascii_format():
    assert _format_fixed_width_table([]) == "\n\n"
    rows = [
        {"Sequence": "SIINFEKL", "IC50": "12.3 nM", "Allele": "HLA-A*02:01"},
        {"Sequence": "GILGFVFTL", "IC50": "1234.5 nM", "Allele": "HLA-B*07:02"},
    ]
    expected = (
        " Sequence        IC50        Allele\n"
        "---------   ---------   -----------\n"
        " SIINFEKL     12.3 nM   HLA-A*02:01\n"
        "GILGFVFTL   1234.5 nM   HLA-B*07:02\n"
    )
    assert _format_fixed_width_table(rows) == expected
    assert _format_fixed_width_table([{"z": 1, "a": 2}, {"m": 3, "z": 4}]) == (
        "z   a   m\n-   -   -\n1   2    \n4       3\n"
    )
