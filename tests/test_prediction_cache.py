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

"""Tests for --prediction-cache (CachedPredictor integration)."""

import pandas as pd

from topiary import CachedPredictor


def _make_cache_df():
    """Minimal DataFrame that satisfies CachedPredictor's schema."""
    return pd.DataFrame({
        "peptide": ["AALSGGGGX", "EPGQALFNG"],
        "allele": ["HLA-A*02:01", "HLA-A*02:01"],
        "peptide_length": [9, 9],
        "kind": ["affinity", "affinity"],
        "value": [123.4, 456.7],
        "affinity": [123.4, 456.7],
        "percentile_rank": [0.5, 1.2],
        "prediction_method_name": ["netMHCpan", "netMHCpan"],
        "predictor_version": ["4.1", "4.1"],
    })


def test_cached_predictor_cache_hit():
    """Peptides in the cache should be returned without calling fallback."""
    df = _make_cache_df()
    cp = CachedPredictor.from_dataframe(df)
    result = cp.predict_peptides_dataframe(peptides=["AALSGGGGX"])
    assert len(result) == 1
    assert result.iloc[0]["affinity"] == 123.4


def test_cached_predictor_roundtrip(tmp_path):
    """Save cache to disk, reload, and verify lookups still work."""
    df = _make_cache_df()
    cache_path = tmp_path / "predictions.tsv"
    cp = CachedPredictor.from_dataframe(df)
    cp.save(str(cache_path))

    reloaded = CachedPredictor.from_topiary_output(str(cache_path))
    result = reloaded.predict_peptides_dataframe(peptides=["EPGQALFNG"])
    assert len(result) == 1
    assert result.iloc[0]["allele"] == "HLA-A*02:01"


def test_cached_predictor_from_topiary_output_tsv(tmp_path):
    """Verify that from_topiary_output correctly loads a TSV cache file
    and that it can be used as a predictor."""
    df = _make_cache_df()
    cache_path = tmp_path / "cache.tsv"
    df.to_csv(str(cache_path), sep="\t", index=False)

    cp = CachedPredictor.from_topiary_output(str(cache_path))
    assert set(cp.alleles) == {"HLA-A*02:01"}
    result = cp.predict_peptides_dataframe(peptides=["AALSGGGGX", "EPGQALFNG"])
    assert len(result) == 2
