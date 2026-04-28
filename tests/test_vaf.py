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

"""Tests for DNA VAF extraction from varcode VCF metadata."""

from types import SimpleNamespace

from vaxrank.vaf import extract_dna_vaf_by_variant


def _stub_collection(metadata_by_source):
    return SimpleNamespace(source_to_metadata_dict=metadata_by_source)


def _entry(sample_info, alt_allele_index=0):
    return {
        'id': '.',
        'qual': None,
        'filter': 'PASS',
        'info': {},
        'sample_info': sample_info,
        'alt_allele_index': alt_allele_index,
    }


def test_dna_vaf_prefers_af_field():
    v = object()
    coll = _stub_collection({'a.vcf': {v: _entry({'TUMOR': {'AF': 0.42, 'AD': [10, 5]}})}})
    out = extract_dna_vaf_by_variant(coll, tumor_sample_name='TUMOR')
    assert out[v] == 0.42


def test_dna_vaf_falls_back_to_ad():
    v = object()
    coll = _stub_collection({'a.vcf': {v: _entry({'TUMOR': {'AD': [70, 30]}})}})
    out = extract_dna_vaf_by_variant(coll, tumor_sample_name='TUMOR')
    assert abs(out[v] - 0.30) < 1e-9


def test_dna_vaf_uses_alt_allele_index():
    v = object()
    sample = {'TUMOR': {'AD': [60, 20, 20]}}
    coll = _stub_collection({'a.vcf': {v: _entry(sample, alt_allele_index=1)}})
    out = extract_dna_vaf_by_variant(coll, tumor_sample_name='TUMOR')
    assert abs(out[v] - 0.20) < 1e-9


def test_dna_vaf_auto_picks_single_sample():
    v = object()
    coll = _stub_collection({'a.vcf': {v: _entry({'ONLY': {'AF': 0.1}})}})
    out = extract_dna_vaf_by_variant(coll)
    assert out[v] == 0.1


def test_dna_vaf_skips_multi_sample_without_name():
    v = object()
    sample = {'TUMOR': {'AF': 0.4}, 'NORMAL': {'AF': 0.0}}
    coll = _stub_collection({'a.vcf': {v: _entry(sample)}})
    out = extract_dna_vaf_by_variant(coll)
    assert v not in out


def test_dna_vaf_skips_when_neither_af_nor_ad():
    v = object()
    coll = _stub_collection({'a.vcf': {v: _entry({'TUMOR': {'GT': '0/1'}})}})
    out = extract_dna_vaf_by_variant(coll, tumor_sample_name='TUMOR')
    assert v not in out


def test_dna_vaf_handles_af_as_list():
    v = object()
    coll = _stub_collection({'a.vcf': {v: _entry({'TUMOR': {'AF': [0.25]}})}})
    out = extract_dna_vaf_by_variant(coll, tumor_sample_name='TUMOR')
    assert out[v] == 0.25


def test_dna_vaf_empty_collection():
    coll = SimpleNamespace(source_to_metadata_dict={})
    assert extract_dna_vaf_by_variant(coll) == {}


def test_dna_vaf_zero_depth_returns_none():
    v = object()
    coll = _stub_collection({'a.vcf': {v: _entry({'TUMOR': {'AD': [0, 0]}})}})
    out = extract_dna_vaf_by_variant(coll, tumor_sample_name='TUMOR')
    assert v not in out
