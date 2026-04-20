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

"""Smoke tests for PatientInfo's dataclass-based serialization path —
catches regressions if the DataclassSerializable mixin or the upstream
`serializable` wire format ever shifts."""

import pickle

from vaxrank.patient_info import PatientInfo


def _sample_patient_info():
    return PatientInfo(
        patient_id="P1",
        vcf_paths=["a.vcf", "b.vcf"],
        bam_path="sample.bam",
        mhc_alleles=["HLA-A*02:01", "HLA-B*07:02"],
        num_somatic_variants=42,
        num_coding_effect_variants=17,
        num_variants_with_rna_support=9,
        num_variants_with_vaccine_peptides=3,
    )


def test_patient_info_json_roundtrip():
    p = _sample_patient_info()
    restored = PatientInfo.from_json(p.to_json())
    assert restored == p


def test_patient_info_pickle_roundtrip():
    # Dataclass gives us __eq__; DataclassSerializable gives us __reduce__
    # so pickling routes through the to_dict / from_dict envelope.
    p = _sample_patient_info()
    assert pickle.loads(pickle.dumps(p)) == p


def test_patient_info_to_dict_matches_init_kwargs():
    p = _sample_patient_info()
    d = p.to_dict()
    # Round-trip through the plain dict form.
    assert PatientInfo.from_dict(d) == p
