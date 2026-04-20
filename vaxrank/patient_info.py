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

from dataclasses import dataclass, field
from typing import Any

from serializable import DataclassSerializable


@dataclass
class PatientInfo(DataclassSerializable):
    patient_id: str
    vcf_paths: list[str] = field(default_factory=list)
    bam_path: Any = None
    mhc_alleles: list[str] = field(default_factory=list)
    num_somatic_variants: int = 0
    num_coding_effect_variants: int = 0
    num_variants_with_rna_support: int = 0
    num_variants_with_vaccine_peptides: int = 0
