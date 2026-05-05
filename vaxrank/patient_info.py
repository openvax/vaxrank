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
from typing import Optional

from serializable import DataclassSerializable


@dataclass
class PatientInfo(DataclassSerializable):
    patient_id: str
    vcf_paths: list[str] = field(default_factory=list)
    bam_path: Optional[str] = None
    mhc_alleles: list[str] = field(default_factory=list)
    num_somatic_variants: int = 0
    num_coding_effect_variants: int = 0
    num_variants_with_rna_support: int = 0
    num_variants_with_vaccine_peptides: int = 0
    # Generic inputs description: ordered ``[(label, path), ...]``
    # rendered verbatim in the report's Inputs block. Lets the
    # external (LENS / pVACseq) path label its file as "LENS report"
    # / "pVACseq report" rather than misleadingly as "VCF (somatic
    # variants)" via ``vcf_paths``. When empty the renderer falls
    # back to the legacy ``vcf_paths`` / ``bam_path`` fields for
    # backward compatibility with previously-saved JSON.
    inputs: list[tuple[str, str]] = field(default_factory=list)
