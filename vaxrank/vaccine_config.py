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

import msgspec

class VaccineConfig(msgspec.Struct):
    """Parameters for assembling epitope predictions into vaccine peptides"""
    vaccine_peptide_length : int = 25
    
    padding_around_mutation : int = 5

    max_vaccine_peptides_per_variant : int = 1
    
    num_mutant_epitopes_to_keep : int = 1000 