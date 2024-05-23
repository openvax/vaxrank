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

VACCINE_PEPTIDE_LENGTH_DEFAULT = 25

class VaccineConfig(msgspec.Struct):
    """Parameters for assembling epitope predictions into vaccine peptides"""
    vaccine_peptide_length : int = VACCINE_PEPTIDE_LENGTH_DEFAULT
    
    
def vaccine_config_from_args(args):
    pass