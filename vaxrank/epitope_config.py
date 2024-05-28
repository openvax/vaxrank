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

DEFAULT_MIN_EPITOPE_SCORE = 0.00001
DEFAULT_BINDING_AFFINITY_CUTOFF = 5000.0

class EpitopeConfig(msgspec.Struct):

    """Parameters for score, filtering, and ranking both epitopes and vaccine peptides"""
    logistic_epitope_score_midpoint : float = 350.0
    logistic_epitope_score_width : float = 150.0
    
    min_epitope_score : float = DEFAULT_MIN_EPITOPE_SCORE
    binding_affinity_cutoff : float = DEFAULT_BINDING_AFFINITY_CUTOFF
