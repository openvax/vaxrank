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


from .core_logic import run_vaxrank
from .epitope_config import EpitopeConfig
from .epitope_logic import predict_epitopes
from .candidate_epitope import CandidateEpitope, Peptide
from .vaccine_config import VaccineConfig
from .vaccine_antigen import (
    AminoAcidInterval,
    SelfReferenceMatch,
    SelfReferenceSource,
    TargetableMask,
    TumorSpecificityAttestation,
    VaccineAntigen,
)
from .vaccine_peptide import VaccinePeptide
from .safety_assessment import (
    ConstructBoundary,
    EmittedSafetyLigand,
    SafetyAssessmentError,
    SafetyPrediction,
    SafetyPredictionCoverage,
    WindowSafetyAssessment,
    assess_vaccine_antigen_window,
    safety_assessment_from_prediction_frame,
)
from .version import __version__

__all__ = [
    "__version__",
    "AminoAcidInterval",
    "CandidateEpitope",
    "ConstructBoundary",
    "EmittedSafetyLigand",
    "EpitopeConfig",
    "Peptide",
    "SafetyAssessmentError",
    "SafetyPrediction",
    "SafetyPredictionCoverage",
    "SelfReferenceMatch",
    "SelfReferenceSource",
    "TargetableMask",
    "TumorSpecificityAttestation",
    "VaccineAntigen",
    "VaccineConfig",
    "VaccinePeptide",
    "WindowSafetyAssessment",
    "assess_vaccine_antigen_window",
    "predict_epitopes",
    "run_vaxrank",
    "safety_assessment_from_prediction_frame",
]
