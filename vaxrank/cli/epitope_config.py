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

import argparse

DEFAULT_MIN_EPITOPE_SCORE = 0.00001
DEFAULT_BINDING_AFFINITY_CUTOFF = 5000.0

class EpitopeConfig(msgspec.Struct):

    """Parameters for score, filtering, and ranking both epitopes and vaccine peptides"""
    logistic_epitope_score_midpoint : float = 350.0
    logistic_epitope_score_width : float = 150.0
    
    min_epitope_score : float = DEFAULT_MIN_EPITOPE_SCORE
    binding_affinity_cutoff : float = DEFAULT_BINDING_AFFINITY_CUTOFF


def add_epitope_prediction_args(arg_parser : argparse.ArgumentParser):
    epitope_prediction_args = arg_parser.add_argument_group("T-cell epitope prediction options")
    epitope_prediction_args.add_argument(
        "--min-epitope-score",
        type=float,
        help=(
            f"Ignore predicted MHC ligands whose normalized binding score "
            "falls below this threshold. (default: {DEFAULT_MIN_EPITOPE_SCORE})"))
    

def epitope_config_from_args(args : argparse.Namespace) -> EpitopeConfig:
    """
    Extract config path and overrides from argument namespace
    """
    epitope_config_kwargs = {}
    if args.epitope_prediction_config:
        with open(args.epitope_prediction_config) as f:
            epitope_config_kwargs.update(msgspec.yaml.decode(f.read()))

    if args.min_epitope_score is not None:
        epitope_config_kwargs["min_epitope_score"] = args.min_epitope_score
    epitope_config = msgspec.convert(epitope_config_kwargs, EpitopeConfig)
    return epitope_config