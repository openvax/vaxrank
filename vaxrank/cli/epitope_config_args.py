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

import argparse

import msgspec

from ..epitope_config import EpitopeConfig


def add_epitope_prediction_args(arg_parser : argparse.ArgumentParser):
    epitope_prediction_args = arg_parser.add_argument_group("T-cell epitope prediction options")
    epitope_prediction_args.add_argument(
        "--min-epitope-score",
        type=float,
        help=(
            "Ignore predicted MHC ligands whose normalized binding score "
            "falls below this threshold."))
    

def epitope_config_from_args(args : argparse.Namespace) -> EpitopeConfig:
    """
    Extract config path and overrides from argument namespace
    """
    epitope_config_kwargs = {}
    if args.config:
        with open(args.config) as f:
            content = f.read()
            if content.strip():  # Only decode if file is not empty
                yaml_config = msgspec.yaml.decode(content, type=dict)
                # Extract epitope-related config if present
                if yaml_config and "epitope_config" in yaml_config:
                    epitope_config_kwargs.update(yaml_config["epitope_config"])

    if args.min_epitope_score is not None:
        epitope_config_kwargs["min_epitope_score"] = args.min_epitope_score
    epitope_config = msgspec.convert(epitope_config_kwargs, EpitopeConfig)
    return epitope_config