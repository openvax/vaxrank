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

from ..vaccine_config import VaccineConfig



def add_vaccine_peptide_args(arg_parser : argparse.ArgumentParser) -> None:
    vaccine_peptide_group = arg_parser.add_argument_group("Vaccine peptide options")
    vaccine_peptide_group.add_argument(
        "--vaccine-peptide-length",
        default=25,
        type=int,
        help="Number of amino acids in the vaccine peptides. (default: %(default)s)")

    vaccine_peptide_group.add_argument(
        "--padding-around-mutation",
        default=5,
        type=int,
        help=(
            "Number of off-center windows around the mutation to consider "
            "as vaccine peptides. (default: %(default)s)"
        ))

    vaccine_peptide_group.add_argument(
        "--max-vaccine-peptides-per-mutation",
        default=1,
        type=int,
        help=(
            "Number of vaccine peptides to generate for each mutation. "
            "(default: %(default)s)"
        ))

    vaccine_peptide_group.add_argument(
        "--num-epitopes-per-vaccine-peptide",
        type=int,
        help=(
            "Maximum number of mutant epitopes to consider when scoring "
            "each vaccine peptide. (default: %(default)s)"))



    

def vaccine_config_from_args(args : argparse.Namespace) -> VaccineConfig:
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