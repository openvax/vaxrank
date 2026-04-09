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

from ..config.loader import extract_vaccine_config_kwargs, load_merged_vaxrank_config
from ..vaccine_config import VaccineConfig



def add_vaccine_peptide_args(arg_parser : argparse.ArgumentParser) -> None:
    default_vaccine_config = VaccineConfig()

    vaccine_peptide_group = arg_parser.add_argument_group("Vaccine peptide options")
    vaccine_peptide_group.add_argument(
        "--vaccine-peptide-length",
        default=None,
        type=int,
        help=(
            "Number of amino acids in the vaccine peptides. "
            f"(default: {default_vaccine_config.vaccine_peptide_length})"
        ))

    vaccine_peptide_group.add_argument(
        "--padding-around-mutation",
        default=None,
        type=int,
        help=(
            "Number of off-center windows around the mutation to consider "
            "as vaccine peptides. "
            f"(default: {default_vaccine_config.padding_around_mutation})"
        ))

    vaccine_peptide_group.add_argument(
        "--max-vaccine-peptides-per-variant",
        "--max-vaccine-peptides-per-mutation",  # legacy alias
        dest="max_vaccine_peptides_per_variant",
        default=None,
        type=int,
        help=(
            "Number of vaccine peptides to generate for each variant. "
            f"(default: {default_vaccine_config.max_vaccine_peptides_per_variant})"
        ))

    vaccine_peptide_group.add_argument(
        "--num-epitopes-per-vaccine-peptide",
        default=None,
        type=int,
        help=(
            "Maximum number of mutant epitopes to consider when scoring "
            "each vaccine peptide. Set to 0 to keep all epitopes. "
            f"(default: {default_vaccine_config.num_mutant_epitopes_to_keep})"
        ))



    

def vaccine_config_from_args(args : argparse.Namespace, merged_config=None) -> VaccineConfig:
    """
    Extract config path and overrides from argument namespace.

    Parameters
    ----------
    args : argparse.Namespace
    merged_config : dict, optional
        Pre-loaded merged config dict.  When supplied the YAML file is
        not re-read, avoiding duplicate I/O when both epitope and
        vaccine configs are built from the same args.
    """
    if merged_config is None:
        merged_config = load_merged_vaxrank_config(args)
    vaccine_config_kwargs = extract_vaccine_config_kwargs(merged_config)

    if args.vaccine_peptide_length is not None:
        vaccine_config_kwargs["vaccine_peptide_length"] = args.vaccine_peptide_length
    if args.padding_around_mutation is not None:
        vaccine_config_kwargs["padding_around_mutation"] = args.padding_around_mutation
    max_vaccine_peptides_per_variant = getattr(
        args, "max_vaccine_peptides_per_variant", None)
    if max_vaccine_peptides_per_variant is None:
        max_vaccine_peptides_per_variant = getattr(
            args, "max_vaccine_peptides_per_mutation", None)
    if max_vaccine_peptides_per_variant is not None:
        vaccine_config_kwargs["max_vaccine_peptides_per_variant"] = (
            max_vaccine_peptides_per_variant
        )
    if args.num_epitopes_per_vaccine_peptide is not None:
        vaccine_config_kwargs["num_mutant_epitopes_to_keep"] = args.num_epitopes_per_vaccine_peptide

    vaccine_config = msgspec.convert(vaccine_config_kwargs, VaccineConfig)
    return vaccine_config
