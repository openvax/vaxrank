# Copyright (c) 2016. Mount Sinai School of Medicine
#
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

from isovar.args import (
    add_somatic_vcf_args,
    add_rna_args,
    add_reference_context_args,
)

# inherit all commandline options from Topiary
arg_parser = argparse.ArgumentParser(
    prog="vaxrank",
    description=(
        "Select personalized vaccine peptides from cancer variants, "
        "expression data, and patient HLA type."),
)
#    parents=[],
#    add_help=False)

add_somatic_vcf_args(arg_parser)
add_rna_args(arg_parser)
add_reference_context_args(arg_parser)

arg_parser.add_argument(
    "--output-csv",
    default="vaccine_peptides.csv",
    help="Name of CSV file which contains predicted sequences")

arg_parser.add_argument(
    "--vaccine-peptide-min-length",
    default=25,
    type=int)

arg_parser.add_argument(
    "--vaccine-peptide-max-length",
    default=25,
    type=int)

