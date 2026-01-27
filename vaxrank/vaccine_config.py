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

"""
Configuration for vaccine peptide assembly.

This module defines the VaccineConfig class which controls how vaccine
peptides are assembled from epitope predictions. Vaccine peptides are
longer sequences (typically 25 amino acids) that contain one or more
predicted MHC-binding epitopes spanning the mutated residue(s).

The vaccine peptide selection considers:
1. Peptide length requirements for synthesis and delivery
2. Position of the mutation within the peptide
3. Number of strong MHC-binding epitopes contained
4. Total predicted immunogenicity score
"""

import msgspec


class VaccineConfig(msgspec.Struct, frozen=True):
    """
    Configuration parameters for vaccine peptide assembly.

    This immutable struct contains all parameters needed to assemble
    epitope predictions into vaccine peptides suitable for synthesis
    and therapeutic use.

    Attributes
    ----------
    vaccine_peptide_length : int
        Target length of vaccine peptides in amino acids. Longer peptides
        allow for more epitope presentation but may be harder to synthesize.
        Default: 25

    padding_around_mutation : int
        Number of off-center window positions to consider when selecting
        vaccine peptides. A value of 5 means considering windows where
        the mutation is anywhere from position 5 to position 20 in a
        25-mer peptide.
        Default: 5

    max_vaccine_peptides_per_variant : int
        Maximum number of vaccine peptides to select for each variant.
        Multiple peptides may be selected if they have sufficiently
        different epitope content.
        Default: 1

    num_mutant_epitopes_to_keep : int
        Maximum number of epitope predictions to retain per variant
        during processing. Higher values increase computational cost
        but may capture more potential vaccine targets.
        Default: 1000

    Examples
    --------
    >>> config = VaccineConfig()
    >>> config.vaccine_peptide_length
    25

    >>> long_peptide_config = VaccineConfig(vaccine_peptide_length=30)
    >>> long_peptide_config.vaccine_peptide_length
    30

    >>> multi_peptide_config = VaccineConfig(max_vaccine_peptides_per_variant=3)
    >>> multi_peptide_config.max_vaccine_peptides_per_variant
    3
    """

    vaccine_peptide_length: int = 25
    padding_around_mutation: int = 5
    max_vaccine_peptides_per_variant: int = 1
    num_mutant_epitopes_to_keep: int = 1000
