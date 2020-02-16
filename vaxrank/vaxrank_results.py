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

from __future__ import absolute_import, print_function, division

from serializable import Serializable

class VaxrankResults(Serializable):
    """
    Data class used to represent all results captured by running  Vaxrank.
    """

    def __init__(
            self,
            variants,
            variant_to_protein_sequences_dict,
            variant_to_vaccine_peptides_dict,
            ranked_vaccine_peptides):
        """
        Parameters
        ----------
        variants : varcode.VariantCollection
            All variants without any filtering

        variant_to_protein_sequences_dict : dict
            Dictionary mapping each variant to a list of candidate
            protein sequences

        variant_to_vaccine_peptides_dict : dict
            Dictionary mapping variant to a list of possible vaccine peptides

        ranked_vaccine_peptides : list of VaccinePeptide
        """
        self.variants = variants
        self.variant_to_protein_sequences_dict = variant_to_protein_sequences_dict
        self.variant_to_vaccine_peptides_dict = variant_to_vaccine_peptides_dict
        self.ranked_vaccine_peptides = ranked_vaccine_peptides

    def variant_counts(self):
        """
        Summarize Vaxrank counts for total variants, variants with coding effects,
        variants with RNA support, and variants with associated vaccine peptides.


        Returns
        -------
        dict

        """
        # dictionary which will contain some overall variant counts for report display
        counts_dict = {
            'num_total_variants': len(self.variants),
            'num_coding_effect_variants': 0,
            'num_variants_with_rna_support': 0,
            'num_variants_with_vaccine_peptides': 0,
        }
        for variant in self.variants:
            if len(variant.effects().drop_silent_and_noncoding()) > 0:
                counts_dict['num_coding_effect_variants'] += 1
            if variant in self.variant_to_protein_sequences_dict:
                counts_dict['num_variants_with_rna_support'] += 1
            if variant in self.variant_to_vaccine_peptides_dict:
                counts_dict['num_variants_with_vaccine_peptides'] += 1
        return counts_dict
