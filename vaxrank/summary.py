# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

def variant_counts(self, isovar_results, ranked_vaccine_peptides):
    """
    Gather statistics about the number of variants at different filtering
    steps.

    Parameters
    ----------
    isovar_results : list of isovar.IsovarResult

    ranked_vaccine_peptides : OrderedDict
        Ordered dictionary mapping variants to list of VaccinePeptide

    Returns
    -------
    Dictionary from keys such as 'num_total_variants' to int
    """

    # dictionary which will contain some overall variant counts for report
    # display
    counts_dict = {
        'num_total_variants': len(isovar_results),
        'num_coding_effect_variants': 0,
        'num_variants_with_rna_support': 0,
        'num_variants_with_vaccine_peptides': 0,
    }
    for isovar_result in isovar_results:
        variant = isovar_result.variant
        effect = isovar_result.predicted_effect
        if effect.modifies_protein_sequence:
            counts_dict['num_coding_effect_variants'] += 1
        if isovar_result.num_alt_fragments > 0:
            counts_dict['num_variants_with_rna_support'] += 1
        if variant in ranked_vaccine_peptides:
            counts_dict['num_variants_with_vaccine_peptides'] += 1
    return counts_dict
