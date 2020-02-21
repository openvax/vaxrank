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

from collections import OrderedDict

from serializable import Serializable

class VaxrankResults(Serializable):
    """
    Data class used to represent all results captured by running  Vaxrank.
    """
    def __init__(
            self,
            isovar_results,
            variant_to_vaccine_peptides_dict,
            ranked_vaccine_peptides):
        """
        Parameters
        ----------
        isovar_results : list of isovar.IsovarResult
            IsovarResult object for each variant without any filtering

        variant_to_vaccine_peptides_dict : dict
            Dictionary mapping variant to a list of possible vaccine peptides

        ranked_vaccine_peptides : list of VaccinePeptide
        """
        self.isovar_results = isovar_results
        self.variant_to_vaccine_peptides_dict = variant_to_vaccine_peptides_dict
        self.ranked_vaccine_peptides = ranked_vaccine_peptides


    @property
    def variants(self):
        """
        Unfiltered list of variants

        Returns
        -------
        list of varcode.Variant
        """
        return [isovar_result.variant for isovar_result in self.isovar_results]


    @property
    def variant_to_protein_sequences_dict(self):
        # TODO: find out if we can safely get rid of this property
        return {
            isovar_result.variant: isovar_result.sorted_protein_sequences[0]
            for isovar_result in self.isovar_results
            if isovar_result.passes_all_filters
               and len(isovar_result.sorted_protein_sequences) > 0
        }

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
            'num_total_variants': len(self.isovar_results),
            'num_coding_effect_variants': 0,
            'num_variants_with_rna_support': 0,
            'num_variants_with_vaccine_peptides': 0,
        }
        for isovar_result in self.isovar_results:
            variant = isovar_result.variant
            if isovar_result.predicted_effect_modifies_protein_sequence:
                counts_dict['num_coding_effect_variants'] += 1
            if isovar_result.has_mutant_protein_sequence_from_rna:
                counts_dict['num_variants_with_rna_support'] += 1
            if variant in self.variant_to_vaccine_peptides_dict:
                counts_dict['num_variants_with_vaccine_peptides'] += 1
        return counts_dict


    def variant_properties(self, gene_pathway_check=None):
        """
        Parameters
        ----------
        gene_pathway_check : GenePathwayCheck (optional)
            Used to look up whether a mutation or its affected gene are in some
            biologically important pathway.

        Returns
        -------
        list of dictionaries containing properties we want to analyze later,
        e.g. whether this variant is part of a pathway of interest,
        is a strong MHC binder, etc.
        """
        variant_properties_list = []
        for variant in self.variants:
            gene_name = ''
            if variant.gene_names:
                gene_name = variant.effects().top_priority_effect().gene_name
            variant_dict = OrderedDict((
                ('contig', variant.contig),
                ('start', variant.start),
                ('ref', variant.ref),
                ('alt', variant.alt),
                ('is_coding_nonsynonymous', False),
                ('rna_support', False),
                ('mhc_binder', False),
                ('gene_name', gene_name),
            ))
            if gene_pathway_check is not None:
                pathway_dict = gene_pathway_check.make_variant_dict(variant)
                variant_dict.update(pathway_dict)
            if len(variant.effects().drop_silent_and_noncoding()) > 0:
                variant_dict['is_coding_nonsynonymous'] = True
            if variant in self.variant_to_protein_sequences_dict:
                variant_dict['rna_support'] = True
            # TODO: compute MHC binder status for variants that don't have RNA support
            if variant in self.variant_to_vaccine_peptides_dict:
                variant_dict['mhc_binder'] = True
            variant_properties_list.append(variant_dict)
        return variant_properties_list
