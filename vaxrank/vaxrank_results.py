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
        return [
            isovar_result.variant
            for isovar_result
            in self.isovar_results
        ]

    def variant_counts(self):
        """
        Summarize Vaxrank counts for total variants, variants with coding effects,
        variants with RNA support, and variants with associated vaccine peptides.

        Returns
        -------
        dict
        """
        variant_properties = self.variant_properties()

        # dictionary which will contain some overall variant counts
        # for report display
        counts_dict = {}
        counts_dict['num_total_variants'] = len(self.isovar_results)
        counts_dict['num_coding_effect_variants'] = \
            sum([v['is_coding_nonsynonymous'] for v in variant_properties])
        counts_dict['num_variants_with_rna_support'] = \
            sum([v['rna_support'] for v in variant_properties])

        counts_dict['num_variants_with_vaccine_peptides'] =  \
            sum([v['has_vaccine_peptide'] for v in variant_properties])
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
        for isovar_result in self.isovar_results:
            variant = isovar_result.variant

            variant_dict = OrderedDict((
                ('gene_name', isovar_result.top_gene_name),
                ('contig', variant.contig),
                ('start', variant.start),
                ('ref', variant.ref),
                ('alt', variant.alt),
                ('is_coding_nonsynonymous',
                    isovar_result.predicted_effect_modifies_protein_sequence),
                ('rna_support',
                    isovar_result.has_mutant_protein_sequence_from_rna),
            ))

            # TODO:
            #  compute MHC binder status for variants that don't have RNA support
            variant_dict['mhc_binder'] = \
                variant_dict["has_vaccine_peptide"] = \
                    variant in self.variant_to_vaccine_peptides_dict

            if gene_pathway_check is not None:
                pathway_dict = gene_pathway_check.make_variant_dict(variant)
                variant_dict.update(pathway_dict)

            variant_properties_list.append(variant_dict)
        return variant_properties_list
