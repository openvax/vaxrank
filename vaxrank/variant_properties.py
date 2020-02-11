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

import pandas as pd

from .gene_pathway_check import GenePathwayCheck

default_gene_pathways = GenePathwayCheck()

def make_variant_properties_dataframe(
        self,
        isovar_results,
        mhc_predictions,
        gene_pathways=default_gene_pathways):
    """
    Parameters
    ----------
    isovar_results : list of isovar.IsovarResult

    gene_pathways : GenePathwayCheck

    Returns
    -------
    List of  dictionary of properties we want to
    analyze later, e.g. whether this variant is part of a pathway of interest,
    is a strong MHC binder, etc.
    """
    variant_properties_dicts = []

    for isovar_result in isovar_results:
        variant = isovar_result.variant
        effect = isovar_result.predicted_effect

        variant_dict = OrderedDict([
            ('contig', variant.contig),
            ('start', variant.start),
            ('ref', variant.ref),
            ('alt', variant.alt),
            ('gene_name', effect.gene_name),
            ("predicted_effect", effect.short_description),
            ("predicted_effect_class", effect.__class__.__name__),
            ('is_coding_nonsynonymous', effect.modifies_protein_sequence),
            ("rna_num_alt_fragments", isovar_result.num_alt_fragments),
            ("rna_num_alt_reads", isovar_result.num_alt_reads),
            ("rna_num_total_fragments", isovar_result.num_total_fragments),
            ("rna_vaf", isovar_result.fraction_alt_fragments),
            ('rna_support', isovar_result.passes_all_filters),
        ])
        if gene_pathways is not None:
            variant_dict.update(gene_pathways.make_variant_dict(variant))

        if variant in self.vaccine_peptides:
            variant_dict['mhc_binder'] = True

        variant_dict[variant] = variant
        variant_properties_dicts.append(variant_dict)

    return pd.DataFrame.from_records(variant_properties_dicts)
