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

def variant_counts(
        variants,
        protein_sequence_dict,
        vaccine_peptides_dict):
    """
    Summarize Vaxrank counts for total variants, variants with coding effects,
    variants with RNA support, and variants with associated vaccine peptides.

    Parameters
    ----------
    variants : list of varcode.Variant

    protein_sequence_dict : dict
        Dictionary from variants to isovar protein sequences

    vaccine_peptides_dict : dict
        Dictionary from variants to VaccinePeptide objects

    Returns
    -------
    dict

    """
    # dictionary which will contain some overall variant counts for report display
    counts_dict = {
        'num_total_variants': len(variants),
        'num_coding_effect_variants': 0,
        'num_variants_with_rna_support': 0,
        'num_variants_with_vaccine_peptides': 0,
    }
    for variant in variants:
        if len(variant.effects().drop_silent_and_noncoding()) > 0:
            counts_dict['num_coding_effect_variants'] += 1
        if variant in protein_sequence_dict:
            counts_dict['num_variants_with_rna_support'] += 1
        if variant in vaccine_peptides_dict:
            counts_dict['num_variants_with_vaccine_peptides'] += 1
    return counts_dict


def variant_properties(
        variants,
        protein_sequence_dict,
        vaccine_peptide_dict,
        gene_pathway_check=None):
    """

    Parameters
    ----------
    variants : list of varcode.Variant

    protein_sequence_dict : dict
        Dictionary from variants to isovar protein sequences

    vaccine_peptides_dict : dict
        Dictionary from variants to VaccinePeptide objects

    gene_pathway_check : GenePathwayCheck (optional)
        Used to look up whether a mutation or its affected gene are in some
        biologically important pathway.

    Returns
    -------
    list of dictionaries containing properties we want to analyze later, e.g.
    whether this variant is part of a pathway of interest,
    is a strong MHC binder, etc.
    """
    variant_properties_list = []
    for variant in variants:
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
        if variant in protein_sequence_dict:
            variant_dict['rna_support'] = True
        # TODO: compute MHC binder status for variants that don't have RNA support
        if variant in vaccine_peptide_dict:
            variant_dict['mhc_binder'] = True
        variant_properties_list.append(variant_dict)
    return variant_properties_list
