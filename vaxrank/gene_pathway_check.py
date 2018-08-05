# Copyright (c) 2018. Mount Sinai School of Medicine
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

from __future__ import absolute_import, print_function, division
from collections import OrderedDict

import pandas as pd


_ENSEMBL_GENE_ID = 'Ensembl Gene ID'
_MUTATION = 'Mutation'

_IFG_RESPONSE = 'interferon_gamma_response'
_CLASS_I_MHC = 'class1_mhc_presentation_pathway'
_DRIVER_GENE = 'cancer_driver_gene'
_DRIVER_VARIANT = 'cancer_driver_variant'


class GenePathwayCheck:
    """
    This class is meant for use with gene/variant list files from 
    https://github.com/openvax/gene-lists. Other files can be used as well, but need to follow
    a similar column structure. Most logic is based on Ensembl gene IDs.

    Parameters
    ----------
    interferon_gamma_response_csv : str, optional
        Local path to interferon-gamma response CSV file.

    class1_mhc_presentation_pathway_csv : str, optional
        Local path to MHC class I presentation pathway CSV file.

    cancer_driver_genes_csv : str, optional
        Local path to cancer driver genes CSV file.

    cancer_driver_variants_csv : str, optional
        Local path to cancer driver variants CSV file.
    """
    def __init__(
            self,
            interferon_gamma_response_csv=None,
            class1_mhc_presentation_pathway_csv=None,
            cancer_driver_genes_csv=None,
            cancer_driver_variants_csv=None):
        self.interferon_gamma_response = pd.read_csv(interferon_gamma_response_csv) \
            if interferon_gamma_response_csv else pd.DataFrame()
        self.class1_mhc_presentation_pathway = pd.read_csv(class1_mhc_presentation_pathway_csv) \
            if class1_mhc_presentation_pathway_csv else pd.DataFrame()
        self.cancer_driver_genes = pd.read_csv(cancer_driver_genes_csv) \
            if cancer_driver_genes_csv else pd.DataFrame()
        self.cancer_driver_variants = pd.read_csv(cancer_driver_variants_csv) \
            if cancer_driver_variants_csv else pd.DataFrame()
        self.check_expected_columns()

    def check_expected_columns(self):
        if len(self.interferon_gamma_response) > 0:
            if not _ENSEMBL_GENE_ID in self.interferon_gamma_response.columns:
                raise ValueError(
                    "Interferon gamma response file needs column: %s" % _ENSEMBL_GENE_ID)
        if len(self.class1_mhc_presentation_pathway) > 0:
            if not _ENSEMBL_GENE_ID in self.class1_mhc_presentation_pathway.columns:
                raise ValueError(
                    "Class I MHC presentation pathway file needs column: %s" % _ENSEMBL_GENE_ID)
        if len(self.cancer_driver_genes) > 0:
            if not _ENSEMBL_GENE_ID in self.cancer_driver_genes.columns:
                raise ValueError(
                    "Cancer driver gene file needs column: %s" % _ENSEMBL_GENE_ID)
        if len(self.cancer_driver_variants) > 0:
            if not _ENSEMBL_GENE_ID in self.cancer_driver_variants.columns or \
                    not _MUTATION in self.cancer_driver_variants.columns:
                raise ValueError(
                    "Cancer driver variant file needs columns: %s, %s" % (
                        _ENSEMBL_GENE_ID, _MUTATION))

    @staticmethod
    def is_present(df, value, column):
        return len(df.loc[df[column] == value]) > 0

    def make_variant_dict(self, variant):
        """
        Returns a dictionary of boolean values, depending on whether we see this variant in any
        relevant pathway or cancer driver files.

        Parameters
        ----------
        variant : varcode.Variant
            Variant object to evaluate
        """
        variant_dict = OrderedDict([
            (_IFG_RESPONSE, False),
            (_CLASS_I_MHC, False),
            (_DRIVER_GENE, False),
            (_DRIVER_VARIANT, False),
        ])
        effect = variant.effects().top_priority_effect().short_description
        gene_ids = variant.gene_ids

        for gene_id in gene_ids:
            if len(self.interferon_gamma_response) > 0 and \
                    is_present(self.interferon_gamma_response, gene_id, _ENSEMBL_GENE_ID):
                variant_dict[_IFG_RESPONSE] = True
            if len(self.class1_mhc_presentation_pathway) > 0 and \
                    is_present(self.class1_mhc_presentation_pathway, gene_id, _ENSEMBL_GENE_ID):
                variant_dict[_CLASS_I_MHC] = True
            if len(self.cancer_driver_genes) > 0 and \
                    is_present(self.cancer_driver_genes, gene_id, _ENSEMBL_GENE_ID):
                variant_dict[_DRIVER_GENE] = True
            if len(self.cancer_driver_variants) > 0 and \
                    len(self.cancer_driver_variants.loc[
                        (self.cancer_driver_variants[_ENSEMBL_GENE_ID] == gene_id) & 
                        (self.cancer_driver_variants[_MUTATION] == effect)]) > 0:
                variant_dict[_DRIVER_VARIANT] = True

        return variant_dict