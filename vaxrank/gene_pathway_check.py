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
from os.path import join, dirname

import pandas as pd


_ENSEMBL_GENE_ID_COLUMN_NAME = 'Ensembl Gene ID'
_MUTATION_COLUMN_NAME = 'Mutation'

_IFNG_RESPONSE_COLUMN_NAME = 'interferon_gamma_response'
_CLASS_I_MHC_COLUMN_NAME = 'class1_mhc_presentation_pathway'
_DRIVER_GENE_COLUMN_NAME = 'cancer_driver_gene'
_DRIVER_VARIANT_COLUMN_NAME = 'cancer_driver_variant'

_CURRENT_DIR = dirname(__file__)
_DATA_DIR = join(_CURRENT_DIR, "data")


class GenePathwayCheck(object):
    """
    This class is meant for use with gene/variant list files from
    https://github.com/openvax/gene-lists. Other files can be used as well, but
    need to follow a similar column structure. Most logic is based on Ensembl
    gene IDs.

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

        self.interferon_gamma_response_gene_set = self._load_set_from_csv(
            csv_path=interferon_gamma_response_csv,
            default_filename="interferon-gamma-response.csv",
            description="Interferon gamma response pathway",
            column_names=[_ENSEMBL_GENE_ID_COLUMN_NAME])

        self.class1_mhc_presentation_pathway_gene_set = self._load_set_from_csv(
            csv_path=class1_mhc_presentation_pathway_csv,
            default_filename="class1-mhc-presentation-pathway.csv",
            description="Class I MHC presentation pathway",
            column_names=[_ENSEMBL_GENE_ID_COLUMN_NAME])

        self.cancer_driver_genes_set = self._load_set_from_csv(
            csv_path=cancer_driver_genes_csv,
            default_filename="cancer-driver-genes.csv",
            description="Cancer driver genes",
            column_names=[_ENSEMBL_GENE_ID_COLUMN_NAME])
        # set of gene ID, variant description pairs
        self.cancer_driver_variants_set = self._load_set_from_csv(
            csv_path=cancer_driver_variants_csv,
            default_filename="cancer-driver-variants.csv",
            description="Driver variants",
            column_names=[_ENSEMBL_GENE_ID_COLUMN_NAME, _MUTATION_COLUMN_NAME])

    @classmethod
    def _load_set_from_csv(cls, csv_path, default_filename, description, column_names):
        if not csv_path:
            csv_path = join(_DATA_DIR, default_filename)
        df = pd.read_csv(csv_path)
        columns = []
        for column_name in column_names:
            if column_name not in df.columns:
                raise ValueError("%s file (%s) needs column '%s'" % (
                    description,
                    csv_path,
                    column_name))
            columns.append(df[column_name].values)
        if len(columns) == 1:
            return set(columns[0])
        else:
            return set(zip(*columns))

    def make_variant_dict(self, variant):
        """
        Returns a dictionary of boolean values, depending on whether we see this
        variant in any relevant pathway or cancer driver files.

        Parameters
        ----------
        variant : varcode.Variant
            Variant object to evaluate
        """
        effect_description = variant.effects().top_priority_effect().short_description
        overlapping_gene_ids = variant.gene_ids
        variant_dict = OrderedDict()
        variant_dict[_IFNG_RESPONSE_COLUMN_NAME] = any([
            gene_id in self.interferon_gamma_response_gene_set
            for gene_id in overlapping_gene_ids
        ])
        variant_dict[_CLASS_I_MHC_COLUMN_NAME] = any([
            gene_id in self.class1_mhc_presentation_pathway_gene_set
            for gene_id in overlapping_gene_ids
        ])
        variant_dict[_DRIVER_GENE_COLUMN_NAME] = any([
            gene_id in self.cancer_driver_genes_set
            for gene_id in overlapping_gene_ids
        ])

        variant_dict[_DRIVER_VARIANT_COLUMN_NAME] = any([
            (gene_id, effect_description) in self.cancer_driver_variants_set
            for gene_id in overlapping_gene_ids
        ])
        return variant_dict
