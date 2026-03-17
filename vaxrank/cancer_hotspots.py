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
Cancer hotspot mutation lookup using bundled data from cancerhotspots.org.

Data source: Chang et al. 2016/2017 cancer hotspots, via maftools.
"""

import logging
import os
from functools import lru_cache

import pkg_resources

logger = logging.getLogger(__name__)


def _get_hotspots_path():
    """Get path to bundled cancer hotspots TSV file."""
    return pkg_resources.resource_filename(
        'vaxrank', 'data/cancer_hotspots_v2.tsv'
    )


@lru_cache(maxsize=1)
def _load_hotspots():
    """
    Load cancer hotspots data into a lookup dictionary.

    Returns dict mapping (gene_upper, protein_change) -> hotspot info dict
    """
    hotspots = {}
    path = _get_hotspots_path()

    if not os.path.exists(path):
        logger.warning("Cancer hotspots file not found at %s", path)
        return hotspots

    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) < 7:
                continue

            chrom = parts[0]
            pos = parts[1]
            ref = parts[2]
            alt = parts[3]
            gene = parts[4]
            mutation_type = parts[5]
            protein_change = parts[6]

            # Key by uppercase gene and protein change
            key = (gene.upper(), protein_change.upper())
            hotspots[key] = {
                'gene': gene,
                'protein_change': protein_change,
                'mutation_type': mutation_type,
                'chromosome': chrom,
                'position': pos,
                'ref': ref,
                'alt': alt,
            }

    logger.info("Loaded %d cancer hotspots", len(hotspots))
    return hotspots


def is_hotspot(gene_name, protein_change):
    """
    Check if a mutation is a known cancer hotspot.

    Parameters
    ----------
    gene_name : str
        Gene symbol (e.g., "BRAF", "TP53")
    protein_change : str
        Protein change notation, with or without "p." prefix
        (e.g., "V600E", "p.V600E")

    Returns
    -------
    bool
        True if the mutation is a known cancer hotspot
    """
    hotspots = _load_hotspots()

    # Normalize protein change - strip "p." prefix if present
    if protein_change.startswith('p.'):
        protein_change = protein_change[2:]

    key = (gene_name.upper(), protein_change.upper())
    return key in hotspots


def get_hotspot_info(gene_name, protein_change):
    """
    Get information about a cancer hotspot mutation.

    Parameters
    ----------
    gene_name : str
        Gene symbol (e.g., "BRAF", "TP53")
    protein_change : str
        Protein change notation, with or without "p." prefix
        (e.g., "V600E", "p.V600E")

    Returns
    -------
    dict or None
        Hotspot information dict if found, None otherwise
    """
    hotspots = _load_hotspots()

    # Normalize protein change - strip "p." prefix if present
    if protein_change.startswith('p.'):
        protein_change = protein_change[2:]

    key = (gene_name.upper(), protein_change.upper())
    return hotspots.get(key)


def get_hotspot_url(gene_name, protein_change):
    """
    Get a URL for information about a cancer hotspot.

    Since cancerhotspots.org doesn't have direct variant links,
    this returns a link to the gene search on the site.

    Parameters
    ----------
    gene_name : str
        Gene symbol (e.g., "BRAF", "TP53")
    protein_change : str
        Protein change notation

    Returns
    -------
    str or None
        URL to cancerhotspots.org gene page if mutation is a hotspot,
        None otherwise
    """
    if is_hotspot(gene_name, protein_change):
        # Link to the gene page on cancerhotspots.org
        return f"https://www.cancerhotspots.org/#/gene/{gene_name.upper()}"
    return None
