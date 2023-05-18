# Copyright (c) 2016-2018. Mount Sinai School of Medicine
#
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

from __future__ import absolute_import, print_function, division
import os
import logging

import six
from datacache import get_data_dir
import shellinford


logger = logging.getLogger(__name__)


def fm_index_path(genome):
    """
    Returns a path for cached reference peptides, for the given genome.
    """
    # if $VAXRANK_REF_PEPTIDES_DIR is set, that'll be the location of the cache
    cache_dir = get_data_dir(envkey='VAXRANK_REF_PEPTIDES_DIR')
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)

    return os.path.join(cache_dir, '%s_%d_%d.fm' % (
        genome.species.latin_name, genome.release, 2 if six.PY2 else 3))


def generate_protein_sequences(genome):
    """
    Generator whose elements are protein sequences from the given genome.

    Parameters
    ----------
    genome : pyensembl.EnsemblRelease
        Input genome to load for reference peptides
    """
    for t in genome.transcripts():
        if t.is_protein_coding:
            protein_sequence = t.protein_sequence
            if six.PY2:
                # shellinford on PY2 seems to sometimes fail with
                # unicode strings
                protein_sequence = protein_sequence.encode("ascii")
            yield protein_sequence


def load_reference_peptides_index(genome, force_reload=False):
    """
    Loads the FM index containing reference peptides.

    Parameters
    ----------
    genome : pyensembl.EnsemblRelease
        Input genome to load for reference peptides

    force_reload : bool, optional
        If true, will recompute index for this genome even if it already exists.

    Returns
    -------
    fm : shellinford.FMIndex
        Index populated with reference peptides from the genome
    """
    path = fm_index_path(genome)
    if force_reload or not os.path.exists(path):
        logger.info("Building FM index at %s", path)
        fm = shellinford.FMIndex()
        fm.build(generate_protein_sequences(genome), path)
        logger.info("Done building FM index")
        return fm
    return shellinford.FMIndex(filename=path)


class ReferenceProteome(object):
    def __init__(self, genome):
        self.fm_index = load_reference_peptides_index(genome)

    def contains(self, kmer):
        return len(self.fm_index.search(kmer)) > 0
