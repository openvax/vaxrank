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

import os
import logging

from datacache import get_data_dir
import shellinford
import pandas as pd
from tqdm import tqdm
from collections import defaultdict


logger = logging.getLogger(__name__)


def get_cache_dir() -> str:
    # if $VAXRANK_REF_PEPTIDES_DIR is set, that'll be the location of the cache
    cache_dir = get_data_dir(envkey="VAXRANK_REF_PEPTIDES_DIR")
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)
    return cache_dir


def fm_index_path(genome) -> str:
    """
    Returns a path for cached reference peptides, for the given genome.
    """
    return os.path.join(
        get_cache_dir(), "%s_%d_3.fm" % (genome.species.latin_name, genome.release)
    )


def kmer_index_path(genome, min_len=8, max_len=15) -> str:
    return os.path.join(
        get_cache_dir(),
        "%s_%d_kmers_%d_%d.csv"
        % (genome.species.latin_name, genome.release, min_len, max_len),
    )


def generate_protein_sequences_dict(genome) -> dict[str, str]:
    """
    Generator whose elements are protein sequences from the given genome.

    Parameters
    ----------
    genome : pyensembl.EnsemblRelease
        Input genome to load for reference peptides
    """
    return {
        t.id: t.protein_sequence for t in genome.transcripts() if t.is_protein_coding
    }


def generate_protein_sequences(genome) -> list[str]:
    return list(generate_protein_sequences_dict(genome).values())


def build_reference_peptides_fm_index(genome, path):
    """
    Parameters
    ----------
    genome : pyensembl.EnsemblRelease
        Input genome to load for reference peptides

    """
    logger.info("Building FM index at %s", path)
    fm = shellinford.FMIndex()
    fm.build(generate_protein_sequences(genome), path)
    logger.info("Done building FM index")
    return fm


def build_reference_peptides_kmer_index(
    genome, path, min_len=8, max_len=15
) -> dict[str, list[str, int]]:
    """
    Parameters
    ----------
    genome : pyensembl.EnsemblRelease
        Input genome to load for reference peptides

    path : str
        Path to save the kmer index
    """
    logger.info("Building kmer index at %s", path)
    kmers = defaultdict(list)
    for t in tqdm(genome.transcripts()):
        if not t.is_protein_coding:
            continue
        t_id = t.id
        protein = t.protein_sequence
        for i in range(len(protein) - min_len + 1):

            for k in range(min_len, max_len + 1):
                kmer = protein[i : i + k]
                if len(kmer) == k:
                    kmers[kmer].append((t_id, i))
    cols = {"kmer": [], "transcript": [], "offset": []}
    for kmer, pairs in tqdm(kmers.items()):
        for transcript, offset in pairs:
            cols["kmer"].append(kmer)
            cols["transcript"].append(transcript)
            cols["offset"].append(offset)
    df = pd.DataFrame(cols)
    df.to_csv(path, index=False)
    logger.info("Done building kmer index")
    return kmers


def load_reference_peptides_fm_index(genome, force_reload=False) -> shellinford.FMIndex:
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
        return build_reference_peptides_fm_index(genome, path)
    return shellinford.FMIndex(filename=path)


def load_reference_peptides_kmer_index(
    genome, force_reload=False
) -> dict[str, list[str, int]]:
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
    kmer_index : dict[str, list[str, int]]
        Index mapping kmer peptides to the list of transcript IDs and amino acid positions
    """
    path = kmer_index_path(genome)
    if force_reload or not os.path.exists(path):
        return build_reference_peptides_kmer_index(genome, path)
    df = pd.read_csv(path)
    return {k: list(v) for k, v in df.set_index("kmer").to_dict(orient="index").items()}


class ReferenceProteome(object):
    def __init__(self, genome, index_type="kmer"):
        self.genome = genome
        self.index_type = index_type

        if index_type == "kmer":
            self.fm_index = None
            self.kmer_index = load_reference_peptides_kmer_index(genome)
        elif index_type == "fm":
            self.fm_index = load_reference_peptides_fm_index(genome)
        else:
            raise ValueError("Uknown index type: %s" % index_type)

    def contains(self, kmer):
        if self.index_type == "kmer":
            return kmer in self.kmer_index
        elif self.index_type == "fm":
            return len(self.fm_index.search(kmer)) > 0
