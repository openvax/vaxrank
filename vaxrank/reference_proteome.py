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

"""
Reference proteome indexing for checking if peptides exist in the reference.

Uses a set-based index for O(1) membership testing of peptide kmers.
"""

import os
import logging
import pickle

from datacache import get_data_dir

logger = logging.getLogger(__name__)

# Default kmer length range for epitope peptides
DEFAULT_MIN_KMER_LENGTH = 8
DEFAULT_MAX_KMER_LENGTH = 15


def get_cache_dir() -> str:
    """Get the cache directory for reference peptide indices."""
    cache_dir = get_data_dir(envkey="VAXRANK_REF_PEPTIDES_DIR")
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)
    return cache_dir


def kmer_set_index_path(genome, min_len: int, max_len: int) -> str:
    """Returns path for the cached kmer set index."""
    return os.path.join(
        get_cache_dir(),
        "%s_%d_kmer_set_%d_%d.pkl"
        % (genome.species.latin_name, genome.release, min_len, max_len),
    )


def build_kmer_set_index(
    genome,
    min_len: int = DEFAULT_MIN_KMER_LENGTH,
    max_len: int = DEFAULT_MAX_KMER_LENGTH,
) -> set[str]:
    """
    Build a set of all kmers from protein sequences in the genome.

    Parameters
    ----------
    genome : pyensembl.EnsemblRelease
        Input genome to load for reference peptides

    min_len : int
        Minimum kmer length to index

    max_len : int
        Maximum kmer length to index

    Returns
    -------
    set[str]
        Set of all kmers found in the reference proteome
    """
    logger.info(
        "Building kmer set index for %s release %d (lengths %d-%d)",
        genome.species.latin_name,
        genome.release,
        min_len,
        max_len,
    )
    kmers = set()
    transcripts = genome.transcripts()
    n_transcripts = len(transcripts)

    for i, t in enumerate(transcripts):
        if not t.is_protein_coding:
            continue
        protein = t.protein_sequence
        if protein is None:
            continue

        # Extract all kmers of each length
        for k in range(min_len, max_len + 1):
            for j in range(len(protein) - k + 1):
                kmers.add(protein[j : j + k])

        if (i + 1) % 10000 == 0:
            logger.info("Processed %d/%d transcripts", i + 1, n_transcripts)

    logger.info(
        "Done building kmer set index: %d unique kmers from %d transcripts",
        len(kmers),
        n_transcripts,
    )
    return kmers


def load_kmer_set_index(
    genome,
    min_len: int = DEFAULT_MIN_KMER_LENGTH,
    max_len: int = DEFAULT_MAX_KMER_LENGTH,
    force_reload: bool = False,
) -> set[str]:
    """
    Load or build the kmer set index for the given genome.

    Parameters
    ----------
    genome : pyensembl.EnsemblRelease
        Input genome to load for reference peptides

    min_len : int
        Minimum kmer length to index

    max_len : int
        Maximum kmer length to index

    force_reload : bool
        If true, rebuild index even if cached version exists

    Returns
    -------
    set[str]
        Set of all kmers found in the reference proteome
    """
    path = kmer_set_index_path(genome, min_len, max_len)

    if not force_reload and os.path.exists(path):
        logger.info("Loading cached kmer set index from %s", path)
        with open(path, "rb") as f:
            return pickle.load(f)

    kmers = build_kmer_set_index(genome, min_len, max_len)

    logger.info("Saving kmer set index to %s", path)
    with open(path, "wb") as f:
        pickle.dump(kmers, f)

    return kmers


class ReferenceProteome:
    """
    Index for checking if peptide sequences exist in a reference proteome.

    Uses a set-based index for O(1) membership testing.
    """

    def __init__(
        self,
        genome,
        min_kmer_length: int = DEFAULT_MIN_KMER_LENGTH,
        max_kmer_length: int = DEFAULT_MAX_KMER_LENGTH,
    ):
        """
        Parameters
        ----------
        genome : pyensembl.EnsemblRelease or None
            Input genome for reference peptides. If None, contains() always
            returns False.

        min_kmer_length : int
            Minimum peptide length to index

        max_kmer_length : int
            Maximum peptide length to index
        """
        self.genome = genome
        self.min_kmer_length = min_kmer_length
        self.max_kmer_length = max_kmer_length

        if genome is not None:
            self._kmer_set = load_kmer_set_index(
                genome, min_kmer_length, max_kmer_length
            )
        else:
            self._kmer_set = set()

    def contains(self, peptide: str) -> bool:
        """
        Check if a peptide sequence exists in the reference proteome.

        Parameters
        ----------
        peptide : str
            Peptide sequence to check

        Returns
        -------
        bool
            True if the peptide exists in the reference proteome
        """
        return peptide in self._kmer_set
