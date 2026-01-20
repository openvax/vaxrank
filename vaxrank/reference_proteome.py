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

import gzip
import io
import os
import logging
import pickle

from datacache import get_data_dir
from tqdm import tqdm

logger = logging.getLogger(__name__)

# Default kmer length range for epitope peptides
DEFAULT_MIN_KMER_LENGTH = 8
DEFAULT_MAX_KMER_LENGTH = 15

# In-memory cache for loaded kmer sets to avoid repeated disk reads
# Key: (species, release, min_len, max_len) -> set of kmers
_kmer_set_cache: dict[tuple, set[str]] = {}


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
        "%s_%d_kmer_set_%d_%d.pkl.gz"
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

    for t in tqdm(transcripts, desc="Extracting kmers", unit="transcripts"):
        if not t.is_protein_coding:
            continue
        protein = t.protein_sequence
        if protein is None:
            continue

        # Extract all kmers of each length
        for k in range(min_len, max_len + 1):
            for j in range(len(protein) - k + 1):
                kmers.add(protein[j : j + k])

    logger.info(
        "Done building kmer set index: %d unique kmers from %d transcripts",
        len(kmers),
        len(transcripts),
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
    # Check in-memory cache first to avoid repeated disk reads
    cache_key = (genome.species.latin_name, genome.release, min_len, max_len)
    if not force_reload and cache_key in _kmer_set_cache:
        logger.debug("Using in-memory cached kmer set for %s", cache_key)
        return _kmer_set_cache[cache_key]

    path = kmer_set_index_path(genome, min_len, max_len)
    # Also check for legacy uncompressed path
    legacy_path = path.replace(".pkl.gz", ".pkl")

    if not force_reload and os.path.exists(path):
        file_size = os.path.getsize(path)
        with tqdm(
            total=file_size, unit="B", unit_scale=True, desc="Loading reference proteome"
        ) as pbar:
            with open(path, "rb") as raw_f:
                # Wrap raw file to track compressed bytes read
                class ProgressReader:
                    def __init__(self, f, pbar):
                        self._f = f
                        self._pbar = pbar
                    def read(self, n=-1):
                        data = self._f.read(n)
                        self._pbar.update(len(data))
                        return data

                with gzip.GzipFile(fileobj=ProgressReader(raw_f, pbar)) as f:
                    kmers = pickle.load(f)
        _kmer_set_cache[cache_key] = kmers
        return kmers

    # Check for legacy uncompressed file
    if not force_reload and os.path.exists(legacy_path):
        file_size = os.path.getsize(legacy_path)
        with tqdm(
            total=file_size, unit="B", unit_scale=True, desc="Loading reference proteome"
        ) as pbar:
            with open(legacy_path, "rb") as f:
                data = io.BytesIO()
                while chunk := f.read(1024 * 1024):
                    data.write(chunk)
                    pbar.update(len(chunk))
                data.seek(0)
                kmers = pickle.load(data)
        _kmer_set_cache[cache_key] = kmers
        # Save as compressed for next time
        logger.info("Converting to compressed format: %s", path)
        with gzip.open(path, "wb", compresslevel=6) as f:
            pickle.dump(kmers, f, protocol=pickle.HIGHEST_PROTOCOL)
        return kmers

    kmers = build_kmer_set_index(genome, min_len, max_len)

    logger.info("Saving kmer set index to %s", path)
    with gzip.open(path, "wb", compresslevel=6) as f:
        pickle.dump(kmers, f, protocol=pickle.HIGHEST_PROTOCOL)

    _kmer_set_cache[cache_key] = kmers
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
