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
import threading

from platformdirs import user_cache_dir
from tqdm import tqdm

from .config.defaults import DEFAULT_MAX_KMER_LENGTH, DEFAULT_MIN_KMER_LENGTH

logger = logging.getLogger(__name__)

# In-memory cache for loaded kmer sets to avoid repeated disk reads
# Key: (species, release, min_len, max_len) -> set of kmers
_kmer_set_cache: dict[tuple, set[str]] = {}
_kmer_set_cache_lock = threading.Lock()


def get_cache_dir() -> str:
    """Get the cache directory for reference peptide indices."""
    env_dir = os.environ.get("VAXRANK_REF_PEPTIDES_DIR")
    if env_dir:
        cache_dir = env_dir
    else:
        cache_dir = user_cache_dir("vaxrank")
    os.makedirs(cache_dir, exist_ok=True)
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

    # First, collect unique proteins - many transcripts share the same protein
    # This avoids redundant kmer extraction (250k transcripts -> ~20k proteins)
    unique_proteins = set()
    transcripts = genome.transcripts()
    for t in tqdm(transcripts, desc="Collecting unique proteins", unit="transcripts"):
        if t.is_protein_coding and t.protein_sequence:
            unique_proteins.add(t.protein_sequence)

    logger.info(
        "Found %d unique proteins from %d transcripts",
        len(unique_proteins),
        len(transcripts),
    )

    # Extract kmers from unique proteins only, using batch updates
    kmers = set()
    for protein in tqdm(unique_proteins, desc="Extracting kmers", unit="proteins"):
        kmers.update(
            protein[j : j + k]
            for k in range(min_len, max_len + 1)
            for j in range(len(protein) - k + 1)
        )

    logger.info(
        "Done building kmer set index: %d unique kmers from %d proteins",
        len(kmers),
        len(unique_proteins),
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
    cache_key = (genome.species.latin_name, genome.release, min_len, max_len)

    with _kmer_set_cache_lock:
        # Check in-memory cache first to avoid repeated disk reads
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


def extract_kmers(sequences, min_len, max_len):
    """Extract all kmers from an iterable of protein sequences."""
    kmers = set()
    for seq in sequences:
        kmers.update(
            seq[j: j + k]
            for k in range(min_len, max_len + 1)
            for j in range(len(seq) - k + 1)
        )
    return kmers


def genome_protein_dict(genome, exclude_gene_ids=None):
    """
    Build a dict of transcript_id -> protein_sequence from a pyensembl genome,
    optionally excluding proteins from specific genes.

    Parameters
    ----------
    genome : pyensembl.EnsemblRelease
    exclude_gene_ids : set of str, optional
        Ensembl gene IDs whose proteins should be excluded.

    Returns
    -------
    dict[str, str]
    """
    if exclude_gene_ids:
        exclude_tids = set()
        for gene_id in exclude_gene_ids:
            try:
                exclude_tids.update(genome.transcript_ids_of_gene_id(gene_id))
            except ValueError:
                logger.warning("Gene ID %s not found in genome, skipping", gene_id)
        logger.info(
            "Excluding %d transcripts from %d genes",
            len(exclude_tids), len(exclude_gene_ids))
    else:
        exclude_tids = set()

    proteins = {}
    for t in genome.transcripts():
        if t.is_protein_coding and t.protein_sequence and t.transcript_id not in exclude_tids:
            proteins[t.transcript_id] = t.protein_sequence
    return proteins


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
        Backwards-compatible constructor: load from a pyensembl genome.

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

    @classmethod
    def from_sequences(cls, proteins, min_kmer_length=None, max_kmer_length=None,
                       epitope_lengths=None):
        """
        Build a ReferenceProteome from a dict of protein sequences.

        Parameters
        ----------
        proteins : dict[str, str]
            Mapping of protein name/ID to amino acid sequence.

        min_kmer_length : int, optional
            Minimum kmer length. If None, derived from epitope_lengths.

        max_kmer_length : int, optional
            Maximum kmer length. If None, derived from epitope_lengths.

        epitope_lengths : collection of int, optional
            Epitope lengths being predicted. Used to derive kmer range
            when min/max not explicitly set.
        """
        if epitope_lengths and (min_kmer_length is None or max_kmer_length is None):
            lengths = sorted(epitope_lengths)
            if min_kmer_length is None:
                min_kmer_length = lengths[0]
            if max_kmer_length is None:
                max_kmer_length = lengths[-1]
        if min_kmer_length is None:
            min_kmer_length = DEFAULT_MIN_KMER_LENGTH
        if max_kmer_length is None:
            max_kmer_length = DEFAULT_MAX_KMER_LENGTH

        unique_seqs = set(proteins.values())
        logger.info(
            "Building kmer set from %d proteins (%d unique), lengths %d-%d",
            len(proteins), len(unique_seqs), min_kmer_length, max_kmer_length)
        kmer_set = extract_kmers(unique_seqs, min_kmer_length, max_kmer_length)
        logger.info("Built kmer set with %d unique kmers", len(kmer_set))

        obj = cls.__new__(cls)
        obj.genome = None
        obj.min_kmer_length = min_kmer_length
        obj.max_kmer_length = max_kmer_length
        obj._kmer_set = kmer_set
        return obj

    @classmethod
    def from_genome(cls, genome, exclude_gene_ids=None,
                    exclude_pirlygenes_cta=False, exclude_fasta=None,
                    epitope_lengths=None,
                    min_kmer_length=None, max_kmer_length=None):
        """
        Build a ReferenceProteome from a pyensembl genome with optional
        exclusions.

        Parameters
        ----------
        genome : pyensembl.EnsemblRelease
        exclude_gene_ids : set of str, optional
            Ensembl gene IDs to exclude.
        exclude_pirlygenes_cta : bool
            If True, exclude cancer-testis antigen genes from pirlygenes.
            Requires pirlygenes to be installed.
        exclude_fasta : str or Path, optional
            Path to a FASTA file whose sequences should be excluded.
        epitope_lengths : collection of int, optional
            Used to derive kmer range.
        min_kmer_length : int, optional
        max_kmer_length : int, optional
        """
        all_exclude_ids = set(exclude_gene_ids or [])

        if exclude_pirlygenes_cta:
            try:
                from pirlygenes.gene_sets_cancer import CTA_gene_ids
            except ImportError:
                raise ImportError(
                    "Config has exclude_pirlygenes_cta: true but pirlygenes "
                    "is not installed. Install with: pip install pirlygenes"
                ) from None
            cta_ids = CTA_gene_ids()
            logger.info("Excluding %d CTA gene IDs from pirlygenes", len(cta_ids))
            all_exclude_ids.update(cta_ids)

        proteins = genome_protein_dict(genome, exclude_gene_ids=all_exclude_ids)

        if exclude_fasta:
            exclude_seqs = set(_read_fasta(exclude_fasta).values())
            before = len(proteins)
            proteins = {
                tid: seq for tid, seq in proteins.items()
                if seq not in exclude_seqs
            }
            logger.info(
                "Excluded %d proteins matching FASTA %s",
                before - len(proteins), exclude_fasta)

        obj = cls.from_sequences(
            proteins,
            min_kmer_length=min_kmer_length,
            max_kmer_length=max_kmer_length,
            epitope_lengths=epitope_lengths,
        )
        obj.genome = genome
        return obj

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


def _read_fasta(path):
    """Read a FASTA file into a dict of name -> sequence."""
    proteins = {}
    current_name = None
    current_seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name is not None:
                    proteins[current_name] = "".join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            elif current_name is not None:
                current_seq.append(line)
    if current_name is not None:
        proteins[current_name] = "".join(current_seq)
    return proteins
