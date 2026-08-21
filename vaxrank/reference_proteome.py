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
import hashlib
import io
import logging
import os
import pickle
import threading
from functools import lru_cache

from platformdirs import user_cache_dir
from tqdm import tqdm

from .config.defaults import DEFAULT_MAX_KMER_LENGTH, DEFAULT_MIN_KMER_LENGTH

logger = logging.getLogger(__name__)


@lru_cache(maxsize=1)
def oncoref_cta_source_gene_ids() -> frozenset[str]:
    """Full human CTA candidate universe owned by oncoref.

    The broad, unfiltered universe is intentional for self-reference source
    exclusions: a valid CTA peptide must not cancel itself merely because the
    same sequence occurs in a held-out CTA paralog. Target admission uses the
    narrower ``oncoref.cta.cta_gene_ids()`` set in the antigen layer.
    """
    from oncoref.cta import cta_unfiltered_gene_ids

    return frozenset(str(gene_id).split(".")[0]
                     for gene_id in cta_unfiltered_gene_ids())


def _is_human_genome(genome) -> bool:
    if genome is None:
        return False
    species = getattr(genome, "species", None)
    latin_name = str(getattr(species, "latin_name", ""))
    normalized = latin_name.strip().casefold().replace("_", " ")
    return normalized in {"homo sapiens", "human"}


def cta_source_gene_ids_for_genome(genome) -> frozenset[str]:
    """Human CTA source exclusions, or an empty set for other species."""
    if not _is_human_genome(genome):
        return frozenset()
    return oncoref_cta_source_gene_ids()

# In-memory cache for loaded kmer sets to avoid repeated disk reads
# Key: (species, release, min_len, max_len) -> set of kmers
_kmer_set_cache: dict[tuple, set[str]] = {}
_kmer_set_cache_lock = threading.Lock()

# Filtered indexes are shared within a run. Their keys describe the retained
# protein content and filtering policy, so equivalent inputs can share work
# while distinct genome objects or releases cannot contaminate one another.
_filtered_kmer_set_cache: dict[tuple, set[str]] = {}
_filtered_kmer_set_cache_lock = threading.Lock()


def _protein_content_digest(proteins: dict[str, str]) -> str:
    """Return a deterministic digest of transcript IDs and protein sequences."""
    digest = hashlib.sha256()
    for transcript_id, sequence in sorted(proteins.items()):
        for value in (transcript_id, sequence):
            encoded = value.encode("utf-8")
            digest.update(len(encoded).to_bytes(8, "big"))
            digest.update(encoded)
    return digest.hexdigest()


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


def _resolve_kmer_lengths(min_kmer_length, max_kmer_length, epitope_lengths):
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
    return min_kmer_length, max_kmer_length


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
    excluded_gene_ids = {
        str(gene_id).split(".")[0] for gene_id in (exclude_gene_ids or [])
    }

    proteins = {}
    num_excluded_transcripts = 0
    for t in genome.transcripts():
        gene_id = str(t.gene_id).split(".")[0]
        if gene_id in excluded_gene_ids:
            num_excluded_transcripts += 1
        elif t.is_protein_coding and t.protein_sequence:
            proteins[t.transcript_id] = t.protein_sequence
    if excluded_gene_ids:
        logger.info(
            "Excluded %d transcripts from %d source genes",
            num_excluded_transcripts,
            len(excluded_gene_ids),
        )
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
        self.excluded_gene_ids = frozenset()

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
            CandidateEpitope lengths being predicted. Used to derive kmer range
            when min/max not explicitly set.
        """
        min_kmer_length, max_kmer_length = _resolve_kmer_lengths(
            min_kmer_length, max_kmer_length, epitope_lengths)

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
        obj.excluded_gene_ids = frozenset()
        return obj

    @classmethod
    def from_genome(cls, genome, exclude_gene_ids=None,
                    exclude_cta_genes=False, exclude_fasta=None,
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
        exclude_cta_genes : bool
            If True, exclude cancer-testis antigen (CTA) genes.
        exclude_fasta : str or Path, optional
            Path to a FASTA file whose sequences should be excluded.
        epitope_lengths : collection of int, optional
            Used to derive kmer range.
        min_kmer_length : int, optional
        max_kmer_length : int, optional
        """
        min_kmer_length, max_kmer_length = _resolve_kmer_lengths(
            min_kmer_length, max_kmer_length, epitope_lengths)

        if genome is None:
            return cls(
                None,
                min_kmer_length=min_kmer_length,
                max_kmer_length=max_kmer_length,
            )

        all_exclude_ids = {
            str(gene_id).split(".")[0] for gene_id in (exclude_gene_ids or [])
        }

        if exclude_cta_genes:
            cta_gene_ids = cta_source_gene_ids_for_genome(genome)
            logger.info(
                "Excluding %d oncoref CTA source gene IDs", len(cta_gene_ids))
            all_exclude_ids.update(cta_gene_ids)

        # With no exclusions this is exactly the ordinary reference proteome;
        # reuse its disk/in-memory cache rather than rebuilding from transcripts.
        if exclude_cta_genes and not all_exclude_ids and not exclude_fasta:
            return cls(
                genome,
                min_kmer_length=min_kmer_length,
                max_kmer_length=max_kmer_length,
            )

        proteins = genome_protein_dict(
            genome, exclude_gene_ids=all_exclude_ids)
        cache_key = None
        kmer_set = None

        if not exclude_fasta:
            cache_key = (
                _protein_content_digest(proteins),
                frozenset(all_exclude_ids),
                min_kmer_length,
                max_kmer_length,
            )
            with _filtered_kmer_set_cache_lock:
                kmer_set = _filtered_kmer_set_cache.get(cache_key)

        if kmer_set is None:
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

            kmer_set = extract_kmers(
                set(proteins.values()), min_kmer_length, max_kmer_length)
            if not exclude_fasta:
                with _filtered_kmer_set_cache_lock:
                    _filtered_kmer_set_cache[cache_key] = kmer_set

        obj = cls.__new__(cls)
        obj.genome = genome
        obj.min_kmer_length = min_kmer_length
        obj.max_kmer_length = max_kmer_length
        obj._kmer_set = kmer_set
        obj.excluded_gene_ids = frozenset(all_exclude_ids)
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
