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
import json
import logging
import os
import pickle
import threading
from functools import lru_cache
from typing import Iterable, Optional

from platformdirs import user_cache_dir
from pyensembl import Genome
from tqdm import tqdm

from .config.defaults import DEFAULT_MAX_KMER_LENGTH, DEFAULT_MIN_KMER_LENGTH
from .identifiers import normalize_ensembl_gene_id
from .vaccine_antigen import (
    SelfReferenceMatch,
    SelfReferenceSource,
    VaccineAntigen,
)

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

    return frozenset(
        normalize_ensembl_gene_id(gene_id)
        for gene_id in cta_unfiltered_gene_ids()
    )


def cta_source_gene_ids_for_genome(genome) -> frozenset[str]:
    """Human CTA source exclusions, or an empty set for other species."""
    species = getattr(genome, "species", None)
    latin_name = str(getattr(species, "latin_name", ""))
    normalized = latin_name.strip().casefold().replace("_", " ")
    if normalized not in {"homo sapiens", "human"}:
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

# Protein/source snapshots avoid reconstructing every pyensembl Transcript for
# each antigen. Only concrete pyensembl genomes are cached, using their public
# serialized definition and resolved local source-file identities.
_protein_source_snapshot_cache: dict[str, tuple] = {}
_protein_source_snapshot_cache_lock = threading.Lock()

# Source-file digests make dataset identities independent of installation
# paths. File metadata avoids re-reading multi-gigabyte Ensembl inputs when
# they have not changed.
_ensembl_source_digest_cache: dict[tuple, str] = {}
_ensembl_source_digest_cache_lock = threading.Lock()


def _protein_content_digest(proteins: dict[str, str]) -> str:
    """Return a deterministic digest of transcript IDs and protein sequences."""
    digest = hashlib.sha256()
    for transcript_id, sequence in sorted(proteins.items()):
        for value in (transcript_id, sequence):
            encoded = value.encode("utf-8")
            digest.update(len(encoded).to_bytes(8, "big"))
            digest.update(encoded)
    return digest.hexdigest()


def ensembl_dataset_cache_identity(genome) -> Optional[str]:
    """Return a content-based identity for an installed Ensembl dataset.

    The identity combines the non-location parts of the public pyensembl
    genome definition with SHA-256 digests of every resolved annotation,
    transcript, and protein source file. Equivalent datasets installed at
    different paths therefore share an identity, while a changed source file
    invalidates cached protein provenance.

    ``None`` means the object is not a concrete pyensembl ``Genome`` or its
    required local files cannot be resolved.
    """
    if not isinstance(genome, Genome):
        return None
    files = []
    try:
        definition = genome.to_dict()
        source_paths = [definition.get("gtf_path_or_url")]
        source_paths.extend(
            definition.get("transcript_fasta_paths_or_urls") or []
        )
        source_paths.extend(
            definition.get("protein_fasta_paths_or_urls") or []
        )
        source_paths = [path for path in source_paths if path]
        cached_paths = genome.required_local_files()
        if len(source_paths) != len(cached_paths):
            return None
        for source_path, cached_path in zip(source_paths, cached_paths):
            use_source_directly = (
                "://" not in source_path
                and not definition.get("copy_local_files_to_cache")
                and not definition.get("decompress_on_download")
            )
            path = source_path if use_source_directly else cached_path
            resolved_path = os.path.realpath(path)
            stat = os.stat(resolved_path)
            digest_key = (
                resolved_path,
                stat.st_dev,
                stat.st_ino,
                stat.st_size,
                stat.st_mtime_ns,
            )
            with _ensembl_source_digest_cache_lock:
                source_digest = _ensembl_source_digest_cache.get(digest_key)
            if source_digest is None:
                digest = hashlib.sha256()
                with open(resolved_path, "rb") as source_file:
                    while chunk := source_file.read(1024 * 1024):
                        digest.update(chunk)
                source_digest = digest.hexdigest()
                with _ensembl_source_digest_cache_lock:
                    _ensembl_source_digest_cache[digest_key] = source_digest
            files.append((stat.st_size, source_digest))
    except OSError:
        return None
    location_keys = {
        "cache_directory_path",
        "copy_local_files_to_cache",
        "decompress_on_download",
        "gtf_path_or_url",
        "protein_fasta_paths_or_urls",
        "transcript_fasta_paths_or_urls",
    }
    content_definition = {
        key: value for key, value in definition.items()
        if key not in location_keys
    }
    payload = json.dumps(
        {
            "genome": content_definition,
            "files": files,
        },
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def clear_reference_proteome_caches() -> None:
    """Clear every in-memory reference-proteome and source-data cache."""
    for lock, cache in (
        (_kmer_set_cache_lock, _kmer_set_cache),
        (_filtered_kmer_set_cache_lock, _filtered_kmer_set_cache),
        (_protein_source_snapshot_cache_lock, _protein_source_snapshot_cache),
        (_ensembl_source_digest_cache_lock, _ensembl_source_digest_cache),
    ):
        with lock:
            cache.clear()
    oncoref_cta_source_gene_ids.cache_clear()


def _protein_source_snapshot(genome):
    """Return distinct protein sequences with every Ensembl source."""
    cache_key = ensembl_dataset_cache_identity(genome)
    if cache_key is not None:
        with _protein_source_snapshot_cache_lock:
            cached = _protein_source_snapshot_cache.get(cache_key)
        if cached is not None:
            return cached

    species = str(
        getattr(getattr(genome, "species", None), "latin_name", "") or ""
    )
    sources_by_sequence: dict[str, set[SelfReferenceSource]] = {}
    for transcript in genome.transcripts():
        if not transcript.is_protein_coding or not transcript.protein_sequence:
            continue
        gene_id = normalize_ensembl_gene_id(
            getattr(transcript, "gene_id", "") or ""
        )
        if not gene_id:
            raise ValueError(
                "Protein-coding self-reference transcript is missing a gene ID"
            )
        source = SelfReferenceSource(
            gene_id=gene_id,
            transcript_id=str(
                getattr(transcript, "transcript_id", "") or ""
            ),
            protein_id=str(getattr(transcript, "protein_id", "") or ""),
            gene_name=str(getattr(transcript, "gene_name", "") or ""),
            species=species,
        )
        sources_by_sequence.setdefault(
            transcript.protein_sequence, set()
        ).add(source)

    snapshot = tuple(
        (
            sequence,
            tuple(sorted(
                sources,
                key=lambda source: (
                    source.gene_id,
                    source.transcript_id,
                    source.protein_id,
                    source.gene_name,
                    source.species,
                ),
            )),
        )
        for sequence, sources in sorted(sources_by_sequence.items())
    )
    if cache_key is not None:
        with _protein_source_snapshot_cache_lock:
            _protein_source_snapshot_cache[cache_key] = snapshot
    return snapshot


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


def resolve_reference_kmer_lengths(
    min_kmer_length, max_kmer_length, epitope_lengths
):
    """Resolve an index range from explicit bounds, epitope lengths, or defaults."""
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
        normalize_ensembl_gene_id(gene_id)
        for gene_id in (exclude_gene_ids or [])
    }

    proteins = {}
    num_excluded_transcripts = 0
    for t in genome.transcripts():
        gene_id = normalize_ensembl_gene_id(t.gene_id)
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


def self_reference_matches(
        peptides: Iterable[str],
        antigen: VaccineAntigen,
        genome) -> dict[str, SelfReferenceMatch]:
    """Find every retained Ensembl source of each exact peptide.

    The genome is traversed once for the batch. Source-gene exclusions come
    exclusively from ``antigen.self_reference_excluded_gene_ids``; excluding a
    CTA source therefore cannot erase a match encoded by a retained non-CTA
    gene. With no genome, results remain explicitly provenance-incomplete.
    """
    unique_peptides = tuple(dict.fromkeys(peptides))
    if not unique_peptides:
        return {}
    genome_release = str(getattr(genome, "release", "") or "")
    if genome is None:
        return {
            peptide: antigen.self_reference_match(
                peptide, False, genome_release=genome_release)
            for peptide in unique_peptides
        }

    excluded_gene_ids = set(antigen.self_reference_excluded_gene_ids)
    sources_by_peptide: dict[str, set[SelfReferenceSource]] = {
        peptide: set() for peptide in unique_peptides
    }
    for sequence, sequence_sources in _protein_source_snapshot(genome):
        retained_sources = tuple(
            source for source in sequence_sources
            if source.gene_id not in excluded_gene_ids
        )
        if not retained_sources:
            continue
        matching_peptides = [
            peptide for peptide in unique_peptides if peptide in sequence
        ]
        if not matching_peptides:
            continue
        for peptide in matching_peptides:
            sources_by_peptide[peptide].update(retained_sources)

    result = {}
    for peptide, sources in sources_by_peptide.items():
        ordered_sources = tuple(sorted(
            sources,
            key=lambda source: (
                source.gene_id,
                source.transcript_id,
                source.protein_id,
                source.gene_name,
                source.species,
            ),
        ))
        result[peptide] = SelfReferenceMatch(
            peptide=peptide,
            occurs=bool(ordered_sources),
            antigen_kind=antigen.kind,
            excluded_gene_ids=antigen.self_reference_excluded_gene_ids,
            sources=ordered_sources,
            source_provenance_complete=True,
            genome_release=genome_release,
        )
    return result


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
        min_kmer_length, max_kmer_length = resolve_reference_kmer_lengths(
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
        min_kmer_length, max_kmer_length = resolve_reference_kmer_lengths(
            min_kmer_length, max_kmer_length, epitope_lengths)

        if genome is None:
            return cls(
                None,
                min_kmer_length=min_kmer_length,
                max_kmer_length=max_kmer_length,
            )

        all_exclude_ids = {
            normalize_ensembl_gene_id(gene_id)
            for gene_id in (exclude_gene_ids or [])
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
