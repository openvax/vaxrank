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

# Cancer-testis antigen (CTA) Ensembl gene IDs, derived from pirlygenes
# CTA_gene_ids() (filtered & expressed). Inlined to avoid the heavy dependency.
CTA_GENE_IDS = frozenset({
    "ENSG00000006047", "ENSG00000007350", "ENSG00000010318", "ENSG00000042813",
    "ENSG00000046774", "ENSG00000054796", "ENSG00000068985", "ENSG00000077800",
    "ENSG00000077935", "ENSG00000092345", "ENSG00000095627", "ENSG00000099399",
    "ENSG00000102021", "ENSG00000103023", "ENSG00000104755", "ENSG00000104804",
    "ENSG00000104901", "ENSG00000106013", "ENSG00000106304", "ENSG00000114487",
    "ENSG00000117148", "ENSG00000118434", "ENSG00000118491", "ENSG00000121101",
    "ENSG00000122304", "ENSG00000123576", "ENSG00000123584", "ENSG00000124092",
    "ENSG00000124227", "ENSG00000124260", "ENSG00000124467", "ENSG00000124490",
    "ENSG00000125207", "ENSG00000126467", "ENSG00000126752", "ENSG00000126856",
    "ENSG00000129862", "ENSG00000129864", "ENSG00000131126", "ENSG00000132446",
    "ENSG00000135248", "ENSG00000137090", "ENSG00000137948", "ENSG00000139351",
    "ENSG00000140478", "ENSG00000141096", "ENSG00000141255", "ENSG00000141316",
    "ENSG00000141371", "ENSG00000141437", "ENSG00000142025", "ENSG00000142698",
    "ENSG00000143006", "ENSG00000143452", "ENSG00000144962", "ENSG00000145309",
    "ENSG00000146047", "ENSG00000147081", "ENSG00000147183", "ENSG00000147378",
    "ENSG00000147381", "ENSG00000148156", "ENSG00000150628", "ENSG00000151962",
    "ENSG00000152430", "ENSG00000152670", "ENSG00000153060", "ENSG00000153779",
    "ENSG00000155087", "ENSG00000155495", "ENSG00000156009", "ENSG00000156269",
    "ENSG00000158639", "ENSG00000159289", "ENSG00000159708", "ENSG00000161860",
    "ENSG00000161973", "ENSG00000162039", "ENSG00000162641", "ENSG00000162771",
    "ENSG00000162843", "ENSG00000163114", "ENSG00000163530", "ENSG00000164113",
    "ENSG00000164304", "ENSG00000165583", "ENSG00000165584", "ENSG00000166049",
    "ENSG00000166069", "ENSG00000166118", "ENSG00000166796", "ENSG00000168594",
    "ENSG00000168757", "ENSG00000169059", "ENSG00000169551", "ENSG00000169800",
    "ENSG00000170516", "ENSG00000170748", "ENSG00000170950", "ENSG00000170965",
    "ENSG00000171402", "ENSG00000171794", "ENSG00000172073", "ENSG00000172717",
    "ENSG00000174015", "ENSG00000174016", "ENSG00000174898", "ENSG00000175646",
    "ENSG00000175809", "ENSG00000175820", "ENSG00000175877", "ENSG00000176256",
    "ENSG00000176566", "ENSG00000176635", "ENSG00000176774", "ENSG00000176988",
    "ENSG00000177294", "ENSG00000177673", "ENSG00000177689", "ENSG00000177947",
    "ENSG00000177992", "ENSG00000178021", "ENSG00000178279", "ENSG00000178997",
    "ENSG00000179046", "ENSG00000179088", "ENSG00000179407", "ENSG00000180138",
    "ENSG00000180336", "ENSG00000181323", "ENSG00000181433", "ENSG00000182308",
    "ENSG00000182459", "ENSG00000182583", "ENSG00000183206", "ENSG00000183434",
    "ENSG00000183668", "ENSG00000184033", "ENSG00000184361", "ENSG00000184507",
    "ENSG00000184650", "ENSG00000184735", "ENSG00000185264", "ENSG00000185686",
    "ENSG00000185863", "ENSG00000186075", "ENSG00000186118", "ENSG00000186451",
    "ENSG00000186788", "ENSG00000187003", "ENSG00000187191", "ENSG00000187475",
    "ENSG00000187772", "ENSG00000188120", "ENSG00000188425", "ENSG00000188782",
    "ENSG00000189023", "ENSG00000189064", "ENSG00000189186", "ENSG00000189252",
    "ENSG00000189326", "ENSG00000189357", "ENSG00000196406", "ENSG00000196553",
    "ENSG00000196862", "ENSG00000197172", "ENSG00000198021", "ENSG00000198033",
    "ENSG00000198573", "ENSG00000198681", "ENSG00000198759", "ENSG00000198765",
    "ENSG00000198870", "ENSG00000203907", "ENSG00000203909", "ENSG00000203926",
    "ENSG00000204279", "ENSG00000204363", "ENSG00000204379", "ENSG00000204450",
    "ENSG00000204849", "ENSG00000204941", "ENSG00000205359", "ENSG00000205642",
    "ENSG00000205777", "ENSG00000212710", "ENSG00000213030", "ENSG00000213218",
    "ENSG00000213401", "ENSG00000214107", "ENSG00000215113", "ENSG00000215115",
    "ENSG00000215269", "ENSG00000215529", "ENSG00000215817", "ENSG00000216649",
    "ENSG00000221826", "ENSG00000221867", "ENSG00000224659", "ENSG00000224902",
    "ENSG00000226023", "ENSG00000226650", "ENSG00000226685", "ENSG00000226929",
    "ENSG00000227234", "ENSG00000228517", "ENSG00000228927", "ENSG00000229549",
    "ENSG00000230594", "ENSG00000230601", "ENSG00000231924", "ENSG00000233803",
    "ENSG00000234068", "ENSG00000234414", "ENSG00000235631", "ENSG00000235699",
    "ENSG00000236126", "ENSG00000236362", "ENSG00000236371", "ENSG00000236424",
    "ENSG00000236446", "ENSG00000236761", "ENSG00000237671", "ENSG00000237957",
    "ENSG00000238269", "ENSG00000240021", "ENSG00000241476", "ENSG00000242362",
    "ENSG00000242389", "ENSG00000242875", "ENSG00000243130", "ENSG00000244395",
    "ENSG00000244588", "ENSG00000250741", "ENSG00000258713", "ENSG00000258992",
    "ENSG00000261649", "ENSG00000267978", "ENSG00000268009", "ENSG00000268447",
    "ENSG00000268606", "ENSG00000268629", "ENSG00000268651", "ENSG00000268696",
    "ENSG00000268988", "ENSG00000269058", "ENSG00000269586", "ENSG00000270806",
    "ENSG00000271321", "ENSG00000271449", "ENSG00000273696", "ENSG00000274274",
    "ENSG00000274391", "ENSG00000275793", "ENSG00000275969", "ENSG00000276040",
    "ENSG00000277322",
})

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
            CandidateEpitope lengths being predicted. Used to derive kmer range
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
        all_exclude_ids = set(exclude_gene_ids or [])

        if exclude_cta_genes:
            logger.info("Excluding %d CTA gene IDs", len(CTA_GENE_IDS))
            all_exclude_ids.update(CTA_GENE_IDS)

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
