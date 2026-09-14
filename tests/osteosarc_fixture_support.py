"""Shared pinned-fixture verification and reference indexing.

Both osteosarc fixture families check committed bytes against a manifest of
checksums and build an isolated pyensembl ``Genome`` over a partial GRCh38
annotation. Those two steps were reimplemented per helper, so a fix to
checksum reporting or cache isolation had to be found and applied in each
copy, and the copies could drift apart unnoticed.

Verification raises rather than asserting. A bare ``assert`` disappears
under ``python -O`` / ``PYTHONOPTIMIZE``, which would load tampered or
truncated fixture bytes with no check at all -- defeating the only
tamper-detection these pinned-data designs have.
"""

import gzip
from hashlib import sha256


def verify_digest(path, expected, *, label=None):
    """Raise unless ``path``'s bytes hash to ``expected``. Returns the bytes."""
    data = path.read_bytes()
    digest = sha256(data).hexdigest()
    if digest != expected:
        raise ValueError(
            "Fixture checksum mismatch for %s: manifest expects %s, found %s"
            % (label or path.name, expected, digest))
    return data


def verify_gzip_digests(path, compressed, uncompressed, *, label=None):
    """Verify a gzipped fixture as stored and as decompressed.

    Returns the decompressed bytes. Checking both catches a re-compression
    that preserves content as well as a content change that preserves size.
    """
    data = verify_digest(path, compressed, label=label)
    plain = gzip.decompress(data)
    digest = sha256(plain).hexdigest()
    if digest != uncompressed:
        raise ValueError(
            "Uncompressed fixture checksum mismatch for %s: manifest expects "
            "%s, found %s" % (label or path.name, uncompressed, digest))
    return plain


def verify_manifest_files(root, files, *, label=None):
    """Verify a ``{relative path: sha256}`` mapping under ``root``."""
    for filename, digest in files.items():
        verify_digest(root / filename, digest,
                      label="%s/%s" % (label, filename) if label else filename)


def indexed_genome(
        reference_directory,
        *,
        reference_name,
        annotation_name,
        annotation_version,
        cache_directory,
        gtf="reference.gtf.gz",
        cdna="reference.cdna.fa.gz",
        peptides="reference.pep.fa.gz"):
    """Build and index a partial-annotation Genome with its own cache identity.

    ``reference_name`` must be distinct per fixture family: pyensembl and
    varcode key reference-derived data by that identity, so two partial
    annotations sharing a name would serve each other's cached transcripts.
    A partial annotation must never populate the full-genome GRCh38 cache.
    """
    from pyensembl import Genome

    genome = Genome(
        reference_name=reference_name,
        annotation_name=annotation_name,
        annotation_version=annotation_version,
        gtf_path_or_url=str(reference_directory / gtf),
        transcript_fasta_paths_or_urls=[str(reference_directory / cdna)],
        protein_fasta_paths_or_urls=[str(reference_directory / peptides)],
        copy_local_files_to_cache=True,
        cache_directory_path=str(cache_directory),
    )
    genome.index()
    return genome
