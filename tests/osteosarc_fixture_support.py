"""Shared Osteosarc integrity checks and consumer-owned reference construction."""

from osteosarc import (
    verify_digest as verify_digest,
    verify_gzip_digests as verify_gzip_digests,
    verify_manifest_files as verify_manifest_files,
)


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
