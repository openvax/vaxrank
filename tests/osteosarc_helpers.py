"""Offline original-RNA fixtures shared by downstream integration regressions."""

import gzip
from hashlib import sha256
import json
from pathlib import Path

import pysam
from pyensembl import Genome
from varcode import Variant


DATA = Path(__file__).parent / "data" / "osteosarc"
SELECTION = json.loads((DATA / "selection.json").read_text())
MANIFEST = json.loads((DATA / "manifest.json").read_text())


def load_osteosarc(directory):
    """Check pinned data and build isolated reference and alignment indices."""
    reference = DATA / "protein_reference"
    metadata = json.loads((reference / "protein_reference_manifest.json").read_text())
    for filename, checksums in metadata["files"].items():
        data = (reference / filename).read_bytes()
        assert sha256(data).hexdigest() == checksums["subset_sha256"]
        assert sha256(gzip.decompress(data)).hexdigest() == checksums["uncompressed_sha256"]
    genome = Genome(
        reference_name="GRCh38-osteosarc-six-transcript-subset",
        annotation_name="osteosarc-ensembl-subset", annotation_version=87,
        gtf_path_or_url=str(reference / "reference.gtf.gz"),
        transcript_fasta_paths_or_urls=[str(reference / "reference.cdna.fa.gz")],
        protein_fasta_paths_or_urls=[str(reference / "reference.pep.fa.gz")],
        copy_local_files_to_cache=True,
        cache_directory_path=str(directory / "reference"),
    )
    genome.index()
    variants = {
        record["gene"]: Variant(
            record["chrom"].removeprefix("chr"), int(record["pos"]),
            record["ref"], record["alt"], ensembl=genome)
        for record in SELECTION["variants"]
    }
    bams = {}
    for name, dataset in MANIFEST["datasets"].items():
        sam_data = gzip.decompress((DATA / dataset["file"]).read_bytes())
        assert sha256(sam_data).hexdigest() == dataset["sam_sha256"]
        sam = directory / (name + ".sam")
        sam.write_bytes(sam_data)
        bam = directory / (name + ".bam")
        pysam.sort("--no-PG", "-o", str(bam), str(sam))
        pysam.index(str(bam))
        bams[name] = str(bam)
    # Independent expectations: original reference proteins plus the reported
    # protein edits, not output snapshots or a second call to the annotator.
    with pysam.FastxFile(str(reference / "reference.pep.fa.gz")) as fasta:
        proteins = {
            token.split(":", 1)[1].split(".")[0]: record.sequence
            for record in fasta for token in record.comment.split()
            if token.startswith("transcript:")
        }
    edits = {"DYNC1H1": (313, "V", "I"), "EXOC4": (33, "S", "I"),
             "H1-2": (196, "AAKPK", ""), "GTF3C5": (502, "EEEE", "")}
    expected = {}
    for gene, (start, ref, alt) in edits.items():
        protein = proteins[metadata["transcripts"][gene]]
        assert protein[start:start + len(ref)] == ref
        expected[gene] = protein[:start] + alt + protein[start + len(ref):]
    return variants, bams, expected
