"""Offline original-RNA fixtures shared by downstream integration regressions."""

import gzip
from hashlib import sha256
import json
from pathlib import Path

import pysam
from varcode import Variant

from .osteosarc_fixture_support import indexed_genome, verify_gzip_digests


DATA = Path(__file__).parent / "data" / "osteosarc"
SELECTION = json.loads((DATA / "selection.json").read_text())
MANIFEST = json.loads((DATA / "manifest.json").read_text())


def load_osteosarc(directory):
    """Check pinned data and build isolated reference and alignment indices."""
    reference = DATA / "protein_reference"
    metadata = json.loads((reference / "protein_reference_manifest.json").read_text())
    for filename, checksums in metadata["files"].items():
        verify_gzip_digests(
            reference / filename, checksums["subset_sha256"],
            checksums["uncompressed_sha256"], label="protein_reference/" + filename)
    genome = indexed_genome(
        reference,
        reference_name="GRCh38-osteosarc-six-transcript-subset",
        annotation_name="osteosarc-ensembl-subset", annotation_version=87,
        cache_directory=directory / "reference")
    variants = {
        record["gene"]: Variant(
            record["chrom"].removeprefix("chr"), int(record["pos"]),
            record["ref"], record["alt"], ensembl=genome)
        for record in SELECTION["variants"]
    }
    bams = {}
    for name, dataset in MANIFEST["datasets"].items():
        sam_data = gzip.decompress((DATA / dataset["file"]).read_bytes())
        digest = sha256(sam_data).hexdigest()
        if digest != dataset["sam_sha256"]:
            raise ValueError(
                "Fixture checksum mismatch for %s: manifest expects %s, found %s"
                % (dataset["file"], dataset["sam_sha256"], digest))
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
        if protein[start:start + len(ref)] != ref:
            raise ValueError(
                "Reference protein for %s does not carry %r at %d; the pinned "
                "annotation and the documented edit disagree" % (gene, ref, start))
        expected[gene] = protein[:start] + alt + protein[start + len(ref):]
    return variants, bams, expected
