"""Maintainer-only original-read import for the independently documented pilot.

Usage: python import_rna.py ISOVAR_REPO ISOVAR_REGIONAL_CACHE
Never called by tests. Reads immutable released Git objects and checks cached
regional BAMs against the published receipts. No network or whole-BAM download.
Selected reads are not random samples, full-source counts or VAF estimates.
"""

import argparse
from collections import Counter
import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import pysam


COMMIT = "9297bb7bf19def6c33dd8295030de2b040d7dee3"
PREFIX = "tests/data/osteosarc/expansion/"
VARIANT = "DYNC1H1-chr14-101980529"
# Timepoints/modalities are explicit curated source-label/path claims, not
# inferred donor or cell identities. The verbatim upstream evidence is retained.
SOURCES = {
    "0066232879babe83": ("PacBio single-cell", "T1"),
    "53f498a544883d51": ("ONT single-cell", "T1"),
    "1120a096937e29e1": ("ONT single-cell", "T2"),
    "58c4d68e5f7e55b5": ("ONT single-cell", "T3"),
    "71f33f4e5748f7d1": ("short-read single-cell", "T1"),
    "bc138ca69b6b790e": ("short-read single-cell", "T2"),
    "8b8a9a02cbcbcf4b": ("short-read single-cell CD45neg", "T3"),
    "1b66c15da594a3ef": ("bulk RNA", "T0"),
    "50da1d13e05059fd": ("bulk RNA", "T2"),
}


def digest(path):
    checksum = sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(chunk)
    return checksum.hexdigest()


def original_records(source_bam, output, preferred_names):
    """Keep original records by deterministic read-name sampling, not sequence.

    Retain up to 32 names supporting the upstream default protein, plus 48
    hash-ordered regional names regardless of allele/quality. Preserve every
    regional alignment for each selected name, all tags and original qualities.
    No selection step consults the historical vaccine strings.
    """
    region = ("chr14", 101980527, 101980530)
    with pysam.AlignmentFile(source_bam) as source:
        names = {r.query_name for r in source.fetch(*region)}
        preferred = set(preferred_names[:32])
        if not preferred <= names:
            raise ValueError("Published supporting names missing from regional source")
        selected = preferred | set(sorted(names, key=lambda n: (sha256(n.encode()).hexdigest(), n))[:48])
        records = []
        with pysam.AlignmentFile(output, "wb", header=source.header) as target:
            for ordinal, read in enumerate(source.fetch(*region)):
                if read.query_name in selected:
                    target.write(read)
                    records.append(dict(regional_ordinal=ordinal, name=read.query_name,
                                        sam_sha256=sha256(read.to_string().encode("ascii")).hexdigest()))
        header_sha256 = sha256(str(source.header).encode()).hexdigest()
    pysam.index(str(output))
    return dict(region=list(region), coordinate_system="0-based half-open",
                header_sha256=header_sha256, selected_name_count=len(selected), records=records)


def main(repo, regional_cache, output):
    def read(path):
        return subprocess.check_output(["git", "-C", str(repo), "show", COMMIT + ":" + PREFIX + path])

    published = json.loads(read("audit/manifest.json"))
    matrix_bytes = read("audit/matrix.json.gz")
    index_bytes = read("audit/counts-index.json.gz")
    # Validate the named upstream artifacts as well as pinning the Git commit.
    for filename, data in [("matrix.json.gz", matrix_bytes), ("counts-index.json.gz", index_bytes)]:
        if sha256(data).hexdigest() != published["files"][filename]:
            raise ValueError("Released audit checksum mismatch: " + filename)
    matrix = json.loads(gzip.decompress(matrix_bytes))
    index = json.loads(gzip.decompress(index_bytes))
    rows = {r["source_id"]: r for r in matrix["rows"] if r["variant_id"] == VARIANT}
    sources = {s["source_id"]: s for s in matrix["sources"]}
    variant = next(v for v in matrix["variants"] if v["variant_id"] == VARIANT)
    output.mkdir(parents=True, exist_ok=False)
    cases, files = [], {}
    for sid, (modality, timepoint) in SOURCES.items():
        source, row = sources[sid], rows[sid]
        receipts = [r for r in source["acquisition_receipts"]
                    if r["status"] == "ok" and r["bam_sha256"] == row["source_bam_sha256"]]
        if len(receipts) != 1:
            raise ValueError("No unique released regional receipt: " + sid)
        receipt, = receipts
        bam = regional_cache / "alignments" / sid / "regions-GRCh38.bam"
        if digest(bam) != receipt["bam_sha256"]:
            raise ValueError("Cached BAM differs from released source: " + sid)
        protein = row["defaults"]["proteins"][0]
        names = [row["supporting_read_name_table"][i] for i in protein["supporting_read_name_indices"]]
        filename = sid + ".bam"
        selection = original_records(bam, output / filename, names)
        for name in (filename, filename + ".bai"):
            files[name] = digest(output / name)
        cases.append(dict(source_id=sid, modality=modality, timepoint=timepoint,
                          attribution=source["attribution"], source=source,
                          variant=variant, bam=filename, selection=selection,
                          upstream_default_protein=protein,
                          upstream_region_counts=row["defaults"]["counts"],
                          upstream_independent_primary=row["independent_primary"]))
        print(sid, modality, timepoint, len(selection["records"]), "original records", flush=True)
    # Every RNA product remains represented, including missing/failed inputs.
    inventory = dict(
        variant_ids=index["variant_ids"], variant_summaries=index["summary"],
        variants=matrix["variants"], source_ids=index["source_ids"],
        sources=[sources[sid] for sid in index["source_ids"]],
        rows=[r for r in index["rows"] if r["variant_id"] == VARIANT])
    data = gzip.compress((json.dumps(inventory, sort_keys=True, separators=(",", ":")) + "\n").encode(), mtime=0)
    (output / "inventory.json.gz").write_bytes(data)
    files["inventory.json.gz"] = sha256(data).hexdigest()
    manifest = dict(
        upstream_repo="https://github.com/openvax/isovar", upstream_commit=COMMIT,
        upstream_release="1.8.1", upstream_matrix_sha256=sha256(matrix_bytes).hexdigest(),
        upstream_counts_index_sha256=sha256(index_bytes).hexdigest(),
        variant_id=VARIANT, files=files, cases=cases,
        selection="32 published default-protein supporting names plus 48 SHA-256-ordered regional names; all original alignments for selected names",
        scope="Selected fixtures are not unbiased VAF, cell fractions or biological replicates. Timepoints/modalities retain source-label/path attribution; cell barcodes do not independently establish tumor-cell identity.",
        inventory_status_counts=dict(Counter(r["status"] for r in inventory["rows"])))
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repo", type=Path)
    parser.add_argument("regional_cache", type=Path)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "rna")
    args = parser.parse_args()
    main(args.repo, args.regional_cache, args.output)
