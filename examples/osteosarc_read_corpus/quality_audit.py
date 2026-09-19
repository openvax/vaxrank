"""Audit real PacBio QUAL availability and the installed Isovar read gate."""

import argparse
from collections import Counter
import importlib.metadata
import inspect
import json
import logging
from pathlib import Path

import pysam
from isovar import ReadCollector
from pyensembl import EnsemblRelease
from varcode import Variant

from build import digest, now, write_json


def audit(root, label="", require_qualities=False):
    logging.disable(logging.WARNING)
    dataset = root / "dataset"
    path = dataset / "alignments/0066232879babe83/reads.bam"
    counts = Counter()
    rows = []
    collector = ReadCollector(**({"use_reads_without_base_qualities": False} if require_qualities else {}))
    genome = EnsemblRelease(87)
    with pysam.AlignmentFile(path) as bam:
        for read in bam:
            counts["records"] += 1
            counts["missing_quality"] += read.query_qualities is None
            counts["missing_sequence"] += read.query_sequence is None
            counts["primary"] += not (read.is_secondary or read.is_supplementary or read.is_unmapped)
        for event in json.loads((dataset / "events.json").read_text()):
            if "original_44_loci" not in event.get("reasons", []):
                continue
            contig = "MT" if event["chrom"] == "chrM" else event["chrom"].removeprefix("chr")
            variant = Variant(contig, event["pos"], event["ref"], event["alt"], ensembl=genome)
            evidence = collector.read_evidence_for_variant(variant, bam)
            rows.append(dict(event_id=event["event_id"], ref=len(evidence.ref_reads),
                             alt=len(evidence.alt_reads), other=len(evidence.other_reads)))
    output = dict(created_utc=now(), source_bam=str(path), source_bam_sha256=digest(path),
        counts=dict(counts), installed_collector_counts=rows,
        installed_versions={p: importlib.metadata.version(p) for p in ("isovar", "varcode", "pysam")},
        collector_code_sha256=digest(inspect.getfile(ReadCollector)),
        collector_defaults=str(inspect.signature(ReadCollector)),
        accepts_missing_base_qualities=getattr(collector, "use_reads_without_base_qualities", False),
        interpretation="Missing QUAL is not low quality and has not been imputed. Counts describe the installed collector, not a platform capability or full-ORF benchmark.",
        upstream_issue="https://github.com/openvax/isovar/issues/294")
    write_json(root / ("pacbio_full/quality-audit%s.json" % ("-" + label if label else "")), output)
    print(json.dumps(dict(counts=counts, loci_with_retained_alt=sum(r["alt"] > 0 for r in rows))), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--label", default="")
    parser.add_argument("--require-base-qualities", action="store_true")
    args = parser.parse_args()
    audit(args.run, args.label, args.require_base_qualities)
