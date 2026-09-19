"""Consolidate acquisition parts without double-counting overlapping records."""

import argparse
from collections import Counter
import json
from pathlib import Path
import shutil
import subprocess

import pysam

from build import canonical_allele, digest, now, record_identity, summarize, write_json


def snapshot(root):
    """Keep required local input snapshots, not just ephemeral external paths."""
    provenance = json.loads((root / "provenance.json").read_text())
    copies = []
    for key in ("matrix", "listing"):
        source = Path(provenance[key])
        expected = provenance[key + "_sha256"]
        target = root / "sources" / source.name
        if not target.exists():
            if digest(source) != expected:
                raise ValueError("Preparation source changed: " + str(source))
            shutil.copyfile(source, target)
        if digest(target) != expected:
            raise ValueError("Pinned preparation input changed")
        copies.append(dict(path=str(target.relative_to(root)), sha256=expected, original_path=str(source)))
    old_corpus_path = Path(provenance["matrix"]).parent.parent / "corpus/manifest.json"
    corpus_snapshot = root / "sources/original-49-case-manifest.json"
    if not corpus_snapshot.exists():
        shutil.copyfile(old_corpus_path, corpus_snapshot)
    cases = json.loads(corpus_snapshot.read_text())["cases"]
    events = json.loads((root / "events.json").read_text())
    identities = {canonical_allele(e["chrom"], e["pos"], e["ref"], e["alt"]): e["event_id"]
                  for e in events if "ref" in e}
    crosswalk = []
    for case in cases:
        allele = case["variant"]
        original = allele.get("original_identity", allele)
        key = canonical_allele(original["chrom"], original["pos"], original["ref"], original["alt"])
        crosswalk.append(dict(case_id=case["case_id"], source_id=case["source_id"], event_id=identities[key],
                              native_allele={k: allele[k] for k in ("assembly", "chrom", "pos", "ref", "alt")}))
    if len(crosswalk) != 49 or len({c["event_id"] for c in crosswalk}) != 44:
        raise ValueError("Original case/locus denominator changed; review explicitly")
    write_json(root / "sources/original-case-crosswalk.json", crosswalk)
    generator = root / "generator"
    generator.mkdir(exist_ok=True)
    for source in sorted(Path(__file__).parent.iterdir()):
        if source.suffix in (".py", ".json", ".md"):
            target = generator / source.name
            shutil.copyfile(source, target)
            copies.append(dict(path=str(target.relative_to(root)), sha256=digest(target)))
    write_json(root / "reproduction.json", dict(created_utc=now(), files=copies,
        case_manifest_sha256=digest(corpus_snapshot), case_crosswalk_sha256=digest(root / "sources/original-case-crosswalk.json"),
        original_cases=len(cases), original_unique_loci=len({c["event_id"] for c in crosswalk}),
        samtools=subprocess.check_output(["samtools", "--version"]).decode("utf-8", errors="replace").splitlines()[:2],
        reproduction="Exact commands and queried BED intervals are retained per acquisition; downloaded subset BAMs and indexes have independent hashes"))


def combine_bams(first, second, output):
    """Multiset union: preserve true source multiplicity, not query overlap."""
    with pysam.AlignmentFile(second) as bam:
        header = bam.header.to_dict()
        remaining, originals = Counter(), {}
        for read in bam:
            key = record_identity(read)
            remaining[key] += 1
            originals[key] = read
    temporary = output.with_name("combined-unsorted.bam")
    with pysam.AlignmentFile(first) as bam:
        if bam.header.to_dict().get("SQ") != header.get("SQ"):
            raise ValueError("Cannot combine different reference dictionaries")
        if bam.header.to_dict().get("RG") != header.get("RG"):
            raise ValueError("Cannot combine different read groups")
        with pysam.AlignmentFile(temporary, "wb", header=bam.header) as target:
            for read in bam:
                target.write(read)
                key = record_identity(read)
                if remaining[key]:
                    remaining[key] -= 1
            for key, count in remaining.items():
                for _ in range(count):
                    target.write(originals[key])
    pysam.sort("--no-PG", "-o", str(output), str(temporary))
    pysam.index(str(output))
    temporary.unlink()


def finalize(root):
    snapshot(root)
    parts = [root, root / "phase-and-sv-context"]
    target = root / "dataset"
    sources = json.loads((root / "source_inventory.json").read_text())
    events, native = [], []
    for part in parts:
        events.extend(json.loads((part / "events.json").read_text()))
        native.extend(json.loads((part / "events-GRCh37.json").read_text()))
    if len({e["event_id"] for e in events}) != len(events):
        raise ValueError("Duplicate event identities between acquisition parts")
    write_json(target / "events.json", events)
    write_json(target / "events-GRCh37.json", native)
    write_json(target / "source_inventory.json", sources)
    selection_path = root / "source-selection.json"
    if selection_path.exists():
        if json.loads(selection_path.read_text()) != json.loads((parts[1] / "source-selection.json").read_text()):
            raise ValueError("Acquisition parts must use the same source selection")
        shutil.copyfile(selection_path, target / "source-selection.json")
    write_json(target / "provenance.json", dict(created_utc=now(), parts=[str(p.resolve()) for p in parts],
        generator_sha256=digest(__file__), base_provenance=json.loads((root / "provenance.json").read_text()),
        combination="Per-source multiset union of identical SAM records: maximum multiplicity across parts, preserving duplicates within each input",
        read_count_cap=None, species="Homo sapiens", taxon_id=9606))
    for source in sources:
        source_id = source["source_id"]
        receipts = []
        for part in parts:
            path = part / "alignments" / source_id / "receipt.json"
            receipts.append(json.loads(path.read_text()) if path.exists() else {"status": "not_acquired"})
        directory = target / "alignments" / source_id
        directory.mkdir(parents=True, exist_ok=True)
        result = dict(source_id=source_id, component_receipts=receipts)
        if not all(r["status"] == "ok" for r in receipts):
            statuses = {r["status"] for r in receipts}
            result["status"] = receipts[0]["status"] if len(statuses) == 1 else "incomplete_components"
        else:
            if len({r["assembly"] for r in receipts}) != 1:
                raise ValueError("Source assembly changed between queries")
            if len({r["header_sha256"] for r in receipts}) != 1 or len({r["source_index"]["sha256"] for r in receipts}) != 1:
                raise ValueError("Source header/index changed between queries")
            inputs = [p / "alignments" / source_id / "reads.bam" for p in parts]
            hashes = [digest(p) for p in inputs]
            if hashes != [r["bam_sha256"] for r in receipts]:
                raise ValueError("Input BAM hash changed")
            output = directory / "reads.bam"
            previous = directory / "receipt.json"
            old = json.loads(previous.read_text()) if previous.exists() else {}
            if not (old.get("input_hashes") == hashes and old.get("combination_schema") == 2 and output.exists() and digest(output) == old.get("bam_sha256")):
                combine_bams(*inputs, output)
            with pysam.AlignmentFile(output) as bam:
                count = sum(1 for _ in bam)
            result.update(status="ok", assembly=receipts[0]["assembly"], input_hashes=hashes, combination_schema=2,
                bam_sha256=digest(output), index_sha256=digest(str(output) + ".bai"),
                bam_bytes=output.stat().st_size, record_counts={"records": count},
                unmapped_event_ids=sorted({i for r in receipts for i in r.get("unmapped_event_ids", [])}),
                native_event_overrides={k: v for r in receipts for k, v in r.get("native_event_overrides", {}).items()})
        write_json(directory / "receipt.json", result)
    summarize(target)
    report(root, target)


def report(root, target):
    manifest = json.loads((target / "manifest.json").read_text())
    sources = manifest["sources_receipts"]
    events = json.loads((target / "events.json").read_text())
    coverage = json.loads((target / "coverage.json").read_text())
    kinds = Counter(e["kind"] for e in events)
    lines = ["# Expanded osteosarc read dataset", "", "Generated: " + now(), "",
        "Status: " + ("all acquisition attempts resolved; inspect explicit unavailable products below."
                       if manifest["acquisition_complete"] else "**INCOMPLETE — downloads are pending or require retry. Do not use this as a finished platform comparison.**"), "",
        "The earlier **49 cases represent 44 unique vaccine loci**, not 49 distinct mutations. "
        "All 44 remain in this dataset. Candidates and linked alleles are not a somatic truth set.", "",
        "## Inventory", "", f"{len(events)} event hypotheses; {len(sources)} processing products. "
        "Several products reprocess the same biological libraries and must not be pooled as independent samples.", "",
        "Kinds: " + ", ".join(f"{n} {k}" for k, n in sorted(kinds.items())) + ".", "",
        "Acquisition prioritizes one reviewed representative per library/preparation family, plus tagged long-read evidence and DNA controls. "
        "Alternatives remain inventoried and already acquired data are retained. Grouping is not a claim of independent biological replicates. "
        "Deferred products are not failed assays or zero coverage; unassigned CellRanger partitions remain outside the primary assigned-sample set.", "",
        "| Assay | Platform | Product timepoint claim | Acquired products | Selected products | Inventoried products |",
        "|---|---|---|---:|---:|---:|"]
    groups = sorted({(s["assay"], s["platform"], ", ".join(s["timepoints"])) for s in sources})
    for assay, platform, timepoint in groups:
        group = [s for s in sources if (s["assay"], s["platform"], ", ".join(s["timepoints"])) == (assay, platform, timepoint)]
        lines.append(f"| {assay} | {platform} | {timepoint} | {sum(s['status'] == 'ok' for s in group)} | {sum(s.get('acquisition_selection', {}).get('acquire', True) for s in group)} | {len(group)} |")
    lines.extend(["", "PacBio processing stages are one T1 library, not nine samples. "
        "No T0/T2/T3 PacBio or T0 ONT RNA source was identified in the inspected catalogue. "
        "No missing assay is counted as a failed ORF.", "",
        "## Acquisition", "", "Source dispositions: `" + json.dumps(manifest["source_statuses"], sort_keys=True) + "`.", "",
        f"Consolidated subset records: {sum(s.get('record_counts', {}).get('records', 0) for s in sources):,}. "
        "These are alignment records, not independent molecules or ALT-support counts.", "",
        "Acquisition requests query each selected source at every target in its native coordinate system: "
        "±2 kb around literal alleles and each structural-event endpoint, with paired-mate retrieval. "
        "Complete **stored** sequences/tags/flags and source record multiplicity are retained. "
        "There is no read-count, ALT, MAPQ, base-quality, duplicate or supplementary filter. "
        "This is regional acquisition, not whole-transcript or whole-genome sequencing: "
        "unplaced unmapped mates are not guaranteed; unmapped reads outside the PacBio audit, hard-clipped bases, and unqueried supplementary partners are not comprehensively recovered.", "",
        "GRCh37 queries use uniquely mapped ungapped chain intervals and independently verified reference alleles. "
        "MUC2 and the symbolic MUC3A anchor lack accepted GRCh37 mappings. Their native-GRCh37 coverage is unavailable, not zero. "
        "Both conflicting catalogue and VAF-table timepoint labels are retained; metadata is not silently reconciled.", "",
        "## PacBio input audit", ""])
    accounting_path = root / "pacbio_full/molecule-accounting.json"
    if accounting_path.exists():
        counts = json.loads(accounting_path.read_text())["counts"]
        lines.append(f"The complete deduplicated input contains {counts['input_molecules']:,} molecules; "
            f"{counts['input_molecules_represented_in_mapping']:,} appear in the published genomic BAM; "
            f"{counts['input_molecules_absent_from_mapping']:,} are absent. Exact molecule identities were checked, not just file counts.")
        lines.extend(["", "The absent molecules were aligned against the full GRCh38 reference with minimap2 splice:hq. "
            "Original cell/UMI/molecule tags and qualities were restored by exact identity. "
            "See `../pacbio_full/realignment.json` and `rescued-genome.bam`; these are derived alignments of the same T1 library, not new sequencing."])
    quality_path = root / "pacbio_full/quality-audit.json"
    if quality_path.exists():
        quality = json.loads(quality_path.read_text())
        counts = quality["counts"]
        lines.extend(["", f"**Missing QUAL is a separate software-compatibility limitation:** "
            f"{counts['missing_quality']:,}/{counts['records']:,} consolidated PacBio records lack base qualities. "
            f"The audit records Isovar {quality['installed_versions']['isovar']} and its exact collector-code hash. "
            f"Its recorded missing-QUAL acceptance policy is {quality.get('accepts_missing_base_qualities', False)}. "
            "A strict collector does not measure all available RNA evidence. "
            "See [Isovar #294](https://github.com/openvax/isovar/issues/294) and `../pacbio_full/quality-audit.json`. "
            "Do not interpret the earlier 3/44 translated-window result as a PacBio platform failure or attribute it solely to missing genomic coverage. "
            "No artificial quality scores have been assigned."])
        newer_path = root / "pacbio_full/quality-audit-isovar-1.18.1.json"
        if newer_path.exists():
            newer = json.loads(newer_path.read_text())
            if newer["source_bam_sha256"] != quality["source_bam_sha256"]:
                raise ValueError("Quality-policy comparison requires identical input reads")
            old_alt = sum(r["alt"] > 0 for r in quality["installed_collector_counts"])
            new_alt = sum(r["alt"] > 0 for r in newer["installed_collector_counts"])
            lines.extend(["", f"On the **identical** consolidated BAM, Isovar {quality['installed_versions']['isovar']} retains ALT reads at {old_alt}/44 original loci; "
                f"the isolated released Isovar {newer['installed_versions']['isovar']} retains ALT reads at {new_alt}/44. "
                "This measures allele-read collection, **not** established proteins/ORFs. "
                "[Isovar PR #298](https://github.com/openvax/isovar/pull/298) introduced missing-QUAL support in 1.18.0. "
                "The shared Python environment was not upgraded; both audits and exact code hashes are retained."])
            strict_path = root / "pacbio_full/quality-audit-isovar-1.18.1-strict.json"
            if strict_path.exists():
                strict = json.loads(strict_path.read_text())
                if strict["source_bam_sha256"] != newer["source_bam_sha256"] or strict["collector_code_sha256"] != newer["collector_code_sha256"]:
                    raise ValueError("Strict/adaptive comparison needs identical data and code")
                strict_alt = sum(r["alt"] > 0 for r in strict["installed_collector_counts"])
                lines.extend(["", f"Within Isovar 1.18.1 itself, requiring measured base qualities yields {strict_alt}/44 ALT-positive loci, "
                    f"versus {new_alt}/44 with missing qualities retained as unknown. This isolates the quality-availability policy from other release changes."])
    lines.extend(["", "## Expanded event coverage", "",
        "Counts below are the number of RNA **products** with at least one primary read having an aligned block "
        "overlapping an event coordinate. They are not ALT counts, established ORFs, exon-path confirmation or unique libraries. "
        "For indels this includes the reference anchor; it does not prove an indel read crosses both flanks.", "",
        "Each cell is covered/assessed products; **— means no usable acquisition**, not zero expression.", "",
        "| Event | Kind / origin status | ILMN | ONT | PacBio |", "|---|---|---:|---:|---:|"])
    metadata = {s["source_id"]: s for s in sources}
    for event in sorted(events, key=lambda e: e["event_id"]):
        assessed = [r for r in coverage if r["event_id"] == event["event_id"] and r["status"] == "ok"
                    and metadata[r["source_id"]]["assay"] == "RNA"]
        totals = Counter(metadata[r["source_id"]]["platform"] for r in assessed)
        counts = Counter(metadata[r["source_id"]]["platform"] for r in assessed if r["primary_fragment_names"] > 0)
        cells = [f"{counts[p]}/{totals[p]}" if totals[p] else "—" for p in ("ILMN", "ONT", "PacBio")]
        lines.append(f"| {event['event_id']} | {event['kind']} / {event['somatic_status']} | "
                     + " | ".join(cells) + " |")
    lines.extend(["", "## Files and downstream use", "",
        "- `events.json`: literal alleles or explicit unresolved breakpoint hypotheses, with selection reasons and source claims.",
        "- `source_inventory.json`: URLs, assay/platform/timepoint/tissue claims, processing identities and available indexes.",
        "- `source-selection.json`: explicit representative/evidence choices, deferred alternatives and grouping rationale.",
        "- `alignments/<source_id>/reads.bam` + `.bai`: consolidated per-source reads. Use these, **not** both acquisition parts concatenated.",
        "- `coverage.json`: measured coverage, per-endpoint counts, and explicit unavailable states.",
        "- `manifest.json`: hashes and per-source acquisition/combination receipts.",
        "- Parent `sources/` and `pacbio_full/`: pinned online tables and complete PacBio inputs/audits.", "",
        "Use Isovar for RNA-supported assembly/translation and Varcode for annotation; do not treat this inventory as a new vaccine ranking. "
        "The matched-normal DNA subsets enable checking linked-allele origin, but genotype confidence, somatic status, direct phasing, "
        "coding frame and exact translated sequence remain separate evidence gates.", ""])
    (target / "REPORT.md").write_text("\n".join(lines))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path)
    finalize(parser.parse_args().run)
