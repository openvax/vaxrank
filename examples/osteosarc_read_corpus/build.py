"""Acquire uncapped regional osteosarc RNA/DNA records with explicit provenance.

This is an analysis corpus, not a somatic truth set or a small test fixture.
Run prepare, acquire and summarize in order; successful downloads are resumable.
"""

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import importlib.util
import shutil
import subprocess
import struct
from urllib.parse import quote, urlsplit

from osteosarc import (
    Asset, Cache, IntegrityError, OsteosarcError, Region, __version__ as osteosarc_version,
    digest, extract_reads, inspect_alignment,
)
import pysam


BUCKET = "https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/"
HERE = Path(__file__).resolve().parent


def now():
    return datetime.now(timezone.utc).isoformat()


def record_identity(read):
    """SAM identity augmented with float bits lost in text formatting."""
    floats = []
    for tag, value, kind in read.get_tags(with_value_type=True):
        if kind in ("f", "d"):
            floats.append((tag, kind, struct.pack("<" + kind, value)))
        elif kind == "B" and value.typecode in ("f", "d"):
            floats.append((tag, kind, value.typecode, value.tobytes()))
    return read.to_string(), tuple(floats)


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def download(url, path, timeout=900, *, cache=None):
    """Export an osteosarc cache object, importing verified older downloads."""
    cache = cache if cache is not None else Cache(timeout=timeout)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path = path.with_name(path.name + ".receipt.json")
    if receipt_path.exists():
        receipt = json.loads(receipt_path.read_text())
        if receipt["url"] != url or digest(path) != receipt["sha256"]:
            raise ValueError("Cached download identity mismatch: %s" % path)
        cache.import_file(path, url, sha256=receipt["sha256"], size=receipt["bytes"])
        return receipt
    downloaded = cache.fetch(url)
    partial = path.with_name(path.name + ".download")
    shutil.copyfile(cache.path(downloaded), partial)
    partial.replace(path)
    receipt = dict(url=url, retrieved_utc=downloaded.retrieved_at, bytes=downloaded.size,
                   sha256=downloaded.sha256, osteosarc_version=osteosarc_version,
                   osteosarc_receipt=downloaded.to_dict())
    write_json(receipt_path, receipt)
    return receipt


def canonical_allele(chrom, pos, ref, alt):
    """Trim shared sequence for identity, without conflating nearby haplotypes."""
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while ref and alt and ref[0] == alt[0]:
        pos, ref, alt = pos + 1, ref[1:], alt[1:]
    return chrom, pos, ref, alt


def build_events(matrix, vafs_path, calls_path, reference_path):
    events = {}

    def add(row, reason, evidence):
        chrom, pos, ref, alt = row["chrom"], int(row["pos"]), row["ref"], row["alt"]
        if not ref or not alt or not re.fullmatch("[ACGT]+", ref + alt):
            return
        key = canonical_allele(chrom, pos, ref, alt)
        if key not in events:
            gene = row["gene"]
            events[key] = dict(event_id="%s-%s-%s-%s-%s" % (gene, chrom, pos, ref, alt),
                               gene=gene, chrom=chrom, pos=pos, ref=ref, alt=alt,
                               kind="SNV" if len(ref) == len(alt) == 1 else
                               "MNV" if len(ref) == len(alt) else "indel",
                               assembly="GRCh38", reasons=[], source_claims=[],
                               somatic_status="source_reported_candidate")
        event = events[key]
        if reason not in event["reasons"]:
            event["reasons"].append(reason)
        if evidence not in event["source_claims"]:
            event["source_claims"].append(evidence)

    for row in matrix["variants"]:
        add(row, "original_44_loci", row)
    # All site-listed literal alleles are cheap to retain once the same libraries
    # are being queried. This avoids selecting extra candidates by RNA outcome.
    seen = set()
    with vafs_path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["on_variants_page"] != "true" or row["variant_id"] in seen:
                continue
            seen.add(row["variant_id"])
            add(row, "site_variant_catalogue", {k: row[k] for k in (
                "variant_id", "gene", "chrom", "pos", "ref", "alt", "protein_change")})
    with calls_path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if "organoid" in row["pair"].lower():
                continue
            if row["coding_effect"] not in ("MISSENSE", "NONSENSE_OR_FRAMESHIFT", "SPLICE"):
                continue
            try:
                eligible = (float(row["tumor_dp"]) >= 10 and
                            float(row["normal_dp"]) >= 10 and
                            float(row["normal_vaf"]) <= 0.02 and
                            float(row["tumor_vaf"]) * float(row["tumor_dp"]) >= 3)
            except ValueError:
                eligible = False
            if eligible:
                add(row, "coding_tumor_DNA_candidate", row)
    # The adjacent NTF3 change and the compound allele are distinct queries.
    add(dict(gene="NTF3", chrom="chr12", pos=5494382, ref="G", alt="T"),
        "compound_haplotype_component", {"url": "https://osteosarc.com/variant/NTF3-chr12-5494381/"})
    add(dict(gene="NTF3", chrom="chr12", pos=5494381, ref="AG", alt="GT"),
        "compound_haplotype", {"url": "https://osteosarc.com/variant/NTF3-chr12-5494381/"})
    result = list(events.values())
    with pysam.FastaFile(str(reference_path)) as reference:
        for event in result:
            observed = reference.fetch(event["chrom"], event["pos"] - 1,
                                       event["pos"] - 1 + len(event["ref"])).upper()
            if observed != event["ref"]:
                raise ValueError("Reference mismatch for %s: %s" % (event["event_id"], observed))
            event["reference_verified"] = True
    audit = json.loads((HERE.parent / "osteosarc_complex_results/source/additional_candidate_audit.json").read_text())
    for row in audit["larger_unresolved_leads"]:
        match = re.fullmatch(r"(chr\w+):(\d+)-(\d+) deletion", row["event"])
        chrom, start, end = match.groups()
        result.append(dict(event_id=row["id"], gene=row["id"].split("-")[0],
                           kind="large_deletion", assembly="GRCh38",
                           breakpoints=[[chrom, int(start)], [chrom, int(end)]],
                           source_claims=[row], reasons=["previous_SV_investigation"],
                           somatic_status="source_reported_candidate"))
    for row in audit["structural_variants"]:
        result.append(dict(event_id=row["id"], gene=row["id"], kind="rearrangement",
                           assembly="GRCh38", breakpoints=[
                               [p.split(":")[0], int(p.split(":")[1])] for p in row["breakpoints"]],
                           source_claims=[row], reasons=["previous_SV_investigation"],
                           somatic_status="source_reported_candidate"))
    return sorted(result, key=lambda x: x["event_id"])


def metadata(source, vaf_metadata):
    key, name = source["key"], source["name"]
    basename = Path(key).name
    claims = vaf_metadata.get(basename, [])
    if not claims and basename.endswith(".out.bam"):
        claims = vaf_metadata.get(basename.removesuffix(".bam") + ".md.bam", [])
    # Vendor/reprocessed aliases retain their library's published metadata.
    if not claims:
        identifiers = [token for token in ("BG003082", "BG009368", "SARC0277", "TL-24-ALMY2X4KMV", "TL-24-KCVBE1UI1P") if token in key]
        for alias, values in vaf_metadata.items():
            if any(token in alias for token in identifiers):
                claims = [*claims, *[v for v in values if v not in claims]]
    timepoints = sorted({r["timepoint"] for r in claims if r["timepoint"]})
    path_points = sorted(set(re.findall(r"(?:^|[/_ ])(T[0-3])(?:[/_ .]|$)", key + " " + name)))
    if not timepoints:
        timepoints = path_points
    date_claims = []
    for prefix, point in (("genomics/genomics-bulk/2022.12.16/", "T0"),
                          ("genomics/genomics-bulk/2024.06.06/", "T1"),
                          ("genomics/genomics-bulk/2025.01.06/", "T2"),
                          ("vendor/cegat/P116686_2_S000048/", "T0")):
        if key.startswith(prefix):
            date_claims.append(point)
    if not timepoints:
        timepoints = date_claims
    library_group = None
    for token in ("BG003082", "BG009368", "SARC0277", "TL-24-ALMY2X4KMV", "TL-24-KCVBE1UI1P"):
        if token in key:
            library_group = token
    if "sclrs_ONT" in key:
        library_group = "UCSF-ONT-" + "-".join(timepoints)
    elif "pacbio" in key.lower():
        library_group = "UCSF-PacBio-T1"
    text = (key + " " + name).lower()
    platform = "PacBio" if "pacbio" in text else "ONT" if "/ont/" in ("/" + text) or "sclrs_ont" in text else "ILMN"
    assay = source.get("category", "RNA")
    if assay not in ("WGS", "WES"):
        assay = "RNA"
    return dict(platform=platform, assay=assay, timepoints=timepoints or ["unresolved"],
                catalogue_timepoint_claims=path_points, path_date_timepoint_claims=date_claims, sample_metadata_claims=claims,
                timepoint_conflict=bool(timepoints and path_points and {p.split("-")[0] for p in timepoints} != set(path_points)),
                tissue=source.get("tissue", "unresolved"),
                source_library_group=library_group,
                library_group_basis="explicit path sample identifier" if library_group else "unresolved; retain separate processing products",
                independence="processing_product; do not count as independent biological replicate")


def prepare(root, isovar_repo, listing_path, reference_path, *, cache=None):
    snapshots = root / "sources"
    snapshots.mkdir(parents=True, exist_ok=True)
    for filename, url in {
        "bams.json": "https://osteosarc.com/bams/bams.json",
        "variant_vafs_long.tsv": "https://osteosarc.com/variants/variant_vafs_long.tsv",
        "snv_top.tsv": "https://osteosarc.com/oncoanalyser/tables/snv_top.tsv",
    }.items():
        download(url, snapshots / filename, cache=cache)
    matrix_path = isovar_repo / "tests/data/osteosarc/expansion/audit/matrix.json.gz"
    with gzip.open(matrix_path, "rt") as handle:
        matrix = json.load(handle)
    listing = json.loads(listing_path.read_text())
    files = {r[0]: r[1] for r in listing["files"]}
    catalogue = json.loads((snapshots / "bams.json").read_text())
    vaf_metadata = defaultdict(list)
    with (snapshots / "variant_vafs_long.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            claim = {k: row[k] for k in ("sample_label", "data_source", "pipeline", "assay_type", "timepoint", "tissue", "sample_date")}
            if claim not in vaf_metadata[row["bam_file"]]:
                vaf_metadata[row["bam_file"]].append(claim)
    selected = {s["key"]: dict(s) for s in matrix["sources"]
                if s["scope"] == "rna_candidate" and s["tissue"] != "blood"}
    for category in catalogue["categories"]:
        if category["name"] == "Blood scRNA":
            continue
        for row in category["bams"]:
            key = row["url"].removeprefix(BUCKET)
            source = selected.setdefault(key, dict(key=key, source_id=hashlib.sha256((BUCKET + key).encode()).hexdigest()[:16],
                                                   url=BUCKET + quote(key, safe="/"), name=row["name"], tissue=row.get("tissue", "unresolved")))
            source["category"] = category["name"]
    for source in selected.values():
        source["listed_indexes"] = [k for k in (source["key"] + ".bai", source["key"][:-4] + ".bai", source["key"] + ".csi") if k in files]
        source["metadata"] = metadata(source, vaf_metadata)
        source["bucket_bytes"] = files.get(source["key"])
    events = build_events(matrix, snapshots / "variant_vafs_long.tsv", snapshots / "snv_top.tsv", reference_path)
    write_json(root / "events.json", events)
    write_json(root / "source_inventory.json", sorted(selected.values(), key=lambda s: (s["metadata"]["assay"], s["metadata"]["platform"], s["key"])))
    write_json(root / "provenance.json", dict(created_utc=now(), matrix=str(matrix_path), matrix_sha256=digest(matrix_path),
        listing=str(listing_path), listing_sha256=digest(listing_path), listing_generated_at=listing.get("generated_at"),
        reference=str(reference_path), reference_sha256=digest(reference_path), generator_sha256=digest(__file__),
        baseline_unique_events=44, baseline_source_variant_cases=49,
        scope="Tumor RNA across all catalogued clinical time points, plus DNA controls; pooled blood RNA excluded because donor identity is unresolved",
        selection="All site-listed literal alleles; additional coding patient calls with tumor/normal depth >=10, normal VAF <=0.02 and estimated tumor ALT reads >=3; known SV leads. No claim that all are proven somatic."))
    print("Prepared", len(events), "events and", len(selected), "source products", flush=True)


def intervals(events, header, padding):
    lengths = {s["SN"]: s["LN"] for s in header["SQ"]}
    regions = []
    missing = []
    for event in events:
        points = event.get("breakpoints", [[event.get("chrom"), event.get("pos")]])
        for chrom, pos in points:
            names = [n for n in lengths if n.removeprefix("chr").replace("MT", "M") == chrom.removeprefix("chr").replace("MT", "M")]
            if len(names) != 1:
                missing.append(event["event_id"])
                continue
            end = pos - 1 + max(1, len(event.get("ref", "")))
            regions.append((names[0], max(0, pos - 1 - padding), min(lengths[names[0]], end + padding)))
    return sorted(set(regions)), sorted(set(missing))


def lift_events(root, isovar_repo, reference_path, *, cache=None):
    """Retain independently validated GRCh37 alleles for native vendor reads."""
    module_path = isovar_repo / "tests/data/osteosarc/expansion/liftover.py"
    spec = importlib.util.spec_from_file_location("osteosarc_coordinate_mapping", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    chain = root / "sources/hg38ToHg19.over.chain.gz"
    receipt = download(module.CHAIN_URL, chain, cache=cache)
    events = json.loads((root / "events.json").read_text())
    with gzip.open(chain, "rt") as handle:
        blocks = list(module.chain_blocks(handle))
    by_contig = defaultdict(list)
    for block in blocks:
        by_contig[block["source_contig"]].append(block)
    mapped, failed = [], []
    with pysam.FastaFile(str(reference_path)) as reference:
        for event in events:
            try:
                if event.get("chrom") == "chrM":
                    # human_g1k_v37 uses rCRS MT, as do native Ensembl/Tempus.
                    lifted = dict(event, assembly="GRCh37", coordinate_mapping="rCRS MT retained")
                elif "breakpoints" in event:
                    points, mappings = [], []
                    for chrom, pos in event["breakpoints"]:
                        item = module.map_interval(by_contig[chrom], chrom, pos - 1, pos)
                        points.append([item["chrom"], item["start"] + 1])
                        mappings.append(item)
                    lifted = dict(event, assembly="GRCh37", breakpoints=points, coordinate_mapping=mappings)
                else:
                    lifted = module.lift_variant(event, by_contig[event["chrom"]])
                if "ref" in lifted:
                    chrom = lifted["chrom"].removeprefix("chr")
                    chrom = "MT" if chrom == "M" else chrom
                    observed = reference.fetch(chrom, lifted["pos"] - 1,
                                               lifted["pos"] - 1 + len(lifted["ref"])).upper()
                    if observed != lifted["ref"]:
                        raise ValueError("GRCh37 reference mismatch: " + observed)
                    lifted["native_reference_verified"] = True
                mapped.append(lifted)
            except (ValueError, KeyError) as error:
                failed.append(dict(event_id=event["event_id"], reason=str(error)))
    write_json(root / "events-GRCh37.json", mapped)
    write_json(root / "coordinate-mapping.json", dict(chain=receipt, reference=str(reference_path),
        reference_sha256=digest(reference_path), mapping_code=str(module_path), mapping_code_sha256=digest(module_path),
        events_sha256=digest(root / "events.json"), mapped=len(mapped), failures=failed))
    print("GRCh37 coordinate mapping:", len(mapped), "verified events;", len(failed), "unresolved", flush=True)


def hg19_mito(root, isovar_repo, *, cache=None):
    """Pin the distinct UCSC hg19 mitochondrial sequence, not rCRS MT."""
    module_path = isovar_repo / "tests/data/osteosarc/expansion/liftover.py"
    spec = importlib.util.spec_from_file_location("osteosarc_mito_mapping", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    source = root / "sources/hg19-chrM.json"
    receipt = download("https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chrM", source, cache=cache)
    sequence = json.loads(source.read_text())["dna"].upper()
    if len(sequence) != 16571:
        raise ValueError("Unexpected UCSC hg19 mitochondrial reference length")
    with gzip.open(root / "sources/hg38ToHg19.over.chain.gz", "rt") as handle:
        blocks = list(module.chain_blocks(handle, contigs={"chrM"}))
    mapped = []
    for event in json.loads((root / "events.json").read_text()):
        if event.get("chrom") != "chrM":
            continue
        lifted = module.lift_variant(event, blocks)
        if sequence[lifted["pos"] - 1:lifted["pos"] - 1 + len(lifted["ref"])] != lifted["ref"]:
            raise ValueError("UCSC hg19 mitochondrial allele mismatch")
        lifted.update(native_reference_verified=True, mitochondrial_reference="UCSC hg19 chrM (16571 bp)",
                      reference_receipt=receipt)
        mapped.append(lifted)
    write_json(root / "events-hg19-MT.json", mapped)


def acquire_one(source, events, root, padding, *, cache=None):
    """Acquire through osteosarc while retaining this run's reviewed event panel."""
    cache = cache if cache is not None else Cache()
    directory = root / "alignments" / source["source_id"]
    directory.mkdir(parents=True, exist_ok=True)
    receipt_path = directory / "receipt.json"
    index_urls = tuple(key if urlsplit(key).scheme or Path(key).is_absolute()
                       else BUCKET + quote(key, safe="/") for key in source["listed_indexes"])
    asset = Asset(source["source_id"], source.get("key", source["url"]), source["url"],
                  "alignment", Path(urlsplit(source["url"]).path).suffix.lstrip("."),
                  size=source.get("bucket_bytes"), index_urls=index_urls)
    snapshot = root / "provenance.json"
    snapshot_id = digest(snapshot if snapshot.exists() else root / "events.json")
    request = dict(schema_version=2, osteosarc_version=osteosarc_version,
                   source_url=source["url"], source_size=asset.size, index_urls=list(index_urls),
                   events_sha256=digest(root / "events.json"), padding_bp=padding,
                   snapshot_id=snapshot_id, fetch_pairs=True,
                   platform=source.get("metadata", {}).get("platform"),
                   coordinate_files={name: digest(root / name) for name in (
                       "events-GRCh37.json", "events-hg19-MT.json", "coordinate-mapping.json")
                       if (root / name).exists()})
    old = None
    if receipt_path.exists():
        old = json.loads(receipt_path.read_text())
        if old["request"] != request:
            raise ValueError("Existing acquisition request changed; choose a new run directory")
        if old["status"] == "integrity_error":
            return old  # preserve damaged evidence; a retry must not silently repair it
    result = dict(source_id=source["source_id"], request=request, started_utc=now())

    def save_result():
        if receipt_path.exists():
            archive = directory / ("receipt-attempt-%s.json" % now().replace(":", ""))
            receipt_path.rename(archive)
        write_json(receipt_path, result)
        return result

    try:
        info = inspect_alignment(asset, cache=cache, snapshot_id=snapshot_id, timeout=120)
        header = info.header
        header_path = directory / "header.sam"
        result.update(assembly=info.assembly, header_sha256=digest(info.path),
                      header_receipt=info.receipt, programs=header.get("PG", []))
        if result["assembly"] == "GRCh37" and (root / "events-GRCh37.json").exists():
            native_events = json.loads((root / "events-GRCh37.json").read_text())
            result["unmapped_event_ids"] = sorted({e["event_id"] for e in events} - {e["event_id"] for e in native_events})
            events = native_events
            result["coordinate_mapping_sha256"] = digest(root / "coordinate-mapping.json")
            mito_lengths = {s["LN"] for s in header["SQ"] if s["SN"] in ("M", "MT", "chrM", "chrMT")}
            if mito_lengths == {16571} and (root / "events-hg19-MT.json").exists():
                overrides = {e["event_id"]: e for e in json.loads((root / "events-hg19-MT.json").read_text())}
                events = [overrides.get(e["event_id"], e) for e in events]
                result["native_event_overrides"] = overrides
            elif mito_lengths and mito_lengths != {16569}:
                result["unmapped_event_ids"].extend(e["event_id"] for e in events if e.get("chrom") == "chrM")
                events = [e for e in events if e.get("chrom") != "chrM"]
        elif result["assembly"] != "GRCh38":
            result["status"] = "requires_coordinate_mapping" if result["assembly"] == "GRCh37" else "no_genomic_coordinates"
            return save_result()
        if not source["listed_indexes"]:
            result["status"] = "missing_genomic_index"
            return save_result()
        regions, missing = intervals(events, header, padding)
        if not regions:
            result.update(status="no_mapped_regions", missing_contigs_for=missing)
            return save_result()
        lengths = {row["SN"]: row["LN"] for row in header["SQ"]}
        requested = [Region(contig, start, end, info.assembly, lengths[contig])
                     for contig, start, end in regions]
        unpaired_platform = source.get("metadata", {}).get("platform") in ("ONT", "PacBio")
        subset = extract_reads(asset, requested, cache=cache, fetch_pairs=not unpaired_platform,
                               snapshot_id=snapshot_id, timeout=7200)
        if unpaired_platform:
            with subset.open() as bam:
                unexpectedly_paired = any(read.is_paired for read in bam)
            if unexpectedly_paired:
                subset = extract_reads(asset, requested, cache=cache, fetch_pairs=True,
                                       snapshot_id=snapshot_id, timeout=7200)
        if old and old["status"] == "ok":
            for name, key in (("reads.bam", "bam_sha256"), ("reads.bam.bai", "index_sha256"),
                              ("header.sam", "header_sha256")):
                path = directory / name
                if not path.is_file() or digest(path) != old[key]:
                    raise IntegrityError("Existing acquisition was modified: " + name)
            if subset.receipt != old["osteosarc_receipt"]:
                raise IntegrityError("Source extraction changed; choose a new run directory")
            return old
        output = directory / "reads.bam"
        for original, target in ((subset.path, output),
                                 (subset.index_path, directory / "reads.bam.bai"),
                                 (info.path, header_path),
                                 (subset.path.parent / "regions.bed", directory / "regions.bed")):
            temporary = target.with_name(target.name + ".partial")
            shutil.copyfile(original, temporary)
            temporary.replace(target)
        counts = Counter(records=0)
        with pysam.AlignmentFile(output) as bam:
            for read in bam:
                counts["records"] += 1
                counts["unmapped"] += read.is_unmapped
                counts["secondary"] += read.is_secondary
                counts["supplementary"] += read.is_supplementary
                counts["missing_sequence"] += read.query_sequence is None
                counts["missing_quality"] += read.query_qualities is None
        result.update(status="ok", completed_utc=now(), regions=regions, missing_contigs_for=missing,
                      record_counts=dict(counts), bam_sha256=digest(output), bam_bytes=output.stat().st_size,
                      index_sha256=digest(str(output) + ".bai"), read_count_cap=None,
                      osteosarc_receipt=subset.receipt, command=subset.receipt["command"],
                      source_index=subset.receipt["index_receipt"] or dict(
                          path=subset.receipt["request"]["index"],
                          sha256=subset.receipt["request"]["index_sha256"]),
                      preserved="All retrieved SAM fields; complete stored sequences; no allele, quality, duplicate or secondary/supplementary filter")
    except (OSError, ValueError, OsteosarcError, subprocess.SubprocessError) as error:
        result.update(status="integrity_error" if isinstance(error, IntegrityError) else "acquisition_error",
                      error=str(error))
        if isinstance(error, subprocess.SubprocessError):
            stderr = getattr(error, "stderr", None)
            if stderr:
                (directory / "acquisition.stderr.log").write_bytes(
                    stderr if isinstance(stderr, bytes) else stderr.encode())
    return save_result()


def source_selection(root, sources):
    """Load an explicit acquisition policy without dropping inventoried sources."""
    path = root / "source-selection.json"
    if not path.exists():
        return {}
    policy = json.loads(path.read_text())
    entries = policy["sources"]
    ids = [s["source_id"] for s in entries]
    if len(ids) != len(set(ids)) or set(ids) != {s["source_id"] for s in sources}:
        raise ValueError("Source selection must account for every inventoried product exactly once")
    if any(type(s["acquire"]) is not bool or not s.get("reason") for s in entries):
        raise ValueError("Source selection needs explicit boolean decisions and reasons")
    return {s["source_id"]: s for s in entries}


def acquisition_sources(root, sources, source_ids=None):
    selection = source_selection(root, sources)
    if source_ids:
        unknown = set(source_ids) - {s["source_id"] for s in sources}
        if unknown:
            raise ValueError("Unknown source IDs: " + ", ".join(sorted(unknown)))
        sources = [s for s in sources if s["source_id"] in source_ids]
    return [s for s in sources if selection.get(s["source_id"], {}).get("acquire", True)]


def acquire(root, workers, padding, source_ids=None, *, cache=None):
    sources = acquisition_sources(root, json.loads((root / "source_inventory.json").read_text()), source_ids)
    events = json.loads((root / "events.json").read_text())
    # Interleave assay/platform/timepoint/tissue strata, so alternative products
    # for one stratum do not delay first coverage of another clinical timepoint.
    groups = defaultdict(list)
    for source in sources:
        m = source["metadata"]
        groups[(m["assay"], m["platform"], tuple(m["timepoints"]), m["tissue"])].append(source)
    ordered = []
    for key, group in sorted(groups.items()):
        group.sort(key=lambda s: (s["name"] == s["key"], "unassigned" in s["key"], s["key"]))
        ordered.extend((i, key, source) for i, source in enumerate(group))
    sources = [source for _, _, source in sorted(ordered, key=lambda x: (x[0], x[1]))]
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(acquire_one, s, events, root, padding, cache=cache): s for s in sources}
        for future in as_completed(futures):
            s = futures[future]
            result = future.result()
            print(s["source_id"], result["status"], result.get("record_counts", {}).get("records"), s["name"], flush=True)


def supplement(root, reference_path):
    """Pin additional phased companions and unresolved structural hypotheses."""
    child = root / "phase-and-sv-context"
    source = HERE / "context_events.json"
    events = json.loads(source.read_text())
    with pysam.FastaFile(str(reference_path)) as reference:
        for event in events:
            if "ref" in event:
                observed = reference.fetch(event["chrom"], event["pos"] - 1,
                                           event["pos"] - 1 + len(event["ref"])).upper()
                if observed != event["ref"]:
                    raise ValueError("Companion reference mismatch: " + event["event_id"])
                event["reference_verified"] = True
            event["assembly"] = "GRCh38"
    write_json(child / "events.json", events)
    write_json(child / "source_inventory.json", json.loads((root / "source_inventory.json").read_text()))
    write_json(child / "provenance.json", dict(created_utc=now(), parent=str(root.resolve()),
        source=str(source), source_sha256=digest(source),
        interpretation="RNA-linked companion alleles have unresolved germline/somatic origin; breakpoint hypotheses do not establish a translated fusion"))
    print("Prepared context supplement:", len(events), "events", flush=True)


def pacbio_full(root, *, cache=None):
    sources = json.loads((root / "source_inventory.json").read_text())
    selected = [s for s in sources if s["key"] == "pacbio/IPISRC044_T1_sclrs_live_pbmm2_mapped.bam" or s["key"].endswith(".corr.sort.dedup.bam")]
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = {executor.submit(download, s["url"], root / "pacbio_full" / Path(s["key"]).name, 3600, cache=cache): s for s in selected}
        for future in as_completed(futures):
            receipt = future.result()
            print("PacBio full file", futures[future]["source_id"], receipt["bytes"], flush=True)


def summarize(root):
    sources = json.loads((root / "source_inventory.json").read_text())
    selection = source_selection(root, sources)
    events = json.loads((root / "events.json").read_text())
    rows, receipts = [], []
    for source in sources:
        directory = root / "alignments" / source["source_id"]
        receipt_path = directory / "receipt.json"
        receipt = json.loads(receipt_path.read_text()) if receipt_path.exists() else {"status": "not_acquired"}
        decision = selection.get(source["source_id"])
        if decision:
            receipt = dict(receipt, acquisition_selection=decision)
            if not decision["acquire"] and receipt["status"] != "ok":
                receipt = dict(receipt, observed_acquisition_status=receipt["status"], status="deferred_by_policy")
        receipts.append({"source_id": source["source_id"], **source["metadata"], **receipt})
        if receipt["status"] != "ok":
            for event in events:
                rows.append(dict(source_id=source["source_id"], event_id=event["event_id"], status=receipt["status"], records=None))
            continue
        with pysam.AlignmentFile(directory / "reads.bam") as bam:
            native_by_id = {e["event_id"]: e for e in json.loads((root / "events-GRCh37.json").read_text())} if receipt["assembly"] == "GRCh37" else {}
            for event in events:
                if receipt["assembly"] == "GRCh37":
                    if event["event_id"] in receipt.get("unmapped_event_ids", []):
                        rows.append(dict(source_id=source["source_id"], event_id=event["event_id"], status="unresolved_coordinate_mapping", records=None))
                        continue
                    event = receipt.get("native_event_overrides", {}).get(event["event_id"], native_by_id[event["event_id"]])
                regs, missing = intervals([event], bam.header.to_dict(), 0)
                record_counts, aligned_counts = Counter(), Counter()
                fragments = set()
                endpoints = []
                for chrom, start, end in regs:
                    endpoint_records, endpoint_aligned = Counter(), Counter()
                    for read in bam.fetch(chrom, start, end):
                        identity = record_identity(read)
                        endpoint_records[identity] += 1
                        # Count aligned bases across the locus, not intron-spanning
                        # bounding boxes. Each breakpoint is reported separately.
                        aligned = any(a < end and b > start for a, b in read.get_blocks())
                        if aligned:
                            endpoint_aligned[identity] += 1
                        if aligned and not read.is_secondary and not read.is_supplementary:
                            fragments.add((read.get_tag("RG") if read.has_tag("RG") else "", read.query_name))
                    # Counter union takes the maximum multiplicity, retaining
                    # original duplicate records but not double counting a read
                    # retrieved at both structural-event endpoints.
                    record_counts |= endpoint_records
                    aligned_counts |= endpoint_aligned
                    endpoints.append(dict(contig=chrom, start0=start, end0=end,
                                          records=sum(endpoint_records.values()),
                                          aligned_block_overlap=sum(endpoint_aligned.values())))
                rows.append(dict(source_id=source["source_id"], event_id=event["event_id"], status="missing_contig" if missing else "ok",
                                 records=sum(record_counts.values()), aligned_block_overlap=sum(aligned_counts.values()),
                                 primary_fragment_names=len(fragments), endpoints=endpoints))
    write_json(root / "coverage.json", rows)
    write_json(root / "manifest.json", dict(created_utc=now(), events=len(events), sources=len(sources),
        acquisition_complete=all(r["status"] in ("ok", "missing_genomic_index", "no_genomic_coordinates", "deferred_by_policy") for r in receipts),
        selected_acquisition_complete=all(r["status"] == "ok" for r in receipts if r.get("acquisition_selection", {}).get("acquire", True)) if selection else None,
        source_statuses=dict(Counter(r["status"] for r in receipts)), sources_receipts=receipts,
        files={name: digest(root / name) for name in ("events.json", "events-GRCh37.json", "events-hg19-MT.json", "source_inventory.json", "source-selection.json", "provenance.json", "coverage.json") if (root / name).exists()},
        interpretation="Coverage is not ALT support or a proven ORF; duplicate processing products are not independent samples"))
    print(json.dumps(Counter(r["status"] for r in receipts)), flush=True)


def refresh_metadata(root):
    metadata_by_name = defaultdict(list)
    with (root / "sources/variant_vafs_long.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            claim = {k: row[k] for k in ("sample_label", "data_source", "pipeline", "assay_type", "timepoint", "tissue", "sample_date")}
            if claim not in metadata_by_name[row["bam_file"]]:
                metadata_by_name[row["bam_file"]].append(claim)
    sources = json.loads((root / "source_inventory.json").read_text())
    for source in sources:
        source["metadata"] = metadata(source, metadata_by_name)
    write_json(root / "source_inventory.json", sources)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", "lift", "hg19-mito", "acquire", "pacbio-full", "summarize", "metadata", "supplement"))
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--isovar-repo", type=Path)
    parser.add_argument("--listing", type=Path)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--workers", type=int, default=4, choices=range(1, 9))
    parser.add_argument("--padding", type=int, default=2000)
    parser.add_argument("--source-id", action="append", help="Limit acquisition to explicit source ID(s)")
    parser.add_argument("--cache-root", type=Path, help="Shared OpenVax cache (defaults to osteosarc's cache environment)")
    parser.add_argument("--offline", action="store_true", help="Use cached remote data and local inputs only")
    args = parser.parse_args()
    required = {"prepare": ("isovar_repo", "listing", "reference"),
                "lift": ("isovar_repo", "reference"), "hg19-mito": ("isovar_repo",),
                "supplement": ("reference",)}
    for name in required.get(args.mode, ()):
        if getattr(args, name) is None:
            parser.error("%s requires --%s" % (args.mode, name.replace("_", "-")))
    cache = Cache(args.cache_root, offline=args.offline, timeout=3600 if args.mode == "pacbio-full" else 900)
    if args.mode == "prepare":
        prepare(args.run, args.isovar_repo, args.listing, args.reference, cache=cache)
    elif args.mode == "lift":
        lift_events(args.run, args.isovar_repo, args.reference, cache=cache)
    elif args.mode == "hg19-mito":
        hg19_mito(args.run, args.isovar_repo, cache=cache)
    elif args.mode == "acquire":
        acquire(args.run, args.workers, args.padding, args.source_id, cache=cache)
    elif args.mode == "pacbio-full":
        pacbio_full(args.run, cache=cache)
    elif args.mode == "metadata":
        refresh_metadata(args.run)
    elif args.mode == "supplement":
        supplement(args.run, args.reference)
    else:
        summarize(args.run)


if __name__ == "__main__":
    main()
