"""The packaged read allowlist, original record integrity and offline loading."""

from collections import Counter
import gzip
import importlib.util
import io
import json
from pathlib import Path
import socket
import sys
import zipfile

from osteosarc import OfflineError, ReadSubset, digest
import pysam
import pytest

from vaxrank.sid_test_data import extract_bundle, sid_test_data, sid_reads, sid_variants


ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("sid_builder", ROOT / "examples/osteosarc_test_data/build.py")
builder = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(builder)


def test_legacy_regeneration_cli_respects_empty_offline_cache(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Offline regeneration attempted network access")
    monkeypatch.setattr(socket.socket, "connect", forbidden)
    output = tmp_path / "cohorts.zip"
    monkeypatch.setattr(sys, "argv", [str(builder.__file__), "--cache", str(tmp_path / "cache"),
                                    "--output", str(output), "--offline"])
    with pytest.raises(OfflineError):
        builder.main()
    assert not output.exists()


def test_bundle_has_only_the_explicit_required_read_records(tmp_path):
    root = sid_test_data()
    provenance = json.loads((root / "provenance.json").read_text())
    from packaging.requirements import Requirement
    generator_requirements = (builder.RECIPE.parent / "requirements.txt").read_text()
    osteosarc = next(Requirement(line) for line in generator_requirements.splitlines()
                     if line.startswith("osteosarc"))
    assert osteosarc.specifier.contains(provenance["osteosarc_version"])
    plan = json.loads(gzip.decompress((builder.RECIPE / "selection.json.gz").read_bytes()))
    assert provenance["snapshot_id"] == plan["snapshot_id"]
    assert provenance["corrections"] is False
    assert provenance["recipe_files"] == {p.relative_to(builder.RECIPE).as_posix(): digest(p)
                                           for p in builder.RECIPE.rglob("*") if p.is_file()}
    assert len(provenance["sources"]) == 14
    assert len(provenance["cohorts"]) == len(plan["cohorts"]) == 58
    assert sum(len(c["records"]) for c in plan["cohorts"]) == 1148
    read_files = {c["path"] for c in plan["cohorts"]}
    read_files |= {c["path"] + ".bai" for c in plan["cohorts"] if c["format"] == "bam"}
    observed_files = {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file()
                      and p.name.endswith((".bam", ".bai", ".cram", ".crai", ".sam", ".sam.gz",
                                           ".fastq", ".fastq.gz", ".fq", ".fq.gz", ".input.json.gz"))}
    assert observed_files == read_files
    for cohort in plan["cohorts"]:
        path = root / cohort["path"]
        if cohort["format"] == "bam":
            with sid_reads(cohort["path"]).open() as bam:
                assert bam.check_index()
                actual = [builder.record_digest(r) for r in bam]
        elif cohort["format"] == "sam.gz":
            sam = tmp_path / "input.sam"
            sam.write_bytes(gzip.decompress(path.read_bytes()))
            with pysam.AlignmentFile(sam) as reads:
                actual = [builder.record_digest(r, text_only=True) for r in reads]
        else:
            data = json.loads(gzip.decompress(path.read_bytes()))
            names = sorted({r["sam"].split("\t")[2] for r in data["original_records"]})
            header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": n, "LN": 300000000} for n in names]})
            actual = [builder.record_digest(pysam.AlignedSegment.fromstring(r["sam"], header), text_only=True)
                      for r in data["original_records"]]
        assert actual == cohort["records"], cohort["path"]
    # Only indexed-retrieval smoke cases were reduced to one template. The
    # full reviewed context/ranking cohorts remain separately selected.
    shared = json.loads((root / "osteosarc/shared-v1/manifest.json").read_text())
    assert len(shared["cases"]) == 49 and shared["original_vaccine_loci"] == 44
    for case in shared["cases"]:
        with sid_reads("osteosarc/shared-v1/" + case["bam"]).open() as bam:
            assert len({(r.get_tag("RG") if r.has_tag("RG") else "", r.query_name) for r in bam}) == 1


def test_loader_is_offline_and_retains_native_variant_identity(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Network access while opening packaged Sid data")
    monkeypatch.setattr(socket.socket, "connect", forbidden)
    monkeypatch.setattr(socket, "create_connection", forbidden)
    variant, = sid_variants(["MAP2-chr2-209694768"])
    assert variant.allele == ("chr2", 209694768, "CCTGGGCTACTGTGTGTTCAATA", "C")
    assert sid_test_data().is_dir()


def test_exact_selection_preserves_duplicates_and_float_bits(tmp_path):
    header = pysam.AlignmentHeader.from_dict({"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]})
    a = pysam.AlignedSegment.fromstring("selected\t0\tchr1\t11\t60\t4M\t*\t0\t0\tACGT\tIIII", header)
    a.set_tag("rq", 0.123456789, value_type="f")
    b = pysam.AlignedSegment.fromstring("unneeded\t0\tchr1\t12\t60\t4M\t*\t0\t0\tACGT\tIIII", header)
    path = tmp_path / "source.bam"
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        bam.write(a)
        bam.write(a)
        bam.write(b)
    pysam.index(str(path))
    subset = ReadSubset(path, Path(str(path) + ".bai"), {})
    key = builder.record_digest(a)
    cohort = dict(path="selected.bam", format="bam", records=[key, key])
    chosen = builder.select_records(subset, cohort)
    assert Counter(builder.record_digest(r) for r in chosen) == Counter({key: 2})
    # SAM rounds the tag, so a text round trip must not pass a native BAM pin.
    rounded = pysam.AlignedSegment.fromstring(a.to_string(), header)
    assert rounded.to_string() == a.to_string()
    assert builder.record_digest(rounded) != key
    for wanted in ([key], [key, key, key], [builder.record_digest(rounded)]):
        with pytest.raises(ValueError, match="multiplicity changed"):
            builder.select_records(subset, dict(cohort, records=wanted))


@pytest.mark.parametrize("damage", ["unexpected", "missing", "checksum", "traversal"])
def test_bundle_rejects_damaged_or_unlisted_payload(tmp_path, damage):
    with zipfile.ZipFile(ROOT / "vaxrank/data/sid-test-data.zip") as original:
        members = {name: original.read(name) for name in original.namelist()}
    selected = next(name for name in members if name.endswith(".bam"))
    if damage == "unexpected":
        members["unneeded.bam"] = b"extra read payload"
    elif damage == "missing":
        del members[selected]
    elif damage == "checksum":
        data = bytearray(members[selected])
        data[-1] ^= 1
        members[selected] = bytes(data)
    else:
        manifest = json.loads(members["bundle.json"])
        manifest["files"]["../escape"] = manifest["files"].pop(selected)
        members["../escape"] = members.pop(selected)
        members["bundle.json"] = json.dumps(manifest).encode()
    data = io.BytesIO()
    with zipfile.ZipFile(data, "w") as archive:
        for name, content in members.items():
            archive.writestr(name, content)
    data.seek(0)
    destination = tmp_path / "output"
    destination.mkdir()
    with pytest.raises(ValueError):
        extract_bundle(data, destination)
    assert not (tmp_path / "escape").exists()


def test_retrieval_points_cover_required_alignments_without_per_read_queries():
    header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 1000}]})
    records = []
    for start, end in ((10, 30), (20, 40), (25, 45), (60, 80), (70, 90)):
        read = pysam.AlignedSegment(header)
        read.query_name = "read-%d" % start
        read.reference_id = 0
        read.reference_start = start
        read.query_sequence = "A" * (end - start)
        read.cigarstring = "%dM" % (end - start)
        records.append(read)
    regions = builder.retrieval_regions(records, "GRCh38")
    assert regions == [["chr1", 29, 30, "GRCh38", 1000], ["chr1", 79, 80, "GRCh38", 1000]]
    assert all(any(r.reference_start <= region[1] < r.reference_end for region in regions)
               for r in records)
