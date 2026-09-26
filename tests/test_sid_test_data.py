"""The reviewed read allowlist, exact records from openvax-v1, and offline builds."""

from collections import Counter
import gzip
import json
from pathlib import Path
import socket

import osteosarc
from osteosarc import ReadSubset, digest
from osteosarc.cohort_bundle import record_digest, select_records
from osteosarc.shared import published
from packaging.requirements import Requirement
import pysam
import pytest

from vaxrank.sid_test_data import (
    READS, build_sid_test_data, recipe_directory, sid_test_data, sid_reads, sid_variants)


ROOT = Path(__file__).resolve().parents[1]
RECIPE = recipe_directory()


def test_test_data_has_only_the_explicit_required_read_records(tmp_path):
    root = sid_test_data()
    provenance = json.loads((root / "provenance.json").read_text())
    runtime = next(Requirement(line) for line in (ROOT / "requirements.txt").read_text().splitlines()
                   if line.startswith("osteosarc"))
    assert provenance["osteosarc_version"] == osteosarc.__version__
    assert runtime.specifier.contains(provenance["osteosarc_version"])
    assert provenance["reads"] == dict(bundle=READS, manifest_sha256=published(READS)["manifest_sha256"])
    plan = json.loads(gzip.decompress((RECIPE / "selection.json.gz").read_bytes()))
    assert provenance["snapshot_id"] == plan["snapshot_id"]
    assert provenance["corrections"] is False
    assert provenance["recipe_files"] == {p.relative_to(RECIPE).as_posix(): digest(p)
                                           for p in RECIPE.rglob("*") if p.is_file()}
    assert len(provenance["sources"]) == 14
    assert all(s["identity"]["url"] == url for url, s in provenance["sources"].items())
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
                actual = [record_digest(r) for r in bam]
        elif cohort["format"] == "sam.gz":
            sam = tmp_path / "input.sam"
            sam.write_bytes(gzip.decompress(path.read_bytes()))
            with pysam.AlignmentFile(sam) as reads:
                actual = [record_digest(r, text_only=True) for r in reads]
        else:
            data = json.loads(gzip.decompress(path.read_bytes()))
            names = sorted({r["sam"].split("\t")[2] for r in data["original_records"]})
            header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": n, "LN": 300000000} for n in names]})
            actual = [record_digest(pysam.AlignedSegment.fromstring(r["sam"], header), text_only=True)
                      for r in data["original_records"]]
        assert actual == cohort["records"], cohort["path"]
    # Only indexed-retrieval smoke cases were reduced to one template. The
    # full reviewed context/ranking cohorts remain separately selected.
    shared = json.loads((root / "osteosarc/shared-v1/manifest.json").read_text())
    assert len(shared["cases"]) == 49 and shared["original_vaccine_loci"] == 44
    for case in shared["cases"]:
        with sid_reads("osteosarc/shared-v1/" + case["bam"]).open() as bam:
            assert len({(r.get_tag("RG") if r.has_tag("RG") else "", r.query_name) for r in bam}) == 1


def test_build_is_offline_once_the_reads_are_cached(tmp_path, monkeypatch):
    sid_test_data()  # downloads openvax-v1 into the osteosarc cache if needed
    def forbidden(*args, **kwargs):
        raise AssertionError("Network access while building Sid test data")
    monkeypatch.setattr(socket.socket, "connect", forbidden)
    monkeypatch.setattr(socket, "create_connection", forbidden)
    root = build_sid_test_data(tmp_path)
    files = sorted(p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file())
    assert files == sorted(p.relative_to(sid_test_data()).as_posix()
                           for p in sid_test_data().rglob("*") if p.is_file())
    for name in files:
        if name != "provenance.json":
            assert digest(root / name) == digest(sid_test_data() / name), name
    variant, = sid_variants(["MAP2-chr2-209694768"])
    assert variant.allele == ("chr2", 209694768, "CCTGGGCTACTGTGTGTTCAATA", "C")


def test_build_refuses_a_nonempty_destination(tmp_path):
    (tmp_path / "existing").write_text("keep")
    with pytest.raises(ValueError, match="must be empty"):
        build_sid_test_data(tmp_path)
    assert [p.name for p in tmp_path.iterdir()] == ["existing"]


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
    key = record_digest(a)
    cohort = dict(path="selected.bam", format="bam", records=[key, key])
    chosen = select_records(subset, cohort)
    assert Counter(record_digest(r) for r in chosen) == Counter({key: 2})
    # SAM rounds the tag, so a text round trip must not pass a native BAM pin.
    rounded = pysam.AlignedSegment.fromstring(a.to_string(), header)
    assert rounded.to_string() == a.to_string()
    assert record_digest(rounded) != key
    for wanted in ([key], [key, key, key], [record_digest(rounded)]):
        with pytest.raises(ValueError, match="multiplicity changed"):
            select_records(subset, dict(cohort, records=wanted))
