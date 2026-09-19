"""Scientific acquisition contracts, independent of report wording/layout."""

import importlib.util
from copy import deepcopy
import json
from pathlib import Path
import socket
import sys

from osteosarc import Cache
import pysam
import pytest

from .osteosarc_selection_helpers import DATA, load_selection_inputs, reconstruct_selection


PATH = Path(__file__).resolve().parents[1] / "examples/osteosarc_read_corpus/build.py"
SPEC = importlib.util.spec_from_file_location("osteosarc_read_corpus", PATH)
corpus = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(corpus)


@pytest.fixture(autouse=True)
def offline_corpus_tests(monkeypatch, tmp_path):
    def forbidden(*args, **kwargs):
        raise AssertionError("Corpus regressions must not access the network")

    monkeypatch.setattr(socket.socket, "connect", forbidden)
    monkeypatch.setattr(socket, "create_connection", forbidden)
    monkeypatch.setenv("OSTEOSARC_CACHE", str(tmp_path / "cache"))
    monkeypatch.setenv("VAXRANK_REF_PEPTIDES_DIR", str(tmp_path / "kmers"))


def test_library_first_policy_filters_retries_without_dropping_inventory(tmp_path):
    sources = [dict(source_id=s) for s in ("primary", "alternative", "tagged")]
    policy = {"sources": [dict(source_id=s["source_id"], acquire=s["source_id"] != "alternative", reason="reviewed") for s in sources]}
    corpus.write_json(tmp_path / "source-selection.json", policy)
    selected = corpus.acquisition_sources(tmp_path, sources)
    assert [s["source_id"] for s in selected] == ["primary", "tagged"]
    assert corpus.acquisition_sources(tmp_path, sources, ["alternative"]) == []
    assert len(sources) == 3
    with pytest.raises(ValueError):
        corpus.acquisition_sources(tmp_path, sources, ["unknown"])
    policy["sources"].pop()
    corpus.write_json(tmp_path / "source-selection.json", policy)
    with pytest.raises(ValueError):
        corpus.acquisition_sources(tmp_path, sources)


def test_deferred_products_are_not_failures_or_zero_coverage(tmp_path):
    sources = [dict(source_id="alternative", metadata={})]
    corpus.write_json(tmp_path / "source_inventory.json", sources)
    corpus.write_json(tmp_path / "source-selection.json", {"sources": [dict(source_id="alternative", acquire=False, reason="reviewed alternative")]})
    corpus.write_json(tmp_path / "events.json", [dict(event_id="event")])
    corpus.write_json(tmp_path / "provenance.json", {})
    corpus.summarize(tmp_path)
    coverage = json.loads((tmp_path / "coverage.json").read_text())
    assert coverage[0]["status"] == "deferred_by_policy"
    assert coverage[0]["records"] is None
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert manifest["sources"] == 1
    assert manifest["sources_receipts"][0]["observed_acquisition_status"] == "not_acquired"


def test_minimal_allele_identity_preserves_compound_events():
    assert corpus.canonical_allele("chr1", 10, "AC", "AT") == ("chr1", 11, "C", "T")
    assert corpus.canonical_allele("chr1", 10, "AG", "GT") != corpus.canonical_allele("chr1", 10, "A", "G")


def test_breakpoints_and_contig_boundaries():
    header = {"SQ": [{"SN": "1", "LN": 1000}, {"SN": "MT", "LN": 16569}]}
    events = [dict(event_id="deletion", breakpoints=[["chr1", 2], ["chr1", 998]]),
              dict(event_id="mito", chrom="chrM", pos=100, ref="A"),
              dict(event_id="absent", chrom="chrX", pos=5, ref="C")]
    regions, missing = corpus.intervals(events, header, 20)
    assert regions == [("1", 0, 22), ("1", 977, 1000), ("MT", 79, 120)]
    assert missing == ["absent"]


def test_conflicting_timepoint_claims_are_retained():
    source = dict(key="rna-seq/BG009368.bam", name="T0 BostonGene RNA", tissue="tumor")
    metadata = corpus.metadata(source, {"BG009368.bam": [{"timepoint": "T1"}]})
    assert metadata["timepoint_conflict"]
    assert metadata["timepoints"] == ["T1"]
    assert metadata["catalogue_timepoint_claims"] == ["T0"]


@pytest.mark.parametrize("sequences, status", [
    ([{"SN": "chr1", "LN": 1000}], "no_genomic_coordinates"),
    ([{"SN": "chr1", "LN": 249250621}, {"SN": "chr2", "LN": 243199373}],
     "requires_coordinate_mapping"),
    ([{"SN": "chr1", "LN": 248956422}, {"SN": "chr2", "LN": 242193529}],
     "missing_genomic_index"),
])
def test_unacquired_sources_are_not_reported_as_zero_reads(tmp_path, sequences, status):
    source_bam = tmp_path / "source.bam"
    with pysam.AlignmentFile(source_bam, "wb", header={"HD": {"SO": "coordinate"}, "SQ": sequences}):
        pass
    source = dict(source_id="fixture", url=str(source_bam), listed_indexes=[])
    events = [dict(event_id="one", chrom="chr1", pos=121, ref="A")]
    corpus.write_json(tmp_path / "events.json", events)
    result = corpus.acquire_one(source, events, tmp_path, 5, cache=Cache(tmp_path / "cache", offline=True))
    assert result["status"] == status
    assert "record_counts" not in result


@pytest.mark.parametrize("platform", ["ILMN", "ONT", "PacBio"])
def test_acquisition_preserves_multiplicity_and_fetches_off_region_mate(tmp_path, platform):
    """Fetch a complete pair, overlap regions once, and retain duplicate records."""
    header = {"HD": {"SO": "coordinate"}, "SQ": [
        {"SN": "chr1", "LN": 248956422}, {"SN": "chr2", "LN": 242193529}]}
    source_bam = tmp_path / "source.bam"
    records = []
    for name, start, flag, mate_start in [
        ("pair", 100, 99, 10000), ("duplicate", 110, 1024, -1),
        ("duplicate", 110, 1024, -1), ("secondary", 120, 256, -1),
        ("supplementary", 125, 2048, -1), ("noquality", 126, 0, -1),
        ("pair", 10000, 147, 100),
    ]:
        read = pysam.AlignedSegment()
        read.query_name, read.query_sequence, read.flag = name, "A" * 50, flag
        read.reference_id, read.reference_start = 0, start
        read.mapping_quality, read.cigarstring = 60 if name == "pair" else 0, "50M"
        read.query_qualities = None if name == "noquality" else pysam.qualitystring_to_array("I" * 50)
        read.next_reference_id = 0 if mate_start >= 0 else -1
        read.next_reference_start = mate_start
        read.template_length = 9950 if flag == 99 else -9950 if flag == 147 else 0
        read.set_tag("ZZ", "retained")
        read.set_tag("ZF", 0.12345671, value_type="f")
        records.append(read)
    with pysam.AlignmentFile(source_bam, "wb", header=header) as bam:
        for read in records:
            bam.write(read)
    pysam.index(str(source_bam))
    source = dict(source_id="fixture", url=str(source_bam),
                  listed_indexes=[str(source_bam) + ".bai"], metadata={"platform": platform})
    events = [dict(event_id="one", chrom="chr1", pos=121, ref="A"),
              dict(event_id="two", chrom="chr1", pos=125, ref="A")]
    corpus.write_json(tmp_path / "events.json", events)

    cache = Cache(tmp_path / "cache", offline=True)
    result = corpus.acquire_one(source, events, tmp_path, 5, cache=cache)
    assert result["status"] == "ok", result
    with pysam.AlignmentFile(tmp_path / "alignments/fixture/reads.bam") as bam:
        observed = list(bam)
    assert len(observed) == 7
    with pysam.AlignmentFile(source_bam) as original:
        assert [corpus.record_identity(r) for r in observed] == [corpus.record_identity(r) for r in original]
    assert result["read_count_cap"] is None
    receipt = json.loads((tmp_path / "alignments/fixture/receipt.json").read_text())
    assert receipt["record_counts"]["records"] == 7
    assert receipt["record_counts"]["secondary"] == 1
    assert receipt["record_counts"]["supplementary"] == 1
    assert receipt["record_counts"]["missing_quality"] == 1
    assert receipt["osteosarc_receipt"]["scope"] == "regional_records_and_paired_mates"
    assert receipt["osteosarc_receipt"]["request"]["fetch_pairs"] is True
    assert receipt["source_index"]["sha256"] == corpus.digest(str(source_bam) + ".bai")
    assert corpus.acquire_one(source, events, tmp_path, 5, cache=cache) == receipt

    # An intact BAM cannot hide a damaged exported index during a resume.
    (tmp_path / "alignments/fixture/reads.bam.bai").write_bytes(b"damaged")
    failed = corpus.acquire_one(source, events, tmp_path, 5, cache=cache)
    assert failed["status"] == "integrity_error"
    assert "reads.bam.bai" in failed["error"]
    assert corpus.acquire_one(source, events, tmp_path, 5, cache=cache) == failed
    assert (tmp_path / "alignments/fixture/reads.bam.bai").read_bytes() == b"damaged"


def test_combining_queries_preserves_source_duplicates_not_query_duplicates(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "build", corpus)
    spec = importlib.util.spec_from_file_location("corpus_finalize", PATH.with_name("finalize.py"))
    finalizer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(finalizer)
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    paths = [tmp_path / name for name in ("first.bam", "second.bam", "union.bam")]
    for path, names in zip(paths, [("shared", "shared", "first"), ("shared", "shared", "second")]):
        with pysam.AlignmentFile(path, "wb", header=header) as bam:
            for name in names:
                read = pysam.AlignedSegment()
                read.query_name, read.query_sequence = name, "A" * 20
                read.reference_id, read.reference_start = 0, 100
                read.cigarstring, read.mapping_quality = "20M", 0
                read.flag = 1024
                read.set_tag("ZF", 0.12345671, value_type="f")
                bam.write(read)
    finalizer.combine_bams(*paths)
    with pysam.AlignmentFile(paths[2]) as bam:
        observed = list(bam)
    assert sorted(r.query_name for r in observed) == ["first", "second", "shared", "shared"]
    assert all(r.is_duplicate and r.mapping_quality == 0 for r in observed)
    with pysam.AlignmentFile(paths[0]) as original:
        expected_float = next(original).get_tag("ZF")
    assert all(r.get_tag("ZF") == expected_float for r in observed)


def test_record_identity_retains_sub_sam_precision_float_tags():
    header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 1000}]})
    reads = []
    for value in (0.12345671, 0.12345674):
        read = pysam.AlignedSegment(header)
        read.query_name, read.query_sequence = "same_name", "A"
        read.reference_id, read.reference_start, read.cigarstring = 0, 0, "1M"
        read.set_tag("ZF", value, value_type="f")
        reads.append(read)
    assert reads[0].to_string() == reads[1].to_string()
    assert corpus.record_identity(reads[0]) != corpus.record_identity(reads[1])


def test_structural_coverage_does_not_double_count_one_read_at_both_endpoints(tmp_path):
    source = {"source_id": "one", "metadata": {"assay": "RNA", "platform": "ONT"}}
    event = {"event_id": "sv", "breakpoints": [["chr1", 110], ["chr1", 140]]}
    for filename, value in (("source_inventory.json", [source]), ("events.json", [event]), ("provenance.json", {})):
        corpus.write_json(tmp_path / filename, value)
    directory = tmp_path / "alignments/one"
    directory.mkdir(parents=True)
    corpus.write_json(directory / "receipt.json", {"status": "ok", "assembly": "GRCh38"})
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(directory / "reads.bam", "wb", header=header) as bam:
        read = pysam.AlignedSegment()
        read.query_name, read.query_sequence = "same", "A" * 50
        read.reference_id, read.reference_start = 0, 100
        read.cigarstring = "50M"
        bam.write(read)
        bam.write(read)
    pysam.index(str(directory / "reads.bam"))
    corpus.summarize(tmp_path)
    row = json.loads((tmp_path / "coverage.json").read_text())[0]
    assert row["records"] == row["aligned_block_overlap"] == 2
    assert row["primary_fragment_names"] == 1
    assert [e["records"] for e in row["endpoints"]] == [2, 2]


def test_missing_acquisition_is_not_reported_as_zero_coverage(tmp_path):
    for filename, value in (("source_inventory.json", [{"source_id": "unavailable", "metadata": {}}]),
                            ("events.json", [{"event_id": "locus"}]), ("provenance.json", {})):
        corpus.write_json(tmp_path / filename, value)
    corpus.summarize(tmp_path)
    row = json.loads((tmp_path / "coverage.json").read_text())[0]
    assert row["status"] == "not_acquired"
    assert row["records"] is None


def test_hg19_mitochondrial_coordinate_override_is_used_for_retrieval_and_coverage(tmp_path):
    header = {"HD": {"SO": "coordinate"}, "SQ": [
        {"SN": "chr1", "LN": 249250621}, {"SN": "chr2", "LN": 243199373},
        {"SN": "chrM", "LN": 16571}]}
    source_bam = tmp_path / "hg19.bam"
    with pysam.AlignmentFile(source_bam, "wb", header=header) as bam:
        read = pysam.AlignedSegment()
        read.query_name, read.query_sequence = "native_hg19", "G"
        read.reference_id, read.reference_start = 2, 12994
        read.cigarstring = "1M"
        bam.write(read)
    pysam.index(str(source_bam))
    event = dict(event_id="MT-ND5", chrom="chrM", pos=12994, ref="G", alt="A")
    native = dict(event, pos=12995)
    source = dict(source_id="hg19", url=str(source_bam),
                  listed_indexes=[str(source_bam) + ".bai"], metadata={})
    for name, value in (("events.json", [event]), ("events-GRCh37.json", [event]),
                        ("events-hg19-MT.json", [native]), ("coordinate-mapping.json", {}),
                        ("source_inventory.json", [source]), ("provenance.json", {})):
        corpus.write_json(tmp_path / name, value)

    result = corpus.acquire_one(source, [event], tmp_path, 0, cache=Cache(tmp_path / "cache", offline=True))
    assert result["status"] == "ok", result
    assert result["record_counts"]["records"] == 1
    assert result["native_event_overrides"]["MT-ND5"]["pos"] == 12995
    corpus.summarize(tmp_path)
    row = json.loads((tmp_path / "coverage.json").read_text())[0]
    assert row["records"] == row["primary_fragment_names"] == 1
    request = result["osteosarc_receipt"]["request"]
    assert request["regions"] == [dict(contig="chrM", start=12994, end=12995,
                                       assembly="GRCh37", reference_length=16571)]
    corpus.write_json(tmp_path / "events-hg19-MT.json", [dict(native, pos=12996)])
    with pytest.raises(ValueError, match="request changed"):
        corpus.acquire_one(source, [event], tmp_path, 0)


def test_download_imports_verified_prior_bytes_into_the_shared_cache(tmp_path):
    local = tmp_path / "source.txt"
    local.write_text("published source bytes")
    url = "https://example.test/source.txt"
    old = dict(url=url, sha256=corpus.digest(local), bytes=local.stat().st_size,
               retrieved_utc="2026-09-18T00:00:00Z")
    corpus.write_json(local.with_name("source.txt.receipt.json"), old)
    cache = Cache(tmp_path / "cache", offline=True)
    assert corpus.download(url, local, cache=cache) == old
    exported = tmp_path / "second.txt"
    receipt = corpus.download(url, exported, cache=cache)
    assert exported.read_bytes() == local.read_bytes()
    assert receipt["osteosarc_version"] == "0.1.0"
    cached = cache.path(cache.fetch(url))
    assert cached == cache.root / "objects" / "sha256" / (old["sha256"] + ".txt")
    cached.write_text("tampered")
    with pytest.raises(corpus.OsteosarcError, match="modified"):
        corpus.download(url, tmp_path / "third.txt", cache=Cache(cache.root, offline=True))


@pytest.fixture(scope="module")
def acquired_selection_inputs(tmp_path_factory):
    """Pass all five real pinned BAMs through the adopted acquisition path."""
    root = tmp_path_factory.mktemp("osteosarc-selection-acquisition")
    original = load_selection_inputs(root / "original")
    genome, cases, metadata = original
    acquired_cases = deepcopy(cases)
    cache = Cache(root / "cache", offline=True)
    for variant_id, case in acquired_cases.items():
        source_bam = DATA / "isovar" / case["bam"]
        run = root / variant_id
        run.mkdir()
        events = [dict(case["variant"], event_id=variant_id)]
        corpus.write_json(run / "events.json", events)
        source = dict(source_id=case["source_id"], url=str(source_bam),
                      listed_indexes=[str(source_bam) + ".bai"])
        receipt = corpus.acquire_one(source, events, run, 2000, cache=cache)
        assert receipt["status"] == "ok", receipt
        extracted = run / "alignments" / case["source_id"] / "reads.bam"
        with pysam.AlignmentFile(source_bam) as expected, pysam.AlignmentFile(extracted) as actual:
            assert actual.header.to_dict() == expected.header.to_dict()
            assert [corpus.record_identity(r) for r in actual] == [corpus.record_identity(r) for r in expected]
        case["bam"] = str(extracted)
    return original, (genome, acquired_cases, metadata)


DOCUMENTED = json.loads((DATA / "documented.json").read_text())


@pytest.mark.parametrize("record", DOCUMENTED["records"], ids=lambda record: record["id"])
def test_osteosarc_acquisition_preserves_real_rna_final_selections(acquired_selection_inputs, record):
    from topiary import CachedPredictor, TopiaryPredictor
    from vaxrank.core_logic import vaccine_peptides_for_variant

    baseline, acquired = acquired_selection_inputs
    original, original_config = reconstruct_selection(
        baseline, record["variant_id"], len(record["native_sequence"]))
    result, config = reconstruct_selection(
        acquired, record["variant_id"], len(record["native_sequence"]))
    assert result.variant == original.variant
    assert result.filter_values == original.filter_values
    assert result.num_alt_reads == original.num_alt_reads
    assert result.num_alt_fragments == original.num_alt_fragments
    assert result.top_protein_sequence.amino_acids == original.top_protein_sequence.amino_acids
    predictor = TopiaryPredictor(models=[CachedPredictor.from_topiary_output(str(DATA / "netmhcpan42.tsv"))])
    assert all(model.fallback is None for model in predictor.models)
    original_selected = vaccine_peptides_for_variant(original, predictor, vaccine_config=original_config)
    selected = vaccine_peptides_for_variant(result, predictor, vaccine_config=config)
    sequences = [p.mutant_protein_fragment.amino_acids for p in selected]
    assert sequences == [p.mutant_protein_fragment.amino_acids for p in original_selected]
    if record["id"] in {"dync1h1-mrna-minimal", "dync1h1-cegat", "exoc4-jlf", "exoc4-cegat"}:
        assert sequences == [record["native_sequence"]]
    elif record["id"] == "h1-2-cegat":
        assert sequences == []
    else:
        assert len(sequences) == 1 and sequences[0] != record["native_sequence"]


@pytest.mark.parametrize("ref, alt", [
    ("CCTGGGCTACTGTGTGTTCAATA", "C"),
    ("CCTGGGCTACTGTGTGTTCAATAAGTACACAGT", "CAGGG"),
], ids=["published-deletion", "corrected-complex-replacement"])
def test_osteosarc_acquisition_preserves_map2_rna_gate(acquired_selection_inputs, ref, alt):
    from unittest.mock import Mock
    from varcode import Variant
    from vaxrank.core_logic import vaccine_peptides_for_variant

    _, acquired = acquired_selection_inputs
    variant = Variant("2", 209694768, ref, alt, ensembl=acquired[0])
    result, config = reconstruct_selection(acquired, "MAP2-chr2-209694768", variant=variant)
    assert result.num_alt_reads == result.num_alt_fragments == 1
    assert result.top_protein_sequence is None
    predictor = Mock(side_effect=AssertionError("No protein should reach MHC prediction"))
    assert vaccine_peptides_for_variant(result, predictor, vaccine_config=config) == []
    assert not predictor.mock_calls
