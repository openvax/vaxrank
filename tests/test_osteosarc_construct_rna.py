"""Original multimodal RNA versus independent Sid vaccine sequence definitions.

These tests deliberately distinguish reconstructability from Vaxrank's RNA
selection gates and from final historical ranking agreement (owned by #414).
"""

import gzip
from hashlib import sha256
import json
import socket
from unittest.mock import Mock

import pysam
import pytest

from vaxrank.construct_sequence import ConstructSequence
from vaxrank.core_logic import vaccine_peptides_for_variant
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.native_serialization import from_native_json, to_native_json

from .osteosarc_construct_helpers import (
    DATA, DOCUMENTED, RNA, construct_from_result, reconstruct_source, verify_rna_inputs,
)
from .osteosarc_selection_helpers import load_selection_inputs


@pytest.fixture(autouse=True)
def no_network(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Sid original-RNA regressions must stay offline")
    monkeypatch.setattr(socket.socket, "connect", forbidden)
    monkeypatch.setattr(socket, "create_connection", forbidden)


@pytest.fixture(scope="module")
def inputs(tmp_path_factory):
    verify_rna_inputs()
    return load_selection_inputs(tmp_path_factory.mktemp("sid-construct-reference"))


@pytest.fixture(scope="module")
def results(inputs):
    return {case["source_id"]: reconstruct_source(inputs, case) for case in RNA["cases"]}


def test_documented_sequences_are_independent_and_modalities_explicit():
    assert DOCUMENTED["historical_ranking_settings"] is None
    assert DOCUMENTED["complete_historical_mrna_product"] is None
    assert DOCUMENTED["assay_only_comparators"] == ["ISFDTDTGL", "RFHATISF"]
    for record in DOCUMENTED["records"]:
        assert record["native_sequence"] in DOCUMENTED["source_protein_context"]
        assert record["sequence"] == (record["n_terminal_addition"] + record["native_sequence"]
                                      + record["c_terminal_addition"])
        start, end = record["intended_ligand_interval"]
        assert record["sequence"][start:end] == DOCUMENTED["intended_ligand"]
        start, end = record["native_mutation_interval"]
        assert record["native_sequence"][start:end] == "I"
        assert record["chemical_modifications"] is None
        assert record["modality"] == ("mrna" if record["provider"] == "mRNA" else "peptide")


@pytest.mark.parametrize("case", RNA["cases"], ids=lambda c: c["source_id"])
def test_original_records_preserve_qualities_flags_tags_and_multiplicity(case):
    with pysam.AlignmentFile(DATA / "rna" / case["bam"]) as bam:
        assert sha256(str(bam.header).encode()).hexdigest() == case["selection"]["header_sha256"]
        records = list(bam)
    assert len(records) == len(case["selection"]["records"])
    for read, original in zip(records, case["selection"]["records"]):
        assert sha256(read.to_string().encode("ascii")).hexdigest() == original["sam_sha256"]
        assert read.query_name == original["name"]
    if "single-cell" in case["modality"]:
        assert any(r.has_tag("CB") for r in records)
        if case["modality"] != "PacBio single-cell":
            assert any(r.has_tag("UB") for r in records)
    if case["modality"] == "PacBio single-cell":
        assert any(r.query_qualities is None for r in records)
        assert any(r.query_qualities is not None for r in records)


@pytest.mark.parametrize("case", RNA["cases"], ids=lambda c: c["source_id"])
def test_reconstruction_matches_independent_reference_plus_reported_edit(inputs, results, case):
    result = results[case["source_id"]]
    protein = result.top_protein_sequence
    assert protein is not None
    fragment = MutantProteinFragment.from_isovar_result(result)
    reference = inputs[0].transcript_by_id("ENST00000360184").protein_sequence
    assert reference[313] == "V"
    expected = reference[:313] + "I" + reference[314:]
    start = 313 - fragment.mutant_amino_acid_start_offset
    assert fragment.amino_acids == expected[start:start + len(fragment.amino_acids)]
    assert fragment.amino_acids == protein.amino_acids
    assert fragment.n_alt_reads == result.num_alt_reads
    assert fragment.n_alt_fragments == result.num_alt_fragments
    assert fragment.supporting_reference_transcripts
    assert DOCUMENTED["intended_ligand"] in fragment.amino_acids


@pytest.mark.parametrize("record", DOCUMENTED["records"], ids=lambda r: r["id"])
def test_all_documented_constructs_map_without_giving_tails_rna_support(results, record):
    for sid, result in results.items():
        if sid == "8b8a9a02cbcbcf4b" and record["id"] == "dync1h1-mrna-long":
            with pytest.raises(ValueError, match="absent or non-unique"):
                construct_from_result(result, record, sid)
            continue
        construct = construct_from_result(result, record, sid)
        assert construct.sequence == record["sequence"]
        assert construct.native_sequence == record["native_sequence"]
        assert construct.native_residue_offsets.count(None) == len(record["c_terminal_addition"])
        assert sid in construct.native_antigen.source_identifier
        assert construct.native_antigen.species == "Homo sapiens"
        assert from_native_json(to_native_json(construct), ConstructSequence) == construct
        if record["c_terminal_addition"]:
            assert construct.native_residue_offsets[-2:] == (None, None)
            assert construct.sequence_edits[0].rationale.evidence_level == "reported"


def test_pacbio_reconstructability_does_not_bypass_independent_selection_gate(results):
    result = results["0066232879babe83"]
    assert result.top_protein_sequence is not None
    assert result.num_alt_reads == result.num_alt_fragments == 2
    assert {k for k, passed in result.filter_values.items() if not passed} == {"min_num_alt_reads"}
    model = Mock(side_effect=AssertionError("Failed RNA gate must prevent ranking predictions"))
    assert vaccine_peptides_for_variant(result, model) == []
    assert not model.mock_calls
    upstream = next(c for c in RNA["cases"] if c["source_id"] == "0066232879babe83")
    assert upstream["upstream_independent_primary"]["alt_missing_quality"] == 187
    assert upstream["upstream_independent_primary"]["counts"]["alt"] == 189
    assert upstream["upstream_region_counts"]["reads"]["alt"] == 2


def test_bulk_read_and_fragment_counts_remain_distinct(results):
    for sid in ("1b66c15da594a3ef", "50da1d13e05059fd"):
        result = results[sid]
        assert result.num_alt_reads > result.num_alt_fragments >= 2


def test_inventory_retains_all_sources_and_unassessed_cohort():
    inventory = json.loads(gzip.decompress((DATA / "rna" / "inventory.json.gz").read_bytes()))
    assert len(inventory["variant_ids"]) == 44
    assert len(inventory["source_ids"]) == len(inventory["rows"]) == 164
    assert {r["source_id"] for r in inventory["rows"]} == set(inventory["source_ids"])
    assert {c["timepoint"] for c in RNA["cases"]} == {"T0", "T1", "T2", "T3"}
    assert len({c["source_id"] for c in RNA["cases"]}) == 9
    assert set(RNA["inventory_status_counts"]) != {"audited"}
    assert "not unbiased VAF" in RNA["scope"]
