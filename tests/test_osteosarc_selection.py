"""Real RNA -> real offline predictions -> final documented peptide comparison.

No synthetic scores or network/model fallback. The historical provider settings
are unavailable, so matching sequences are not claims of protocol reproduction.
"""

from hashlib import sha256
import json
import logging
import socket
from unittest.mock import Mock

import msgspec
import pytest
from topiary import CachedPredictor, TopiaryPredictor

from vaxrank.construct_sequence import ConstructEvidence, ConstructSequence, ConstructSequenceEdit
from vaxrank.core_logic import vaccine_peptides_for_variant
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.native_serialization import from_native_json, to_native_json
from vaxrank.vaccine_antigen import VaccineAntigen
from vaxrank.vaccine_config import VaccineConfig

from .osteosarc_selection_helpers import DATA, load_selection_inputs, reconstruct_selection


DOCUMENTED = json.loads((DATA / "documented.json").read_text())
RECORDS = DOCUMENTED["records"]


@pytest.fixture(autouse=True)
def offline_selection(monkeypatch, tmp_path):
    def forbidden(*args, **kwargs):
        raise AssertionError("Sid selection regressions must never access the network")
    monkeypatch.setattr(socket.socket, "connect", forbidden)
    monkeypatch.setattr(socket, "create_connection", forbidden)
    monkeypatch.setenv("VAXRANK_REF_PEPTIDES_DIR", str(tmp_path / "kmers"))


@pytest.fixture(scope="module")
def inputs(tmp_path_factory):
    return load_selection_inputs(tmp_path_factory.mktemp("sid-selection"))


def cached_predictor():
    return TopiaryPredictor(models=[CachedPredictor.from_topiary_output(
        str(DATA / "netmhcpan42.tsv"))])


def test_regeneration_refuses_to_overwrite_existing_evidence(tmp_path, monkeypatch):
    from tests.data.osteosarc.selection_validation import generate_predictions

    reconstruction = Mock(side_effect=AssertionError("Existing output must fail before reconstruction"))
    monkeypatch.setattr(generate_predictions, "load_selection_inputs", reconstruction)
    with pytest.raises(ValueError, match="Choose a new output directory"):
        generate_predictions.main(tmp_path)
    reconstruction.assert_not_called()


def test_pinned_predictions_are_complete_real_model_outputs():
    metadata = json.loads((DATA / "predictions_manifest.json").read_text())
    for filename, key in [("documented.json", "documented_sha256"),
                          ("isovar/manifest.json", "input_manifest_sha256"),
                          ("netmhcpan42.tsv", "output_sha256")]:
        assert sha256((DATA / filename).read_bytes()).hexdigest() == metadata[key]
    assert metadata["package_versions"]["isovar"] == "1.8.1"
    assert metadata["predictor_name"] == "netMHCpan"
    assert metadata["predictor_version"] == "4.2c"
    assert metadata["vaccine_config"] == msgspec.to_builtins(VaccineConfig())
    assert metadata["epitope_config"] == msgspec.to_builtins(EpitopeConfig())
    cache = CachedPredictor.from_topiary_output(str(DATA / "netmhcpan42.tsv"))
    assert cache.fallback is None
    frame = cache.predict_peptides_dataframe(metadata["requested_peptides"])
    expected = {(p, a, k) for p in metadata["requested_peptides"]
                for a in DOCUMENTED["prediction_alleles"]
                for k in ("pMHC_affinity", "pMHC_presentation")}
    assert set(frame[["peptide", "allele", "kind"]].itertuples(index=False, name=None)) == expected
    assert len(frame) == len(expected)
    assert DOCUMENTED["historical_ranking_settings"] is None
    assert "HLA-A*01:11N" in DOCUMENTED["clinical_class_i"]
    assert "HLA-A*01:11N" not in cache.alleles


@pytest.mark.parametrize("record", RECORDS, ids=lambda r: r["id"])
def test_documented_native_and_manufactured_sequences_stay_distinct(inputs, record):
    result, _ = reconstruct_selection(inputs, record["variant_id"], len(record["native_sequence"]))
    fragment = MutantProteinFragment.from_isovar_result(result)
    native = record["native_sequence"]
    offset = fragment.amino_acids.index(native)
    assert fragment.amino_acids.count(native) == 1
    assert [fragment.mutant_amino_acid_start_offset - offset,
            fragment.mutant_amino_acid_end_offset - offset] == record["native_mutation_interval"]
    assert fragment.n_alt_reads == result.num_alt_reads
    assert fragment.n_alt_fragments == result.num_alt_fragments
    assert fragment.supporting_reference_transcripts
    candidates = [f.amino_acids for _, f in fragment.generate_subsequences(len(native))]
    assert native in candidates and len(set(candidates)) > 1
    evidence = ConstructEvidence(record["source"], "documented", "Published provider sequence")
    edits = []
    for boundary, extra in [(offset, record["n_terminal_addition"]),
                            (offset + len(native), record["c_terminal_addition"])]:
        if extra:
            edits.append(ConstructSequenceEdit(boundary, boundary, extra, evidence))
    construct = ConstructSequence(
        record["id"], record["sequence"], record["modality"], evidence,
        native_antigen=VaccineAntigen.from_mutant_protein_fragment(fragment),
        native_start=offset, native_end=offset + len(native), sequence_edits=tuple(edits))
    assert construct.native_sequence == native
    expected_modality = (
        "mrna" if record["id"] in {"dync1h1-mrna-minimal", "dync1h1-mrna-long"}
        else "peptide")
    assert construct.modality == expected_modality
    assert construct.native_residue_offsets.count(None) == len(record["sequence"]) - len(native)
    assert from_native_json(to_native_json(construct), ConstructSequence) == construct
    assert fragment.amino_acids == result.top_protein_sequence.amino_acids


@pytest.mark.parametrize("record", RECORDS, ids=lambda r: r["id"])
def test_real_cached_final_selections_against_documented_sequences(inputs, record, caplog):
    result, config = reconstruct_selection(inputs, record["variant_id"], len(record["native_sequence"]))
    if record["id"] == "h1-2-cegat":
        # Upstream fixtures select read names for useful branch coverage, not
        # unbiased VAF. Preserve this actual gate, not a claim about binding.
        assert not result.filter_values["min_ratio_alt_to_other_fragments"]
        predictor = Mock(side_effect=AssertionError("RNA gate must precede MHC prediction"))
        assert vaccine_peptides_for_variant(result, predictor, vaccine_config=config) == []
        assert not predictor.mock_calls
        return
    assert result.passes_all_filters
    predictor = cached_predictor()
    frame = predictor.predict_from_named_sequences({"rna": result.top_protein_sequence.amino_acids})
    assert frame.columns.is_unique, "Upstream Topiary #296: duplicate cache occurrence coordinates"
    with caplog.at_level(logging.ERROR):
        selected = vaccine_peptides_for_variant(result, predictor, vaccine_config=config)
    assert not [r for r in caplog.records if r.levelno >= logging.ERROR]
    sequences = [p.mutant_protein_fragment.amino_acids for p in selected]
    # Only independently published strings are exact-match goldens. The other
    # cases preserve honest disagreement; current output is not relabelled as
    # the historical expectation, and no ranking parameters are tuned to it.
    matches = {"dync1h1-mrna-minimal", "dync1h1-cegat", "exoc4-jlf", "exoc4-cegat"}
    if record["id"] in matches:
        assert sequences == [record["native_sequence"]]
        fragment = selected[0].mutant_protein_fragment
        assert [fragment.mutant_amino_acid_start_offset,
                fragment.mutant_amino_acid_end_offset] == record["native_mutation_interval"]
        assert fragment.variant == result.variant
    else:
        assert len(sequences) == 1
        assert sequences[0] != record["native_sequence"]
        assert sequences[0] in result.top_protein_sequence.amino_acids


def test_expanded_map2_evidence_does_not_enable_dna_fallback(inputs):
    result, config = reconstruct_selection(inputs, "MAP2-chr2-209694768")
    assert result.num_alt_reads == 1
    assert result.top_protein_sequence is None
    predictor = Mock(side_effect=AssertionError("No protein should reach prediction"))
    assert vaccine_peptides_for_variant(result, predictor, vaccine_config=config) == []
    assert not predictor.mock_calls


def test_two_dync1h1_loci_are_not_collapsed_by_gene_name(inputs):
    first, _ = reconstruct_selection(inputs, "DYNC1H1-chr14-101980529")
    second, _ = reconstruct_selection(inputs, "DYNC1H1-chr14-102030200")
    assert first.variant != second.variant
    assert first.top_protein_sequence.amino_acids != second.top_protein_sequence.amino_acids
    assert first.num_alt_reads > 1 and second.num_alt_reads > 1
