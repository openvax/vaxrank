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
from varcode import Variant

from vaxrank.construct_sequence import ConstructEvidence, ConstructSequence, ConstructSequenceEdit
from vaxrank.core_logic import vaccine_peptides_for_variant
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.native_serialization import from_native_json, to_native_json
from vaxrank.vaccine_antigen import VaccineAntigen
from vaxrank.vaccine_config import VaccineConfig

from .osteosarc_selection_helpers import DATA, load_selection_inputs, reconstruct_selection
from .test_declared_dependencies import requirement_named


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


@pytest.fixture(scope="module")
def reconstruct(inputs):
    """Reconstruct each (variant, length) once for the whole module.

    Two parametrized tests need the same fragment for all seven documented
    records, and every reconstruction re-reads real BAM alignments and
    reassembles RNA context. Without sharing, the expensive half of this
    suite runs twice over for no additional coverage. Results are read-only
    here; no test mutates the isovar result or the config it returns.
    """
    cache = {}

    def reconstruct_once(variant_id, length=25):
        key = (variant_id, length)
        if key not in cache:
            cache[key] = reconstruct_selection(inputs, variant_id, length)
        return cache[key]

    return reconstruct_once


def cached_predictor():
    predictor = TopiaryPredictor(models=[CachedPredictor.from_topiary_output(
        str(DATA / "netmhcpan42.tsv"))])
    # Asserted on the instance the selection tests actually score with, not
    # only on a separate cache built inside the manifest test. A fallback
    # here would quietly reintroduce live prediction into a suite whose
    # whole design claim is that it never leaves the cache.
    assert all(model.fallback is None for model in predictor.models)
    return predictor


def test_regeneration_refuses_to_overwrite_existing_evidence(tmp_path, monkeypatch):
    from examples.osteosarc_test_data import generate_predictions

    reconstruction = Mock(side_effect=AssertionError("Existing output must fail before reconstruction"))
    monkeypatch.setattr(generate_predictions, "load_selection_inputs", reconstruction)
    with pytest.raises(ValueError, match="Choose a new output directory"):
        generate_predictions.main(tmp_path)
    reconstruction.assert_not_called()


def test_pinned_predictions_are_complete_real_model_outputs():
    metadata = json.loads((DATA / "predictions_manifest.json").read_text())
    for filename, key in [("documented.json", "documented_sha256"),
                          ("isovar/manifest.json", "bundled_input_manifest_sha256"),
                          ("netmhcpan42.tsv", "output_sha256")]:
        assert sha256((DATA / filename).read_bytes()).hexdigest() == metadata[key]
    # The cache records the versions that actually produced it; compare those
    # against the declared floors instead of one hardcoded value. Pinning
    # isovar alone both went stale and never covered mhctools or topiary, so
    # a cache generated below the topiary floor carrying the #296 cached-scan
    # fix, or below the mhctools floor, read as verified provenance.
    recorded = metadata["package_versions"]
    for name in ("isovar", "mhctools", "topiary", "varcode"):
        requirement = requirement_named(name)
        assert requirement.specifier.contains(recorded[name], prereleases=True), (
            "prediction cache was generated with %s %s, which violates the "
            "declared requirement %s; regenerate it against the declared "
            "versions rather than trusting its recorded provenance"
            % (name, recorded[name], requirement))
    assert metadata["predictor_name"] == "netMHCpan-4.2"
    assert metadata["predictor_version"] == "4.2c"
    # Compare the manifest to the JSON representation the generator writes.
    # Immutable config tuples necessarily round-trip through JSON as lists.
    expected_vaccine_config = json.loads(json.dumps(
        msgspec.to_builtins(VaccineConfig())
    ))
    assert metadata["vaccine_config"] == expected_vaccine_config
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
def test_documented_native_and_manufactured_sequences_stay_distinct(reconstruct, record):
    result, _ = reconstruct(record["variant_id"], len(record["native_sequence"]))
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
def test_real_cached_final_selections_against_documented_sequences(reconstruct, record, caplog):
    result, config = reconstruct(record["variant_id"], len(record["native_sequence"]))
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
    context = result.top_protein_sequence.amino_acids
    frame = predictor.predict_from_named_sequences({"rna": context})
    assert frame.columns.is_unique, "Upstream Topiary #296: duplicate cache occurrence coordinates"
    # Completeness for the peptides this test actually scores, not only for
    # the manifest's fixed list. A partially regenerated cache missing rows
    # for a window the scan produces would otherwise rank on incomplete MHC
    # evidence with every assertion below still passing.
    scanned = {context[offset:offset + length]
               for length in DOCUMENTED["prediction_peptide_lengths"]
               for offset in range(len(context) - length + 1)}
    covered = set(frame[["peptide", "allele", "kind"]].itertuples(index=False, name=None))
    missing = {(p, a, k) for p in scanned
               for a in DOCUMENTED["prediction_alleles"]
               for k in ("pMHC_affinity", "pMHC_presentation")} - covered
    assert not missing, (
        "cache is missing %d (peptide, allele, kind) rows for windows this "
        "scan produced, e.g. %s" % (len(missing), sorted(missing)[:3]))
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


@pytest.mark.parametrize("corrected", [False, True],
                         ids=["published-deletion", "corrected-complex-replacement"])
def test_map2_alleles_do_not_enable_dna_fallback(inputs, corrected):
    # Keep the published Isovar 1.8.1 allele/bytes intact. Osteosarc 0.1.0's
    # correction is a separate comparison on the same selected-read subset,
    # not evidence of absent expression in the full RNA sample. Primary vendor
    # sources and coordinate conventions are recorded in the fixture README.
    ref, alt = ("CCTGGGCTACTGTGTGTTCAATAAGTACACAGT", "CAGGG") if corrected else (
        "CCTGGGCTACTGTGTGTTCAATA", "C")
    expected_variant = Variant("2", 209694768, ref, alt, ensembl=inputs[0])
    result, config = reconstruct_selection(
        inputs, "MAP2-chr2-209694768",
        variant=expected_variant if corrected else None)
    assert result.variant == expected_variant
    assert (result.num_alt_reads, result.num_ref_reads, result.num_other_reads) == (1, 0, 0)
    assert (result.num_alt_fragments, result.num_ref_fragments, result.num_other_fragments) == (1, 0, 0)
    assert not result.passes_all_filters
    assert result.top_protein_sequence is None
    predictor = Mock(side_effect=AssertionError("No protein should reach prediction"))
    assert vaccine_peptides_for_variant(result, predictor, vaccine_config=config) == []
    assert not predictor.mock_calls
    assert "MAP2-chr2-209694768" not in {record["variant_id"] for record in RECORDS}


def test_two_dync1h1_loci_are_not_collapsed_by_gene_name(reconstruct):
    first, _ = reconstruct("DYNC1H1-chr14-101980529")
    second, _ = reconstruct("DYNC1H1-chr14-102030200")
    assert first.variant != second.variant
    assert first.top_protein_sequence.amino_acids != second.top_protein_sequence.amino_acids
    assert first.num_alt_reads > 1 and second.num_alt_reads > 1


@pytest.mark.parametrize("record", RECORDS, ids=lambda r: r["id"])
def test_wildtype_comparators_are_reference_aligned_or_absent(reconstruct, record):
    """No cached WT comparator may be a peptide the reference cannot support.

    The generator used to slice the reference at the mutant's own offsets.
    That holds for a substitution but not for an indel, which shifts every
    downstream reference position, so the slice returned a real but
    unrelated peptide and cached it as wild type. H1-2 is a real in-frame
    deletion and produced exactly that.
    """
    result, _ = reconstruct(record["variant_id"], len(record["native_sequence"]))
    if result.top_protein_sequence is None:
        pytest.skip("no RNA-backed protein for this record")
    fragment = MutantProteinFragment.from_isovar_result(result)
    effect = fragment.predicted_effect()
    reference = effect.original_protein_sequence
    start = fragment.global_start_pos()
    mutation_start = fragment.mutant_amino_acid_start_offset
    length_preserving = len(effect.aa_ref) == len(effect.aa_alt)
    for length in DOCUMENTED["prediction_peptide_lengths"]:
        for offset in range(len(fragment) - length + 1):
            comparator = fragment.wildtype_peptide_at(offset, length)
            if comparator is None:
                # Only an alignment-breaking variant may withhold one, and
                # only for windows reaching the mutation.
                assert not length_preserving
                assert offset + length > mutation_start
                continue
            assert len(comparator) == length
            # An emitted comparator is exactly the reference at the same
            # protein coordinates, never a shifted or re-searched window.
            assert comparator == reference[start + offset:start + offset + length]
            # Upstream of the mutation the fragment is unmutated reference,
            # so those residues must agree. Downstream residues may differ
            # at more than the annotated position when co-occurring variants
            # appear in the same reconstructed RNA context.
            upstream = min(length, max(0, mutation_start - offset))
            assert fragment.amino_acids[offset:offset + upstream] == comparator[:upstream]


def test_indel_context_emits_no_unverifiable_wildtype_comparator(reconstruct):
    """The H1-2 deletion must not reproduce the AAKPKVVKP class of error."""
    result, _ = reconstruct("H1_2-chr6-26055824", 16)
    fragment = MutantProteinFragment.from_isovar_result(result)
    reference = fragment.predicted_effect().original_protein_sequence
    start = fragment.global_start_pos()
    skipped = 0
    mutation_start = fragment.mutant_amino_acid_start_offset
    for length in DOCUMENTED["prediction_peptide_lengths"]:
        for offset in range(len(fragment) - length + 1):
            comparator = fragment.wildtype_peptide_at(offset, length)
            naive = reference[start + offset:start + offset + length]
            if comparator is None:
                skipped += 1
                # Withheld only where the deletion has shifted coordinates.
                assert offset + length > mutation_start
                continue
            assert comparator == naive
            # Everything emitted here sits wholly upstream of the deletion.
            assert offset + length <= mutation_start
    # The deletion genuinely breaks alignment for a large share of windows;
    # a fix that silently emitted everything would leave this at zero.
    assert skipped > 0
    # AAKPKVVKP is the reference straddling the deleted AAKPK, which the old
    # slice emitted as a wild-type comparator for this context even though
    # it covers residues no mutant window occupies.
    emitted = {fragment.wildtype_peptide_at(offset, length)
               for length in DOCUMENTED["prediction_peptide_lengths"]
               for offset in range(len(fragment) - length + 1)}
    emitted.discard(None)
    assert "AAKPKVVKP" not in emitted
    assert all(peptide in reference for peptide in emitted)
