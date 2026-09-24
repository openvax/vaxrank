"""Isovar reconstruction -> real fusion antigen, with no invented variant."""

from dataclasses import replace
import gzip
from hashlib import sha256
import json
from vaxrank.sid_test_data import sid_test_data

import pytest
from isovar import FusionBlock, FusionBreakpoint, FusionRead, FusionReference, FusionTranscript, reconstruct_fusion
from isovar.fusion import FUSION_INPUT_KEYS, fusion_from_dict
from mhctools import RandomBindingPredictor

from vaxrank import EpitopeConfig, IsovarFusionAntigens, VaccinePeptide, fusion_antigens_from_isovar, predict_epitopes
from vaxrank.vaccine_antigen import TumorSpecificityAttestation


def example(cut=24, insert="", partial=0):
    donor = "ATG" + "GCT" * 12 + "TAA"
    acceptor = "ATG" + "GCC" * 14 + "TAA"
    sequence = donor[partial:cut] + insert + acceptor[3:]
    j0, j1 = cut - partial, cut - partial + len(insert)
    blocks = (FusionBlock(0, j0, "d", 100 + partial, 100 + cut, "+"),
              FusionBlock(j1, len(sequence), "a", 203, 200 + len(acceptor), "+"))
    fusion = FusionTranscript("donor--acceptor", "test", sequence, j0, j1,
        FusionBreakpoint("d", 100 + cut, "+"), FusionBreakpoint("a", 203, "+"), blocks,
        dict(sample_id="sample", method="test assembly", version="1", parameters={}, source="synthetic", contig_id="contig"))
    references = tuple(FusionReference(name, "test", "test annotation", contig, "+",
        ((offset, offset + len(seq)),), seq, 0, len(seq))
        for name, contig, offset, seq in (("donor", "d", 100, donor), ("acceptor", "a", 200, acceptor)))
    reads = tuple(FusionRead("sample", "library", "f" + str(i), "r" + str(i), "synthetic", 0, 0, sequence, blocks)
                  for i in range(2))
    return fusion, references, reads


def admitted():
    return TumorSpecificityAttestation(status="admitted", evidence_kind="independent_somatic_evidence",
        evidence_source="synthetic test only", patient_specific=True, rationale_code="test")


@pytest.mark.parametrize("cut", [24, 25, 26])
@pytest.mark.parametrize("insert", ["", "G", "GG", "GCT" * 8])
def test_exact_junction_masks_including_split_codons_and_insert_only_peptides(cut, insert):
    result = reconstruct_fusion(*example(cut, insert), peptide_lengths=range(1, 12))
    adapted = fusion_antigens_from_isovar(result, tumor_specificity=admitted(), gene_name="D::A", species="test")
    antigen, = adapted.admitted_antigens
    protein, = result["paths"][0]["translations"]
    assert antigen.transcript_ids == ("acceptor", "donor")
    assert antigen.self_reference_excluded_gene_ids == ()
    for length in range(1, 12):
        for start in range(len(antigen.amino_acids) - length + 1):
            end = start + length
            expected = any(3 * start < b < 3 * end for b in (cut, cut + len(insert)))
            assert antigen.interval_is_targetable(start, end) == expected
    assert {(p["protein_interval"][0], p["protein_interval"][1]) for p in protein["candidate_peptides"]} == {
        (start, start + length) for length in range(1, 12)
        for start in range(len(antigen.amino_acids) - length + 1)
        if antigen.interval_is_targetable(start, start + length)}
    assert json.loads(dict(antigen.source_metadata)["isovar_fusion_result"]) == json.loads(json.dumps(result))
    assert IsovarFusionAntigens.from_json(adapted.to_json()) == adapted
    result["provenance"]["source"] = "changed after adaptation"
    assert adapted.reconstruction["provenance"]["source"] == "synthetic"


def test_partial_and_ambiguous_hypotheses_stay_held_out_even_with_attestation():
    fusion, refs, reads = example(partial=3)
    alternative = FusionReference("other_frame", "test", "test annotation", "d", "+",
        ((0, 4), (103, 124), (300, 305)), "ATGA" + refs[0].sequence[3:24] + "CCTAA", 0, 30)
    for models, status, count in [(refs, "translated", 1), (refs + (alternative,), "ambiguous", 2)]:
        result = reconstruct_fusion(fusion, models, reads)
        assert result["status"] == status
        adapted = fusion_antigens_from_isovar(result, tumor_specificity=admitted())
        assert len(adapted.antigens) == count and not adapted.admitted_antigens
        for antigen in adapted.antigens:
            assert antigen.tumor_specificity.requires_review
            assert json.loads(dict(antigen.source_metadata)["requested_tumor_specificity"])["status"] == "admitted"
            with pytest.raises(ValueError, match="held-out antigen"):
                VaccinePeptide(antigen=antigen)
    fusion, refs, reads = example()
    noncoding = replace(refs[0], transcript_id="noncoding", cds_start=None, cds_end=None)
    result = reconstruct_fusion(fusion, refs + (noncoding,), reads)
    adapted = fusion_antigens_from_isovar(result, tumor_specificity=admitted())
    assert adapted.reconstruction["reasons"] == ["donor_CDS_unavailable:noncoding"]
    assert len(adapted.antigens) == 1 and not adapted.admitted_antigens


def test_no_attestation_or_insufficient_support_never_admits():
    fusion, refs, reads = example()
    adapted = fusion_antigens_from_isovar(reconstruct_fusion(fusion, refs, reads))
    assert len(adapted.antigens) == 1 and not adapted.admitted_antigens
    result = reconstruct_fusion(fusion, refs, reads[:1])
    adapted = fusion_antigens_from_isovar(result, tumor_specificity=admitted())
    assert adapted.antigens == ()
    assert adapted.reconstruction["status"] == "insufficient_support"
    assert adapted.reconstruction["reasons"] == result["reasons"]


@pytest.mark.parametrize("name,checksum,status,count", [
    ("ATP5MG--KMT2A", "6c4974ae83898705f647f15f2dcdc9e86a913b9baf5577d2326720620df938a1", "ambiguous", 1),
    ("TPST1--CRCP-T1", "aae3eca08578f1ee36f1e74043acfee2d885d13cf7594b83de5e18b1d3739e74", "unresolved", 0),
])
def test_original_sid_rna_preserves_coding_uncertainty_and_provenance(name, checksum, status, count):
    raw = (sid_test_data() / "isovar_fusions" / (name + ".input.json.gz")).read_bytes()
    data = json.loads(gzip.decompress(raw))
    # Pin original content independently of gzip/Python transport details.
    assert sha256(json.dumps(data, sort_keys=True, separators=(",", ":")).encode()).hexdigest() == checksum
    result = reconstruct_fusion(*fusion_from_dict({k: data[k] for k in FUSION_INPUT_KEYS if k in data}))
    adapted = fusion_antigens_from_isovar(result, tumor_specificity=admitted(), species="Homo sapiens")
    assert adapted.reconstruction == json.loads(json.dumps(result)) and result["status"] == status
    assert len(adapted.antigens) == count and not adapted.admitted_antigens
    assert IsovarFusionAntigens.from_json(adapted.to_json()) == adapted
    if count:
        antigen, = adapted.antigens
        assert antigen.amino_acids == "MAQFVRNLVEKTPALVNG"
        assert antigen.interval_is_targetable(17, 18)  # Mixed junction codon, not a rounded boundary.
        assert not antigen.interval_is_targetable(0, 17)
        assert result["evidence"]["direct_fragments"] == 2
        assert set(result["paths"][0]["compatible_transcripts"]["acceptor"]) <= set(antigen.transcript_ids)


@pytest.mark.parametrize("damage,match", [
    (lambda r: r.update(schema="isovar.fusion_rna.v1"), "schema"),
    (lambda r: r["paths"][0].update(sequence_sha256="wrong"), "checksum"),
    (lambda r: r.update(reasons=["unresolved"]), "unambiguous"),
    (lambda r: r["paths"][0]["translations"][0].update(cds_start=None), "CDS start"),
    (lambda r: r["paths"][0]["translations"][0].update(amino_acids="WRONG"), "protein disagrees"),
    (lambda r: r["paths"][0]["translations"][0].update(junction_in_translated_cds=[1, 1]), "junction"),
])
def test_inconsistent_result_cannot_be_admitted(damage, match):
    result = reconstruct_fusion(*example())
    damage(result)
    with pytest.raises(ValueError, match=match):
        fusion_antigens_from_isovar(result, tumor_specificity=admitted())


def test_downstream_prediction_and_construct_keep_fusion_identity_and_evidence():
    adapted = fusion_antigens_from_isovar(reconstruct_fusion(*example(cut=25)), tumor_specificity=admitted())
    antigen, = adapted.admitted_antigens
    epitopes = predict_epitopes(mhc_predictor=RandomBindingPredictor(["HLA-A*02:01"]),
        epitope_config=EpitopeConfig(min_epitope_score=0), antigen=antigen)
    assert epitopes and any(e.overlaps_targetable for e in epitopes)
    assert all(not e.overlaps_mutation and "wt" not in e.comparators for e in epitopes)
    assert all(e.overlaps_targetable == antigen.interval_is_targetable(e.offset, e.offset + len(e.sequence))
               for e in epitopes)
    peptide = VaccinePeptide(antigen=antigen, epitopes=epitopes,
        combined_score_expr="target_epitope_score", ranking_rules=("target_epitope_score",))
    restored = VaccinePeptide.from_json(peptide.to_json())
    assert restored.antigen == antigen and restored.mutant_protein_fragment is None
