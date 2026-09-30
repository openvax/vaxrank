"""Annotated reference membership is distinct from expression or biotype."""

import gzip
import json
import pickle

import pytest

from vaxrank.native_serialization import from_native_json, to_native_json
from vaxrank.reference_proteome import (
    ReferenceProteome, clear_reference_proteome_caches, ensembl_dataset_cache_identity,
    genome_protein_dict, load_kmer_set_index, self_reference_matches,
)
from vaxrank.vaccine_antigen import (
    AminoAcidInterval, SelfReferenceSource, TargetableMask,
    TumorSpecificityAttestation, VaccineAntigen,
)
from .test_reference_proteome import (
    _standard_release_with_local_sources, create_mock_genome, create_mock_transcript,
)


def antigen(excluded=()):
    return VaccineAntigen(
        kind="CTA", amino_acids="ACDEFGHIKL",
        targetable_mask=TargetableMask((AminoAcidInterval(0, 10),)),
        tumor_specificity=TumorSpecificityAttestation(
            status="admitted", evidence_kind="test", evidence_source="offline fixture"),
        self_reference_excluded_gene_ids=excluded)


@pytest.mark.parametrize("biotype", [
    "nonsense_mediated_decay", "non_stop_decay", "IG_V_gene", "TR_C_gene",
    "protein_coding_LoF", "translated_processed_pseudogene",
])
def test_noncanonical_reference_source_retained_with_its_biotype(tmp_path, monkeypatch, biotype):
    monkeypatch.setenv("VAXRANK_REF_PEPTIDES_DIR", str(tmp_path))
    clear_reference_proteome_caches()
    source = create_mock_transcript(
        "T1", "ACDEFGHIKL", gene_id="G1", is_protein_coding=False, biotype=biotype)
    genome = create_mock_genome([source])
    assert ReferenceProteome(genome, 8, 8).contains("ACDEFGHI")
    assert ReferenceProteome.from_genome(
        genome, exclude_gene_ids={"OTHER"}, min_kmer_length=8, max_kmer_length=8
    ).contains("ACDEFGHI")
    match = self_reference_matches(["ACDEFGHI"], antigen(), genome)["ACDEFGHI"]
    assert match.occurs and match.source_provenance_complete
    assert match.sources[0].transcript_biotype == biotype
    assert from_native_json(to_native_json(match), type(match)) == match


def test_cta_exclusion_does_not_erase_noncanonical_non_cta_source():
    genome = create_mock_genome([
        create_mock_transcript("CTA", "ACDEFGHIKL", gene_id="GCTA"),
        create_mock_transcript("NMD", "ACDEFGHIKL", gene_id="GNORMAL", is_protein_coding=False),
        create_mock_transcript("NO_PROTEIN", None, gene_id="GEMPTY", is_protein_coding=False),
    ])
    match = self_reference_matches(["ACDEFGHI"], antigen(("GCTA",)), genome)["ACDEFGHI"]
    assert match.occurs
    assert [(s.gene_id, s.transcript_id) for s in match.sources] == [("GNORMAL", "NMD")]
    assert genome_protein_dict(genome, {"GCTA"}) == {"NMD": "ACDEFGHIKL"}


def test_old_protein_coding_only_disk_cache_cannot_hide_a_translation(tmp_path, monkeypatch):
    monkeypatch.setenv("VAXRANK_REF_PEPTIDES_DIR", str(tmp_path))
    clear_reference_proteome_caches()
    genome, _ = _standard_release_with_local_sources(tmp_path / "reference")
    genome.transcripts = lambda: [create_mock_transcript("NMD", "ACDEFGHIKL", is_protein_coding=False)]
    identity = ensembl_dataset_cache_identity(genome)
    with gzip.open(tmp_path / f"content_{identity}_kmer_set_8_8.pkl.gz", "wb") as output:
        pickle.dump(set(), output)
    assert "ACDEFGHI" in load_kmer_set_index(genome, 8, 8)


def test_legacy_source_without_biotype_remains_unknown():
    payload = json.loads(to_native_json(SelfReferenceSource("G1")))
    payload.pop("transcript_biotype")
    source = from_native_json(json.dumps(payload), SelfReferenceSource)
    assert source.transcript_biotype == ""
