"""Vaxrank exports retained source products without rerunning selection."""

import importlib

from isovar import __version__ as isovar_version
from isovar import from_protein_hypotheses
from orfdata import ProducerInfo, ProteinSequence, Proteome, ReferenceContext
import pytest
from varcode import Variant

from vaxrank import to_proteome
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.vaxrank_results import VaxrankResults
from tests.test_isovar_support_contract import (
    MUTANT,
    reconstruct,
    reference as support_reference,
)


REFERENCE = ReferenceContext(9606, "human", "GRCh38", "Ensembl", "115")


@pytest.fixture(params=["+", "-"])
def assembly_reference(tmp_path, request):
    return support_reference.__wrapped__(tmp_path, request)


def _protein(hypothesis_id, rank, amino_acids):
    return {
        "hypothesis_id": hypothesis_id,
        "isovar_rank": rank,
        "amino_acids": amino_acids,
        "n_terminus": "annotated_start_codon" if rank == 1 else "local_window",
        "c_terminus": "stop_codon" if rank == 1 else "no_stop_in_context",
        "gene_ids": ["ENSG1"],
        "gene_names": ["GENE A/1"],
        "transcript_ids": [f"ENST{rank}"],
        "rna_support": {},
        "translations": [],
    }


def _isovar_proteome():
    export = {
        "schema": "isovar.protein_hypotheses.v2",
        "isovar_version": isovar_version,
        "sample_id": "tumor 1",
        "source": "rna / run 1",
        "evidence_scope": ["tumor 1", "rna / run 1"],
        "interval_convention": "zero_based_half_open",
        "events": [{
            "event_id": "variant-1",
            "variant": {
                "reference_name": "GRCh38",
                "contig": "1",
                "start": 100,
                "ref": "A",
                "alt": "T",
            },
            "allele_support": {},
            "protein_hypotheses_complete": True,
            "protein_hypotheses": [
                _protein("rna-top", 1, "MPEPTIDE"),
                _protein("rna-alt-1", 2, "MALTERNA"),
                _protein("rna-alt-2", 3, "MSECONDQ"),
            ],
        }],
        "evidence_sets": {},
    }
    return from_protein_hypotheses(export, reference=REFERENCE)


def _fragment(
    amino_acids="MPEPTIDE",
    *,
    variant=None,
    sequence_source="isovar_assembly",
    protein_start_offset=0,
    dna_effect_selection=None,
):
    return MutantProteinFragment(
        variant=variant or Variant("1", 100, "A", "T", "GRCh38"),
        gene_name="GENE A/1",
        amino_acids=amino_acids,
        mutant_amino_acid_start_offset=1,
        mutant_amino_acid_end_offset=2,
        supporting_reference_transcripts=[],
        n_overlapping_reads=10,
        n_alt_reads=6,
        n_ref_reads=4,
        n_alt_reads_supporting_protein_sequence=5,
        sequence_source=sequence_source,
        sequence_source_version=(
            "10.19.0" if sequence_source == "varcode_translation"
            else isovar_version if sequence_source == "isovar_assembly"
            else ""
        ),
        protein_start_offset=protein_start_offset,
        dna_effect_selection=dna_effect_selection,
    )


def _results(*fragments):
    # The adapter under test delegates this marker to the patched, real Isovar
    # proteome below. The full native integration is covered in Isovar itself.
    return VaxrankResults(
        isovar_results=[{"event": "variant-1"}],
        variant_to_vaccine_peptides_dict={},
        ranked_vaccine_peptides=[],
        retained_protein_fragments=list(fragments),
    )


def _convert(monkeypatch, results):
    upstream = _isovar_proteome()
    module = importlib.import_module("vaxrank.proteome")
    monkeypatch.setattr(module, "isovar_to_proteome", lambda *args, **kwargs: upstream)
    return upstream, to_proteome(
        results,
        sample_id="tumor 1",
        source="rna / run 1",
        reference=REFERENCE,
    )


def test_joint_export_keeps_unranked_rna_alternatives_and_dna_window(monkeypatch, tmp_path):
    rna = _fragment()
    dna = _fragment(
        "MDNAWINQ",
        variant=Variant("2", 200, "G", "C", "GRCh38"),
        sequence_source="varcode_translation",
        protein_start_offset=42,
        dna_effect_selection={"transcript_id": "ENST-DNA", "mutation_start": 43},
    )
    _, patient = _convert(monkeypatch, _results(rna, dna))

    reference = Proteome(REFERENCE, [
        ProteinSequence("REF1", "MPEPTIDE", "annotated", 9606, "Ensembl:115"),
    ])
    ncorf = ProteinSequence.from_producer(
        "MNCORFQQ",
        producer=ProducerInfo(
            origin="nc_orf",
            method="catalogue",
            tool="catalogue",
            native_id="ORF1",
            reference=REFERENCE,
        ),
    )
    joint = reference.merge(Proteome.from_proteins(ncorf, reference=REFERENCE), patient)

    assert [protein.sequence for protein in patient.proteins[:3]] == [
        "MPEPTIDE", "MALTERNA", "MSECONDQ"]
    assert len(patient.proteins) == 4
    assert {protein.protein_source for protein in joint.proteins} == {
        "reference_cds", "nc_orf"
    }
    dna_protein = next(protein for protein in patient.proteins if protein.sequence == "MDNAWINQ")
    assert dna_protein.producer.tool == "varcode"
    assert dna_protein.producer.exporters == ("vaxrank",)
    fragment_record = next(
        record for record in patient.records
        if record.namespace == "vaxrank" and record.kind == "fragment"
        and record.payload["sequence_source"] == "varcode_translation"
    )
    assert fragment_record.payload["protein_interval"] == (42, 50)
    assert fragment_record.payload["mutation_interval"] == (1, 2)

    path = tmp_path / "joint.proteome.json.gz"
    joint.save(path)
    loaded = Proteome.load(path)
    assert loaded.identity == joint.identity
    fasta = tmp_path / "joint.fasta"
    loaded.write_fasta(fasta)
    header = next(line for line in fasta.read_text().splitlines() if "varcode" in line)
    assert "version=10.19.0" in header
    assert "start_complete=false" in header
    assert "stop_complete=unknown" in header
    assert "GENE%20A%2F1" in header.split()[0]
    deduplicated = joint.fasta_index(deduplicate=True)
    assert any(
        len(group) == 2 and group[0].sequence == "MPEPTIDE"
        for group in deduplicated.values()
    )


def test_direct_and_vaxrank_wrapped_isovar_share_source_identity(monkeypatch):
    direct, wrapped = _convert(monkeypatch, _results(_fragment()))
    merged = direct.merge(wrapped)

    assert len(direct.proteins) == len(wrapped.proteins) == len(merged.proteins)
    top = next(protein for protein in merged.proteins if protein.sequence == "MPEPTIDE")
    assert top.producer.exporters == ("isovar", "vaxrank")
    hypothesis = next(
        record for record in merged.records
        if record.namespace == "isovar" and record.kind == "hypothesis"
        and record.payload["isovar_rank"] == 1
    )
    assert hypothesis.exporters == ("isovar", "vaxrank")
    assert wrapped.merge(wrapped).identity == wrapped.identity
    retained = [
        link for link in wrapped.links
        if link.relation == "retains_fragment" and link.target == top.accession
    ]
    assert len(retained) == 1


def test_real_isovar_results_merge_direct_and_wrapped(
    tmp_path, assembly_reference
):
    result = reconstruct(tmp_path, assembly_reference, [
        (MUTANT, [(len(MUTANT), "M")], 0),
        (MUTANT, [(len(MUTANT), "M")], 0),
    ])
    fragment = MutantProteinFragment.from_isovar_result(result)
    reference = ReferenceContext(
        9606,
        "human",
        result.variant.reference_name,
        "test",
        "1",
    )
    results = VaxrankResults([result], {result.variant: []}, [(result.variant, [])], [fragment])
    direct = importlib.import_module("isovar").to_proteome(
        [result],
        sample_id="sample",
        source="synthetic-bam",
        reference=reference,
    )
    wrapped = to_proteome(
        results,
        sample_id="sample",
        source="synthetic-bam",
        reference=reference,
    )
    merged = direct.merge(wrapped)
    assert len(merged.proteins) == len(direct.proteins)
    assert all(
        protein.producer.exporters == ("isovar", "vaxrank")
        for protein in merged.proteins
    )
    path = tmp_path / "wrapped.proteome.json.gz"
    wrapped.save(path)
    assert Proteome.load(path).identity == wrapped.identity


def test_saved_rna_results_replay_all_isovar_evidence_offline(
    tmp_path, assembly_reference, monkeypatch
):
    result = reconstruct(tmp_path, assembly_reference, [
        (MUTANT, [(len(MUTANT), "M")], 0),
        (MUTANT, [(len(MUTANT), "M")], 0),
    ])
    fragment = MutantProteinFragment.from_isovar_result(result)
    reference = ReferenceContext(9606, "human", result.variant.reference_name, "test", "1")
    results = VaxrankResults([result], {result.variant: []}, [(result.variant, [])], [fragment])
    native = to_proteome(results, sample_id="sample", source="synthetic-bam", reference=reference)

    with pytest.raises(ValueError, match="sample_id and source"):
        results.to_json()
    saved = results.to_json(sample_id="sample", source="synthetic-bam")
    restored = VaxrankResults.from_json(saved)
    assert restored.isovar_results == []
    assert restored.isovar_protein_hypotheses["schema"] == "isovar.protein_hypotheses.v2"
    with pytest.raises(ValueError, match="live variant"):
        _ = restored.variants
    with pytest.raises(ValueError, match="live variant"):
        restored.variant_counts()
    monkeypatch.setattr(
        "vaxrank.proteome.isovar_to_proteome",
        lambda *args, **kwargs: pytest.fail("replay must use saved Isovar hypotheses"),
    )
    monkeypatch.setattr(
        "pyensembl.Genome.transcripts",
        lambda *args, **kwargs: pytest.fail("replay must not load genome transcripts"),
    )
    replay = to_proteome(restored, sample_id="sample", source="synthetic-bam", reference=reference)
    assert replay.identity == native.identity
    assert replay.proteins == native.proteins
    assert replay.records == native.records
    assert replay.links == native.links
    with pytest.raises(ValueError, match="sample_id and source"):
        to_proteome(restored, sample_id="other", source="synthetic-bam", reference=reference)


def test_saved_results_convert_identically_without_prediction(monkeypatch):
    fragment = _fragment(
        "MDNAWINQ",
        variant=Variant("2", 200, "G", "C", "GRCh38"),
        sequence_source="varcode_translation",
        protein_start_offset=42,
        dna_effect_selection={"transcript_id": "ENST-DNA", "mutation_start": 43},
    )
    results = VaxrankResults(
        [],
        {fragment.variant: []},
        [(fragment.variant, [])],
        [fragment],
    )
    restored = VaxrankResults.from_json(results.to_json())
    monkeypatch.setattr(
        "vaxrank.core_logic.predict_epitopes",
        lambda *args, **kwargs: pytest.fail("export must not run prediction"),
    )
    native = to_proteome(
        results, sample_id="tumor-1", source="run-1", reference=REFERENCE
    )
    offline = to_proteome(
        restored, sample_id="tumor-1", source="run-1", reference=REFERENCE
    )
    assert offline.identity == native.identity
    export = next(
        record for record in native.records
        if record.namespace == "vaxrank" and record.kind == "export"
    )
    assert export.payload["schema"] == "vaxrank.results.v1"
    assert export.payload["fragment_record_ids"]
    assert export.payload["ranked_variant_count"] == 1
    assert export.payload["selection"]["ranked_vaccine_peptides"]


def test_unknown_fragment_stays_unknown_and_requires_explicit_policy():
    fragment = _fragment(
        "MUNKNOWN",
        # A configured host gives Varcode's placeholder a reference name, but
        # the supplied fragment itself made no source-reference claim.
        variant=Variant("3", 300, "C", "A", "GRCh38"),
        sequence_source="caller_supplied",
        protein_start_offset=None,
    )
    patient = to_proteome(
        VaxrankResults([], {}, [], [fragment]),
        sample_id="tumor-1",
        source="supplied",
        reference=REFERENCE,
    )
    protein, = patient.proteins
    assert protein.taxon_id is None
    assert protein.producer.reference is None
    assert patient.reference_policy == "allow_unknown"
    assert "taxon=unknown" in protein.fasta_header()

    strict = Proteome(REFERENCE, ())
    with pytest.raises(ValueError, match="explicit compatibility"):
        strict.merge(Proteome.from_proteins(protein, compatibility="mixed"))
    assert strict.merge(
        Proteome.from_proteins(protein, compatibility="mixed"),
        compatibility="allow_unknown",
    ).proteins == (protein,)


def test_conflicting_fragment_content_under_one_source_identity_fails():
    first = _fragment(
        "MDNAWINQ",
        variant=Variant("2", 200, "G", "C", "GRCh38"),
        sequence_source="varcode_translation",
        protein_start_offset=42,
        dna_effect_selection={"transcript_id": "ENST-DNA", "mutation_start": 43},
    )
    changed = _fragment(
        "MCHANGED",
        variant=Variant("2", 200, "G", "C", "GRCh38"),
        sequence_source="varcode_translation",
        protein_start_offset=42,
        dna_effect_selection={"transcript_id": "ENST-DNA", "mutation_start": 43},
    )
    left = to_proteome(
        VaxrankResults([], {}, [], [first]),
        sample_id="tumor-1",
        source="run-1",
        reference=REFERENCE,
    )
    right = to_proteome(
        VaxrankResults([], {}, [], [changed]),
        sample_id="tumor-1",
        source="run-1",
        reference=REFERENCE,
    )
    with pytest.raises(ValueError, match="Conflicting protein identity"):
        left.merge(right)


def test_adapter_rejects_a_known_wrong_fragment_reference():
    fragment = _fragment(
        "MDNAWINQ",
        variant=Variant("2", 200, "G", "C", "GRCh37"),
        sequence_source="varcode_translation",
        protein_start_offset=42,
        dna_effect_selection={"transcript_id": "ENST-DNA", "mutation_start": 43},
    )
    with pytest.raises(ValueError, match="conflicts with"):
        to_proteome(
            VaxrankResults([], {}, [], [fragment]),
            sample_id="tumor-1",
            source="run-1",
            reference=REFERENCE,
        )
