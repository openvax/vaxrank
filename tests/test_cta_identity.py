"""Consumer regressions using exact full Ensembl PRAME/PRAMEF10 records."""

from dataclasses import replace
import json
from pathlib import Path
from types import SimpleNamespace

import oncoref
import pytest

from vaxrank.cta_admission import CTAAdmissionError, resolve_cta_reference_evidence
from vaxrank.cta_expression import CTAExpressionResult
from vaxrank.cta_identity import snapshot_cta_annotation
from vaxrank.reference_proteome import ReferenceProteome, cta_source_gene_ids_for_genome, self_reference_matches
from .test_cta_expression import run

PRAME = "ENSG00000185686"
ALTERNATE = "ENSG00000275013"
PRAMEF10 = "ENSG00000187545"


class Annotation:
    """Exact fixture loci/proteins; the fixture retains upstream file hashes."""

    species = SimpleNamespace(latin_name="homo_sapiens")
    reference_name = "GRCh38"
    annotation_name = "ensembl"

    def __init__(self, release, primary_only=False, fixture_name=None):
        self.release = self.annotation_version = release
        fixture_name = fixture_name or f"prame-ensembl-{release}.json"
        self.source = json.loads((Path(__file__).parent / "data" / "cta_identity"
                                  / fixture_name).read_text())
        self.genes, self.proteins = {}, {}
        for row in self.source["genes"]:
            if primary_only and row["gene_id"] == ALTERNATE:
                continue
            transcripts = [SimpleNamespace(**t, gene_id=row["gene_id"],
                           gene_name=row["gene_name"], transcript_version=None,
                           biotype="protein_coding") for t in row["transcripts"]]
            self.proteins.update({t.transcript_id: t for t in transcripts})
            self.genes[row["gene_id"]] = SimpleNamespace(**dict(row, transcripts=transcripts))

    def gene_ids(self):
        return list(self.genes)

    def gene_by_id(self, gid):
        if gid not in self.genes:
            raise ValueError(gid)
        return self.genes[gid]

    def transcript_ids(self):
        return list(self.proteins)

    def transcript_by_id(self, tid):
        return self.proteins[tid]

    def transcripts(self):
        return list(self.proteins.values())


@pytest.mark.parametrize("release", [93, 112])
@pytest.mark.parametrize("primary_only", [True, False])
def test_gene_admission_keeps_actual_sources_and_distinct_membership_hashes(tmp_path, release, primary_only):
    genome = Annotation(release, primary_only)
    result = run(tmp_path, genome, [(PRAME, 3), (ALTERNATE + ".2", 7)])
    assert [d.status for d in result.decisions] == ["admitted", "unverified_identity" if primary_only else "admitted"]
    reference = result.annotation_reference
    assert reference.identity_contract_version == 1
    assert len(reference.canonical_candidate_gene_ids) == 2532
    assert set(reference.self_reference_excluded_gene_ids) == ({PRAME} if primary_only else {PRAME, ALTERNATE})
    assert reference.canonical_candidate_gene_ids_sha256 != reference.self_reference_excluded_gene_ids_sha256
    assert PRAMEF10 not in reference.self_reference_excluded_gene_ids
    assert cta_source_gene_ids_for_genome(genome) == frozenset(reference.self_reference_excluded_gene_ids)
    if primary_only:
        assert result.decisions[1].gene_identity["status"] == "absent"
        return
    decision = result.decisions[1]
    antigen = decision.assessment.antigen
    assert antigen.gene_id == ALTERNATE
    assert decision.assessment.tumor_expression.input_identifier == ALTERNATE + ".2"
    assert decision.assessment.reference_evidence.gene_id == PRAME
    assert decision.reference_resolution.gene_identity["input_gene_id"] == ALTERNATE + ".2"
    assert dict(antigen.source_metadata)["sequence_selection"] == "annotation_longest_protein"
    assert antigen.transcript_ids == ("ENST00000539862",)
    assert antigen.protein_ids == ("ENSP00000445097",)
    assert antigen.species == "Homo sapiens"
    assert [d.value for d in result.decisions] == [3, 7]
    assert result.admitted_antigens[0].source_identifier != antigen.source_identifier
    evidence = dict(antigen.tumor_specificity.evidence_records[0].details)
    assert json.loads(evidence["gene_identity"])["source_gene_id"] == ALTERNATE
    assert evidence["canonical_gene_id"] == PRAME
    assert evidence["source_exclusions_sha256"] == reference.self_reference_excluded_gene_ids_sha256


@pytest.mark.parametrize("release", [93, 112])
@pytest.mark.parametrize("tid,pid", [("ENST00000539862", "ENSP00000445097"),
                                    ("ENST00000617728", "ENSP00000484066")])
def test_measured_alternate_transcript_keeps_exact_translation(tmp_path, release, tid, pid):
    result = run(tmp_path, Annotation(release), [(tid, 6)], level="transcript")
    antigen = result.admitted_antigens[0]
    assert antigen.gene_id == ALTERNATE
    assert antigen.transcript_ids == (tid,)
    assert antigen.protein_ids == (pid,)
    assert result.decisions[0].assessment.tumor_expression.transcript_id == tid
    assert dict(antigen.source_metadata)["sequence_selection"] == "measured_transcript"


@pytest.mark.parametrize("release", [93, 112])
def test_exact_prame_exclusion_retains_real_non_cta_pramef_sources(tmp_path, release):
    genome = Annotation(release)
    antigen = run(tmp_path, genome, [(PRAME, 3)]).admitted_antigens[0]
    peptides = {antigen.amino_acids[i:i + 9] for i in range(len(antigen.amino_acids) - 8)}
    assert len(peptides) == 501
    matches = self_reference_matches(peptides, antigen, genome)
    # This pinned subset includes PRAMEF10 only. The full proteome contains 14
    # genuine shared 9-mers across additional PRAMEF genes, checked by smoke.
    assert sum(m.occurs for m in matches.values()) == 7
    assert {s.gene_id for m in matches.values() for s in m.sources} == {PRAMEF10}
    assert {s.transcript_id for m in matches.values() for s in m.sources} == {"ENST00000235347"}
    assert all(m.source_provenance_complete for m in matches.values())
    reference = ReferenceProteome.from_genome(genome, exclude_cta_genes=True,
                                            min_kmer_length=9, max_kmer_length=9)
    assert sum(reference.contains(p) for p in peptides) == 7


@pytest.mark.parametrize("status", ["annotation_conflict", "ambiguous", "unverified"])
def test_rejected_explicit_alias_stays_in_self_background(tmp_path, monkeypatch, status):
    genome = Annotation(112)
    if status == "annotation_conflict":
        genome.genes[ALTERNATE].contig = "22"
    else:
        # Exercise the consumer boundary with an unresolved public API record;
        # OncoRef owns validation of competing targets and unsupported methods.
        actual = oncoref.gene_identity.resolve_gene_identity
        def resolve(gid, *, genome):
            record = actual(gid, genome=genome)
            return replace(record, canonical_gene_id=None, status=status) if record.source_gene_id == ALTERNATE else record
        monkeypatch.setattr(oncoref.gene_identity, "resolve_gene_identity", resolve)
        monkeypatch.setattr(oncoref, "resolve_gene_identity", resolve)
    result = run(tmp_path, genome, [(ALTERNATE, 10), (PRAME, 3)])
    assert result.decisions[0].status == "unverified_identity"
    assert result.decisions[0].gene_identity["status"] == status
    assert result.annotation_reference.self_reference_excluded_gene_ids == (PRAME,)
    antigen = result.admitted_antigens[0]
    matches = self_reference_matches([antigen.amino_acids[:9]], antigen, genome)
    assert ALTERNATE in {s.gene_id for s in next(iter(matches.values())).sources}
    with pytest.raises(CTAAdmissionError, match="Unverified"):
        resolve_cta_reference_evidence(ALTERNATE, genome=genome)


def test_native_identity_replay_never_queries_current_references(tmp_path, monkeypatch):
    genome = Annotation(112)
    result = run(tmp_path, genome, [(ALTERNATE + ".2", 10)])
    path = tmp_path / "admission.json"
    result.save(path)
    (tmp_path / "expression.tsv").unlink()
    def forbidden(*args, **kwargs):
        raise AssertionError("Replay queried current identity evidence")
    monkeypatch.setattr(oncoref, "resolve_gene_identity", forbidden)
    monkeypatch.setattr(oncoref, "cta_annotation_gene_identities", forbidden)
    monkeypatch.setattr(oncoref.cta, "cta_unfiltered_gene_ids", forbidden)
    monkeypatch.setattr(genome, "gene_by_id", forbidden)
    monkeypatch.setattr(genome, "transcripts", forbidden)
    assert CTAExpressionResult.load(path) == result
    with pytest.raises(ValueError, match="disagree"):
        replace(result.annotation_reference, self_reference_excluded_gene_ids=(PRAME,))
    with pytest.raises(ValueError, match="disagree"):
        replace(result, annotation_reference=replace(result.annotation_reference, annotation_sha256="0" * 64))


def test_identity_contract_changes_reference_cache_identity(monkeypatch):
    from vaxrank.reference_proteome import _kmer_dataset_identity
    genome = Annotation(112)
    before = _kmer_dataset_identity(genome)
    monkeypatch.setattr(oncoref.gene_identity, "GENE_IDENTITY_CONTRACT_VERSION", 2)
    assert _kmer_dataset_identity(genome) != before


def test_same_symbol_identical_protein_never_creates_source_alias(tmp_path):
    genome = Annotation(112)
    unknown = "ENSG99999999999"
    source = genome.genes[ALTERNATE]
    translation = SimpleNamespace(**dict(vars(source.transcripts[0]),
                                  gene_id=unknown, transcript_id="ENST99999999999"))
    genome.genes[unknown] = SimpleNamespace(**dict(vars(source), gene_id=unknown, transcripts=[translation]))
    genome.proteins[translation.transcript_id] = translation
    result = run(tmp_path, genome, [(unknown, 10), (PRAME, 3)])
    assert result.decisions[0].status == "not_cta"
    assert result.decisions[0].gene_identity["status"] == "unmapped"
    antigen = result.admitted_antigens[0]
    assert unknown not in antigen.self_reference_excluded_gene_ids
    assert unknown in {s.gene_id for m in self_reference_matches(
        [antigen.amino_acids[:9]], antigen, genome).values() for s in m.sources}


def test_different_annotation_content_cannot_reuse_cta_mapping_snapshot():
    genome = Annotation(112)
    snapshot = snapshot_cta_annotation(genome)
    genome.genes[ALTERNATE].start += 1
    with pytest.raises(CTAAdmissionError, match="mapping changed"):
        resolve_cta_reference_evidence(ALTERNATE, genome=genome, annotation_reference=snapshot)


def test_explicit_annotation_required():
    with pytest.raises(TypeError, match="genome"):
        resolve_cta_reference_evidence(PRAME)
    with pytest.raises(ValueError, match="explicit"):
        snapshot_cta_annotation(SimpleNamespace())


def test_old_native_admission_remains_explicitly_legacy(tmp_path, monkeypatch):
    from vaxrank.cta_admission import CTAReferenceResolution
    from vaxrank.native_serialization import to_native_json
    result = run(tmp_path, Annotation(112, primary_only=True), [(PRAME, 3)])
    decision = result.decisions[0]
    candidates = result.annotation_reference.canonical_candidate_gene_ids
    legacy_resolution = CTAReferenceResolution(decision.reference_resolution.evidence, candidates)
    antigen = decision.assessment.antigen
    membership = antigen.tumor_specificity.evidence_records[0]
    old_details = tuple((k, v) for k, v in membership.details if k in {
        "specificity_status", "specificity_action", "restriction", "restriction_confidence", "source_row_sha256"})
    antigen = replace(antigen, self_reference_excluded_gene_ids=candidates,
                      tumor_specificity=replace(antigen.tumor_specificity, evidence_records=(
                          replace(membership, details=old_details), antigen.tumor_specificity.evidence_records[1])))
    assessment = replace(decision.assessment, antigen=antigen, reference_resolution=None)
    legacy = replace(result, annotation_reference=None, decisions=(replace(
        decision, assessment=assessment, reference_resolution=legacy_resolution, gene_identity=None),))
    payload = json.loads(to_native_json(legacy))
    added = {
        "CTAExpressionResult": ["annotation_reference"],
        "CTAExpressionDecision": ["gene_identity"],
        "CTAAdmissionAssessment": ["reference_resolution"],
        "CTAReferenceResolution": ["identity_contract_version", "gene_identity", "annotation_reference_sha256",
                                   "annotation_sha256", "self_reference_excluded_gene_ids_sha256"],
    }
    def old_schema(value):
        if isinstance(value, list):
            for item in value:
                old_schema(item)
        elif isinstance(value, dict):
            for key in added.get(value.get("__class__", {}).get("__name__"), []):
                value.pop(key, None)
            for item in value.values():
                old_schema(item)
    old_schema(payload)
    path = tmp_path / "old-admission.json"
    path.write_text(json.dumps(payload))
    monkeypatch.setattr(oncoref, "resolve_gene_identity", lambda *a, **k: pytest.fail("Legacy replay queried OncoRef"))
    restored = CTAExpressionResult.load(path)
    assert restored.annotation_reference is None
    assert restored.decisions[0].reference_resolution.identity_contract_version == 0
    assert restored.decisions[0].assessment.reference_resolution is None
    assert "gene_identity" not in dict(restored.admitted_antigens[0].tumor_specificity.evidence_records[0].details)
    assert restored.admitted_antigens[0].self_reference_excluded_gene_ids == candidates
