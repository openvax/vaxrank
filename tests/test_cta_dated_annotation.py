"""Real dated PyEnsembl indexes exercise public OncoRef and Vaxrank admission."""

import gzip
from io import BytesIO

import oncoref
import pytest
import requests
from pyensembl import EnsemblRelease

from vaxrank.cta_expression import CTAExpressionResult
from .test_cta_expression import run

PRAME = "ENSG00000185686"
TRANSCRIPT = "ENST99999999999"
PROTEIN = "ENSP99999999999"


def snapshot(root):
    records = {}
    for entry in root.rglob("*"):
        state = entry.stat()
        records[str(entry.relative_to(root))] = (
            state.st_ino, state.st_mode, state.st_size,
            state.st_mtime_ns, state.st_ctime_ns,
            entry.read_bytes() if entry.is_file() else None)
    return records


@pytest.mark.parametrize("release,options,assembly", [
    ("2013_09", {"assembly": "GRCh37", "reference_name": "historical study"}, "GRCh37"),
    ("2026_04", {"assembly_accession": "GCA_000001405.29"}, "GRCh38"),
])
def test_public_dated_admission_and_replay_preserve_installed_data(
    tmp_path, monkeypatch, release, options, assembly,
):
    cache = tmp_path / "cache"
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(cache))
    genome = EnsemblRelease(release, species="human",
                            server="https://fixture.invalid/ensembl", **options)
    # Synthetic coordinates/sequences exercise the actual source/index path;
    # they do not stand in for an upstream PRAME sequence or dated annotation.
    attributes = (f'gene_id "{PRAME}"; gene_name "PRAME"; gene_biotype "protein_coding"; '
                  f'transcript_id "{TRANSCRIPT}"; transcript_biotype "protein_coding"; '
                  f'protein_id "{PROTEIN}"; exon_id "synthetic-exon"; exon_number "1";')
    rows = "".join(f"22\ttest\t{feature}\t{start}\t{end}\t.\t+\t{frame}\t{attributes}\n"
                   for feature, start, end, frame in [
                       ("gene", 1, 33, "."), ("transcript", 1, 33, "."),
                       ("exon", 1, 33, "."), ("CDS", 1, 30, "0"),
                       ("start_codon", 1, 3, "0"), ("stop_codon", 31, 33, "0")])
    cdna = "ATG" + "GCT" * 9 + "TAA"
    peptide = "M" + "A" * 9
    bodies = {genome.gtf_url: gzip.compress(rows.encode())}
    bodies.update((url, gzip.compress(f">{TRANSCRIPT}\n{cdna}\n".encode()))
                  for url in genome.transcript_fasta_urls)
    bodies.update((url, gzip.compress(f">{PROTEIN}\n{peptide}\n".encode()))
                  for url in genome.protein_fasta_urls)

    def download(url, **kwargs):
        response = requests.Response()
        response.status_code, response.url = 200, url
        response.raw = BytesIO(bodies[url])
        response.headers["Content-Length"] = str(len(bodies[url]))
        return response

    monkeypatch.setattr(requests, "get", download)
    genome.download()
    genome.index()
    encoded = genome.to_json()
    genome.close()
    before = snapshot(cache)
    with EnsemblRelease.from_json(encoded) as reopened:
        result = run(tmp_path, reopened, [(TRANSCRIPT, 3)], level="transcript")
        assert result.reference_assembly == assembly
        assert result.annotation_version == release
        reference = result.annotation_reference
        assert reference.annotation_assembly == assembly
        assert reference.annotation_release == release
        assert type(reference.annotation_release) is str
        assert reference.annotation_hash_kind == "pyensembl-source-files"
        assert reference.self_reference_excluded_gene_ids == (PRAME,)
        assert [decision.status for decision in result.decisions] == ["admitted"]
        antigen = result.admitted_antigens[0]
        assert antigen.gene_id == PRAME
        assert antigen.transcript_ids == (TRANSCRIPT,)
        assert antigen.protein_ids == (PROTEIN,)
        assert antigen.amino_acids == peptide
        assert reopened.reference_name == genome.reference_name
    assert snapshot(cache) == before
    saved = tmp_path / "dated-admission.json"
    result.save(saved)

    def forbidden(*args, **kwargs):
        pytest.fail("Saved admission replay queried current OncoRef data")

    monkeypatch.setattr(oncoref, "resolve_gene_identity", forbidden)
    monkeypatch.setattr(oncoref, "cta_annotation_gene_identities", forbidden)
    monkeypatch.setattr(oncoref.gene_identity, "resolve_gene_identity", forbidden)
    restored = CTAExpressionResult.load(saved)
    assert restored == result
    assert restored.annotation_reference.fingerprint == reference.fingerprint
    assert snapshot(cache) == before
