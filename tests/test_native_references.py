"""Actual mutation evidence and annotations survive relocation independently."""

from dataclasses import replace
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pandas as pd
import pytest
from pyensembl.download_cache import MissingLocalFile

from vaxrank.cli.entry_point import run_cli
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_dataset import EpitopeDataset, DATASET_METADATA, read_table
from vaxrank.epitope_dsl import attach_per_allele_scores
from vaxrank.native_references import NativeReferences, REFERENCE_FIELD
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.native_serialization import from_native_json, to_native_json
from vaxrank.vaccine_antigen import VaccineAntigen
from .test_core_logic_config import _make_epitope
from .test_epitope_dataset import forbidden
from .test_isovar_support_contract import MUTANT, reconstruct, reference as _support_reference


@pytest.fixture(params=["+", "-"])
def reference(tmp_path, request):
    return _support_reference.__wrapped__(tmp_path, request)


def native_dataset(tmp_path, reference):
    result = reconstruct(tmp_path, reference, [(MUTANT, [(len(MUTANT), "M")], 0)] * 2)
    fragment = MutantProteinFragment.from_isovar_result(result)
    epitope = replace(_make_epitope(fragment.amino_acids[1:10], ic50=100., wt_ic50=1000.,
        source_sequence=fragment.amino_acids, offset=1), prediction_id="occurrence")
    config = EpitopeConfig(score_expr="1 / affinity.value", min_epitope_score=0)
    dataset = EpitopeDataset.from_predictions([epitope], config=config)
    dataset.epitopes = tuple(attach_per_allele_scores(dataset.epitopes, config,
                                                     topiary_df=dataset.scoring_frame()))
    dataset.mutation_fragments["occurrence"] = fragment
    dataset.antigens["occurrence"] = VaccineAntigen.from_mutant_protein_fragment(fragment)
    path = tmp_path / "archive" / "native.tsv"
    dataset.save(path)
    return dataset, path


def delete_original_reference(tmp_path):
    for name in ("ref.gtf", "cdna.fa", "protein.fa"):
        (tmp_path / name).unlink()
    shutil.rmtree(tmp_path / "cache")


def test_mutation_native_bundle_relocates_without_old_sources_or_predictors(tmp_path, reference, monkeypatch):
    import mhctools.cli
    import pyensembl.download_cache
    before, path = native_dataset(tmp_path, reference)
    shutil.move(str(path.parent), str(tmp_path / "relocated"))
    path = tmp_path / "relocated" / "native.tsv"
    delete_original_reference(tmp_path)
    monkeypatch.setattr(mhctools.cli, "predictors_from_args", forbidden)
    monkeypatch.setattr(pyensembl.download_cache.DownloadCache, "_fetch", forbidden)
    loaded = EpitopeDataset.load(path)
    assert all(item["available"] for item in loaded.native_references.availability)
    fragment = loaded.mutation_fragments["occurrence"]
    transcript, = fragment.supporting_reference_transcripts
    assert transcript.id == "t" and transcript.gene_id == "g"
    assert transcript.genome.reference_name == before.mutation_fragments["occurrence"].variant.reference_name
    assert (transcript.genome.annotation_name, transcript.genome.annotation_version) == ("test", 1)
    pd.testing.assert_frame_equal(loaded.result.df[before.result.df.columns], before.result.df)
    rescored = attach_per_allele_scores(loaded.epitopes, loaded.config, topiary_df=loaded.scoring_frame())
    assert [e.per_allele_scores for e in rescored] == [e.per_allele_scores for e in before.epitopes]
    # Explicitly prepare local annotation indexes, never a remote/default release.
    loaded.index_native_references()
    assert transcript.gene.name == "G" and transcript.protein_sequence == "MAAAAKGGGG"
    assert fragment.predicted_effect() is not None
    again = tmp_path / "second" / "again.tsv"
    loaded.save(again)
    second = EpitopeDataset.load(again)
    assert second.antigens == loaded.antigens
    assert second.mutation_fragments["occurrence"].supporting_reference_transcripts[0].id == "t"
    # A fresh process must resolve the moved relative resources independently.
    code = """
import sys
from pyensembl.download_cache import DownloadCache
from vaxrank.epitope_dataset import EpitopeDataset
DownloadCache._fetch = lambda *a, **k: (_ for _ in ()).throw(AssertionError('no downloads'))
d = EpitopeDataset.load(sys.argv[1])
assert d.mutation_fragments['occurrence'].supporting_reference_transcripts[0].id == 't'
d.index_native_references()
assert d.mutation_fragments['occurrence'].supporting_reference_transcripts[0].gene.name == 'G'
"""
    subprocess.run([sys.executable, "-c", code, str(again)], check=True)


def test_missing_bundle_allows_evidence_replay_and_reports_unavailable_annotation(tmp_path, reference):
    _, path = native_dataset(tmp_path, reference)
    delete_original_reference(tmp_path)
    shutil.rmtree(Path(str(path) + ".references"))
    loaded = EpitopeDataset.load(path)
    assert not any(item["available"] for item in loaded.native_references.availability)
    assert loaded.epitopes[0].per_allele_scores == {"HLA-A*02:01": .01}
    assert loaded.mutation_fragments["occurrence"].supporting_reference_transcripts[0].id == "t"
    with pytest.raises(ValueError, match="resources are unavailable"):
        loaded.index_native_references()
    with pytest.raises(MissingLocalFile, match="references"):
        loaded.mutation_fragments["occurrence"].supporting_reference_transcripts[0].gene


def test_missing_reference_identity_survives_cli_resave(tmp_path, reference, monkeypatch):
    import mhctools.cli
    import pyensembl.download_cache
    _, path = native_dataset(tmp_path, reference)
    expected = EpitopeDataset.load(path).native_references.availability
    delete_original_reference(tmp_path)
    shutil.rmtree(Path(str(path) + ".references"))
    monkeypatch.setattr(mhctools.cli, "predictors_from_args", forbidden)
    monkeypatch.setattr(pyensembl.download_cache.DownloadCache, "_fetch", forbidden)
    second = tmp_path / "replay" / "native.tsv"
    run_cli(["--input-epitopes", str(path), "--output-epitopes", str(second),
             "--vaccine-peptide-length", "10", "--no-processing-aware-annotation"])
    loaded = EpitopeDataset.load(second)
    def identity(resources):
        return {(r['field'], r['index'], r['sha256'], r['size']) for r in resources}
    assert identity(loaded.native_references.availability) == identity(expected)
    assert not any(r['available'] for r in loaded.native_references.availability)
    assert loaded.antigens
    with pytest.raises(ValueError, match="resources are unavailable"):
        loaded.index_native_references()


def test_legacy_archive_reuses_recorded_protein_ids_without_annotation(tmp_path, reference, monkeypatch):
    import pyensembl.download_cache
    _, path = native_dataset(tmp_path, reference)
    table = read_table(path)
    payload = table.extra[DATASET_METADATA]
    def legacy(value):
        if isinstance(value, list):
            return [legacy(v) for v in value]
        if isinstance(value, dict):
            return {k: legacy(v) for k, v in value.items()
                    if k not in (REFERENCE_FIELD, 'supporting_reference_protein_ids')}
        return value
    for key, value in payload['mutation_fragments'].items():
        payload['mutation_fragments'][key] = json.dumps(legacy(json.loads(value)))
    # Write through the public table codec, preserving its typed evidence.
    table.to_tsv(path)
    delete_original_reference(tmp_path)
    shutil.rmtree(Path(str(path) + ".references"))
    monkeypatch.setattr(pyensembl.download_cache.DownloadCache, "_fetch", forbidden)
    loaded = EpitopeDataset.load(path)
    assert loaded.mutation_fragments['occurrence'].supporting_reference_protein_ids == ('p',)
    second = tmp_path / 'legacy-replay.tsv'
    run_cli(['--input-epitopes', str(path), '--output-epitopes', str(second),
             '--vaccine-peptide-length', '10', '--no-processing-aware-annotation'])
    assert EpitopeDataset.load(second).antigens


def test_changed_bundle_is_rejected_even_if_original_reference_still_exists(tmp_path, reference):
    _, path = native_dataset(tmp_path, reference)
    source = next(Path(str(path) + ".references").glob("*.gtf"))
    source.write_text(source.read_text().replace("gene_name", "gene_fame"))  # Same byte length.
    with pytest.raises(ValueError, match="checksum mismatch"):
        EpitopeDataset.load(path)


def test_reference_change_after_loading_is_rejected_on_resave(tmp_path, reference):
    _, path = native_dataset(tmp_path, reference)
    loaded = EpitopeDataset.load(path)
    source = next(Path(str(path) + ".references").glob("*.gtf"))
    source.write_text(source.read_text().replace("gene_name", "gene_fame"))
    with pytest.raises(ValueError, match="checksum mismatch"):
        loaded.save(tmp_path / "changed.tsv")


def test_ensembl_release_preserves_identity_and_bundles_attached_dna(tmp_path, monkeypatch):
    from pyensembl import EnsemblRelease
    from varcode import Variant
    import pyensembl.download_cache
    monkeypatch.setattr(pyensembl.download_cache.DownloadCache, "_fetch", forbidden)
    dna = tmp_path / "dna.fa"
    dna.write_text(">1\nACGTACGTACGT\n")
    genome = EnsemblRelease(104, species="mouse", genome_fasta=str(dna))
    dataset = EpitopeDataset.from_predictions([_make_epitope("SIINFEKL")])
    dataset.direct_sources = [{'variant': to_native_json(Variant('1', 2, 'C', 'T', ensembl=genome)),
                               'properties': {}}]
    path = tmp_path / "archive" / "native.tsv"
    dataset.save(path)
    shutil.move(str(path.parent), str(tmp_path / 'moved'))
    dna.unlink()
    loaded = EpitopeDataset.load(tmp_path / 'moved' / 'native.tsv')
    variant = from_native_json(loaded.direct_sources[0]['variant'], Variant)
    assert variant.genome.release == 104 and variant.genome.species.latin_name == 'mus_musculus'
    assert variant.reference_name == genome.reference_name
    assert variant.genome.sequence('1', 1, 4) == 'ACGT'


@pytest.mark.parametrize("damage,match", [
    (lambda resource: resource.update(path="../outside.gtf"), "inside"),
    (lambda resource: resource.update(field="annotation_version"), "resource field"),
    (lambda resource: resource.update(index=True), "scalar index"),
    (lambda resource: resource.update(sha256="unknown"), "checksum"),
])
def test_invalid_reference_manifest_rejected(tmp_path, reference, damage, match):
    _, path = native_dataset(tmp_path, reference)
    table = read_table(path)
    payload = table.extra[DATASET_METADATA]
    fragment = json.loads(payload["mutation_fragments"]["occurrence"])
    genome = fragment["variant"]["ensembl"] if "ensembl" in fragment["variant"] else fragment["variant"]["genome"]
    damage(genome[REFERENCE_FIELD]["resources"][0])
    payload["mutation_fragments"]["occurrence"] = json.dumps(fragment)
    with pytest.raises(ValueError, match=match):
        NativeReferences(path).transform(payload, decode=True)


def test_native_cli_scores_and_constructs_replay_after_relocation(tmp_path, reference, monkeypatch):
    import mhctools.cli
    import pyensembl.download_cache
    _, path = native_dataset(tmp_path, reference)
    first = tmp_path / "first.tsv"
    options = ["--vaccine-peptide-length", "10", "--no-processing-aware-annotation"]
    run_cli(["--input-epitopes", str(path), "--output-epitopes", str(first)] + options)
    # Move the original archive, then remove both raw annotations and indexes.
    moved = tmp_path / "moved"
    shutil.move(str(path.parent), str(moved))
    path = moved / "native.tsv"
    delete_original_reference(tmp_path)
    monkeypatch.setattr(mhctools.cli, "predictors_from_args", forbidden)
    monkeypatch.setattr(pyensembl.download_cache.DownloadCache, "_fetch", forbidden)
    output = tmp_path / "report"
    run_cli(["--input-epitopes", str(path), "--index-native-references",
             "--output-dir", str(output), "--output-epitopes", str(tmp_path / "second.tsv")] + options)
    before, after = EpitopeDataset.load(first), EpitopeDataset.load(tmp_path / "second.tsv")
    assert [e.per_allele_scores for e in before.epitopes] == [e.per_allele_scores for e in after.epitopes]
    assert list(output.rglob("*.fasta")) and list(output.rglob("*.pdf"))


def test_index_option_requires_native_input(tmp_path):
    from .test_epitope_dataset import evidence
    table = tmp_path / "table.tsv"
    evidence().to_tsv(table)
    with pytest.raises(ValueError, match="requires a native"):
        run_cli(["--input-topiary", str(table), "--index-native-references",
                 "--output-csv", str(tmp_path / "report.csv")])
