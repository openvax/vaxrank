"""Explicit host declarations gate mixed inputs without changing source species."""

from dataclasses import replace
import hashlib
import json
from types import SimpleNamespace

import pandas as pd
from pyensembl import EnsemblRelease
import pytest

from vaxrank import AntigenDataset
from vaxrank.cli.entry_point import run_cli
from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.external_rescoring import prepare_reports
from vaxrank.input_scope import (
    InputProvenance, InputScope, read_input_manifest, report_declarations,
    scope_dict, scope_from_mapping, validate_input_scopes,
)
from vaxrank.native_serialization import from_native_json, to_native_json

from . import test_antigen_input
from .test_antigen_manifest import design_options, write_manifest
from .test_epitope_dataset import forbidden

inputs = test_antigen_input.inputs
dependencies = test_antigen_input.dependencies


@pytest.mark.parametrize('value', [True, False, 0, -1, 9606.0, '9606', [], {}])
def test_host_scope_requires_a_positive_integer(tmp_path, value):
    path = write_manifest(tmp_path, [('antigens', 'new.json')], host_taxon_id=value)
    with pytest.raises(ValueError, match='host_taxon_id must be a positive integer'):
        read_input_manifest(path)
    with pytest.raises(ValueError, match='host_taxon_id must be a positive integer'):
        InputScope(host_taxon_id=value)


def test_host_scope_conflicts_between_manifest_and_producer(tmp_path):
    path = write_manifest(tmp_path, [('antigens', 'new.json')], host_taxon_id=9606)
    document = json.loads(path.read_text())
    document['inputs'][0]['host_taxon_id'] = 9615
    path.write_text(json.dumps(document))
    with pytest.raises(ValueError, match='conflicting host_taxon_id'):
        read_input_manifest(path)
    frame = pd.DataFrame({'host_taxon_id': ['9606', '9615']})
    with pytest.raises(ValueError, match='conflicting host_taxon_id'):
        report_declarations('producer.tsv', 'topiary', frame=frame)
    scope, _ = report_declarations('producer.tsv', 'topiary', frame=frame.iloc[:1])
    assert scope.host_taxon_id == 9606


@pytest.mark.parametrize('species,taxon,allele', [
    ('human', 9606, 'HLA-A*02:01'),
    ('mouse', 10090, 'H2-Kb'),
    ('dog', 9615, 'DLA-88*508:01'),
])
def test_host_scope_validates_real_annotation_metadata_without_sequence_access(species, taxon, allele):
    genome = EnsemblRelease(93, species=species)
    scope = scope_from_mapping(dict(host_taxon_id=taxon, mhc_alleles=[allele]), 'test')
    provenance = InputProvenance('test', 'antigens', 'test.json', '', scope)
    assert validate_input_scopes([provenance], genome=genome) == scope
    wrong = replace(provenance, scope=replace(scope, host_taxon_id=taxon + 1, mhc_alleles=None))
    with pytest.raises(ValueError, match='conflicting host_taxon_id'):
        validate_input_scopes([wrong], genome=genome)


def test_grouped_allele_cannot_bypass_host_annotation():
    scope = scope_from_mapping(dict(host_taxon_id=9606, mhc_alleles=['DLA-88*508:01']), 'test')
    source = InputProvenance('test', 'antigens', 'test.json', '', scope)
    with pytest.raises(ValueError, match='allele .* conflicts with host_taxon_id'):
        validate_input_scopes([source], genome=EnsemblRelease(93))


def test_shared_host_applies_to_legacy_sources_without_a_taxon(inputs, dependencies, tmp_path):
    _, genome = dependencies
    genome.species = 'homo_sapiens'
    genome.annotation_name = 'ensembl'
    genome.annotation_version = 93
    inputs.save(tmp_path / 'legacy.json')
    replace(inputs, scope=replace(inputs.scope, host_taxon_id=9615)).save(tmp_path / 'dog.json')
    with pytest.raises(ValueError, match='host_taxon_id'):
        prepare_reports([('antigens', tmp_path / 'legacy.json'), ('antigens', tmp_path / 'dog.json')],
                        genome=genome, antigen_predictor_factory=forbidden)


def test_host_conflict_fails_before_predictor_or_output(inputs, dependencies, tmp_path, monkeypatch):
    inputs.save(tmp_path / 'new.json')
    path = write_manifest(tmp_path, [('antigens', 'new.json')], host_taxon_id=9615)
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    output = tmp_path / 'out'
    with pytest.raises(ValueError, match='host_taxon_id'):
        run_cli(['--input-manifest', str(path), '--ensembl-release', '93',
                 '--mhc-predictor', 'random', '--output-dir', str(output)])
    assert not output.exists()


@pytest.mark.parametrize('modality', ['peptide', 'mrna'])
def test_viral_source_keeps_human_host_through_design_and_offline_replay(
        inputs, dependencies, tmp_path, monkeypatch, modality):
    _, genome = dependencies
    genome.species = 'homo_sapiens'
    genome.annotation_name = 'ensembl'
    genome.annotation_version = 93
    viral = replace(inputs, antigens={'viral': inputs.antigens['viral']})
    viral.save(tmp_path / 'virus.json')
    path = write_manifest(tmp_path, [('antigens', 'virus.json')], host_taxon_id=9606)
    output = tmp_path / 'design'
    run_cli(['--input-manifest', str(path), '--ensembl-release', '93',
             '--mhc-predictor', 'random'] + design_options(output, modality))
    saved = EpitopeDataset.load(output / 'candidate_predictions.tsv')
    assert saved.provenance[0].scope.host_taxon_id == 9606
    exported = AntigenDataset.from_epitope_dataset(saved)
    assert exported.scope.host_taxon_id == 9606
    assert {a.species for a in exported.antigens.values()} == {'fixture virus'}
    exported.save(tmp_path / 'export.json')
    assert AntigenDataset.load(tmp_path / 'export.json') == exported
    (tmp_path / 'virus.json').unlink()
    path.unlink()
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    monkeypatch.setattr('vaxrank.antigen_input.self_reference_matches', forbidden)
    monkeypatch.setattr('vaxrank.cli.entry_point.resolve_ensembl_release',
                        lambda args: setattr(args, 'genome', None))
    monkeypatch.setattr('pyensembl.Genome.index', forbidden)
    monkeypatch.setattr('pyensembl.Genome.download', forbidden)
    monkeypatch.setattr('pyensembl.Genome.transcripts', forbidden)
    replay = tmp_path / 'replay'
    run_cli(['--input-epitopes', str(output / 'candidate_predictions.tsv'),
             '--output-dir', str(replay), '--no-processing-aware-annotation', '--pdf-backend', 'weasyprint'])
    for original in output.glob('*.fasta'):
        assert original.read_bytes() == (replay / original.name).read_bytes()
    assert (output / 'manifest.json').read_bytes() == (replay / 'manifest.json').read_bytes()
    assert EpitopeDataset.load(replay / 'candidate_predictions.tsv').provenance[0].scope.host_taxon_id == 9606


def test_old_scope_and_antigen_identities_are_unchanged(inputs):
    # Recorded from released Vaxrank 3.44.1 (c8822ed), before adding host_taxon_id.
    scope = InputScope(patient_id='patient', reference_assembly='GRCh38', annotation='ensembl:93',
                       mhc_alleles=('HLA-A*02:01',))
    assert hashlib.sha256(json.dumps(scope_dict(scope), sort_keys=True).encode()).hexdigest() == (
        'bafe799ff8798dea6a5cb8df26674633c557838077a4be1354e9067430470ad3')
    provenance = inputs.prediction_provenance(genome=SimpleNamespace(reference_name='GRCh38', release=93))
    assert provenance.input_id == 'antigens:7a69f4db79bca347e75581bd918012195000e622781f591f7066f0c9b98b2bb9'
    raw = json.loads(to_native_json(scope))
    del raw['host_taxon_id']
    assert from_native_json(json.dumps(raw), InputScope) == scope


def test_declared_host_changes_native_occurrence_identity(inputs, dependencies):
    _, genome = dependencies
    genome.species = 'homo_sapiens'
    genome.annotation_name = 'ensembl'
    genome.annotation_version = 93
    declared = replace(inputs, scope=replace(inputs.scope, host_taxon_id=9606))
    assert declared.prediction_provenance(genome=genome).input_id != inputs.prediction_provenance(genome=genome).input_id
