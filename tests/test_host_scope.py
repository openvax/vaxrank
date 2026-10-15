"""Explicit host declarations gate mixed inputs without changing source species."""

from dataclasses import replace
import hashlib
import json
from types import SimpleNamespace

import pandas as pd
from pyensembl import EnsemblRelease
import pytest

from vaxrank import AntigenDataset
from vaxrank.cli.entry_point import run_cli, resolve_ensembl_release
from vaxrank.epitope_dataset import EpitopeDataset, read_table
from vaxrank.external_rescoring import prepare_reports, read_reports, external_input_specs
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
@pytest.mark.parametrize('host,taxon,assembly,mhc_species', [
    ('homo_sapiens', 9606, 'GRCh38', None),
    ('mus_musculus', 10090, 'GRCm38', 'human'),
])
def test_viral_source_keeps_human_host_through_design_and_offline_replay(
        inputs, dependencies, tmp_path, monkeypatch, modality, host, taxon, assembly, mhc_species):
    _, genome = dependencies
    genome.species = host
    genome.reference_name = assembly
    genome.annotation_name = 'ensembl'
    genome.annotation_version = 93
    viral = replace(inputs, antigens={'viral': inputs.antigens['viral']},
                    scope=replace(inputs.scope, reference_assembly=assembly))
    viral.save(tmp_path / 'virus.json')
    path = write_manifest(tmp_path, [('antigens', 'virus.json')], host_species=host, mhc_species=mhc_species)
    output = tmp_path / 'design'
    run_cli(['--input-manifest', str(path), '--ensembl-release', '93',
             '--mhc-predictor', 'random'] + design_options(output, modality))
    saved = EpitopeDataset.load(output / 'candidate_predictions.tsv')
    assert saved.provenance[0].scope.host_taxon_id == taxon
    assert saved.provenance[0].scope.mhc_species == (('Homo sapiens',) if mhc_species else None)
    exported = AntigenDataset.from_epitope_dataset(saved)
    assert exported.scope.host_taxon_id == taxon
    assert exported.scope.mhc_species == saved.provenance[0].scope.mhc_species
    assert exported.scope.reference_assembly == assembly
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
    assert EpitopeDataset.load(replay / 'candidate_predictions.tsv').provenance[0].scope == saved.provenance[0].scope


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
    del raw['mhc_species']
    assert from_native_json(json.dumps(raw), InputScope) == scope


def test_declared_host_changes_native_occurrence_identity(inputs, dependencies):
    _, genome = dependencies
    genome.species = 'homo_sapiens'
    genome.annotation_name = 'ensembl'
    genome.annotation_version = 93
    declared = replace(inputs, scope=replace(inputs.scope, host_taxon_id=9606))
    assert declared.prediction_provenance(genome=genome).input_id != inputs.prediction_provenance(genome=genome).input_id


def write_species_table(tmp_path, *, suffix='tsv', hosts=('9606', ''),
                        alleles=('HLA-A*02:01', 'HLA-A*03:01')):
    path = tmp_path / ('species.' + suffix)
    pd.DataFrame(dict(peptide=['SIINFEKL'] * len(alleles), allele=alleles,
                      kind=['pMHC_affinity'] * len(alleles), value=[10] * len(alleles),
                      host_taxon_id=hosts)).to_csv(path, sep='\t' if suffix == 'tsv' else ',', index=False)
    return path


@pytest.mark.parametrize('suffix', ['tsv', 'csv'])
def test_sparse_host_columns_are_normalized_at_ingestion(tmp_path, suffix, monkeypatch):
    path = write_species_table(tmp_path, suffix=suffix)
    result = read_table(path)
    assert result.df.host_taxon_id.tolist() == [9606, None]
    assert type(result.df.host_taxon_id.iloc[0]) is int
    adapted = EpitopeDataset.from_topiary

    def checked(result, **kwargs):
        assert result.df.host_taxon_id.tolist() == [9606, None]
        assert type(result.df.host_taxon_id.iloc[0]) is int
        return adapted(result, **kwargs)

    monkeypatch.setattr(EpitopeDataset, 'from_topiary', checked)
    report = read_reports([('topiary', path)])[0]
    assert report.input_provenance.scope.host_taxon_id == 9606
    assert len(report.epitopes) == 2


@pytest.mark.parametrize('value', ['true', 'false', '0', '-1', '9606.5', 'Infinity', '"human"'])
def test_invalid_table_host_stops_before_candidate_adaptation(tmp_path, monkeypatch, value):
    path = write_species_table(tmp_path, hosts=(value, ''))
    monkeypatch.setattr(EpitopeDataset, 'from_topiary', forbidden)
    with pytest.raises(ValueError, match='host_taxon_id'):
        read_reports([('topiary', path)])


@pytest.mark.parametrize('allele', ['DLA-88*508:01', 'Mamu-A*01:01', 'H2-Kb'])
def test_annotation_free_wrong_species_is_rejected(tmp_path, allele):
    path = write_species_table(tmp_path, hosts=('9606',), alleles=(allele,))
    with pytest.raises(ValueError, match='conflicts with host_taxon_id'):
        read_reports([('topiary', path)])


def test_annotation_free_dog_group_is_compatible(tmp_path):
    path = write_species_table(tmp_path, hosts=('9615',), alleles=('DLA-88*508:01',))
    assert read_reports([('topiary', path)])[0].input_provenance.scope.host_taxon_id == 9615


@pytest.mark.parametrize('alleles,species', [
    (['HLA-A*02:01'], ['human']),
    (['HLA-A*02:01', 'H2-Kb'], ['human', 'mouse']),
])
def test_hla_mouse_explicit_repertoire_preserves_mouse_reference(tmp_path, alleles, species):
    path = write_species_table(tmp_path, hosts=['10090'] * len(alleles), alleles=alleles)
    scope = dict(mhc_species=species, reference_assembly='GRCm38')
    report = read_reports([('topiary', path)], scopes=[scope], genome=EnsemblRelease(93, species='mouse'))[0]
    assert report.input_provenance.scope.host_taxon_id == 10090
    assert report.input_provenance.scope.reference_assembly == 'GRCm38'
    with pytest.raises(ValueError, match='conflicting host_taxon_id'):
        read_reports([('topiary', path)], scopes=[scope], genome=EnsemblRelease(93))
    with pytest.raises(ValueError, match='conflicts with mhc_species'):
        validate_input_scopes([report.input_provenance], prediction_alleles=['DLA-88*508:01'])


@pytest.mark.parametrize('species', [[], '', ['unknown species'], ['Canis sp.'], True, 9606.5, [None]])
def test_mhc_species_must_be_concrete_species(species):
    with pytest.raises(ValueError, match='mhc_species'):
        InputScope(mhc_species=species)


def test_manifest_mhc_species_and_cli_assertions(tmp_path):
    path = write_manifest(tmp_path, [('topiary', 'input.tsv')], host_taxon_id=10090,
                          mhc_species=['human', 'mouse'])
    args = SimpleNamespace(input_manifest=path, mhc_species=['Mus musculus', 'Homo sapiens'])
    specs = external_input_specs(args)
    assert specs[0][2]['mhc_species'] == ['Homo sapiens', 'Mus musculus']
    args.mhc_species = ['human']
    with pytest.raises(ValueError, match='conflicting mhc_species'):
        external_input_specs(args)


def test_actual_annotation_species_checked_without_taxon_declaration():
    scope = InputScope(mhc_alleles=('HLA-A*02:01',))
    source = InputProvenance('test', 'antigens', 'test.json', '', scope)
    with pytest.raises(ValueError, match='conflicts with host species'):
        validate_input_scopes([source], genome=EnsemblRelease(93, species='mouse'))
    assert validate_input_scopes([replace(source, scope=replace(scope, mhc_species=['human']))],
                                 genome=EnsemblRelease(93, species='mouse')).mhc_species == ('Homo sapiens',)


@pytest.mark.parametrize('location', ['prediction', 'comparator', 'table'])
def test_native_replay_checks_actual_prediction_species(tmp_path, location):
    from mhctools.pred import Prediction
    from vaxrank.candidate_epitope import CandidateEpitope, Peptide
    scope = InputScope(host_taxon_id=9606)
    source = InputProvenance('test', 'topiary', 'original.tsv', '', scope)
    prediction = Prediction(kind='pMHC_affinity', peptide='SIINFEKL', allele='DLA-88*508:01',
                            value=10, score=None, predictor_name='fixture')
    candidate = CandidateEpitope(sequence='SIINFEKL', prediction_id='one', input_provenance=source,
        predictions=(prediction if location == 'prediction' else replace(prediction, allele='HLA-A*02:01'),))
    if location == 'comparator':
        candidate = replace(candidate, comparators={'wt': Peptide(sequence='SIINFEKL', predictions=(prediction,))})
    path = tmp_path / 'native.tsv'
    dataset = EpitopeDataset.from_predictions([candidate])
    if location == 'table':
        dataset.result.df['allele'] = 'DLA-88*508:01'
    dataset.save(path)
    with pytest.raises(ValueError, match='conflicts with host_taxon_id'):
        read_reports([('epitopes', path)])


def test_direct_hla_mouse_checks_species_before_model_acquisition(monkeypatch):
    from varcode import Variant, VariantCollection
    from vaxrank.cli import entry_point as entry
    from vaxrank.cli.arg_parser import parse_vaxrank_args
    genome = EnsemblRelease(93, species='mouse')
    variants = VariantCollection([Variant('1', 100, 'A', 'T', ensembl=genome)])
    monkeypatch.setattr(entry, 'variant_collection_from_args', lambda args: variants)
    monkeypatch.setattr(entry, 'filter_unannotatable_variants', lambda variants: variants)
    monkeypatch.setattr(entry, 'predictors_from_args', forbidden)
    args = parse_vaxrank_args(['--vcf', 'mouse.vcf', '--bam', 'mouse.bam',
                               '--mhc-predictor', 'random', '--mhc-alleles', 'HLA-A*02:01',
                               '--host-taxon-id', '10090'])
    args.genome = genome
    with pytest.raises(ValueError, match='conflicts with host_taxon_id'):
        entry.run_vaxrank_from_parsed_args(args)
    args.mhc_species = ['human']
    # The same query now reaches model acquisition; no annotation/RNA I/O needed.
    with pytest.raises(AssertionError, match='must not instantiate a predictor'):
        entry.run_vaxrank_from_parsed_args(args)


def test_chimeric_scope_changes_identity_and_survives_native_export(inputs, dependencies):
    _, genome = dependencies
    original = inputs.prediction_provenance(genome=genome)
    declared = replace(inputs, scope=replace(inputs.scope, mhc_species=['human', 'mouse']))
    assert declared.prediction_provenance(genome=genome).input_id != original.input_id
    assert from_native_json(to_native_json(declared.scope), InputScope) == declared.scope


def test_hla_mouse_antigen_cli_uses_real_mouse_reference_resolution(inputs, dependencies, tmp_path, monkeypatch):
    monkeypatch.setattr('vaxrank.cli.entry_point.resolve_ensembl_release', resolve_ensembl_release)
    monkeypatch.setattr('pyensembl.Genome.download', forbidden)
    monkeypatch.setattr('pyensembl.Genome.index', forbidden)
    viral = replace(inputs, antigens={'viral': inputs.antigens['viral']}, scope=InputScope.from_dict({
        **scope_dict(inputs.scope), 'reference_assembly': 'GRCm38',
        'host_species': 'mouse', 'mhc_species': 'human'}))
    path = tmp_path / 'virus.json'
    viral.save(path)
    saved = tmp_path / 'predictions.tsv'
    run_cli(['--input-antigens', str(path), '--mhc-predictor', 'random',
             '--genome', 'GRCm38', '--ensembl-release', '93',
             '--host-species', 'Mus musculus', '--mhc-species', '9606',
             '--no-processing-aware-annotation', '--output-epitopes', str(saved)])
    scope = EpitopeDataset.load(saved).provenance[0].scope
    assert scope.reference_assembly == 'GRCm38'
    assert scope.annotation == 'ensembl:93'
    assert scope.host_taxon_id == 10090
    assert scope.mhc_species == ('Homo sapiens',)


@pytest.mark.parametrize('host', ['mouse', 'Mus musculus', 'mus_musculus', 10090, '10090'])
@pytest.mark.parametrize('mhc', ['human', 'Homo sapiens', 9606, '9606', ['human'], ['human', 9606]])
def test_species_aliases_share_canonical_scope_and_identity(host, mhc):
    expected = InputScope(host_taxon_id=10090, mhc_species=['Homo sapiens'])
    actual = InputScope.from_dict(dict(host_species=host, mhc_species=mhc))
    assert actual == expected
    assert scope_from_mapping(dict(host_species=host, mhc_species=mhc), 'input') == expected
    assert scope_dict(actual) == scope_dict(expected)
    assert to_native_json(actual) == to_native_json(expected)
    assert from_native_json(to_native_json(actual), InputScope) == expected


@pytest.mark.parametrize('host', ['', 'unknown species', True, 10090.0, ['mouse']])
def test_invalid_host_species_is_rejected(host):
    with pytest.raises(ValueError, match='host_species'):
        InputScope.from_dict(dict(host_species=host))


def test_host_alias_conflicts_are_rejected_in_api_manifest_and_cli(tmp_path):
    with pytest.raises(ValueError, match='conflicting host_species and host_taxon_id'):
        InputScope.from_dict(dict(host_species='mouse', host_taxon_id=9606))
    path = write_manifest(tmp_path, [('topiary', 'input.tsv')], host_species='mouse')
    document = json.loads(path.read_text())
    document['inputs'][0]['host_taxon_id'] = 9606
    path.write_text(json.dumps(document))
    with pytest.raises(ValueError, match='conflicting host_taxon_id'):
        read_input_manifest(path)
    path = write_manifest(tmp_path, [('topiary', 'input.tsv')], host_species='mouse')
    with pytest.raises(ValueError, match='conflicting host_taxon_id'):
        external_input_specs(SimpleNamespace(input_manifest=path, host_species='human'))


@pytest.mark.parametrize('host,mhc', [('mouse', 'human'), (10090, 9606)])
def test_sparse_species_alias_columns_load_and_validate_alleles(tmp_path, host, mhc):
    path = write_species_table(tmp_path)
    frame = pd.read_csv(path, sep='\t').drop(columns='host_taxon_id')
    frame['host_species'] = [host, None]
    frame['mhc_species'] = [mhc, None]
    frame.to_csv(path, sep='\t', index=False)
    report = read_reports([('topiary', path)])[0]
    scope = report.input_provenance.scope
    assert scope == InputScope(host_taxon_id=10090, mhc_species=['human'])
    frame['allele'] = ['H2-Kb', 'H2-Kb']
    frame.to_csv(path, sep='\t', index=False)
    with pytest.raises(ValueError, match='conflicts with mhc_species'):
        read_reports([('topiary', path)])
