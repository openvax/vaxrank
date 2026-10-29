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
    inputs = replace(inputs, antigens={k: a for k, a in inputs.antigens.items() if a.kind not in ('SV', 'isoform')})
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
    import mhctools.cli
    from varcode import Variant, VariantCollection
    from vaxrank.cli import entry_point as entry
    from vaxrank.cli.arg_parser import parse_vaxrank_args
    genome = EnsemblRelease(93, species='mouse')
    variants = VariantCollection([Variant('1', 100, 'A', 'T', ensembl=genome)])
    monkeypatch.setattr(entry, 'variant_collection_from_args', lambda args: variants)
    monkeypatch.setattr(entry, 'filter_unannotatable_variants', lambda variants: variants)
    monkeypatch.setattr(mhctools.cli, 'predictors_from_args', forbidden)
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
    assert report_declarations(path, 'topiary')[0] == scope
    frame['allele'] = ['H2-Kb', 'H2-Kb']
    frame.to_csv(path, sep='\t', index=False)
    with pytest.raises(ValueError, match='conflicts with mhc_species'):
        read_reports([('topiary', path)])


@pytest.mark.parametrize('name,scientific,taxon', [
    ('human', 'Homo sapiens', 9606), ('mouse', 'Mus musculus', 10090),
    ('dog', 'Canis lupus familiaris', 9615), ('cat', 'Felis catus', 9685),
    ('rat', 'Rattus norvegicus', 10116), ('cow', 'Bos taurus', 9913),
    ('cattle', 'Bos taurus', 9913),
])
def test_host_and_mhc_aliases_share_taxonomy(name, scientific, taxon):
    scopes = [InputScope.from_dict(dict(host_species=value, mhc_species=value))
              for value in (name, scientific, taxon, str(taxon))]
    assert all(scope == scopes[0] for scope in scopes)
    assert len({to_native_json(scope) for scope in scopes}) == 1
    mixed = InputScope.from_dict(dict(host_species=name, mhc_species=['human', name]))
    assert mixed == InputScope.from_dict(dict(host_species=taxon, mhc_species=[9606, taxon]))


@pytest.mark.parametrize('host', ['cow', 'cattle', 'Bos taurus', 9913])
def test_cattle_table_needs_no_annotation(tmp_path, host, monkeypatch):
    from pyensembl.species import Species
    monkeypatch.setattr(Species, 'all_registered_latin_names', forbidden)
    monkeypatch.setattr('pyensembl.Genome.index', forbidden)
    monkeypatch.setattr('pyensembl.Genome.download', forbidden)
    path = write_species_table(tmp_path, hosts=['9913'], alleles=['BoLA-1*009:01'])
    scope = dict(host_species=host)
    report = read_reports([('topiary', path)], scopes=[scope])[0]
    assert report.input_provenance.scope.host_taxon_id == 9913
    report.dataset.save(tmp_path / 'native.tsv')
    assert read_reports([('epitopes', tmp_path / 'native.tsv')])[0].input_provenance.scope.host_taxon_id == 9913
    with pytest.raises(ValueError, match='conflicts with host_taxon_id'):
        validate_input_scopes([report.input_provenance], prediction_alleles=['HLA-A*02:01'])


def write_legacy_species_table(tmp_path, suffix, payload, *, scope=None, **columns):
    from mhctools.pred import Prediction
    from vaxrank.candidate_epitope import CandidateEpitope
    from vaxrank.epitope_io import NATIVE_EPITOPE_COLUMN, save_predictions
    provenance = InputProvenance('original', 'topiary', 'original.tsv', '', scope) if scope else None
    prediction = Prediction(kind='pMHC_affinity', peptide='SIINFEKL', allele='HLA-A*02:01',
                            value=10, score=1.0, predictor_name='fixture')
    candidates = [CandidateEpitope(sequence='SIINFEKL', prediction_id=str(index),
                                  predictions=(prediction,), input_provenance=provenance)
                  for index in range(2)]
    path = tmp_path / ('legacy.' + suffix)
    save_predictions(candidates, path)
    sep = '\t' if suffix == 'tsv' else ','
    frame = pd.read_csv(path, sep=sep)
    if not payload:
        frame = frame.drop(columns=NATIVE_EPITOPE_COLUMN)
    for name, value in columns.items():
        frame[name] = [value, None]
    frame.to_csv(path, sep=sep, index=False)
    return path


@pytest.mark.parametrize('suffix', ['tsv', 'csv'])
@pytest.mark.parametrize('payload', [False, True])
def test_legacy_native_species_columns_validate_before_replay(tmp_path, suffix, payload):
    path = write_legacy_species_table(tmp_path, suffix, payload, host_species='mouse')
    with pytest.raises(ValueError, match='conflicts with host_taxon_id'):
        read_reports([('epitopes', path)])
    # An explicit HLA repertoire is accepted; adding a conflicting CLI/manifest
    # assertion must not silently replace the table's host.
    path = write_legacy_species_table(tmp_path, suffix, payload,
                                     patient_id='patient#1', host_species='mouse', mhc_species='human')
    with pytest.raises(ValueError, match='conflicting host_taxon_id'):
        read_reports([('epitopes', path)], scopes=[dict(host_species='human')])
    report = read_reports([('epitopes', path)])[0]
    expected = InputScope.from_dict(dict(patient_id='patient#1', host_species='mouse', mhc_species='human'))
    assert report.input_provenance.scope == expected
    assert report.input_provenance.report_declarations['host_species'] == ['mouse']
    assert all(e.input_provenance.scope == expected for e in report.epitopes)
    saved = tmp_path / ('replay.' + suffix)
    report.dataset.save(saved)
    replay = read_reports([('epitopes', saved)])[0]
    assert replay.input_provenance.scope == expected
    assert [e.prediction_group_key for e in replay.epitopes] == [e.prediction_group_key for e in report.epitopes]


@pytest.mark.parametrize('suffix', ['tsv', 'csv'])
def test_legacy_native_sparse_numeric_aliases_through_cli(tmp_path, suffix):
    path = write_legacy_species_table(tmp_path, suffix, True, host_species=10090, mhc_species=9606)
    saved = tmp_path / 'cli.tsv'
    run_cli(['--input-epitopes', str(path), '--output-epitopes', str(saved),
             '--no-processing-aware-annotation'])
    scope = EpitopeDataset.load(saved).provenance[0].scope
    assert scope == InputScope(host_taxon_id=10090, mhc_species='human')
    replay = tmp_path / 'replay.tsv'
    run_cli(['--input-epitopes', str(saved), '--output-epitopes', str(replay),
             '--no-processing-aware-annotation'])
    assert EpitopeDataset.load(replay).provenance[0].scope == scope


@pytest.mark.parametrize('column,value', [('host_taxon_id', 10090.5), ('host_species', True),
                                         ('mhc_species', 9606.5)])
def test_invalid_legacy_native_species_stops_before_adaptation(tmp_path, monkeypatch, column, value):
    from vaxrank.external_rescoring import READERS
    path = write_legacy_species_table(tmp_path, 'tsv', True, **{column: value})
    monkeypatch.setitem(READERS, 'epitopes', forbidden)
    with pytest.raises(ValueError, match=column):
        read_reports([('epitopes', path)])


def test_legacy_native_table_scope_reconciles_with_saved_payload(tmp_path):
    path = write_legacy_species_table(tmp_path, 'tsv', True,
                                     scope=InputScope(host_taxon_id=9606), host_species='mouse')
    with pytest.raises(ValueError, match='conflicting host_taxon_id'):
        read_reports([('epitopes', path)])
    path = write_legacy_species_table(tmp_path, 'tsv', True,
                                     scope=InputScope(host_taxon_id=10090), mhc_species='human')
    report = read_reports([('epitopes', path)])[0]
    assert report.input_provenance.scope == InputScope(host_taxon_id=10090, mhc_species='human')


def pooled_hla_mouse_inputs(tmp_path, native_indices):
    common = dict(patient_id='patient', reference_assembly='GRCm38',
                  mhc_alleles=['HLA-A*02:01'], host_species='mouse')
    scopes = [common, {**common, 'mhc_species': 'human'}]
    inputs = []
    for index in range(2):
        path = tmp_path / f'source-{index}.tsv'
        pd.DataFrame(dict(peptide=[('SIINFEKL', 'GILGFVFTL')[index]], allele=['HLA-A*02:01'],
                          kind=['pMHC_affinity'], value=[10 + index])).to_csv(path, sep='\t', index=False)
        inputs.append(('topiary', path))
    reports = read_reports(inputs, scopes=scopes)
    for index in native_indices:
        # Older exports retain only the declarations supplied on that source.
        report = reports[index]
        provenance = replace(report.input_provenance, scope=InputScope.from_dict(scopes[index]))
        dataset = replace(report.dataset, provenance=(provenance,), epitopes=tuple(
            replace(e, input_provenance=provenance) for e in report.epitopes))
        path = tmp_path / f'native-{index}.tsv'
        dataset.save(path)
        inputs[index] = ('epitopes', path)
        scopes[index] = {}
    return inputs, scopes


@pytest.mark.parametrize('native_indices', [(), (0,), (1,), (0, 1)])
@pytest.mark.parametrize('reverse', [False, True])
def test_pooled_species_scope_survives_individual_export(tmp_path, native_indices, reverse):
    inputs, scopes = pooled_hla_mouse_inputs(tmp_path, native_indices)
    if reverse:
        inputs, scopes = inputs[::-1], scopes[::-1]
    reports = read_reports(inputs, scopes=scopes)
    for index, report in enumerate(reports):
        assert report.input_provenance.scope.mhc_species == ('Homo sapiens',)
        assert all(p.scope.mhc_species == ('Homo sapiens',) for p in report.dataset.provenance)
        assert all(e.input_provenance.scope.mhc_species == ('Homo sapiens',) for e in report.epitopes)
        saved = tmp_path / f'export-{index}.tsv'
        report.dataset.save(saved)
        replay, = read_reports([('epitopes', saved)])
        assert replay.input_provenance.scope == report.input_provenance.scope
        assert [e.prediction_group_key for e in replay.epitopes] == [
            e.prediction_group_key for e in report.epitopes]


@pytest.mark.parametrize('field,value', [('host_species', 'human'), ('mhc_species', 'mouse')])
def test_pooled_native_species_conflicts_fail_before_prediction(tmp_path, field, value):
    inputs, scopes = pooled_hla_mouse_inputs(tmp_path, (0, 1))
    scopes[0][field] = value
    with pytest.raises(ValueError, match='conflicting (host_taxon_id|mhc_species)'):
        prepare_reports(inputs, scopes=scopes, mode='additive', prediction_prefix='new',
                        model_factory=forbidden)


@pytest.mark.parametrize('invalid', ['missing_patient', 'outside_genotype'])
def test_pooled_native_species_keeps_individual_source_checks(tmp_path, invalid):
    inputs, _ = pooled_hla_mouse_inputs(tmp_path, (0, 1))
    first, second = [EpitopeDataset.load(path) for _, path in inputs]
    provenance, = first.provenance
    if invalid == 'missing_patient':
        provenance = replace(provenance, scope=replace(provenance.scope, patient_id=None))
        message = 'combining reports requires declared patient_id'
    else:
        provenance = replace(provenance, observed_mhc_alleles=('HLA-A*03:01',))
        message = 'reported alleles outside declared mhc_alleles'
    combined = EpitopeDataset.from_predictions([
        *(replace(e, input_provenance=provenance) for e in first.epitopes), *second.epitopes])
    saved = tmp_path / 'invalid-native.tsv'
    combined.save(saved)
    with pytest.raises(ValueError, match=message):
        read_reports([('epitopes', saved)])


def test_pooled_native_species_cli_export_and_replay(tmp_path):
    inputs, scopes = pooled_hla_mouse_inputs(tmp_path, (0,))
    manifest = tmp_path / 'inputs.json'
    manifest.write_text(json.dumps(dict(schema='vaxrank.input_manifest.v1', inputs=[
        dict(format=fmt, path=str(path), **scope) for (fmt, path), scope in zip(inputs, scopes)])))
    saved = tmp_path / 'combined.tsv'
    run_cli(['--input-manifest', str(manifest), '--output-epitopes', str(saved),
             '--no-processing-aware-annotation'])
    replay = tmp_path / 'replay.tsv'
    run_cli(['--input-epitopes', str(saved), '--output-epitopes', str(replay),
             '--no-processing-aware-annotation'])
    original, restored = EpitopeDataset.load(saved), EpitopeDataset.load(replay)
    assert len(original.provenance) == 2
    assert all(p.scope.mhc_species == ('Homo sapiens',) for p in original.provenance)
    assert restored.provenance == original.provenance
    assert [e.prediction_group_key for e in restored.epitopes] == [
        e.prediction_group_key for e in original.epitopes]
