"""Inventory membership is distinct from registered catalogue availability."""

from dataclasses import replace
import json

from orfdata import Proteome
from orfdata.reference import ReferenceContext
import pytest

from vaxrank import HostReference
from vaxrank.cli.entry_point import run_cli
from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.host_reference import combine_host_references
from vaxrank.input_scope import InputScope, read_input_manifest
from vaxrank.native_references import NativeReferences
from vaxrank.reference_coverage import _digest, coverage_view

from .test_antigen_input import inputs  # noqa: F401
from .test_epitope_dataset import forbidden
from .test_host_reference import reference
from .test_product_reports import render, template


def dog_reference():
    original = reference(extra=False)
    context = ReferenceContext(9615, 'canis_lupus_familiaris', 'ROS_Cfam_1.0', 'ensembl', '93')
    return Proteome(context, [replace(e, taxon_id=9615) for e in original],
                    [dict(layer, reference=context.key) for layer in original.layers])


def mouse_reference(source='ncorf-mammals-mouse:1', *, retained=True, layers=True):
    original = reference(species='mouse', assembly='GRCm39')
    proteins = [replace(e, source=source) if e.kind == 'extra_orf' else e for e in original
                if retained or e.kind != 'extra_orf']
    provenance = [dict(layer, source=source, selected=1, added=1, min_protein_length=8,
                       catalog={'source': source, 'records': 1, 'sha256': 'a' * 64})
                  if layer['kind'] == 'extra_orf' else dict(layer) for layer in original.layers
                  if layers or layer['kind'] != 'extra_orf']
    return Proteome(original.reference, proteins, provenance)


def test_dog_unavailable_human_unsearched_and_mouse_releases_are_species_scoped():
    dog = HostReference(dog_reference(), catalogues=[])
    assert dog.coverage['catalogues'] == []
    view = coverage_view(dog.coverage)
    assert view['included'] == []
    assert 'Annotated-only' in ' '.join(view['summary'])
    assert 'registered ncORF coverage is unavailable' in ' '.join(view['summary'])
    human = coverage_view(HostReference(reference(extra=False)).coverage)
    assert human['catalogues'] and {r['status'] for r in human['catalogues']} == {'Not searched'}
    mouse = coverage_view(HostReference(mouse_reference()).coverage)
    riborf = [r for r in mouse['catalogues'] if r['name'] == 'riborf-mouse']
    assert len(riborf) == 1 and riborf[0]['version'] == '3'
    assert {r['name'] for r in mouse['catalogues']} == {
        'ncorf-mammals-mouse', 'riborf-mouse', 'openprot-mouse'}
    assert next(r for r in mouse['catalogues'] if r['name'] == 'ncorf-mammals-mouse')['status'] == 'Included'


def test_compatible_liftover_is_reported_without_acquiring_or_opening_catalogues(monkeypatch):
    from orfdata.catalogs import CatalogRegistry
    for method in ('download', 'open', 'build'):
        monkeypatch.setattr(CatalogRegistry, method, forbidden)
    host = HostReference(reference(species='mouse', extra=False, assembly='GRCm38'))
    assert host.coverage['catalogues']
    assert all(r['route'] and r['source'].endswith('@GRCm38') for r in host.coverage['catalogues'])
    assert HostReference(reference(species='mouse', extra=False, assembly='custom')).coverage['catalogues'] == []


def test_zero_retained_orfs_preserve_selection_and_source_identity():
    host = HostReference(mouse_reference(retained=False), catalogues=['ncorf-mammals-mouse:1'])
    row, = host.coverage['inventory']['included']
    assert row['retained_entries'] == 0 and row['layers'][0]['added'] == 1
    assert row['layers'][0]['catalog']['sha256'] == 'a' * 64
    view = coverage_view(host.coverage)
    assert 'no additional ORF sequences were retained' in ' '.join(view['summary'])
    assert next(r for r in view['catalogues'] if r['name'] == 'ncorf-mammals-mouse')['status'] == 'Included'
    detached = host.coverage
    detached['inventory']['included'].clear()
    assert host.coverage['inventory']['included']


def test_coverage_detaches_nested_source_metadata_from_the_in_memory_proteome():
    proteome = mouse_reference()
    host = HostReference(proteome)
    before = host.coverage
    proteome.layers[-1]['catalog']['records'] = 999
    assert host.coverage == before


@pytest.mark.parametrize('extension', ['html', 'txt'])
def test_reports_preserve_preselection_when_no_further_evidence_filter_is_applied(extension):
    original = mouse_reference()
    selection = dict(complete_sources=['synthetic-complete-source'], evidence={'kind': 'ribosome'},
                     input={'sha256': 'a' * 64})
    layers = [dict(layer, selection=selection, evidence=None) if layer['kind'] == 'extra_orf'
              else dict(layer) for layer in original.layers]
    host = HostReference(Proteome(original.reference, original.proteins, layers))
    report = render(template(host_reference=host), extension)
    assert 'Source preselection:' in report and 'synthetic-complete-source' in report
    assert 'ribosome' in report and 'Additional evidence filter: none' in report
    assert host.coverage['inventory']['included'][0]['layers'][0]['selection'] == selection


def test_unspecified_host_assembly_reports_unknown_compatibility(monkeypatch):
    original = reference(extra=False)
    context = replace(original.reference, genome_build=None, source_sha256='a' * 64)
    proteome = Proteome(context, original.proteins,
                       [dict(layer, reference=context.key) for layer in original.layers])
    monkeypatch.setattr('orfdata.catalog_status', forbidden)
    host = HostReference(proteome)
    assert host.coverage['catalogues'] is None
    summary = ' '.join(coverage_view(host.coverage)['summary'])
    assert 'availability is unknown' in summary
    assert 'coverage is unavailable' not in summary
    assert HostReference(proteome, coverage=host.coverage).coverage == host.coverage


@pytest.mark.parametrize('source,layers', [
    ('custom', True), ('ncorf-mammals-mouse:1+riborf-mouse:3', True), ('ncorf-mammals-mouse:1', False)])
def test_custom_merged_or_unprovenanced_sources_do_not_invent_catalogue_exclusions(source, layers):
    view = coverage_view(HostReference(mouse_reference(source, layers=layers)).coverage)
    assert view['included'][0]['source'] == source
    assert {r['status'] for r in view['catalogues']} == {'Inclusion not established'}


@pytest.mark.parametrize('catalogues', [['unavailable-dog:1'], ['x', 'x'], [''], [1], 'all', {}])
def test_invalid_or_unavailable_catalogue_declarations_fail(catalogues):
    with pytest.raises(ValueError, match='catalogues'):
        HostReference(dog_reference(), catalogues=catalogues)


def test_manifest_asserts_exact_catalogue_labels_before_predictor_acquisition(tmp_path, monkeypatch):
    host = HostReference(dog_reference())
    host.save(tmp_path / 'host.json.gz')
    declaration = dict(path='host.json.gz', identity=host.identity, catalogues=[])
    manifest = tmp_path / 'inputs.json'
    document = dict(schema='vaxrank.input_manifest.v1', host_reference=declaration,
                    inputs=[dict(format='antigens', path='absent.json')])
    manifest.write_text(json.dumps(document))
    resolved = read_input_manifest(manifest)[0][2]['host_reference']
    assert HostReference.load(resolved).coverage == host.coverage
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    declaration['catalogues'] = ['dog-ncorfs:1']
    manifest.write_text(json.dumps(document))
    with pytest.raises(ValueError, match='Requested host ORF catalogues differ'):
        run_cli(['--input-manifest', str(manifest), '--mhc-predictor', 'random',
                 '--output-epitopes', str(tmp_path / 'predictions.tsv')])
    declaration['catalogues'] = None
    with pytest.raises(ValueError, match='explicit list'):
        HostReference.load(dict(declaration, path=tmp_path / 'host.json.gz'))


@pytest.mark.parametrize('historical', [False, True])
def test_native_coverage_replays_without_current_registry(tmp_path, monkeypatch, historical):
    host = HostReference(mouse_reference())
    transport = NativeReferences(tmp_path / 'native.tsv')
    payload = transport.transform(dict(host_reference=host, epitopes=[], provenance=[]))
    if historical:
        payload['host_reference'].pop('coverage')
    monkeypatch.setattr('orfdata.catalog_status', forbidden)
    restored = transport.transform(payload, decode=True)['host_reference']
    if historical:
        assert restored.coverage['catalogues'] is None
        assert 'not recorded' in ' '.join(coverage_view(restored.coverage)['summary'])
    else:
        assert restored.coverage == host.coverage
    resaved = transport.transform(dict(host_reference=restored, epitopes=[], provenance=[]))
    assert resaved['host_reference']['coverage'] == restored.coverage
    assert restored.identity == host.identity


@pytest.mark.parametrize('change', ['inventory', 'digest', 'registry', 'labels', 'unrecorded'])
def test_native_coverage_rejects_tampering_and_malformed_registry(change):
    host = HostReference(mouse_reference())
    record = host.coverage
    if change == 'inventory':
        record['inventory']['extra_orf_entries'] += 1
    elif change == 'digest':
        record['orfdata_version'] = 'different'
    elif change == 'registry':
        record['catalogues'][0]['route'] = 'not a list'
    elif change == 'labels':
        record['known_source_labels'] = [{}]
    else:
        record['orfdata_version'] = None
    if change != 'digest':
        record.pop('sha256')
        record['sha256'] = _digest(record)
    with pytest.raises(ValueError, match='coverage|availability|registry claims'):
        HostReference(host.proteome, coverage=record)


def test_combining_different_coverage_cannot_silently_replace_history():
    host = HostReference(dog_reference())
    assert combine_host_references([host, HostReference(host.proteome, catalogues=[])]).coverage == host.coverage
    old = HostReference(host.proteome, coverage=None)
    assert old.identity == host.identity
    for references in ([host, old], [old, host]):
        with pytest.raises(ValueError, match='conflicting frozen'):
            combine_host_references(references)
    with pytest.raises(ValueError, match='replace frozen'):
        HostReference.load(host, coverage=old.coverage)


@pytest.mark.parametrize('extension', ['html', 'txt'])
def test_zero_candidate_reports_show_missing_dog_coverage(extension):
    report = render(template(host_reference=HostReference(dog_reference())), extension)
    assert 'No candidates passed ranking' in report
    assert 'Host ORF catalogue coverage' in report
    assert '1 annotated protein records; 0 additional ORF records' in report
    assert 'Annotated-only reference' in report
    assert 'registered ncORF coverage is unavailable' in report
    assert 'does not establish absence from unsearched ORFs' in report
    assert 'was not recorded' in render(template(), extension)


@pytest.mark.parametrize('species', ['dog', 'mouse'])
def test_cli_reports_native_coverage_and_replays_offline(inputs, tmp_path, monkeypatch, species):  # noqa: F811
    monkeypatch.setattr('vaxrank.cli.entry_point.resolve_ensembl_release',
                        lambda args: setattr(args, 'genome', None))
    proteome = dog_reference() if species == 'dog' else mouse_reference()
    catalogues = [] if species == 'dog' else ['ncorf-mammals-mouse:1']
    host = HostReference(proteome, catalogues=catalogues)
    host.save(tmp_path / 'host.json.gz')
    replace(inputs, antigens={'viral': inputs.antigens['viral']}, scope=InputScope(
        patient_id='synthetic-' + species, sample_id='synthetic-tumor',
        host_taxon_id=proteome.reference.taxon_id, reference_assembly=proteome.reference.genome_build,
        annotation='ensembl:93', mhc_alleles=(('DLA-88*001:01',) if species == 'dog' else ('H-2-Kb',))
        )).save(tmp_path / 'antigens.json')
    manifest = tmp_path / 'inputs.json'
    manifest.write_text(json.dumps(dict(schema='vaxrank.input_manifest.v1',
        host_reference=dict(path='host.json.gz', identity=host.identity, catalogues=catalogues),
        inputs=[dict(format='antigens', path='antigens.json')])))
    output = tmp_path / 'original'
    run_cli(['--input-manifest', str(manifest), '--mhc-predictor', 'random',
        '--output-dir', str(output), '--pdf-backend', 'weasyprint',
        '--config-text', 'epitopes.score_expr=1.0', '--no-processing-aware-annotation'])
    native = output / 'candidate_predictions.tsv'
    saved = EpitopeDataset.load(native)
    assert saved.host_reference.coverage == host.coverage
    monkeypatch.setattr('orfdata.catalog_status', forbidden)
    monkeypatch.setattr('mhctools.cli.predictors_from_args', forbidden)
    replay = tmp_path / 'replay'
    run_cli(['--input-epitopes', str(native), '--output-dir', str(replay),
             '--pdf-backend', 'weasyprint', '--no-processing-aware-annotation'])
    assert EpitopeDataset.load(replay / native.name).host_reference.coverage == host.coverage
    for directory in (output, replay):
        for filename in ('run_summary.txt', 'vaccine_report.txt'):
            text = (directory / filename).read_text()
            if species == 'dog':
                assert 'registered ncORF coverage is unavailable' in text
                assert 'ncORF catalogues were not searched' in text
            else:
                assert 'Catalogue ncorf-mammals-mouse 1: Included' in text
                assert 'Catalogue riborf-mouse 3: Not searched' in text
                if filename == 'vaccine_report.txt':
                    assert 'Minimum protein length: 8 aa' in text
