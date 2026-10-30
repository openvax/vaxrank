"""Frozen coverage of a supplied host inventory, without acquiring catalogues."""

from collections import Counter
from copy import deepcopy
import hashlib
import json


DISCOVER_CATALOGUES = object()
SCHEMA = 'vaxrank.host_catalogue_coverage.v1'


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def _inventory(proteome):
    counts = Counter(e.kind for e in proteome)
    sources = Counter(e.source for e in proteome if e.kind == 'extra_orf')
    layers = [deepcopy(dict(layer)) for layer in proteome.layers
              if layer.get('role') == 'host' and layer.get('kind') == 'extra_orf']
    labels = sorted(set(sources) | {layer['source'] for layer in layers})
    return dict(host_reference_sha256=proteome.identity,
        taxon_id=proteome.reference.taxon_id, species=proteome.reference.species,
        genome_build=proteome.reference.genome_build,
        annotated_entries=counts['annotated'], extra_orf_entries=counts['extra_orf'],
        included=[dict(source=source, retained_entries=sources[source],
                       layers=[layer for layer in layers if layer['source'] == source]) for source in labels])


def validate_catalogues(catalogues, proteome):
    """A declaration asserts the exact source labels already in the snapshot."""
    if catalogues is None:
        return
    if (not isinstance(catalogues, (list, tuple))
            or any(not isinstance(x, str) or not x.strip() or x != x.strip() for x in catalogues)
            or len(set(catalogues)) != len(catalogues)):
        raise ValueError('host_reference catalogues must be a list of distinct nonempty source labels')
    inventory = _inventory(proteome)
    actual = {row['source'] for row in inventory['included']}
    if set(catalogues) != actual or any(not row['layers'] for row in inventory['included']):
        raise ValueError('Requested host ORF catalogues differ from the supplied reference; '
                         f'requested {sorted(catalogues)}, recorded {sorted(actual)}. '
                         'Compose the required catalogues with ORFdata before prediction.')


def _catalogues(reference):
    if reference.genome_build is None:
        return None, None
    from orfdata import catalog_status
    from orfdata.catalogs import CatalogRegistry
    from orfdata.reference import species_taxon

    registry = CatalogRegistry()
    rows, known_labels = [], set()
    for status in catalog_status(registry=registry):
        if species_taxon(status.species) != reference.taxon_id:
            continue
        route = next((r for r in status.lifts_to if r['genome_build'] == reference.genome_build), None)
        if status.genome_build != reference.genome_build and route is None:
            continue
        # ORFCatalog.source exposes these exact release labels. Opaque merged
        # or user-supplied labels are retained, never split or guessed here.
        label = (f'orfdata-catalog:{status.version}' if status.name == 'bundled'
                 else f'{status.name}:{status.version}')
        if status.genome_build != reference.genome_build:
            label += '@' + reference.genome_build
        known_labels.add(label)
        if status.name != 'bundled' and status.version != registry.info(status.name).version:
            continue
        rows.append(dict(name=status.name, version=status.version, source=label,
                         genome_build=status.genome_build, records=status.records,
                         route=list(route['chains']) if route is not None else [],
                         acquire=status.acquire, optional=status.optional))
    return rows, sorted(known_labels)


def reference_coverage(proteome, recorded=DISCOVER_CATALOGUES):
    """Capture once, or restore without consulting the current registry."""
    inventory = _inventory(proteome)
    if recorded is DISCOVER_CATALOGUES:
        from orfdata import __version__
        catalogues, known_labels = _catalogues(proteome.reference)
        record = dict(schema=SCHEMA, inventory=inventory, orfdata_version=__version__,
                      catalogues=catalogues, known_source_labels=known_labels)
    elif recorded is None:
        record = dict(schema=SCHEMA, inventory=inventory, orfdata_version=None,
                      catalogues=None, known_source_labels=None)
    else:
        record = deepcopy(recorded)
        if (not isinstance(record, dict) or set(record) != {
                'schema', 'inventory', 'orfdata_version', 'catalogues', 'known_source_labels', 'sha256'}
                or record['schema'] != SCHEMA or record['inventory'] != inventory):
            raise ValueError('Host catalogue coverage differs from its reference inventory')
        digest = record.pop('sha256')
        if digest != _digest(record):
            raise ValueError('Host catalogue coverage checksum differs')
    _validate_availability(record)
    record['sha256'] = _digest(record)
    return record


def _validate_availability(record):
    version, rows, labels = (record[k] for k in ('orfdata_version', 'catalogues', 'known_source_labels'))
    if version is None:
        if rows is not None or labels is not None:
            raise ValueError('Unrecorded catalogue availability cannot contain registry claims')
        return
    def text(value):
        return isinstance(value, str) and bool(value.strip())
    if not text(version):
        raise ValueError('Malformed frozen host catalogue availability')
    if rows is None and labels is None:
        return
    if (not isinstance(rows, list) or not isinstance(labels, list)
            or any(not text(label) for label in labels) or len(set(labels)) != len(labels)):
        raise ValueError('Malformed frozen host catalogue availability')
    names = set()
    for row in rows:
        if (not isinstance(row, dict) or set(row) != {
                'name', 'version', 'source', 'genome_build', 'records', 'route', 'acquire', 'optional'}
                or any(not text(row[k]) for k in ('name', 'version', 'source', 'genome_build'))
                or (row['acquire'] is not None and not text(row['acquire']))
                or row['name'] in names or row['source'] not in labels
                or type(row['optional']) is not bool
                or (row['records'] is not None and (type(row['records']) is not int or row['records'] < 0))
                or not isinstance(row['route'], list) or any(not text(chain) for chain in row['route'])):
            raise ValueError('Malformed frozen host catalogue availability')
        names.add(row['name'])


def _selection(layers):
    if not layers:
        return 'Catalogue selection provenance was not recorded'
    parts = []
    for layer in layers:
        counts = ', '.join(f'{key.replace("_", " ")}: {layer[key]}' for key in (
            'selected', 'added', 'excluded', 'unlifted', 'not_on_source_transcript',
            'identical_to_annotated', 'length_filtered') if key in layer)
        if counts:
            parts.append('At composition: ' + counts)
        if 'selection' in layer:
            selection = layer['selection']
            if isinstance(selection, dict):
                selection = {key: value for key, value in selection.items() if key != 'input'}
            parts.append('Source preselection: ' + json.dumps(selection, sort_keys=True))
        if 'evidence' in layer:
            label = 'Additional evidence filter' if 'selection' in layer else 'Evidence filter'
            parts.append(label + ': ' + (json.dumps(layer['evidence'], sort_keys=True)
                                        if layer['evidence'] is not None else 'none'))
        if 'min_protein_length' in layer:
            parts.append(f"Minimum protein length: {layer['min_protein_length']} aa")
    parts.append('Full source and selection metadata are retained in the native dataset.')
    return '; '.join(parts)


def coverage_view(record):
    """One view for reports and summaries; availability never implies inclusion."""
    if record is None:
        return dict(summary=['Host ORF catalogue coverage was not recorded for this run.'],
                    included=[], catalogues=[], identity=None)
    inventory = record['inventory']
    included = inventory['included']
    known = set(record['known_source_labels'] or ())
    labels = {r['source'] for r in included if r['layers']}
    unresolved = any(not r['layers'] or r['source'] not in known for r in included)
    summary = [f"Host: {inventory['species']} / {inventory['genome_build'] or 'unspecified assembly'}.",
        f"Reference inventory: {inventory['annotated_entries']} annotated protein records; "
        f"{inventory['extra_orf_entries']} additional ORF records."]
    if not included:
        summary.append('Annotated-only reference: ncORF catalogues were not searched.')
    elif not inventory['extra_orf_entries']:
        summary.append('Catalogue layers were selected, but no additional ORF sequences were retained.')
    if record['catalogues'] is None:
        summary.append('Host assembly is unspecified; compatible catalogue availability is unknown.'
                       if inventory['genome_build'] is None else
                       'Catalogue availability was not recorded when this native reference was saved.')
    elif not record['catalogues']:
        summary.append(f"ORFdata {record['orfdata_version']} has no registered ncORF catalogue "
                       'compatible with this host species and assembly; registered ncORF coverage is unavailable.')
    else:
        summary.append(f"Catalogue availability recorded with ORFdata {record['orfdata_version']}; "
                       'registration does not mean installation or inclusion.')
    if unresolved:
        summary.append('Some included sources cannot be attributed to registered catalogue releases; '
                       'coverage of other catalogues is not established.')
    rows = []
    for row in record['catalogues'] or ():
        state = ('Included' if row['source'] in labels else
                 'Inclusion not established' if unresolved else 'Not searched')
        rows.append(dict(row, status=state))
    return dict(summary=summary, identity=inventory['host_reference_sha256'], catalogues=rows,
        included=[dict(source=r['source'], retained_entries=r['retained_entries'],
                       selection=_selection(r['layers'])) for r in included])


def coverage_lines(record):
    view = coverage_view(record)
    lines = list(view['summary'])
    lines.extend(f"Included source: {r['source']}; retained ORF records: {r['retained_entries']}"
                 for r in view['included'])
    lines.extend(f"Catalogue {r['name']} {r['version']}: {r['status']}" for r in view['catalogues'])
    return lines
