"""Source-namespace adapters for the common CTA expression input."""

import re

from .cta_identity import identity_digest


def require_common_cta_data():
    from oncoref import __version__
    from packaging.version import Version
    if Version(__version__) < Version('1.23.1'):
        raise ValueError("Common CTA restriction requires OncoRef >=1.23.1; install 'vaxrank[cta]' in a compatible environment")


def reference_namespace(reference):
    namespaces = reference['identifier_namespaces']
    if namespaces == dict(gene='Ensembl gene', transcript='Ensembl transcript', protein='Ensembl protein'):
        return 'ensembl'
    if namespaces in (
            dict(gene='NCBI Gene', transcript='NCBI RefSeq RNA', protein='NCBI RefSeq protein'),
            dict(gene='NCBI Gene', transcript='RefSeq', protein='RefSeq')):
        return 'ncbi'
    raise ValueError('CTA reference has unsupported identifier namespaces')


def feature_id(value, level, namespace):
    if namespace == 'ncbi':
        pattern = r'[1-9][0-9]*' if level == 'gene' else r'[NX]M_[0-9]+\.[1-9][0-9]*'
        expected = 'NCBI Gene ID or versioned RefSeq transcript'
    elif namespace == 'ensembl':
        pattern = r'ENS[A-Z]*G[0-9]+' if level == 'gene' else r'ENS[A-Z]*T[0-9]+\.[1-9][0-9]*'
        expected = 'Ensembl gene ID or exact versioned Ensembl transcript'
    else:
        raise ValueError('Unsupported CTA expression namespace')
    if not isinstance(value, str) or not re.fullmatch(pattern, value):
        raise ValueError(f'Expected an {expected} for {level}, got {value!r}')
    return value


def _reference_roles(reference):
    namespace = reference_namespace(reference)
    roles = ('protein', 'feature_table') if namespace == 'ncbi' else ('protein', 'gtf', 'rna', 'ncrna')
    if not set(roles) <= set(reference['content_provenance']):
        raise ValueError('CTA evidence requires a complete annotated protein inventory reference; candidate subsets are insufficient')
    return roles


def _reference_context(identity, species):
    """Use provider metadata labels while requiring the owner's exact identity."""
    from orfdata import ReferenceContext

    if reference_namespace(identity) == 'ensembl':
        from pyensembl import EnsemblRelease
        # Metadata only: constructing a release does not download or index it.
        reference = ReferenceContext.from_genome(
            EnsemblRelease(int(identity['annotation_release']), species=species))
    else:
        reference = ReferenceContext(identity['taxon_id'], species,
            identity['assembly_accession'], identity['annotation_source'], identity['annotation_release'],
            assembly_accession=identity['assembly_accession'])
    _validate_reference_identity(reference, identity)
    return reference


def _validate_reference_identity(reference, identity):
    expected = (identity['taxon_id'], identity['assembly_accession'],
                identity['annotation_source'].casefold(), identity['annotation_release'])
    actual = (reference.taxon_id, reference.assembly_accession,
              (reference.annotation_name or '').casefold(), reference.annotation_version)
    if actual != expected:
        raise ValueError('Host reference differs from the species CTA evidence annotation')


def build_cta_host_reference(species, *, evidence_snapshot=None, cache_root=None,
                             local_overrides=None, offline=False, max_bytes=0):
    """Build the snapshot's annotated host using verified OncoRef source assets.

    The transfer budget covers all missing required assets, before acquisition.
    A zero budget permits verified cache reuse. No annotation indexes are built.
    Extra ORFs must be selected and composed separately with ORFdata.
    """
    from oncoref import for_species
    from oncoref.species_assets import inspect_asset
    from orfdata import Proteome
    from .host_reference import HostReference

    require_common_cta_data()
    if type(max_bytes) is not int or max_bytes < 0:
        raise ValueError('max_bytes must be a nonnegative integer transfer budget')
    collection = for_species(species)
    snapshot = collection.antigens.cta_snapshot(path=evidence_snapshot)
    references = snapshot.as_dict()['references']
    if len(references) != 1:
        raise ValueError('CTA host construction requires one explicit mapped reference')
    identity = references[0]['identity']
    if identity['taxon_id'] != collection.species.taxon_id:
        raise ValueError('CTA snapshot reference differs from requested species')
    roles = _reference_roles(identity)
    registered = next((r for r in collection.references.registered_references()
                       if r['identity'] == identity), None)
    if registered is None:
        raise ValueError('CTA reference has no matching registered source assets')
    reference = _reference_context(identity, collection.species.name)
    assets = {a.id: a for a in collection.references.assets()}
    selected = {role: assets[registered['assets'][role]] for role in roles}
    options = dict(cache_root=cache_root, local_overrides=local_overrides)
    needed = sum(a.size_bytes for a in selected.values()
                 if inspect_asset(a, **options)['status'] != 'available')
    if not offline and needed > max_bytes:
        raise ValueError(f'CTA host reference needs {needed} bytes; transfer budget is {max_bytes}')
    paths = {role: collection.references.ensure(asset.id, **options, offline=offline, max_bytes=max_bytes)
             for role, asset in selected.items()}
    pins = {role: identity['content_provenance'][role] for role in roles}
    if reference_namespace(identity) == 'ensembl':
        proteome = Proteome.from_ensembl(paths['protein'], paths['gtf'], paths['rna'], paths['ncrna'],
            reference=reference, expected_sha256=pins)
    else:
        proteome = Proteome.from_refseq(paths['protein'], paths['feature_table'],
            reference=reference, expected_sha256=pins)
    host = HostReference(proteome)
    validate_host_mapping(host, identity)
    return host


def validate_host_mapping(host_reference, mapped_reference):
    """Bind one complete annotated host layer to the source's pinned assets."""
    reference = host_reference.proteome.reference
    _validate_reference_identity(reference, mapped_reference)
    namespace = reference_namespace(mapped_reference)
    roles = _reference_roles(mapped_reference)
    assets = {role: mapped_reference['content_provenance'][role] for role in roles}
    layers = [layer for layer in host_reference.proteome.layers
              if layer.get('kind') == 'annotated' and layer.get('role') == 'host']
    if (len(layers) != 1 or dict(layers[0].get('content_provenance', {})) != assets
            or layers[0].get('identifier_namespaces') != mapped_reference['identifier_namespaces']
            or reference.source_sha256 != identity_digest(assets)):
        raise ValueError('CTA host reference requires the mapped source SHA-256 hashes and identifier namespaces')
    inventory = [(e.accession, e.sequence) for e in host_reference.proteome if e.kind == 'annotated']
    if (reference.sequence_sha256 != identity_digest(sorted(inventory))
            or layers[0].get('added') != len(inventory)):
        raise ValueError('CTA input requires the complete annotated protein inventory')
    return namespace
