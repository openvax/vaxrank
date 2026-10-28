"""Retain actual modality products and their complete assembly inputs."""

from dataclasses import asdict
import json
from pathlib import Path

import msgspec

from .native_serialization import from_native_json, to_native_json


def pipeline_native_dataset(args, variants, epitope_config):
    """Capture original VCF/BAM evidence without requiring a pooling manifest."""
    from pathlib import Path
    from .direct_input import _file_digest
    from .epitope_dataset import EpitopeDataset
    from .input_scope import (
        InputProvenance, InputScope, normalize_alleles, scope_dict,
        combine_scopes, species_scope_from_args,
    )
    from .selection_policy import canonical_digest
    from mhctools.cli import mhc_alleles_from_args

    paths = list(getattr(args, 'vcf', None) or []) + list(getattr(args, 'germline_vcf', None) or [])
    if getattr(args, 'bam', None):
        paths.append(args.bam)
    hashes = {str(Path(path).resolve()): _file_digest(path) for path in paths}
    genomes = {(v.ensembl.reference_name, getattr(v.ensembl, 'release', None)) for v in variants}
    reference, release = next(iter(genomes)) if len(genomes) == 1 else (None, None)
    scope = InputScope(patient_id=getattr(args, 'output_patient_id', None) or None,
                       reference_assembly=reference,
                       annotation='ensembl:' + str(release) if release is not None else None,
                       mhc_alleles=normalize_alleles(mhc_alleles_from_args(args), 'Direct prediction'))
    scope = combine_scopes(scope, species_scope_from_args(args), 'Direct prediction')
    content_hash = canonical_digest(list(hashes.values()))
    provenance = InputProvenance(
        input_id='vcf_bam:' + content_hash + ':' + canonical_digest(scope_dict(scope)),
        source_format='vcf_bam', path='; '.join(paths), content_sha256=content_hash,
        scope=scope, report_declarations=dict(file_sha256=hashes))
    dataset = EpitopeDataset.from_predictions((), config=epitope_config)
    dataset.provenance = (provenance,)
    return dataset


def saved_modalities(dataset):
    from .selection_policy import modality_configuration
    record = dataset.selection.get('run_configuration') or {}
    assembly = dataset.selection.get('construct_assembly', dataset.selection.get('cta_assembly', {}))
    return {
        modality: {**record.get('configuration', {}).get(modality, {}),
                   **modality_configuration(modality, record.get('effective_constructs', {}).get(modality, {})),
                   **assembly.get(modality, {}).get('configuration', {})}
        for modality in ('peptide', 'mrna')}


def resolved_saved_design(args, config, dataset):
    """Compare saved defaults after explicit current YAML and CLI overrides."""
    from .config.loader import extract_construct_kwargs
    from .cli.vaccine_config_args import vaccine_config_from_args
    from .selection_policy import canonical_digest
    cli = getattr(args, '_external_vaccine_args', args)
    explicit = getattr(cli, '_explicit_cli_args', ())
    modalities = saved_modalities(dataset)
    for modality, values in modalities.items():
        values.update(extract_construct_kwargs(config, modality))
        for key in ('antigen_content', 'epitopes_per_antigen'):
            if key in explicit:
                values[key] = getattr(cli, key)
        for attribute in explicit:
            if attribute.startswith(modality + '_'):
                key = attribute[len(modality) + 1:]
                key = {'5p_utr': 'utr_5p', '3p_utr': 'utr_3p'}.get(key, key)
                if key in values:
                    values[key] = getattr(cli, attribute)
    limit = (getattr(cli, 'max_mutations_in_report', None) if 'max_mutations_in_report' in explicit
             else dataset.selection.get('run_configuration', {}).get('effective_design', {}).get('max_ranked_sources'))
    vaccine = msgspec.to_builtins(vaccine_config_from_args(cli, merged_config=config))
    types = (getattr(cli, 'vaccine_type', None) if 'vaccine_type' in explicit else
             dataset.selection.get('run_configuration', {}).get('effective_design', {}).get('vaccine_types'))
    return canonical_digest(dict(vaccine_peptides=vaccine, modalities=modalities, max_ranked_sources=limit,
                                 vaccine_types=types,
                                 metadata={k: v for k, v in config.items()
                                           if k not in ('vaccine_peptides', 'peptide', 'mrna')}))


def initialize_construct_replay(args, dataset, vaccine_config, manufacturability_config=None):
    """Attach one native design to its writers, without acquiring evidence."""
    args._design_dataset = dataset
    from .native_references import NativeReferences
    path = (getattr(args, 'output_epitopes', None)
            or str(Path(getattr(args, 'output_dir', None) or '.') / 'assembly.tsv'))
    dataset._assembly_references = NativeReferences(path, dataset.native_references)
    args._construct_assembly = dataset.selection.setdefault(
        'construct_assembly', dataset.selection.pop('cta_assembly', {}))
    record = dataset.selection.get('run_configuration') or {}
    if ('max_mutations_in_report' not in getattr(args, '_explicit_cli_args', ())
            and getattr(args, 'max_mutations_in_report', None) is None):
        args.max_mutations_in_report = record.get('effective_design', {}).get(
            'max_ranked_sources', getattr(args, '_inherited_source_limit', None))
    inherited = getattr(args, '_inherited_modality_config', {})
    args._saved_modality_config = {modality: {**inherited.get(modality, {}), **values}
                                   for modality, values in saved_modalities(dataset).items()}
    args._construct_policy = json.loads(json.dumps(dict(
        epitopes=msgspec.to_builtins(dataset.config),
        vaccine_peptides=msgspec.to_builtins(vaccine_config),
        manufacturability=msgspec.to_builtins(manufacturability_config)), allow_nan=False))


def _source_graph(ranked, options):
    from .vaccine_library import iter_named_antigens, top_target_epitopes
    names = {id(p): name for name, _, p in iter_named_antigens(ranked, options.candidates_per_slot)}
    def emitted_names(peptide):
        name = names.get(id(peptide))
        if name is None:
            return []
        if options.antigen_content != 'minimal_epitope':
            return [name]
        targets = top_target_epitopes(peptide, n=options.epitopes_per_antigen)
        return [name + ('_epitope' if len(targets) == 1 else '_epitope%d' % (i + 1))
                for i in range(len(targets))]
    return [dict(source=source, candidates=[dict(
        name=names.get(id(p)), antigen=p.antigen, amino_acids=p.amino_acids,
        emitted_antigen_names=emitted_names(p),
        mutation_fragment=p.mutant_protein_fragment,
        epitopes=tuple(p.epitopes), target_epitope_score=p.target_epitope_score,
        window_selection=p.window_selection_audit) for p in peptides])
        for source, peptides in ranked]


def _canonical_native(value):
    if isinstance(value, list):
        return [_canonical_native(item) for item in value]
    if not isinstance(value, dict):
        return value
    result = {key: _canonical_native(item) for key, item in value.items()}
    metadata = result.get('__class__')
    if (metadata == dict(__module__='vaxrank.peptide', __name__='PeptideConstruct')
            and not result.get('antigens')):
        result.pop('antigens', None)
    if metadata == dict(__module__='vaxrank.input_scope', __name__='InputScope'):
        for key in ('host_taxon_id', 'mhc_species', 'host_reference_sha256'):
            if result.get(key) is None:
                result.pop(key, None)
    if metadata == dict(__module__='vaxrank.vaccine_antigen', __name__='SelfReferenceSource'):
        for key in ('orf_id', 'reference_entry', 'reference_occurrence'):
            if result.get(key) in ('', None):
                result.pop(key, None)
    if metadata == dict(__module__='vaxrank.vaccine_antigen', __name__='SelfReferenceMatch'):
        if not result.get('reference_sha256'):
            result.pop('reference_sha256', None)
    prediction_tuple = (metadata == dict(__module__='builtins', __name__='tuple')
                        and result['__value__'] and all(
                            isinstance(item, dict) and item.get('__class__') ==
                            dict(__module__='mhctools.pred', __name__='Prediction')
                            for item in result['__value__']))
    if metadata == dict(__module__='builtins', __name__='set') or prediction_tuple:
        result['__value__'].sort(key=lambda item: json.dumps(item, sort_keys=True, allow_nan=False))
    return result


def design_digest(args, value):
    from .selection_policy import canonical_digest
    references = args._design_dataset._assembly_references
    return canonical_digest(_canonical_native(references.source_identity(to_native_json(value))))


def _assembly_identity(args, ranked, options, context):
    definition = (_source_graph(ranked, options), asdict(options), context, args._construct_policy)
    if args._design_dataset.host_reference is not None:
        definition += (args._design_dataset.host_reference.identity,)
    return design_digest(args, definition)


def _payload_digest(args, products, window_audit, graph):
    from .selection_policy import canonical_digest
    references = args._design_dataset._assembly_references
    return canonical_digest((_canonical_native(references.source_identity(products)), window_audit,
                             _canonical_native(references.source_identity(graph))))


def restore_products(args, ranked, modality, options, context=None):
    """Reuse actual products only for identical sources, policy and settings."""
    assembly = getattr(args, '_construct_assembly', None)
    saved = assembly.get(modality) if assembly is not None else None
    # Earlier CTA snapshots lack the complete identity contract. Their settings
    # remain useful defaults, but cannot establish an identical modern assembly.
    if not saved or saved.get('contract') is None:
        return None
    if saved['contract'] != 1:
        raise ValueError('Unsupported saved construct replay contract')
    if saved['payload_sha256'] != _payload_digest(args,
            saved['products'], saved.get('window_audit'), saved['source_graph']):
        raise ValueError('Saved construct products do not match their checksum')
    if saved['identity'] == _assembly_identity(args, ranked, options, context):
        products = from_native_json(saved['products'], list)
        # Older snapshots can contain the duplicate-slot bug. Rebuild those
        # designs under the corrected policy, retaining their input evidence.
        # Unique historical products remain exactly replayable offline.
        if products and not any(p.antigen_sources for p in products):
            if modality == 'mrna':
                sequences = [a['aa'] for p in products for a in p.antigens]
            elif options.antigens_per_construct <= 1:
                sequences = [p.sequence for p in products]
            else:
                from .assembly_antigens import antigen_occurrences
                selected_names = {name for p in products for name in p.antigen_names}
                sequences = [a.amino_acids for a in antigen_occurrences(ranked, options, modality)
                             if a.name in selected_names]
            if len(sequences) != len(set(sequences)):
                return None
        return products
    return None


def record_products(args, ranked, modality, options, products, context=None, window_audit=None):
    assembly = getattr(args, '_construct_assembly', None)
    if assembly is None:
        return
    from .selection_policy import modality_configuration
    configuration = modality_configuration(modality, asdict(options))
    configuration.update(context or {})
    payload = to_native_json(products)
    window_audit = json.loads(json.dumps(window_audit, allow_nan=False))
    inputs = _source_graph(ranked, options)
    if any(p.antigen_sources for p in products):
        mappings = [(name, source) for product in products
                    for name, sources in product.antigen_sources.items() for source in sources]
        for row in inputs:
            for candidate in row['candidates']:
                antigen = candidate['antigen']
                prediction_ids = {e.prediction_group_source for e in candidate['epitopes']}
                input_ids = {e.input_provenance.input_id for e in candidate['epitopes'] if e.input_provenance}
                # Window refinement can crop a candidate and renumber its
                # alternates. Link by retained source/evidence identities and
                # exact contained sequence, not the mutable display label.
                candidate['emitted_antigen_names'] = sorted({name for name, source in mappings
                    if source['source_identifier'] == antigen.source_identifier
                    and source['gene_name'] == antigen.gene_name
                    and set(source['input_ids']) <= input_ids
                    and set(source['prediction_ids']) <= prediction_ids
                    and source['candidate_sequence'] in candidate['amino_acids']})
    graph = to_native_json(dict(
        inputs=inputs,
        products={p.name: list(p.antigen_names) for p in products}))
    assembly[modality] = dict(
        contract=1, identity=_assembly_identity(args, ranked, options, context),
        configuration=configuration, products=payload, window_audit=window_audit,
        source_graph=graph, payload_sha256=_payload_digest(args, payload, window_audit, graph))
