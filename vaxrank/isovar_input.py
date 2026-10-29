"""Run or replay Isovar sources through the native Vaxrank input pipeline."""

from contextlib import ExitStack
from importlib.metadata import version
from pathlib import Path
from types import SimpleNamespace
import gzip
import json

import msgspec

from .input_scope import InputScope, combine_scopes, scope_from_mapping
from .isovar_antigens import RESULT_FAMILIES, antigens_from_isovar
from .isovar_options import FAMILIES, IsovarSourceOptions, merge_options, producer_options, source_options
from .native_serialization import _json_loads


ISOVAR_INPUT_SCHEMA = 'vaxrank.isovar_input.v1'


def _read_json(path):
    path = Path(path)
    content = gzip.decompress(path.read_bytes()).decode() if path.suffix == '.gz' else path.read_text()
    return _json_loads(content)


def _configured_read_collector(family, settings):
    """Validate controls against the owner and preserve its family defaults."""
    from isovar import ReadCollector

    excluded = ('read_filter', 'read_end_profile')
    baseline = ReadCollector.for_rna_paths() if family in ('viral', 'isoforms') else ReadCollector()
    names = producer_options(ReadCollector, {}, exclude=excluded)
    resolved = producer_options(ReadCollector, settings, exclude=excluded,
                                defaults={name: getattr(baseline, name) for name in names})
    return ReadCollector(**resolved), resolved


def reconstruct_isovar_antigens(family, *, inputs, options=None, scope=None, species='',
                                tumor_specificity=None, candidate_attestations=None):
    """Call a public Isovar producer and adapt its complete result.

    ``inputs`` contains the selected function's typed data arguments: for
    example a Varcode Variant, ReadEvidence and ORFdata TranscriptORFs for
    ``variant_orfs``. A BAM plus ReadCollector controls can replace supplied
    ReadEvidence. Functions, controls and evidence have the same names as the
    CLI source-file contract; no global Isovar state is modified.
    """
    import isovar
    if family not in FAMILIES:
        raise ValueError('Unknown Isovar source family: ' + str(family))
    options = source_options(options)
    data = dict(inputs)
    read_families = ('variant_orfs', 'sv', 'viral', 'isoforms')
    if options.reads and family not in read_families:
        raise ValueError('ReadCollector controls apply only to variant_orfs, sv, viral and isoforms')
    if 'read_collector' in data:
        raise ValueError('Configure ReadCollector through options.reads, or reconstruct in Isovar and replay its result')
    if options.protein and family != 'variant_orfs':
        raise ValueError('ProteinSequenceCreator controls apply only to variant_orfs')
    resolved = dict(admission=msgspec.to_builtins(options.admission))
    if family in read_families:
        collector, reads = _configured_read_collector(family, options.reads)
        resolved['reads'] = reads
        if family != 'variant_orfs':
            data['read_collector'] = collector
        elif 'read_evidence' not in data:
            data['read_evidence'] = collector.read_evidence_for_variant(data['variant'], data.pop('bam'))
        elif 'bam' in data or options.reads:
            raise ValueError('Supplied ReadEvidence cannot also use BAM/ReadCollector controls')
    functions = dict(variant_orfs=isovar.reconstruct_variant_orfs, sv=isovar.reconstruct_sv_rna,
                     fusion=isovar.reconstruct_fusion, viral=isovar.reconstruct_viral_sequences,
                     isoforms=isovar.reconstruct_antigen_isoforms,
                     protein_hypotheses=isovar.export_protein_hypotheses)
    if family == 'variant_orfs':
        protein = producer_options(isovar.ProteinSequenceCreator, options.protein,
                                   defaults={'max_protein_sequences_per_variant': 0},
                                   nullable=('max_protein_sequences_per_variant',))
        if protein['max_protein_sequences_per_variant'] not in (0, None):
            raise ValueError('variant_orfs retains all reference models; protein count caps do not apply')
        data['protein_sequence_creator'] = isovar.ProteinSequenceCreator(**protein)
        resolved['protein'] = protein
    function = functions[family]
    # Data objects and paths have their own namespace, even when Isovar makes
    # them optional. This prevents ignored controls and accidental overrides.
    data_keys = ('bam', 'variant', 'read_evidence', 'reference_orfs', 'fusion', 'references', 'reads',
                 'targets', 'sample_id', 'source', 'genome_fasta', 'contigs', 'reference_identity',
                 'event_id', 'reference_name', 'donor', 'acceptor', 'regions', 'event_provenance',
                 'read_collector', 'protein_sequence_creator', 'isovar_results', 'alignment_header',
                 'library_index', 'read_filter_identity',
                 'cell_labels', 'base_quality_policy', 'rna_support_scopes', 'scratch_dir', 'reference_cache')
    controls = producer_options(function, options.reconstruction, exclude=data_keys,
                                nullable=('reference_profile',) if family == 'viral' else ())
    overlap = set(data).intersection(controls)
    if overlap:
        raise ValueError('Reconstruction controls belong in options, not inputs: ' + ', '.join(sorted(overlap)))
    resolved['reconstruction'] = controls
    result = function(**data, **controls)
    dataset = antigens_from_isovar(result, scope=scope, species=species,
        tumor_specificity=tumor_specificity, candidate_attestations=candidate_attestations,
        options=IsovarSourceOptions(admission=options.admission))
    for evidence in dataset.source_evidence.values():
        evidence['reconstruction'] = dict(controls=resolved,
            versions={name: version(name) for name in ('isovar', 'varcode', 'orfdata')})
    return dataset


def _resolve_path(value, base):
    if not isinstance(value, str) or not value.strip():
        raise ValueError('Isovar input paths must be nonempty strings')
    return str((base / value).expanduser().resolve())


def _reference_orfs(entries, genome, base):
    from orfdata import TranscriptORF, load_catalog
    if genome is None:
        raise ValueError('variant_orfs requests require a selected host annotation')
    references = []
    for entry in entries:
        entry = dict(entry)
        transcript = genome.transcript_by_id(entry.pop('transcript_id'))
        if 'catalog' in entry:
            path = _resolve_path(entry.pop('catalog'), base)
            checksum = entry.pop('expected_sha256')
            if entry:
                raise ValueError('Catalog ORFs accept transcript_id, catalog and expected_sha256 only')
            references.extend(load_catalog(path, expected_sha256=checksum).for_transcript(transcript))
        elif entry.get('annotated') is True:
            if set(entry) != {'annotated'}:
                raise ValueError('Annotated ORF requires only transcript_id and annotated: true')
            references.append(TranscriptORF.from_annotation(transcript))
        else:
            # Source-defined positions are validated and translated by ORFdata;
            # optional expected protein/source evidence must agree with them.
            references.append(TranscriptORF.at(transcript, **entry))
    return references


def _request_inputs(family, inputs, genome, base, stack):
    import isovar
    from isovar.cli.rna_args import alignment_file_from_args
    from isovar.default_parameters import NUM_RNA_DECOMPRESSION_THREADS
    data = dict(inputs)
    if family == 'protein_hypotheses':
        raise ValueError('Replay a protein_hypotheses result file; live exports accept typed IsovarResult objects in Python')
    if family == 'fusion':
        from isovar.fusion import fusion_from_dict
        fusion, references, reads = fusion_from_dict(data)
        return dict(fusion=fusion, references=references, reads=reads)
    bam_path = _resolve_path(data.pop('bam'), base)
    bam = stack.enter_context(alignment_file_from_args(SimpleNamespace(
        bam=bam_path, num_rna_decompression_threads=data.pop(
            'num_rna_decompression_threads', NUM_RNA_DECOMPRESSION_THREADS))))
    data['bam'] = bam
    if family == 'variant_orfs':
        from varcode import Variant
        data['reference_orfs'] = _reference_orfs(data['reference_orfs'], genome, base)
        data['variant'] = Variant(ensembl=genome, **data['variant'])
    elif family == 'sv':
        from isovar.sv_rna import sv_rna_input_from_dict
        event = data.pop('event')
        if isinstance(event, str):
            event = _read_json(_resolve_path(event, base))
        event = sv_rna_input_from_dict(event)
        if set(data).intersection(event):
            raise ValueError('Declare SV event identity once inside event')
        data.update(event)
    elif family == 'viral':
        data['genome_fasta'] = _resolve_path(data['genome_fasta'], base)
    elif family == 'isoforms':
        data['targets'] = [isovar.IsoformTarget(**{
            **row, 'reference': isovar.FusionReference(**row['reference'])}) for row in data['targets']]
    for key in ('scratch_dir', 'reference_cache'):
        if key in data:
            data[key] = _resolve_path(data[key], base)
    return data


def load_isovar_input(path, *, scope=None, genome=None, options=None, configured_sources=None):
    """Read one source request/result, resolving file paths beside its envelope.

    Raw producer JSON is also accepted, with default held-out admission. The
    envelope adds scope, species, independent attestations and options. Saved
    results never run reconstruction or consult the current ORF annotation.
    """
    path = Path(path).resolve()
    document = _read_json(path)
    if not isinstance(document, dict):
        raise ValueError('Isovar input must be an object')
    if document.get('schema') in RESULT_FAMILIES:
        document = dict(schema=ISOVAR_INPUT_SCHEMA, result=document)
    if document.get('schema') != ISOVAR_INPUT_SCHEMA:
        raise ValueError('Expected ' + ISOVAR_INPUT_SCHEMA + ' or a supported Isovar result')
    allowed = {'schema', 'scope', 'species', 'tumor_specificity', 'candidate_attestations',
               'options', 'result', 'reconstruct'}
    if set(document) - allowed:
        raise ValueError('Unknown Isovar input fields: ' + ', '.join(sorted(set(document) - allowed)))
    if ('result' in document) == ('reconstruct' in document):
        raise ValueError('Isovar input requires exactly one of result or reconstruct')
    scope = combine_scopes(scope or InputScope(), scope_from_mapping(document.get('scope', {}), str(path)), str(path))
    adapter = dict(scope=scope, species=document.get('species', ''),
                   tumor_specificity=document.get('tumor_specificity'),
                   candidate_attestations=document.get('candidate_attestations'))
    sources = configured_sources or {}
    if set(sources) - set(FAMILIES):
        raise ValueError('Unknown Isovar configured source families: ' + ', '.join(sorted(set(sources) - set(FAMILIES))))
    if 'result' in document:
        result = document['result']
        if isinstance(result, str):
            result = _read_json(_resolve_path(result, path.parent))
        family = RESULT_FAMILIES.get(result.get('schema'))
        resolved = merge_options(document.get('options'), sources.get(family), options)
        dataset = antigens_from_isovar(result, options=resolved, **adapter)
    else:
        request = document['reconstruct']
        if set(request) != {'family', 'inputs'} or request['family'] not in FAMILIES:
            raise ValueError('Isovar reconstruct requires a known family and inputs')
        family = request['family']
        resolved = merge_options(document.get('options'), sources.get(family), options)
        raw_inputs, base = request['inputs'], path.parent
        if isinstance(raw_inputs, str):
            input_path = Path(_resolve_path(raw_inputs, base))
            raw_inputs, base = _read_json(input_path), input_path.parent
        with ExitStack() as stack:
            data = _request_inputs(family, raw_inputs, genome, base, stack)
            dataset = reconstruct_isovar_antigens(family, inputs=data, options=resolved, **adapter)
        for evidence in dataset.source_evidence.values():
            evidence['reconstruction']['request'] = json.loads(json.dumps(request))
    return dataset
