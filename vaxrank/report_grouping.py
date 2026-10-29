"""Group report displays while preserving every source and measurement."""

from dataclasses import replace
from hashlib import sha256
import json

from .coding_context import isovar_candidates
from .native_serialization import to_native_json


# Move these opaque identities into the source directory, never discard them.
# Biological reference identity, ORF identity and transcript remain in context.
ORIGIN_LABELS = {'Source identifier', 'RNA candidate identity', 'RNA sequence identity', 'RNA path'}


def _key(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def window_evidence_key(peptide, processing):
    """Exact evidence, including flanks/comparators and values before rounding.

    Candidate prediction-group IDs do not define biological equivalence; they
    remain unchanged in native data. Source names move to the source directory.
    Reference-comparator provenance, input
    evidence, support, source sequences and offsets are not guessed equivalent.
    """
    epitopes = [replace(e, prediction_id='', source_name='', predictions=[replace(p, source_sequence_name='')
        for p in e.predictions_flat()]) for e in peptide.epitopes]
    return _key(to_native_json(dict(sequence=peptide.amino_acids,
        target_mask=peptide.antigen.targetable_mask, epitopes=epitopes,
        processing=processing, score=peptide.combined_score,
        target_score=peptide.target_epitope_score, scoring_expression=peptide.combined_score_expr,
        manufacturability=peptide.manufacturability_scores._asdict())))


def prediction_table_key(peptide, processing):
    """Same recorded prediction values; occurrence evidence stays in its group.

    Coordinates and parent sequences are displayed in the linked window table.
    Keep predictor-specific inputs/flanks, comparator provenance, raw values,
    model identities and processing context in the shared prediction key.
    """
    epitopes = [replace(e, prediction_id='', input_provenance=None, input_evidence=(),
        source_sequence='', source_name='', offset=0,
        predictions=[replace(p, source_sequence_name='') for p in e.predictions_flat()])
        for e in peptide.epitopes]
    return _key(to_native_json(dict(epitopes=epitopes, processing=processing)))


def window_measurements_key(peptide):
    return _key(to_native_json(dict(sequence=peptide.amino_acids,
        target_mask=peptide.antigen.targetable_mask, score=peptide.combined_score,
        target_score=peptide.target_epitope_score, expression=peptide.combined_score_expr,
        manufacturability=peptide.manufacturability_scores._asdict())))


def window_contexts(peptide):
    contexts = []
    for e in peptide.epitopes:
        context = next((c for c in contexts if c['sequence'] == e.source_sequence), None)
        if context is None:
            context = dict(sequence=e.source_sequence, occurrences=[])
            contexts.append(context)
        row = dict(peptide=e.sequence, offset=e.offset, n_flank=e.n_flank, c_flank=e.c_flank)
        if row not in context['occurrences']:
            context['occurrences'].append(row)
    return contexts


def prediction_origins(peptides):
    """Retain the occurrence identifiers omitted from shared evidence keys."""
    rows = []
    for peptide in peptides:
        for epitope in peptide.epitopes:
            fields = [('Epitope source name', epitope.source_name)]
            fields.extend(('Predictor source name', p.source_sequence_name)
                          for p in epitope.predictions_flat())
            for label, value in fields:
                if value and (label, value) not in rows:
                    rows.append((label, value))
    return rows


def source_evidence_key(antigen):
    candidates = []
    for candidate in isovar_candidates(antigen):
        # These IDs stay in the source directory. Keep all measured support,
        # coding/reference context and any future producer evidence fields.
        candidates.append({k: v for k, v in candidate.items()
                           if k not in {'model_id', 'rna_candidate_id', 'sequence_id'}})
    context = replace(antigen, source_identifier='', source_metadata=tuple(
        (name, value) for name, value in antigen.source_metadata if name != 'isovar_candidate'))
    return _key(to_native_json(dict(antigen=context, candidates=candidates)))


def grouped_candidate_details(variants):
    """Build compact candidate groups and an exact shared-table inventory."""
    groups, tables, by_key, table_keys = [], [], {}, {}

    def intern(kind, value, source_ids, guard=None):
        key = _key([kind, value, guard])
        if key not in table_keys:
            table_keys[key] = len(tables)
            tables.append(dict(identifier='E%d' % (len(tables) + 1), kind=kind,
                               value=value, sources=[]))
        table = tables[table_keys[key]]
        table['sources'] += [identity for identity in source_ids if identity not in table['sources']]
        return table['identifier']

    for variant in variants:
        source_id = variant['source_id']
        origin = dict(identifier=source_id, full_identifier=variant['short_description'],
                      rank=variant['num'], rows=list(variant['prediction_origins']))

        def split_rows(rows):
            kept = []
            for label, value in rows:
                if label in ORIGIN_LABELS:
                    pair = (label, value)
                    if pair not in origin['rows']:
                        origin['rows'].append(pair)
                else:
                    kept.append((label, value))
            return kept

        sections = [dict(title=variant['entity_label'], rows=split_rows(variant['variant_data'].items()),
                         note=variant['admission_note']),
                    dict(title=variant['effect_label'], rows=split_rows(variant['effect_data'].items()), note='')]
        sections += [dict(s, rows=split_rows(s['rows'])) for s in variant['source_evidence_sections']]
        key = _key([sections, variant['peptides'], variant['databases'], variant['_evidence_key']])
        if key not in by_key:
            by_key[key] = len(groups)
            groups.append(dict(identifier='G%d' % (len(groups) + 1),
                gene=variant['variant_data'].get('Gene name', ''), sources=[],
                sections=sections, peptides=variant['peptides'], databases=variant['databases']))
        groups[by_key[key]]['sources'].append(origin)

    for group in groups:
        source_ids = [source['identifier'] for source in group['sources']]
        group['evidence'] = [intern('section', s, source_ids) for s in group.pop('sections')]
        group['windows'] = []
        group['predictions'] = []
        for peptide in group.pop('peptides'):
            prediction_fields = {'epitopes', 'epitope_displays', 'ascii_epitopes', 'wt_epitopes', 'ascii_wt_epitopes'}
            predictions = {k: v for k, v in peptide.items() if k in prediction_fields}
            prediction_id = intern('predictions', predictions, source_ids, peptide['_prediction_key'])
            if prediction_id not in group['predictions']:
                group['predictions'].append(prediction_id)
            window = {k: v for k, v in peptide.items() if k not in prediction_fields and not k.startswith('_')}
            window['prediction_table'] = prediction_id
            group['windows'].append(intern('window', window, source_ids, peptide['_measurement_key']))
    return groups, tables


def group_report_evidence(data):
    """Connect short source IDs, grouped outcomes and shared evidence tables."""
    directory, identities = [], {}

    def add_source(key, identity, gene='', kind='', context='', origin='Not recorded'):
        if key not in identities:
            entry = dict(identifier='S%d' % (len(directory) + 1), full_identifier=identity,
                         gene=gene, kind=kind, context=context, origin=origin, rows=[], groups=[])
            identities[key] = entry
            directory.append(entry)
        return identities[key]['identifier']

    for row in data['ranked_candidates']:
        row['source_id'] = add_source(row['source_key'], row['identifier'], row['gene'], row['kind'])
    for row in data['input_coverage']['sources']:
        key = row.get('source_key') or _key(row)
        row['source_id'] = add_source(key, row['identity'], row['gene'], row['family'],
                                      row['context'], row['origin'])
        identities[key].update(context=row['context'], origin=row['origin'])
    for variant in data['variants']:
        variant['source_id'] = add_source(variant['source_key'], variant['short_description'],
                                         variant['variant_data'].get('Gene name', ''))
    for group in data['product_groups']:
        for product in group['products']:
            for antigen in product['antigens']:
                antigen['source_ids'] = list(dict.fromkeys(identities[key]['identifier']
                    for key in antigen['source_keys'] if key in identities))
    groups, tables = grouped_candidate_details(data['variants'])
    by_id = {row['identifier']: row for row in directory}
    for group in groups:
        for source in group['sources']:
            entry = by_id[source['identifier']]
            if group['identifier'] not in entry['groups']:
                entry['groups'].append(group['identifier'])
            entry['rows'] += [row for row in source['rows']
                             if row[0] != 'Source identifier' and row[1] != entry['full_identifier']
                             and row not in entry['rows']]

    outcomes, by_key = [], {}
    for row in data['input_coverage']['sources']:
        content = {k: v for k, v in row.items() if k not in {'identity', 'source_key', 'source_id', 'rank'}}
        key = _key(content)
        if key not in by_key:
            by_key[key] = len(outcomes)
            outcomes.append(dict(content, members=[]))
        outcomes[by_key[key]]['members'].append(dict(identifier=row['source_id'], rank=row['rank']))
    data.update(candidate_groups=groups, evidence_tables=tables, source_directory=directory)
    data['input_coverage']['source_groups'] = outcomes
