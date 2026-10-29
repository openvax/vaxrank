"""Recorded coding reference, separate from the source's event kind."""

import json


REFERENCE_FIELDS = ('kind', 'orf_type', 'orf_id', 'source', 'reference_key',
                    'source_transcript_id')


def reference_context(reference):
    """Retain the producer's reference fields without inferring missing values."""
    return {name: reference.get(name) for name in REFERENCE_FIELDS}


def isovar_candidates(antigen):
    candidates = []
    for name, value in antigen.source_metadata:
        if name != 'isovar_candidate':
            continue
        try:
            candidate = json.loads(value)
        except (ValueError, TypeError) as error:
            raise ValueError('Invalid isovar_candidate report metadata: expected a JSON object') from error
        if not isinstance(candidate, dict):
            raise ValueError('Invalid isovar_candidate report metadata: expected a JSON object')
        candidates.append(candidate)
    return candidates


def coding_contexts(antigen, source_evidence=None):
    """Return all recorded contexts; resolve legacy models only when unambiguous."""
    candidates = isovar_candidates(antigen)
    identities = [value for name, value in antigen.source_metadata if name == 'isovar_evidence_id']
    contexts = []
    for candidate in candidates:
        reference = candidate.get('coding_reference')
        if reference is None and len(candidates) == len(identities) == 1:
            result = (source_evidence or {}).get(identities[0], {}).get('isovar_result', {})
            models = [m for m in result.get('models', ())
                      if candidate.get('model_id') and m.get('model_id') == candidate['model_id']]
            if len(models) == 1:
                reference = models[0].get('reference')
        if reference is not None:
            if not isinstance(reference, dict):
                raise ValueError('Invalid Isovar coding_reference: expected a JSON object')
            context = reference_context(reference)
            if context not in contexts:
                contexts.append(context)
    return contexts


def coding_context_label(contexts):
    labels = []
    for context in contexts:
        kind, orf_type = context.get('kind'), context.get('orf_type')
        if kind == 'annotated':
            label = 'annotated CDS'
        elif kind == 'extra_orf':
            label = 'noncanonical ' + (orf_type or 'ORF')
        elif kind:
            label = str(kind) + (' / ' + str(orf_type) if orf_type else '')
        else:
            label = 'coding context not recorded'
        if label not in labels:
            labels.append(label)
    return '; '.join(labels) or 'coding context not recorded'


def antigen_source_label(antigen, source_evidence=None):
    contexts = coding_contexts(antigen, source_evidence)
    if contexts or ('isovar_family', 'variant_orfs') in antigen.source_metadata:
        return antigen.kind + ' / ' + coding_context_label(contexts)
    return antigen.kind
