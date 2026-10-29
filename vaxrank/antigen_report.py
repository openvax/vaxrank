"""Views of recorded antigen evidence, without inferring biological confidence."""

from .coding_context import coding_contexts, isovar_candidates


UNRECORDED = 'Not recorded'
ADMISSION_NOTE = (
    'Admission is a recorded construct-selection decision. It does not by itself '
    'establish tumor specificity or completed scientific review.')
SUPPORT_NOTE = (
    'RNA reconstruction, coding-interval and ORF-occurrence support have different '
    'subjects. These counts are not exact peptide-occurrence support, normalized '
    'expression or direct evidence of translation.')


def _display(value):
    if value is None or value == '':
        return UNRECORDED
    if isinstance(value, bool):
        return 'yes' if value else 'no'
    if isinstance(value, (list, tuple)):
        return '; '.join(str(item) for item in value) or 'none recorded'
    return value


def _section(title, rows, note=''):
    return dict(title=title, rows=[(label, _display(value)) for label, value in rows], note=note)


def _evidence_rows(value, label):
    """Keep catalog citations and their qualifiers as labelled scalar rows."""
    if isinstance(value, dict):
        return [row for key, item in value.items()
                for row in _evidence_rows(item, label + ' / ' + key.replace('_', ' '))]
    if isinstance(value, (list, tuple)) and any(isinstance(item, dict) for item in value):
        return [row for index, item in enumerate(value, 1)
                for row in _evidence_rows(item, label + ' %d' % index)]
    return [(label, value)]


def _support(support, subject):
    support = support or {}
    rows = [(subject + ' — reads', support.get('reads')),
            (subject + ' — fragments', support.get('fragments'))]
    if support.get('scope'):
        rows.append((subject + ' — recorded scope', support['scope']))
    if support.get('evidence_scope'):
        rows.append((subject + ' — sample / input', support['evidence_scope']))
    return rows


def metadata_values(antigen):
    """Retain all source metadata values, including merged/conflicting origins."""
    metadata = {}
    for name, value in antigen.source_metadata:
        metadata.setdefault(name, []).append(value)
    return metadata


def patient_background_sections(peptides, source_evidence):
    """Show exact pair provenance once, with bounded construction/phase states."""
    import json

    def interval(value):
        return '[%s, %s)' % tuple(value) if value is not None else None

    def state(value):
        return value.replace('_', ' ') if value else None

    references = {}
    for peptide in peptides:
        for epitope in peptide.epitopes:
            evidence = epitope.comparison_evidence
            if evidence is None:
                continue
            for alternative in evidence.get('alternatives', ()):
                key = (evidence.get('evidence_id'), evidence['event_id'], evidence['hypothesis_id'],
                       alternative['translation_id'], alternative['transcript_id'], alternative.get('pair_id'))
                references.setdefault(key, []).append((epitope.sequence, alternative))
    sections = []
    for index, ((identity, event_id, hypothesis_id, translation_id, transcript_id, pair_id), comparisons) in enumerate(references.items(), 1):
        source = (source_evidence or {}).get(identity, {}).get('result', {})
        event = next((e for e in source.get('events', ()) if e['event_id'] == event_id), {})
        pair = next((p for p in event.get('pairs', ()) if p['pair_id'] == pair_id), {})
        mutant, baseline = pair.get('mutant', {}), pair.get('baseline') or {}
        rows = [('Sample', source.get('sample_id')), ('RNA read set', source.get('source')),
                ('Source evidence identity', identity), ('Event identity', event_id),
                ('Hypothesis identity', hypothesis_id), ('Pair identity', pair_id),
                ('Translation identity', translation_id),
                ('Transcript identity', transcript_id),
                ('Somatic reversion policy', state(source.get('baseline_policy'))),
                ('Pair state', state(pair.get('state')) or 'Exact pair not recorded'),
                ('Pair reason', state(pair.get('reason'))),
                ('Retained hypothesis set complete', event.get('protein_hypotheses_complete')),
                ('Producer hypothesis cap', 'Unlimited (0)' if event.get('protein_sequence_limit') == 0
                 else event.get('protein_sequence_limit')),
                ('Mutant nucleotide identity', mutant.get('nucleotide_sequence_id')),
                ('Mutant coding interval (nt, 0-based half-open)', interval(mutant.get('coding_interval'))),
                ('Mutant bounded protein', mutant.get('protein_sequence')),
                ('Mutant coding-window phase', state(pair.get('phase', {}).get('state'))),
                ('Baseline construction', state(baseline.get('construction'))),
                ('Baseline observation', state(baseline.get('observation'))),
                ('Baseline nucleotide identity', baseline.get('nucleotide_sequence_id')),
                ('Baseline coding interval (nt, 0-based half-open)', interval(baseline.get('coding_interval'))),
                ('Baseline available interval (nt, 0-based half-open)', interval(baseline.get('available_coding_interval'))),
                ('Baseline bounded protein', baseline.get('protein_sequence')),
                ('Original REF-allele phase', state(baseline.get('reference_allele_phase')))]
        for label, field in (('Removed edit', 'removed_edits'), ('Retained edit', 'retained_edits')):
            edits = pair.get(field)
            if not edits:
                rows.append((label + 's', 'None recorded' if edits == [] else None))
            for number, edit in enumerate(edits or (), 1):
                prefix = '%s %d' % (label, number)
                rows += [(prefix + ' — origin', edit['origin'].replace('_', ' ')),
                         (prefix + ' — reference → alternate', '%s → %s' % (
                             edit['reference_bases'] or '∅', edit['alt_bases'] or '∅')),
                         (prefix + ' — transcript', edit['transcript_id']),
                         (prefix + ' — reference cDNA interval', interval(edit['cdna_interval'])),
                         (prefix + ' — observed cDNA interval', interval(edit['observed_interval'])),
                         (prefix + ' — source event', edit['source_event_id'])]
        rows += _support(mutant.get('rna_support'), 'Mutant RNA reconstruction')
        support = baseline.get('reference_allele_coding_support') or {}
        rows += _support(support.get('single_read_span_support'), 'Separate original REF-allele single-read spans')
        rows += _support(support.get('allele_observation_span_support'), 'Separate original REF-allele observation spans')
        seen = set()
        for sequence, comparison in comparisons:
            value = json.dumps(comparison, sort_keys=True)
            if value not in seen:
                detail = [comparison['state'], 'target aa ' + interval(comparison['mutant_protein_interval'])]
                if comparison.get('baseline_protein_interval') is not None:
                    detail.append('baseline aa ' + interval(comparison['baseline_protein_interval']))
                if comparison.get('sequence') is not None:
                    detail.append('baseline sequence ' + comparison['sequence'])
                if comparison.get('reason'):
                    detail.append(comparison['reason'].replace('_', ' '))
                rows.append(('Peptide comparison: ' + sequence, '; '.join(detail)))
                seen.add(value)
        sections.append(_section('Patient-background comparison %d' % index, rows,
            'A computed bounded counterfactual retains germline and unexplained edits. '
            'Separate REF-allele RNA spans are not evidence of patient-normal tissue or protein expression. '
            'All displayed intervals are zero-based, half-open; peptide predictions use peptide-only context. '
            'Complete source and read-set records remain in the native evidence export.'))
    return sections


def antigen_evidence_sections(antigen, source_evidence=None):
    """Render every recorded attestation and source context without merging counts.

    ``source_evidence`` is the existing native source catalog, keyed by evidence
    identity. Missing historical fields remain unknown. Multiple metadata values
    are retained separately; an ambiguous origin never supplies another origin's
    reference model or sample evidence.
    """
    sections = []
    attestation = antigen.tumor_specificity
    for index, record in enumerate(attestation.evidence_records, 1):
        rows = [('Evidence type', record.evidence_kind), ('Source', record.evidence_source),
                ('Source version', record.evidence_version), ('Subject', record.subject_id),
                ('Sample', record.sample_id), ('Declared patient-specific', record.patient_specific),
                ('Recorded decision', record.passed)]
        if record.numeric_value is not None or record.threshold is not None:
            rows += [('Value', record.numeric_value), ('Units', record.unit),
                     ('Threshold', record.threshold)]
        rows += [(name.replace('_', ' '), value) for name, value in record.details]
        sections.append(_section('Admission evidence %d' % index, rows,
            'Scope and decisions are supplied by the evidence source; they are not independent validation.'))

    metadata = metadata_values(antigen)
    identities = metadata.get('isovar_evidence_id', [])
    sources = [(identity, (source_evidence or {}).get(identity, {})) for identity in identities]
    candidates = isovar_candidates(antigen)
    for index, context in enumerate(coding_contexts(antigen, source_evidence), 1):
        sections.append(_section('Coding reference context %d' % index,
            [(label, context.get(name)) for name, label in (
                ('kind', 'Reference kind'), ('orf_type', 'ORF type'),
                ('orf_id', 'ORF identity'), ('source', 'Reference catalog'),
                ('reference_key', 'Reference identity'),
                ('source_transcript_id', 'Source transcript'))],
            'The coding reference is separate from the source event and sample evidence.'))
    for index, candidate in enumerate(candidates, 1):
        # Older exports already record this fact as sequence_extent. Reading
        # that explicit state does not re-run admission or infer missing RNA.
        partial = candidate.get('partial_coding_candidate', {
            'complete_coding_candidate': False,
            'partial_coding_candidate': True,
        }.get(candidate.get('sequence_extent')))
        rows = [('Coding extent', candidate.get('sequence_extent')),
                ('Partial coding candidate', partial),
                ('Ambiguous coding interpretation', candidate.get('ambiguous_coding_candidate'))]
        if 'occurrence_rna_support' in candidate:
            rows += _support(candidate.get('rna_support'), 'ORF candidate union support')
            rows += _support(candidate.get('occurrence_rna_support'), 'ORF occurrence support')
        elif 'rna_support' in candidate:
            rows += _support(candidate['rna_support'], 'RNA reconstruction support')
        else:
            rows += [('Candidate RNA support', None)]
        coding = candidate.get('coding_support')
        if coding is not None:
            rows += [('Coding interval fully tiled', coding.get('fully_tiled')),
                     ('Minimum coding-interval fragment depth', coding.get('minimum_fragment_depth')),
                     ('Coding-support scope', coding.get('scope'))]
            rows += _support(coding.get('single_read_span_support'), 'Single-read coding-interval span support')
        for name, label in (
                ('rna_candidate_id', 'RNA candidate identity'),
                ('model_id', 'Reference model identity'),
                ('reference_orf_source_gene_id', 'Reference ORF source gene'),
                ('sequence_id', 'RNA sequence identity'),
                ('origin', 'RNA sequence origin'),
                ('all_bases_observed', 'All coding bases observed'),
                ('sequence_status', 'RNA sequence status'),
                ('translation_table', 'Translation table'),
                ('support_scope', 'Admission support scope'),
                ('uncertainty_flags', 'Recorded uncertainties')):
            if name in candidate:
                rows.append((label, candidate[name]))
        sections.append(_section('Candidate RNA and coding evidence %d' % index, rows, SUPPORT_NOTE))
        for reference_index, reference in enumerate(candidate.get('reference_orfs', ()), 1):
            context = reference.get('reference_context', {})
            rows = [(label, reference.get(name)) for name, label in (
                ('gene_name', 'Gene'), ('protein_id', 'Protein / ORF'),
                ('reference_key', 'Reference identity'), ('source', 'Reference catalog'),
                ('protein_matches_reference', 'Protein matches reference'),
                ('start_status', 'Reference start coverage'), ('stop_status', 'Reference stop coverage'))]
            rows += [(label, context.get(name)) for name, label in (
                ('species', 'Source species'), ('strain', 'Strain'),
                ('assembly_accession', 'Assembly accession'), ('genome_build', 'Reference assembly'),
                ('annotation_name', 'Annotation'), ('annotation_version', 'Annotation version'))]
            sections.append(_section('Coding reference %d.%d' % (index, reference_index), rows,
                'Reference attribution is separate from RNA evidence in this sample.'))

    for index, (identity, evidence) in enumerate(sources, 1):
        result = evidence.get('isovar_result', {})
        controls = evidence.get('options', {}).get('admission', {})
        rows = [('Evidence identity', identity), ('Producer schema', result.get('schema')),
                ('Producer family', evidence.get('family')),
                ('Sample', result.get('sample_id') or result.get('provenance', {}).get('sample_id')),
                ('Input source', result.get('source') or result.get('provenance', {}).get('source')),
                ('Allow partial coding candidates', controls.get('allow_partial')),
                ('Allow ambiguous coding candidates', controls.get('allow_ambiguous')),
                ('Minimum admission fragments', controls.get('min_fragments'))]
        unambiguous = len(candidates) == 1 and len(sources) == 1
        note = '' if evidence else 'Linked source record unavailable; only its identity is retained.'
        if not unambiguous:
            note += (' ' if note else '') + (
            'Candidate-to-source pairing is not uniquely recorded. Contexts are retained '
            'separately; reference-model and RNA-path evidence are not assigned by guesswork.')
        sections.append(_section('Recorded Isovar source %d' % index, rows, note))

        # Only one unambiguous candidate/context pair can identify its reference
        # model. In merged historical metadata, retain contexts without guessing.
        if not unambiguous:
            continue
        candidate = candidates[0]
        for model in result.get('models', ()):
            if not candidate.get('model_id') or model.get('model_id') != candidate['model_id']:
                continue
            reference = model['reference']
            rows = [(label, reference.get(name)) for name, label in (
                ('orf_id', 'Reference ORF'), ('orf_type', 'Reference ORF type'),
                ('source', 'Reference catalog'), ('reference_key', 'Reference identity'),
                ('source_transcript_id', 'Reference transcript'), ('source_gene_id', 'Reference gene'))]
            for name, value in (reference.get('evidence') or {}).items():
                rows += _evidence_rows(value, 'Catalog ' + name.replace('_', ' '))
            sections.append(_section('Reference ORF evidence', rows,
                'Reference/catalog evidence is not evidence observed in this sample.'))
        for path in result.get('paths', ()):
            single_fusion_path = evidence.get('family') == 'fusion' and len(result['paths']) == 1
            if not single_fusion_path and (not candidate.get('path_id')
                    or path.get('path_id') != candidate['path_id']):
                continue
            rows = [('RNA path', path.get('path_id')), ('Path frame status', path.get('frame_status'))]
            for partner, transcripts in path.get('compatible_transcripts', {}).items():
                rows.append((partner.title() + ' compatible transcripts', transcripts))
            for junction_index, junction in enumerate(path.get('junctions', ()), 1):
                rows.append(('Junction %d annotated' % junction_index, junction.get('annotated')))
                for name in ('relation', 'query_interval', 'left', 'right'):
                    if name in junction:
                        rows.append(('Junction %d %s' % (junction_index, name.replace('_', ' ')), junction[name]))
                rows += _support(junction.get('direct_support'), 'Junction %d direct support' % junction_index)
            sections.append(_section('Recorded RNA path evidence', rows, SUPPORT_NOTE))
    return sections
