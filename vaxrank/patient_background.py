"""Exact peptide comparisons against Isovar's bounded patient backgrounds.

Isovar constructs the counterfactual; Topiary validates its exchange format.
This consumer only joins exact source identities and projects an already chosen
peptide interval through the recorded cDNA edits. It never aligns proteins.
"""

from copy import deepcopy
from hashlib import sha256
import json

from topiary import read_isovar_hypotheses, read_isovar_patient_background_pairs


def _project_boundary(position, edits, *, right):
    shift = 0
    for edit in edits:
        begin, end = edit['observed_interval']
        delta = len(edit['reference_bases']) - (end - begin)
        if begin < position < end and delta:
            return None
        if end < position or (end == position and (begin < end or right)):
            shift += delta
    return position + shift


def _project_pair(pair, start, end):
    """Project half-open amino-acid bounds through a validated exact pair."""
    result = dict(pair_id=pair['pair_id'], translation_id=pair['translation_id'],
                  transcript_id=pair['transcript_id'], state='unavailable',
                  sequence=None, reason=pair['reason'],
                  mutant_protein_interval=[start, end])
    baseline = pair['baseline']
    if baseline is None:
        return result
    if not baseline['protein_sequence']:
        return dict(result, reason='zero_length_baseline')
    origin = pair['mutant']['coding_interval'][0]
    bounds = [origin + 3 * start, origin + 3 * end]
    result['mutant_nucleotide_interval'] = bounds
    projected = [_project_boundary(p, pair['removed_edits'], right=right)
                 for p, right in zip(bounds, (False, True))]
    if None in projected:
        return dict(result, reason='reverted_edit_crosses_peptide_boundary')
    result['baseline_nucleotide_interval'] = projected
    origin = baseline['coding_interval'][0]
    relative = [p - origin for p in projected]
    if any(p % 3 for p in relative):
        return dict(result, reason='peptide_boundaries_do_not_share_baseline_frame')
    low, high = (p // 3 for p in relative)
    result['baseline_protein_interval'] = [low, high]
    if not 0 <= low <= high <= len(baseline['protein_sequence']):
        return dict(result, reason='peptide_exceeds_bounded_baseline_translation')
    sequence = baseline['protein_sequence'][low:high]
    if not sequence:
        return dict(result, reason='zero_length_peptide_comparator')
    internal_indel = any(
        len(edit['reference_bases']) != edit['observed_interval'][1] - edit['observed_interval'][0]
        and bounds[0] <= edit['observed_interval'][0] <= bounds[1]
        for edit in pair['removed_edits'])
    return dict(result, state='available', sequence=sequence, reason=None,
                same_length=len(sequence) == end - start,
                position_aligned=len(sequence) == end - start and not internal_indel)


class PatientBackgroundPairs:
    """Validated source evidence, shared by all peptide comparisons in a run.

    Parameters
    ----------
    data : mapping or path
        Isovar's complete patient-background pair export, or its JSON/JSON.gz
        path. The full source is retained once; comparison records reference
        it by ``evidence_id`` rather than copying its read sets per peptide.
    """

    def __init__(self, data):
        imported = read_isovar_patient_background_pairs(data)
        self.export = imported.extra['isovar_patient_background_pairs']
        self.evidence_id = 'patient_background:' + sha256(json.dumps(
            self.export, sort_keys=True, separators=(',', ':'), allow_nan=False
        ).encode()).hexdigest()
        # Topiary validates repeated pair identities. Its source export keeps
        # original repetitions, while consumers use each exact context once.
        self.events = {event['event_id']: dict(event, pairs=list({p['pair_id']: p
            for p in event['pairs']}.values())) for event in self.export['events']}

    @property
    def source_evidence(self):
        """A record for the existing native source-evidence store."""
        return {self.evidence_id: dict(kind='patient_background', result=deepcopy(self.export))}

    def match_hypothesis(self, hypotheses, event_id, hypothesis_id):
        """Join every nucleotide/transcript alternative to its exact pair.

        ``hypotheses`` is the independent Isovar protein-hypothesis export
        that supplied the candidate. Missing pairs remain unavailable; a
        matching ID with inconsistent source evidence is rejected.
        """
        imported = read_isovar_hypotheses(hypotheses)
        source = imported.extra['isovar_hypotheses']
        return self._match_hypothesis(source, event_id, hypothesis_id)

    def contexts_for_results(self, results):
        """Bind each live result to the exact protein Vaxrank selects.

        The public export retains the producer's ranked object order. Rank 1
        is ``IsovarResult.top_protein_sequence``; no sequence search or private
        Isovar identity function is used.
        """
        from isovar import export_protein_hypotheses

        source = export_protein_hypotheses(
            results, sample_id=self.export['sample_id'], source=self.export['source'])
        source = read_isovar_hypotheses(source).extra['isovar_hypotheses']
        contexts = {}
        for result, event in zip(results, source['events'], strict=True):
            proteins = event['protein_hypotheses']
            if proteins:
                protein = next(p for p in proteins if p['isovar_rank'] == 1)
                contexts[id(result)] = self._match_hypothesis(
                    source, event['event_id'], protein['hypothesis_id'])
        return contexts

    def _match_hypothesis(self, source, event_id, hypothesis_id):
        if any(source[key] != self.export[key] for key in ('sample_id', 'source', 'evidence_scope')):
            raise ValueError('Patient-background and candidate sample/read-set provenance disagree')
        events = [e for e in source['events'] if e['event_id'] == event_id]
        if len(events) != 1:
            raise ValueError('Selected patient-background event is absent or ambiguous')
        event = events[0]
        proteins = [p for p in event['protein_hypotheses'] if p['hypothesis_id'] == hypothesis_id]
        if len(proteins) != 1:
            raise ValueError('Selected patient-background hypothesis is absent or ambiguous')
        protein = proteins[0]
        paired_event = self.events.get(event_id)
        if paired_event is not None and paired_event['variant'] != event['variant']:
            raise ValueError('Patient-background and candidate event provenance disagree')
        pairs = [p for p in (paired_event or {}).get('pairs', []) if p['hypothesis_id'] == hypothesis_id]
        expected = {}
        for translation in protein['translations']:
            for transcript in translation['reference_context']['transcript_ids']:
                key = (translation['translation_id'], transcript)
                if key in expected and expected[key] != translation:
                    raise ValueError('Conflicting candidate translation provenance')
                expected[key] = translation
        matched = {}
        for pair in pairs:
            key = (pair['translation_id'], pair['transcript_id'])
            if key not in expected or key in matched:
                raise ValueError('Patient-background translation/transcript provenance disagrees')
            translation = expected[key]
            mutant = pair['mutant']
            same = (
                mutant['protein_sequence'] == protein['amino_acids']
                and mutant['nucleotide_sequence'] == translation['nucleotide_sequence']
                and mutant['nucleotide_sequence_id'] == translation['nucleotide_sequence_id']
                and mutant['coding_interval'] == translation['translated_interval']
                and pair['translation_table'] == translation['reference_context']['translation_table']
                and pair['strand'] == translation['reference_context']['strand']
                and all(mutant['rna_support'][k] == translation['rna_support'][k]
                        for k in ('reads', 'fragments')))
            left = mutant['rna_support']['evidence_set_id']
            right = translation['rna_support']['evidence_set_id']
            if not same or left != right or (
                    left is not None and self.export['evidence_sets'][left] != source['evidence_sets'][right]):
                raise ValueError('Patient-background nucleotide/read evidence disagrees with candidate')
            matched[key] = pair
        return MatchedPatientBackground(
            self.evidence_id, self.export['baseline_policy'], event_id,
            hypothesis_id, protein['amino_acids'], expected, matched)


class MatchedPatientBackground:
    """An exact hypothesis join; retained alternatives are never ranked here."""

    def __init__(self, evidence_id, baseline, event_id, hypothesis_id, sequence, expected, matched):
        self.evidence_id, self.baseline = evidence_id, baseline
        self.event_id, self.hypothesis_id = event_id, hypothesis_id
        self.sequence, self.expected, self.matched = sequence, expected, matched

    def compare(self, start, end):
        """Reconcile a peptide interval across every compatible exact pair.

        Different sequences and mixtures of available/unavailable contexts are
        ambiguous. Equality is retained as evidence, including when candidate
        and comparator peptides are identical.
        """
        if (type(start) is not int or type(end) is not int
                or not 0 <= start < end <= len(self.sequence)):
            raise ValueError('Patient-background peptide interval must be nonempty and in bounds')
        alternatives = []
        for key in sorted(self.expected):
            pair = self.matched.get(key)
            alternatives.append(_project_pair(pair, start, end) if pair is not None else dict(
                pair_id=None, translation_id=key[0], transcript_id=key[1],
                state='unavailable', sequence=None, reason='exact_pair_not_recorded',
                mutant_protein_interval=[start, end]))
        sequences = {a['sequence'] for a in alternatives if a['state'] == 'available'}
        all_available = bool(alternatives) and all(a['state'] == 'available' for a in alternatives)
        available = all_available and len(sequences) == 1
        state = 'available' if available else 'ambiguous' if sequences else 'unavailable'
        sequence = next(iter(sequences)) if available else None
        comparison = ('identical' if sequence == self.sequence[start:end] else 'amino_acid_divergent') if available else state
        return dict(kind='patient_background_comparison', evidence_id=self.evidence_id,
                    event_id=self.event_id, hypothesis_id=self.hypothesis_id,
                    baseline_policy=self.baseline, state=state, sequence=sequence,
                    comparison=comparison, same_length=None if sequence is None else len(sequence) == end - start,
                    position_aligned=available and all(a['position_aligned'] for a in alternatives),
                    mutant_protein_interval=[start, end], alternatives=alternatives)


def attach_patient_background(epitopes, predictor, context, *, baseline='focal'):
    """Attach exact comparator measurements through the existing Topiary run.

    Candidate coordinates must still refer to the selected, unsliced hypothesis.
    Per-occurrence evidence retains unavailable and identical comparisons. The
    existing native evidence store retains the full source referenced by IDs.
    """
    from dataclasses import replace
    from .allele_validation import validate_prediction_frame
    from .candidate_epitope import Peptide
    from .epitope_logic import _prediction_from_row
    from .prediction_input import finite_prediction_value

    records, inputs = {}, {}
    lengths = {length for model in predictor.models
               for length in (getattr(model, 'default_peptide_lengths', None) or ())}
    for index, epitope in enumerate(epitopes):
        if not epitope.overlaps_mutation:
            continue
        if context is None:
            record = dict(kind='patient_background_comparison', state='unavailable',
                          sequence=None, comparison='unavailable', alternatives=[],
                          baseline_policy=baseline,
                          reason='exact_hypothesis_context_not_recorded')
        else:
            if epitope.source_sequence != context.sequence:
                raise ValueError('Patient-background hypothesis differs from candidate source window')
            record = context.compare(epitope.offset, epitope.offset + len(epitope.sequence))
        record['prediction_context'] = 'peptide_only'
        record['prediction_state'] = 'comparison_unavailable'
        if record['state'] == 'available':
            if lengths and len(record['sequence']) not in lengths:
                record['prediction_state'] = 'outside_requested_peptide_lengths'
            else:
                record['prediction_state'] = 'no_measurements_recorded'
                inputs[str(index)] = record['sequence']
        records[index] = record
    predictions = {}
    if inputs:
        frame = predictor.predict_from_named_peptides(inputs)
        validate_prediction_frame(frame, 'Patient-background comparator predictions')
        for _, row in frame.iterrows():
            key = str(row['source_sequence_name'])
            if key not in inputs or row['peptide'] != inputs[key]:
                raise ValueError('Patient-background predictor returned an unexpected peptide')
            predictions.setdefault(int(key), []).append(_prediction_from_row(
                row, peptide=inputs[key], score=finite_prediction_value(row.get('score')),
                percentile_rank=finite_prediction_value(row.get('percentile_rank'))))
    output = []
    for index, epitope in enumerate(epitopes):
        if index in records:
            record = records[index]
            comparators = dict(epitope.comparators)
            comparators.pop('patient_background', None)
            if record['state'] == 'available':
                leaves = tuple(predictions.get(index, ()))
                if leaves:
                    record['prediction_state'] = 'measurements_recorded'
                comparators['patient_background'] = Peptide(sequence=record['sequence'], predictions=leaves)
            evidence = tuple(e for e in epitope.input_evidence if e.get('kind') != 'patient_background_comparison')
            epitope = replace(epitope, comparators=comparators, input_evidence=evidence + (record,))
        output.append(epitope)
    return output


def load_patient_backgrounds(inputs, config):
    """Resolve pair inputs once, rejecting conflicting or unused controls."""
    sources = [value if isinstance(value, PatientBackgroundPairs) else PatientBackgroundPairs(value)
               for value in inputs or ()]
    if sources and (config is None or config.comparison_reference != 'patient_background'):
        raise ValueError('Patient-background pair inputs require --comparison-reference patient_background')
    scopes = set()
    for source in sources:
        if source.export['baseline_policy'] != config.patient_background_baseline:
            raise ValueError('Patient-background evidence and configured baseline disagree')
        scope = tuple(source.export['evidence_scope'])
        if scope in scopes:
            raise ValueError('Provide one complete patient-background export per sample/read set')
        scopes.add(scope)
    return sources


def validate_native_comparisons(dataset):
    """Check saved comparison links and coordinates without inference or repair."""
    sources = {}
    store = dataset.selection.get('source_evidence', {})
    for epitope in dataset.epitopes:
        record = epitope.comparison_evidence
        if record is None:
            continue
        if dataset.config is not None and (
                dataset.config.comparison_reference != 'patient_background'
                or dataset.config.patient_background_baseline != record.get('baseline_policy')):
            raise ValueError('Native patient-background comparison disagrees with its saved policy')
        identity = record.get('evidence_id')
        comparator = epitope.comparators.get('patient_background')
        if (comparator is not None) != (record.get('state') == 'available'):
            raise ValueError('Native patient-background comparator disagrees with its availability')
        if identity is None:
            if record.get('state') != 'unavailable' or comparator is not None:
                raise ValueError('Patient-background comparison lacks its exact source context')
            continue
        if identity not in sources:
            raw = store.get(identity, {}).get('result')
            if raw is None:
                raise ValueError('Native patient-background comparison is missing its source evidence')
            source = PatientBackgroundPairs(raw)
            if source.evidence_id != identity:
                raise ValueError('Native patient-background source evidence identity disagrees')
            sources[identity] = source
        source = sources[identity]
        pairs = [p for p in source.events.get(record['event_id'], {}).get('pairs', ())
                 if p['hypothesis_id'] == record['hypothesis_id']]
        if not pairs:
            if record['state'] != 'unavailable' or comparator is not None:
                raise ValueError('Native patient-background hypothesis is unavailable')
            continue
        sequences = {p['mutant']['protein_sequence'] for p in pairs}
        if len(sequences) != 1:
            raise ValueError('Native patient-background hypothesis sequences disagree')
        sequence = next(iter(sequences))
        matched = {(p['translation_id'], p['transcript_id']): p for p in pairs}
        if len(matched) != len(pairs):
            raise ValueError('Native patient-background translation/transcript contexts conflict')
        expected = dict(matched)
        for alternative in record['alternatives']:
            if alternative['pair_id'] is None:
                expected[(alternative['translation_id'], alternative['transcript_id'])] = None
        context = MatchedPatientBackground(identity, source.export['baseline_policy'],
            record['event_id'], record['hypothesis_id'], sequence, expected, matched)
        start, end = record['mutant_protein_interval']
        comparison = context.compare(start, end)
        if epitope.sequence != sequence[start:end] or any(record.get(k) != v for k, v in comparison.items()):
            raise ValueError('Native patient-background peptide comparison disagrees with its exact pair')
        if comparator is not None and (comparator.sequence != comparison['sequence']
                or any(p.peptide != comparator.sequence for p in comparator.predictions_flat())):
            raise ValueError('Native patient-background comparator measurements disagree with its sequence')


def attach_antigen_backgrounds(dataset, predictor, sources, source_evidence):
    """Use exact protein-hypothesis source records on the generic antigen path."""
    from dataclasses import replace

    contexts, groups, candidates = {}, {}, list(dataset.epitopes)
    for index, epitope in enumerate(dataset.epitopes):
        antigen = dataset.antigens[epitope.prediction_id]
        if antigen.kind != 'mutation':
            continue
        key = antigen.source_identifier
        metadata = dict(antigen.source_metadata)
        evidence_id = metadata.get('isovar_evidence_id')
        cache_key = (evidence_id, key)
        if cache_key not in contexts:
            context = None
            source = (source_evidence or {}).get(evidence_id, {}).get('isovar_result', {})
            if source.get('schema') in ('isovar.protein_hypotheses.v1', 'isovar.protein_hypotheses.v2'):
                compatible = [s for s in sources if s.export['evidence_scope'] == source['evidence_scope']]
                if compatible:
                    event_id = json.loads(metadata['isovar_candidate'])['event_id']
                    context = compatible[0].match_hypothesis(source, event_id, key)
            contexts[cache_key] = context
        candidate = replace(epitope, overlaps_mutation=epitope.is_targetable)
        groups.setdefault(cache_key, []).append((index, candidate))
    for key, group in groups.items():
        compared = attach_patient_background(
            [e for _, e in group], predictor, contexts[key],
            baseline=dataset.config.patient_background_baseline)
        for (index, _), epitope in zip(group, compared, strict=True):
            candidates[index] = epitope
    dataset.epitopes = tuple(candidates)
    store = dict(dataset.selection.get('source_evidence', {}))
    for source in sources:
        store.update(source.source_evidence)
    dataset.selection['source_evidence'] = store


def direct_rna_scope(args, background=None):
    """Reconcile declared RNA identities for comparison and native provenance.

    Explicit CLI, direct-manifest and paired-export identities must agree.
    Unknown identities stay unknown; filenames never become sample identities.
    """
    from .input_scope import InputScope, combine_scopes, read_input_manifest, scope_from_mapping
    declarations = {}
    if getattr(args, 'input_manifest', None):
        _, declarations = read_input_manifest(args.input_manifest, include_direct=True)
    scope = scope_from_mapping({name: declarations.get(name) for name in ('sample_id', 'library_id')},
                               'Direct RNA manifest')
    cli = scope_from_mapping(dict(sample_id=getattr(args, 'rna_sample_id', None),
                                  library_id=getattr(args, 'rna_source', None)), 'Direct RNA CLI')
    scope = combine_scopes(scope, cli, 'Direct RNA')
    if background is not None:
        scope = combine_scopes(scope, InputScope(sample_id=background.export['sample_id'],
            library_id=background.export['source']), 'Patient-background RNA evidence')
    return scope


def direct_patient_background(args, config, results=None):
    """Load declared pairs or export from already reconstructed direct RNA."""
    scope = direct_rna_scope(args)
    if config.comparison_reference != 'patient_background':
        if getattr(args, 'patient_background_pairs', None):
            raise ValueError('Patient-background pair inputs require --comparison-reference patient_background')
        return None
    sample, source = scope.sample_id, scope.library_id
    supplied = load_patient_backgrounds(getattr(args, 'patient_background_pairs', None), config)
    if supplied:
        matching = [s for s in supplied if (sample is None or s.export['sample_id'] == sample)
                    and (source is None or s.export['source'] == source)]
        if len(matching) != 1:
            raise ValueError('Direct RNA requires one matching patient-background sample/read set')
        direct_rna_scope(args, matching[0])
        return matching[0]
    if not sample or not source:
        raise ValueError('Live patient-background comparisons require --rna-sample-id and --rna-source, '
                         'direct sample_id/library_id in a manifest, or --patient-background-pairs')
    if results is None:
        return None  # Preflight passed; export only after the already planned RNA reconstruction.
    from isovar import export_patient_background_pairs
    return PatientBackgroundPairs(export_patient_background_pairs(
        results, sample_id=sample, source=source, baseline=config.patient_background_baseline))
