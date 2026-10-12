"""Bridge scored Vaxrank antigen windows to explicit mhctools kinetics."""

from dataclasses import replace
import hashlib

from mhcgnomes import parse
from mhctools import PeptideInput, assess_vaccine_windows

from .core_logic import SOURCE_AGNOSTIC_RANKING_RULES
from .epitope_logic import slice_epitopes
from .native_serialization import to_native_json
from .vaccine_peptide import VaccinePeptide


def enumerate_exposure_windows(dataset, min_length, max_length):
    """Enumerate admitted target-containing windows before ranking/truncation.

    Saved per-allele scores remain authoritative. Explicit antigen context is
    required; table-only and held-out occurrences are returned as coverage gaps.
    Each yielded tuple is (stable id, VaccinePeptide, source start, provenance).
    """
    if (type(min_length) is not int or type(max_length) is not int or
            min_length < 1 or max_length < min_length):
        raise ValueError('Window lengths must be positive integers with min <= max')
    selected = dataset.select_representatives()
    generalized_ids = (set(dataset.scoring_frame().loc[
        dataset.result.long_df.source_observation_id.notna(), 'prediction_id'].dropna())
        if 'source_observation_id' in dataset.result.df else set())
    groups, gaps = {}, []
    for epitope in dataset.epitopes:
        key = epitope.prediction_group_source
        if selected is not None and key in generalized_ids:
            epitope = replace(epitope, per_allele_scores={
                allele: score for allele, score in epitope.per_allele_scores.items()
                if (key, allele) in selected})
        if not epitope.per_allele_scores:
            gaps.append({'prediction_id': key, 'sequence': epitope.sequence, 'reason': 'no_saved_per_allele_scores'})
            continue
        antigen = dataset.antigens.get(key)
        if antigen is None or not antigen.tumor_specificity.admits_construct:
            gaps.append({'prediction_id': key, 'sequence': epitope.sequence,
                         'reason': 'missing_antigen_context' if antigen is None else 'antigen_not_admitted'})
            continue
        if antigen.amino_acids[epitope.offset:epitope.offset + len(epitope.sequence)] != epitope.sequence:
            raise ValueError('Scored target does not match its explicit antigen occurrence: ' + key)
        provenance = to_native_json(epitope.input_provenance) if epitope.input_provenance else ''
        identity = to_native_json(antigen) + '\n' + provenance
        groups.setdefault(identity, (antigen, [], provenance))[1].append(
            replace(epitope, prediction_id=key) if not epitope.prediction_id else epitope)
    windows = []
    for identity, (antigen, epitopes, provenance) in groups.items():
        source_id = hashlib.sha256(identity.encode()).hexdigest()[:20]
        source = {name: getattr(antigen, name) for name in (
            'source_identifier', 'gene_name', 'gene_id', 'transcript_ids',
            'protein_ids', 'species', 'source_metadata')}
        source.update(input_provenance=provenance, source_id=source_id,
                      antigen_kind=antigen.kind, antigen=to_native_json(antigen),
                      antigen_sequence=antigen.amino_acids)
        for length in range(min_length, min(max_length, len(antigen.amino_acids)) + 1):
            for start in range(len(antigen.amino_acids) - length + 1):
                if not antigen.interval_is_targetable(start, start + length):
                    continue
                candidates = slice_epitopes(epitopes, start, start + length)
                vaccine = VaccinePeptide(
                    antigen=antigen.sliced(start, start + length), epitopes=candidates,
                    num_target_epitopes_to_keep=None, combined_score_expr='target_epitope_score',
                    ranking_rules=SOURCE_AGNOSTIC_RANKING_RULES)
                if not vaccine.contains_target_epitopes():
                    continue
                identifier = '%s:%s-%s' % (source_id, start, start + length)
                windows.append((identifier, vaccine, start, source))
    return windows, gaps


def window_exposure_request(windows, profile, *, n_term='free', c_term='free', score_source=None):
    """Convert existing candidate windows without changing scores or selection.

    ``profile`` supplies scenario, times_hours, optional common class kinetics,
    dose and per-candidate overrides keyed by the stable window identifier.
    Overrides may supply kinetics and dose only: they cannot replace core
    coordinates or biological provenance. Unused override IDs are rejected.
    """
    unknown = set(profile) - {'scenario', 'times_hours', 'administration', 'organism',
                              'kinetics', 'dose', 'candidates', 'core_intervals'}
    if unknown:
        raise ValueError('Unknown exposure profile fields: ' + ', '.join(sorted(unknown)))
    request = {key: profile[key] for key in ('scenario', 'times_hours')}
    request.update(schema_version=1, administration=profile.get('administration', 'IM'),
                   organism=profile.get('organism'), kinetics=profile.get('kinetics', {}), candidates=[])
    overrides = profile.get('candidates', {})
    ids = set()
    core_intervals = {}
    for item in profile.get('core_intervals', []):
        key = item['prediction_id'], parse(item['allele']).to_string()
        if key in core_intervals:
            raise ValueError('Duplicate supplied binding-core interval')
        if not isinstance(item.get('source'), str) or not item['source'].strip():
            raise ValueError('Supplied core intervals require a source')
        core_intervals[key] = item
    used_intervals = set()
    for identifier, vaccine, start, source in windows:
        ids.add(identifier)
        peptide = PeptideInput(
            vaccine.amino_acids, n_term=n_term, c_term=c_term, occurrence_id=identifier,
            source_start=start, source_sequence_name=source['source_identifier'] or None,
            source_gene=source['gene_name'] or None, source_species=source['species'] or None)
        cores = []
        for epitope in vaccine.target_epitopes:
            for allele, score in epitope.per_allele_scores.items():
                identity = parse(allele)
                interval_key = epitope.prediction_group_source, identity.to_string()
                interval = core_intervals.get(interval_key)
                core_start, core_end = 0, len(epitope.sequence)
                if interval is not None:
                    core_start, core_end = interval['start'], interval['end']
                    if (type(core_start) is not int or type(core_end) is not int or
                            not 0 <= core_start < core_end <= len(epitope.sequence)):
                        raise ValueError('Supplied binding core must lie within its scored ligand')
                    used_intervals.add(interval_key)
                cores.append({'label': '%s@%s' % (epitope.sequence, epitope.offset),
                              'sequence': epitope.sequence[core_start:core_end], 'start': epitope.offset + core_start,
                              'end': epitope.offset + core_end,
                              'scored_ligand': epitope.sequence,
                              'scored_ligand_start': epitope.offset,
                              'scored_ligand_end': epitope.offset + len(epitope.sequence),
                              'core_source': interval['source'] if interval else None,
                              'allele': identity.to_string(), 'mhc_class': 'I' if identity.is_class1 else 'II',
                              'score_contribution': float(score), 'score_source': score_source,
                              'prediction_id': epitope.prediction_group_source,
                              'target_definition': ('supplied binding-core interval within scored ligand' if interval else
                                                    'scored ligand interval; no shorter binding core inferred')})
        # Several scored ligands can document the same physical binding core.
        # Combine their already-computed score contributions, retaining each
        # component. Repeated evidence for the same original ligand is still
        # an ambiguity for the upstream representative policy to resolve.
        merged = {}
        for core in cores:
            key = core['start'], core['end'], core['allele']
            component = {name: core[name] for name in (
                'prediction_id', 'scored_ligand', 'scored_ligand_start', 'scored_ligand_end',
                'score_contribution', 'score_source', 'core_source')}
            if key not in merged:
                merged[key] = {**core, 'score_components': [component]}
            else:
                existing = merged[key]
                if core['core_source'] is None or existing['core_source'] is None or any(
                        (c['scored_ligand'], c['scored_ligand_start']) ==
                        (core['scored_ligand'], core['scored_ligand_start'])
                        for c in existing['score_components']):
                    raise ValueError('Repeated scored ligand occurrence/allele; resolve duplicate evidence before exposure assessment')
                existing['score_components'].append(component)
                existing['score_contribution'] += core['score_contribution']
                existing['prediction_id'] = None
                existing['scored_ligand'] = None
                existing['scored_ligand_start'] = None
                existing['scored_ligand_end'] = None
                existing['core_source'] = 'See each documented score component'
                existing['label'] = '%s@%s' % (core['sequence'], core['start'])
        cores = list(merged.values())
        override = overrides.get(identifier, {})
        unknown = set(override) - {'kinetics', 'dose'}
        if unknown:
            raise ValueError('Unknown candidate profile fields: ' + ', '.join(sorted(unknown)))
        candidate = {'id': identifier, 'peptide': peptide.to_dict(), 'delivery': 'slp', 'cores': cores,
                     'metadata': {**source, 'source_start': start,
                                  'target_epitope_score': vaccine.target_epitope_score}}
        if 'dose' in profile:
            candidate['dose'] = profile['dose']
        candidate.update(override)
        request['candidates'].append(candidate)
    if set(core_intervals) - used_intervals:
        raise ValueError('Supplied binding-core intervals do not match any enumerated scored occurrence/allele')
    unknown = set(overrides) - ids
    if unknown:
        raise ValueError('Profile references unknown candidate window IDs: ' + ', '.join(sorted(unknown)))
    return request


def assess_dataset_windows(dataset, profile, *, min_length, max_length, n_term='free', c_term='free'):
    """Return the editable request and assessment of every available window."""
    windows, gaps = enumerate_exposure_windows(dataset, min_length, max_length)
    request = window_exposure_request(windows, profile, n_term=n_term, c_term=c_term,
                                     score_source=dataset.config.score_expr if dataset.config else 'saved per_allele_scores; expression unavailable')
    result = assess_vaccine_windows(request)
    result['source_coverage_gaps'] = gaps
    result['enumeration'] = {'min_length': min_length, 'max_length': max_length,
                             'selection': 'all admitted target-containing windows; no rank/count pruning',
                             'target_definition': 'supplied binding-core intervals where documented; otherwise exact scored ligands'}
    return request, result
