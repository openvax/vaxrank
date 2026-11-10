"""Opt-in window heuristics over scored ligands and typed cleavage evidence.

These are design utilities, not calibrated immunogenicity or serum survival
probabilities. See docs/selection-policies.md for assumptions and sources.
"""

from dataclasses import asdict, dataclass, replace
import math

import msgspec
from mhctools import PeptideInput
from mhctools.peptidases import get_cleavage_model

from .cleavage_inference import audit_peptidase_inputs
from .epitope_logic import slice_epitopes
from .ranking import DEFAULT_RANKING_RULES


SERUM_PANEL = ('dpp4-qpisa', 'fap-endo-gp', 'fap-dipeptidyl',
               'ace-dipeptidyl', 'cpn-basic')


@dataclass(frozen=True)
class WindowSelection:
    self_weight: float = 0.0
    serum_weight: float = 0.0
    min_target_fraction: float = .95
    enumerate_lengths: bool = True
    serum_models: tuple[str, ...] = SERUM_PANEL
    dpp4_depletion_threshold: float = 1.0

    def __post_init__(self):
        if type(self.enumerate_lengths) is not bool:
            raise ValueError('enumerate_lengths must be a boolean')
        for name in ('self_weight', 'serum_weight', 'min_target_fraction', 'dpp4_depletion_threshold'):
            value = getattr(self, name)
            if isinstance(value, bool) or not math.isfinite(value) or value < 0:
                raise ValueError('%s must be finite and nonnegative' % name)
        if self.serum_weight > 1 or not 0 < self.min_target_fraction <= 1:
            raise ValueError('serum_weight must be <= 1 and min_target_fraction in (0, 1]')
        if not self.serum_models or len(set(self.serum_models)) != len(self.serum_models):
            raise ValueError('serum_models must be a nonempty, unique panel')
        if set(self.serum_models) - set(SERUM_PANEL):
            raise ValueError('Use reviewed serum models: %s' % ', '.join(SERUM_PANEL))

    @classmethod
    def resolve(cls, value):
        if isinstance(value, cls):
            return value
        unknown = set(value) - set(cls.__dataclass_fields__)
        if unknown:
            raise ValueError('Unknown window selection settings: %s' % sorted(unknown))
        return msgspec.convert(value, cls)


def window_metrics(peptide, policy, *, n_term='free', c_term='free'):
    """Count each sequence/allele once; a surviving copy preserves its content.

    Only a supported cut strictly inside an epitope destroys that epitope.
    Terminal trimming outside it may release it and is not penalized here.
    We evaluate first cuts on the actual substrate, not an assumed cascade.
    """
    profiles = ()
    cuts = set()
    if policy.serum_weight:
        substrate = PeptideInput(peptide.amino_acids, n_term=n_term, c_term=c_term)
        profiles = tuple(audit_peptidase_inputs([substrate], get_cleavage_model(name))[0]
                         for name in policy.serum_models)
        for profile in profiles:
            for site in profile.sites:
                if (site.status == 'matched' or (
                        profile.model.name == 'dpp4-qpisa' and site.status == 'scored'
                        and site.score >= policy.dpp4_depletion_threshold)):
                    cuts.add(site.bond)
    target, protected, self_scores = {}, {}, {}
    target_ms, protected_ms, self_ms = {}, {}, {}
    target_details, self_details = [], []
    for epitope in peptide.target_epitopes:
        internal_cuts = sorted(bond for bond in cuts
                               if epitope.offset < bond < epitope.offset + len(epitope.sequence))
        damaged = bool(internal_cuts)
        increment = epitope.ms_increment
        if not math.isfinite(increment) or increment < 0:
            raise ValueError('Window utility requires finite, nonnegative MS increments')
        target_ms[epitope.sequence] = max(target_ms.get(epitope.sequence, 0.), increment)
        protected_ms[epitope.sequence] = max(protected_ms.get(epitope.sequence, 0.),
                                             0. if damaged else increment)
        for allele, score in epitope.per_allele_scores.items():
            if not math.isfinite(score) or score < 0:
                raise ValueError('Window utility requires finite, nonnegative epitope scores')
            key = epitope.sequence, allele
            target[key] = max(target.get(key, 0.), score)
            protected[key] = max(protected.get(key, 0.), 0. if damaged else score)
        target_details.append(dict(sequence=epitope.sequence, offset=epitope.offset,
                                   **({'ms_increment': increment} if increment else {}),
                                   allele_scores=epitope.per_allele_scores, internal_cuts=internal_cuts))
    for epitope in peptide.epitopes:
        if epitope.occurs_in_non_CTA_reference:
            increment = epitope.ms_increment
            if not math.isfinite(increment) or increment < 0:
                raise ValueError('Window utility requires finite, nonnegative MS increments')
            self_ms[epitope.sequence] = max(self_ms.get(epitope.sequence, 0.), increment)
            for allele, score in epitope.per_allele_scores.items():
                if not math.isfinite(score) or score < 0:
                    raise ValueError('Window utility requires finite, nonnegative self scores')
                key = epitope.sequence, allele
                self_scores[key] = max(self_scores.get(key, 0.), score)
            self_details.append(dict(sequence=epitope.sequence, offset=epitope.offset,
                                     allele_scores=epitope.per_allele_scores,
                                     source_match=(asdict(epitope.self_reference_match)
                                                   if epitope.self_reference_match else None)))
    total = sum(target.values()) + sum(target_ms.values())
    at_risk = total - sum(protected.values()) - sum(protected_ms.values())
    burden = sum(self_scores.values()) + sum(self_ms.values())
    utility = (total - policy.serum_weight * at_risk) / (1 + policy.self_weight * burden)
    return dict(target_score=total, non_cta_self_score=burden, at_risk_target_score=at_risk,
                utility=utility, cuts=sorted(cuts),
                target_details=target_details, non_cta_self_details=self_details,
                self_provenance_complete=all(e.self_reference_match is not None and
                    e.self_reference_match.source_provenance_complete for e in peptide.epitopes),
                cleavage_profiles=[asdict(p) for p in profiles],
                unassessed_models=[p.model.name for p in profiles
                                   if p.status in ('unassessed', 'failed')])


def select_windows(candidates, policy, *, preferred_length, limit,
                   n_term='free', c_term='free', audit_out=None):
    """Retain target content first, then minimize self/degradation burden."""
    policy = WindowSelection.resolve(policy)
    if not candidates:
        return []
    metrics = [window_metrics(vp, policy, n_term=n_term, c_term=c_term) for vp in candidates]
    best_target = max(m['target_score'] for m in metrics)
    alternatives = []
    retained = []
    for index, (vp, metric) in enumerate(zip(candidates, metrics)):
        metric['source'] = {key: getattr(vp.antigen, key) for key in (
            'source_identifier', 'gene_name', 'gene_id', 'transcript_ids',
            'protein_ids', 'species', 'source_metadata')}
        eligible = (best_target > 0 and metric['target_score'] >= best_target * policy.min_target_fraction
                    and not metric['unassessed_models'])
        parent_window = vp.window_selection_audit.get('parent_window')
        vp.window_selection_audit = dict(policy=asdict(policy), **metric)
        ordering_key = ([-metric['utility'], abs(len(vp.amino_acids) - preferred_length),
                         list(vp.lexicographic_sort_key()), index] if eligible else None)
        alternatives.append(dict(sequence=vp.amino_acids, eligible=eligible,
                                 input_order=index, ordering_key=ordering_key,
                                 ranking_rules=list(vp.ranking_rules if vp.ranking_rules is not None
                                                    else DEFAULT_RANKING_RULES),
                                 ineligibility_reason=('unassessed_serum_model' if metric['unassessed_models']
                                                       else 'target_retention' if not eligible else None),
                                 **{k: v for k, v in metric.items() if k not in (
                                     'cleavage_profiles', 'target_details', 'non_cta_self_details')}))
        if parent_window is not None:
            alternatives[-1]['parent_window'] = parent_window
        if eligible:
            retained.append(index)
    retained.sort(key=lambda index: alternatives[index]['ordering_key'])
    ranks = {index: rank for rank, index in enumerate(retained, 1)}
    chosen = set(retained[:limit])
    for index, (vp, alternative) in enumerate(zip(candidates, alternatives)):
        alternative.update(selection_rank=ranks.get(index), selected=index in chosen,
                           decision_reason=(alternative['ineligibility_reason'] or
                                            ('within_window_limit' if index in chosen else 'window_limit')))
        vp.window_selection_audit.update(
            decision={key: alternative[key] for key in (
                'input_order', 'ordering_key', 'ranking_rules', 'selection_rank', 'selected', 'decision_reason')},
            preferred_length=preferred_length, limit=limit)
        if 'parent_window' in alternative:
            vp.window_selection_audit['parent_window'] = alternative['parent_window']
    selected = [candidates[index] for index in retained[:limit]]
    for vp in selected:
        vp.window_selection_audit['alternatives'] = alternatives
    if audit_out is not None:
        audit_out.append(dict(policy=asdict(policy), alternatives=alternatives,
                              preferred_length=preferred_length, limit=limit,
                              selected=[vp.amino_acids for vp in selected]))
    return selected


def record_ranked_window_decisions(candidates, retained, selected, *, score_fraction, fraction_applied, limit):
    """Capture the existing score-fraction/lexicographic selection, without reranking."""
    ranks = {id(peptide): rank for rank, peptide in enumerate(retained, 1)}
    chosen = {id(peptide) for peptide in selected}
    alternatives = []
    for index, peptide in enumerate(candidates):
        eligible = id(peptide) in ranks
        alternatives.append(dict(sequence=peptide.amino_acids, input_order=index,
            eligible=eligible, selected=id(peptide) in chosen, selection_rank=ranks.get(id(peptide)),
            decision_reason=('score_fraction' if not eligible else
                             'within_window_limit' if id(peptide) in chosen else 'window_limit'),
            ordering_key=[list(peptide.lexicographic_sort_key()), index] if eligible else None,
            ranking_rules=list(peptide.ranking_rules if peptide.ranking_rules is not None else DEFAULT_RANKING_RULES),
            combined_score=peptide.combined_score, target_score=peptide.target_epitope_score,
            parent_window=peptide.window_selection_audit.get('parent_window')))
    for peptide in selected:
        index = next(i for i, candidate in enumerate(candidates) if candidate is peptide)
        decision = alternatives[index]
        peptide.window_selection_audit = dict(
            algorithm='score_fraction_then_ranking',
            policy=dict(score_fraction_of_best=score_fraction, fraction_applied=fraction_applied),
            limit=limit, decision={key: decision[key] for key in (
                'input_order', 'ordering_key', 'ranking_rules', 'selection_rank', 'selected', 'decision_reason')},
            parent_window=decision['parent_window'], alternatives=alternatives,
            target_score=peptide.target_epitope_score, combined_score=peptide.combined_score,
            target_details=[dict(sequence=e.sequence, offset=e.offset, allele_scores=e.per_allele_scores)
                            for e in peptide.target_epitopes])


def target_content(peptide, *, positioned=False):
    """Target peptide/allele identities, optionally anchored to the source."""
    start = peptide.window_selection_audit.get('parent_window', {}).get('start', 0)
    return {(epitope.sequence, allele, start + epitope.offset) if positioned else (epitope.sequence, allele)
            for epitope in peptide.target_epitopes for allele in epitope.per_allele_scores}


def overlapping_windows(left, right):
    """Shared source residues, using zero-based half-open parent coordinates."""
    return (left['sequence'] == right['sequence']
            and left['start'] < right['end'] and right['start'] < left['end'])


def select_source_windows(candidates, *, limit=3, score_fraction=.99, policy=None,
                          preferred_length=25, n_term='free', c_term='free', audit_out=None,
                          required_content=None):
    """Greedily select distinct regions; near-best comparisons restart per region.

    Candidates must carry coordinates on one translated source. Scores remain
    the configured DSL scores; this adds a geometric constraint, not a score.
    """
    if type(limit) is not int or limit < 1:
        raise ValueError('Source window limit must be a positive integer')
    candidates = list(candidates)
    if not candidates:
        return []
    parents = [vp.window_selection_audit.get('parent_window') for vp in candidates]
    if any(p is None or p['sequence'][p['start']:p['end']] != vp.amino_acids
           for p, vp in zip(parents, candidates)):
        raise ValueError('Source windows require exact parent sequence coordinates')
    if len({p['sequence'] for p in parents}) != 1:
        raise ValueError('Select windows separately for competing translated source contexts')
    if policy is not None:
        resolved_policy = WindowSelection.resolve(policy)
        for vp in candidates:
            # Populate every score before comparing regions: a DSL expression
            # may use window utility, including self/cleavage penalties.
            vp.window_selection_audit.update(window_metrics(vp, resolved_policy, n_term=n_term, c_term=c_term))
    contents = [target_content(vp) for vp in candidates]
    available = [i for i, vp in enumerate(candidates) if not required_content
                 or required_content <= target_content(vp, positioned=True)]
    selected, used_content = [], set()
    decisions = {i: 'target_preservation' for i in range(len(candidates)) if i not in available}
    while available and len(selected) < limit:
        best = max(candidates[i].combined_score for i in available)
        near_best = [i for i in available if math.isclose(best, 0.) or
                     candidates[i].combined_score >= best - abs(best) * (1 - score_fraction)]
        if policy is None:
            chosen = min(near_best, key=lambda i: (candidates[i].lexicographic_sort_key(), i))
        else:
            chosen_peptides = select_windows([candidates[i] for i in near_best], policy,
                preferred_length=preferred_length, limit=1, n_term=n_term, c_term=c_term)
            if not chosen_peptides:
                for i in near_best:
                    decisions[i] = 'window_policy'
                available = [i for i in available if i not in near_best]
                continue
            chosen = next(i for i in near_best if candidates[i] is chosen_peptides[0])
        selected.append(chosen)
        used_content.update(contents[chosen])
        decisions[chosen] = 'selected_source_window'
        remaining = []
        for i in available:
            if i == chosen:
                continue
            if overlapping_windows(parents[i], parents[chosen]):
                decisions[i] = 'overlaps_selected_window'
            elif contents[i] and contents[i] <= used_content:
                decisions[i] = 'no_additional_target_content'
            else:
                remaining.append(i)
        available = remaining
    for i in available:
        decisions[i] = 'source_window_limit'
    alternatives = [dict(sequence=vp.amino_acids, parent_window=parents[i], input_order=i,
        eligible=decisions[i] not in ('window_policy', 'target_preservation'), selected=i in selected,
        selection_rank=selected.index(i) + 1 if i in selected else None,
        decision_reason=decisions[i], combined_score=vp.combined_score, target_score=vp.target_epitope_score,
        utility=vp.window_selection_audit.get('utility'),
        source={key: getattr(vp.antigen, key) for key in (
            'source_identifier', 'gene_name', 'gene_id', 'transcript_ids', 'protein_ids', 'species', 'source_metadata')},
        ranking_rules=list(vp.ranking_rules if vp.ranking_rules is not None else DEFAULT_RANKING_RULES),
        ordering_key=[list(vp.lexicographic_sort_key()), i]) for i, vp in enumerate(candidates)]
    result = []
    for rank, i in enumerate(selected, 1):
        vp = candidates[i]
        required = target_content(vp, positioned=True)
        manufacturing = [dict(sequence=alt.amino_acids, parent_window=parents[j])
            for j, alt in enumerate(candidates) if j != i and required
            and required <= target_content(alt, positioned=True)
            and not any(overlapping_windows(parents[j], parents[k]) for k in selected if k != i)]
        vp.window_selection_audit.update(
            algorithm='non_overlapping_source_windows',
            policy={**(asdict(WindowSelection.resolve(policy)) if policy is not None else {}),
                    'score_fraction_of_best': score_fraction},
            limit=limit, preferred_length=preferred_length, parent_window=parents[i],
            source_window=dict(index=rank, start=parents[i]['start'], end=parents[i]['end'],
                               required_targets=[list(item) for item in sorted(required)]),
            decision={key: alternatives[i][key] for key in (
                'input_order', 'ordering_key', 'ranking_rules', 'selection_rank', 'selected', 'decision_reason')},
            combined_score=vp.combined_score, target_score=vp.target_epitope_score,
            manufacturing_alternatives=manufacturing, alternatives=alternatives)
        result.append(vp)
    if audit_out is not None:
        audit_out.append(dict(algorithm='non_overlapping_source_windows', limit=limit,
                              alternatives=alternatives, selected=[vp.amino_acids for vp in result]))
    return result


def optimize_peptide_windows(ranked, options, audit_out=None, *, source_limit=None,
                             score_fraction=.99, required_content=None):
    """Refine available windows independently for the peptide modality.

    SLPs only: concatenated constructs and minimal epitope products require a
    different enumeration context. No new predictions are run during selection.
    """
    from .ranking import rank_constructs
    policy = WindowSelection.resolve(options.window_selection or {})
    output = []
    for source, peptides in ranked:
        # A manufacturing refinement operates within each biological region;
        # it cannot replace three independent windows with three nearby shifts.
        if source_limit is None and any(p.window_selection_audit.get('source_window') for p in peptides):
            refined = []
            for peptide in peptides:
                required = set(tuple(x) for x in peptide.window_selection_audit['source_window']['required_targets'])
                items = optimize_peptide_windows([(source, [replace(peptide,
                    window_selection_audit={k: v for k, v in peptide.window_selection_audit.items()
                                            if k != 'source_window'})])],
                    options, audit_out, source_limit=1, score_fraction=score_fraction, required_content=required)
                # Refinement must preserve every selected target occurrence.
                chosen = next((p for _, ps in items for p in ps
                               if required <= target_content(p, positioned=True)), None)
                if chosen is None:
                    continue
                else:
                    chosen.window_selection_audit['source_window'] = peptide.window_selection_audit['source_window']
                refined.append(chosen)
            if refined:
                output.append((source, refined))
            continue
        candidates = []
        for peptide in peptides:
            maximum = min(len(peptide.amino_acids), options.max_antigen_length_aa)
            lengths = (range(options.min_antigen_length_aa, maximum + 1)
                       if policy.enumerate_lengths else [maximum])
            if source_limit is not None and required_content is None and maximum < options.min_antigen_length_aa:
                # Imported contexts can be shorter than a requested window;
                # retain their real sequence (the writer reports the length).
                # Manufacturing refinement still enforces its own minimum.
                lengths = [maximum]
            for length in lengths:
                fragment = peptide.mutant_protein_fragment
                windows = (fragment.sorted_subsequences(length) if fragment is not None else
                           [(i, None) for i in range(len(peptide.amino_acids) - length + 1)])
                for start, part in windows:
                    end = start + length
                    if not peptide.antigen.interval_is_targetable(start, end):
                        continue
                    antigen = peptide.antigen.sliced(start, end)
                    original_parent = peptide.window_selection_audit.get('parent_window') or dict(
                        sequence=peptide.amino_acids, start=0, end=len(peptide.amino_acids))
                    candidates.append(replace(
                        peptide, mutant_protein_fragment=part, antigen=antigen,
                        epitopes=slice_epitopes(peptide.epitopes, start, end),
                        combined_score_expr=(peptide.combined_score_expr if options.window_selection is None or any(
                            name in peptide.combined_score_expr for name in (
                                'window_epitope_score', 'window_selection_factor')) else
                            '(' + peptide.combined_score_expr + ') * window_selection_factor'),
                        window_selection_audit={'parent_window': dict(
                            sequence=original_parent['sequence'], start=original_parent['start'] + start,
                            end=original_parent['start'] + end)}))
        if source_limit is not None and getattr(options, 'require_target_epitopes', True):
            candidates = [p for p in candidates if p.contains_target_epitopes()]
        if source_limit is not None:
            chosen = select_source_windows(candidates, policy=options.window_selection,
                limit=source_limit, score_fraction=score_fraction,
                preferred_length=options.max_antigen_length_aa,
                n_term='acetylated' if options.n_terminal_acetylation else 'free',
                c_term='amidated' if options.c_terminal_amidation else 'free', audit_out=audit_out,
                required_content=required_content)
        else:
            chosen = select_windows(candidates, policy,
                                preferred_length=options.max_antigen_length_aa,
                                limit=options.candidates_per_slot,
                                n_term='acetylated' if options.n_terminal_acetylation else 'free',
                                c_term='amidated' if options.c_terminal_amidation else 'free',
                                audit_out=audit_out)
        if chosen:
            output.append((source, chosen))
    return rank_constructs(output)


def optimize_antigen_windows(ranked, config):
    """Use shared window settings for table-backed antigen contexts as well."""
    from types import SimpleNamespace
    options = SimpleNamespace(
        window_selection=config.window_selection,
        min_antigen_length_aa=config.min_peptide_length,
        max_antigen_length_aa=config.max_peptide_length,
        candidates_per_slot=config.max_vaccine_peptides_per_variant,
        require_target_epitopes=config.require_target_epitopes_in_variant,
        n_terminal_acetylation=False, c_terminal_amidation=False)
    return optimize_peptide_windows(ranked, options, source_limit=config.max_windows_per_antigen_source,
                                    score_fraction=config.score_fraction_of_best)
