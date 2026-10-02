"""Opt-in window heuristics over scored ligands and typed cleavage evidence.

These are design utilities, not calibrated immunogenicity or serum survival
probabilities. See docs/selection-policies.md for assumptions and sources.
"""

from dataclasses import asdict, dataclass, replace
import math

import msgspec
from mhctools.cleavage import CleavageInput
from mhctools.peptidases import get_cleavage_model

from .cleavage_inference import audit_peptidase_inputs
from .epitope_logic import slice_epitopes


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
        substrate = CleavageInput(peptide.amino_acids, n_term=n_term, c_term=c_term)
        profiles = tuple(audit_peptidase_inputs([substrate], get_cleavage_model(name))[0]
                         for name in policy.serum_models)
        for profile in profiles:
            for site in profile.sites:
                if (site.status == 'matched' or (
                        profile.model.name == 'dpp4-qpisa' and site.status == 'scored'
                        and site.score >= policy.dpp4_depletion_threshold)):
                    cuts.add(site.bond)
    target, protected, self_scores = {}, {}, {}
    target_details, self_details = [], []
    for epitope in peptide.target_epitopes:
        internal_cuts = sorted(bond for bond in cuts
                               if epitope.offset < bond < epitope.offset + len(epitope.sequence))
        damaged = bool(internal_cuts)
        for allele, score in epitope.per_allele_scores.items():
            if not math.isfinite(score) or score < 0:
                raise ValueError('Window utility requires finite, nonnegative epitope scores')
            key = epitope.sequence, allele
            target[key] = max(target.get(key, 0.), score)
            protected[key] = max(protected.get(key, 0.), 0. if damaged else score)
        target_details.append(dict(sequence=epitope.sequence, offset=epitope.offset,
                                   allele_scores=epitope.per_allele_scores, internal_cuts=internal_cuts))
    for epitope in peptide.epitopes:
        if epitope.occurs_in_non_CTA_reference:
            for allele, score in epitope.per_allele_scores.items():
                if not math.isfinite(score) or score < 0:
                    raise ValueError('Window utility requires finite, nonnegative self scores')
                key = epitope.sequence, allele
                self_scores[key] = max(self_scores.get(key, 0.), score)
            self_details.append(dict(sequence=epitope.sequence, offset=epitope.offset,
                                     allele_scores=epitope.per_allele_scores,
                                     source_match=(asdict(epitope.self_reference_match)
                                                   if epitope.self_reference_match else None)))
    total = sum(target.values())
    at_risk = total - sum(protected.values())
    burden = sum(self_scores.values())
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
    for vp, metric in zip(candidates, metrics):
        metric['source'] = {key: getattr(vp.antigen, key) for key in (
            'source_identifier', 'gene_name', 'gene_id', 'transcript_ids',
            'protein_ids', 'species', 'source_metadata')}
        eligible = (best_target > 0 and metric['target_score'] >= best_target * policy.min_target_fraction
                    and not metric['unassessed_models'])
        alternatives.append(dict(sequence=vp.amino_acids, eligible=eligible,
                                 ineligibility_reason=('unassessed_serum_model' if metric['unassessed_models']
                                                       else 'target_retention' if not eligible else None),
                                 **{k: v for k, v in metric.items() if k not in (
                                     'cleavage_profiles', 'target_details', 'non_cta_self_details')}))
        vp.window_selection_audit = dict(policy=asdict(policy), **metric)
        if eligible:
            retained.append(vp)
    retained.sort(key=lambda vp: (-vp.window_selection_audit['utility'],
                                 abs(len(vp.amino_acids) - preferred_length),
                                 vp.lexicographic_sort_key()))
    selected = retained[:limit]
    for vp in selected:
        vp.window_selection_audit['alternatives'] = alternatives
    if audit_out is not None:
        audit_out.append(dict(policy=asdict(policy), alternatives=alternatives,
                              selected=[vp.amino_acids for vp in selected]))
    return selected


def optimize_peptide_windows(ranked, options, audit_out=None):
    """Refine available windows independently for the peptide modality.

    SLPs only: concatenated constructs and minimal epitope products require a
    different enumeration context. No new predictions are run during selection.
    """
    from .ranking import rank_constructs
    policy = WindowSelection.resolve(options.window_selection)
    output = []
    for source, peptides in ranked:
        candidates = []
        for peptide in peptides:
            maximum = min(len(peptide.amino_acids), options.max_antigen_length_aa)
            lengths = (range(options.min_antigen_length_aa, maximum + 1)
                       if policy.enumerate_lengths else [maximum])
            for length in lengths:
                fragment = peptide.mutant_protein_fragment
                windows = (fragment.sorted_subsequences(length) if fragment is not None else
                           [(i, None) for i in range(len(peptide.amino_acids) - length + 1)])
                for start, part in windows:
                    end = start + length
                    if not peptide.antigen.interval_is_targetable(start, end):
                        continue
                    antigen = peptide.antigen.sliced(start, end)
                    candidates.append(replace(
                        peptide, mutant_protein_fragment=part, antigen=antigen,
                        epitopes=slice_epitopes(peptide.epitopes, start, end),
                        combined_score_expr=(peptide.combined_score_expr if any(
                            name in peptide.combined_score_expr for name in (
                                'window_epitope_score', 'window_selection_factor')) else
                            '(' + peptide.combined_score_expr + ') * window_selection_factor'),
                        window_selection_audit={}))
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
        n_terminal_acetylation=False, c_terminal_amidation=False)
    return optimize_peptide_windows(ranked, options)
