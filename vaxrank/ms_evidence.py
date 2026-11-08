"""Auditable support from exact, positive MHC-ligand MS observations.

Hitlist owns source curation and lineage. This consumer never derives donors
from sample labels, treats an inferred singleton as monoallelic, or converts
missing predictions to binding measurements. Arithmetic and admission policies
are evaluated by Topiary; this module resolves evidence identities and facets.
"""

from dataclasses import dataclass, field
from functools import lru_cache
import json
import math

from topiary import canonical_allele_name


FACET_WEIGHTS = dict(sequence=.25, predicted_sequence=.5, typed=1.,
                     recorded=1., predicted=3., experimental=6., monoallelic=8.)


@dataclass(frozen=True)
class MSParameters:
    gain: float = 1.
    max_relaxation: float | None = None
    reference_discount: float = .25
    weights: dict = field(default_factory=lambda: dict(FACET_WEIGHTS))
    prediction_ic50_cutoff: float = 50.
    direct_admission: bool = True

    def __post_init__(self):
        if type(self.direct_admission) is not bool:
            raise ValueError('MS direct_admission must be boolean')
        if set(self.weights) != set(FACET_WEIGHTS):
            raise ValueError('MS weights must name exactly: ' + ', '.join(FACET_WEIGHTS))
        for name, value in {**self.weights, 'gain': self.gain}.items():
            if not math.isfinite(value) or value < 0:
                raise ValueError(f'MS {name} must be finite and nonnegative')
        if self.max_relaxation is not None and (
                not math.isfinite(self.max_relaxation) or self.max_relaxation < 1):
            raise ValueError('MS max_relaxation must be null or finite and >= 1')
        if not math.isfinite(self.reference_discount) or not 0 <= self.reference_discount <= 1:
            raise ValueError('MS reference_discount must be in [0, 1]')
        if not math.isfinite(self.prediction_ic50_cutoff) or self.prediction_ic50_cutoff <= 0:
            raise ValueError('MS prediction_ic50_cutoff must be finite and positive')

    @classmethod
    def from_dict(cls, value=None):
        value = dict(value or {})
        if 'weights' in value:
            value['weights'] = {**FACET_WEIGHTS, **value['weights']}
        return cls(**value)


@lru_cache(maxsize=8192)
def exact_allele(value):
    """The named Topiary/mhcgnomes exact-allele rule; coarse names stay unknown."""
    if not value:
        return None
    try:
        return canonical_allele_name(value)
    except (ValueError, TypeError):
        return None


def _texts(value):
    if value is None:
        return ()
    if isinstance(value, str):
        if not value.strip():
            return ()
        if value.lstrip().startswith('['):
            value = json.loads(value)
        else:
            value = value.split(';')
    return tuple(sorted({str(x).strip() for x in value if x is not None and str(x).strip()}))


class IdentityGroups:
    """Connected source identities, formed before evidence qualification."""

    def __init__(self):
        self.parent = {}

    def root(self, key):
        self.parent.setdefault(key, key)
        parent = self.parent[key]
        if parent != key:
            self.parent[key] = self.root(parent)
        return self.parent[key]

    def join(self, keys):
        keys = tuple(keys)
        if not keys:
            return
        roots = sorted({self.root(key) for key in keys})
        for root in roots[1:]:
            self.parent[root] = roots[0]


def support_features(observations, allele, *, prediction_qualifies=False, parameters=None):
    """Resolve R/P/J for one exact peptide and patient allele.

    Each normalized observation declares ``eligible``, quality/context weights,
    references, optional producer-resolved ``unit_id`` and its native/engineered
    ``unit_kind``, exact restriction/basis, independently attributed typing and
    verified monoallelic status. Unknown axes are absent, never inferred here.
    ``reuse_ids`` join known reused experiments across publications. Records
    excluded by quality still connect identity groups.
    """
    parameters = parameters or MSParameters()
    allele = exact_allele(allele)
    if allele is None:
        raise ValueError('MS scoring requires an exact, expressed patient allele')
    observations = list({json.dumps(o, sort_keys=True): o for o in observations}.values())
    sequences = {o['peptide'] for o in observations}
    if len(sequences) > 1:
        raise ValueError('MS support cannot combine different peptide sequences')
    references = IdentityGroups()
    for o in observations:
        keys = ['reference:' + r for r in _texts(o.get('references'))]
        keys += ['reuse:' + r for r in _texts(o.get('reuse_ids'))]
        if o.get('unit_id'):
            keys.append('unit:' + o['unit_id'])
        references.join(keys)

    units, ref_support, details = {}, {}, []
    presence, direct = 0., False
    weights = parameters.weights
    for o in observations:
        q, t = o.get('quality', 1.), o.get('context_weight', 1.)
        for name, value in [('quality', q), ('context_weight', t)]:
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f'MS {name} must be finite and in [0, 1]')
        applicable = ['sequence']
        if prediction_qualifies:
            applicable.append('predicted_sequence')
        typed = allele in {exact_allele(a) for a in _texts(o.get('typed_alleles'))}
        match = allele == exact_allele(o.get('restriction'))
        if typed:
            applicable.append('typed')
        if match:
            applicable.append('recorded')
        if (typed and prediction_qualifies) or (
                match and o.get('source_prediction_qualified') is True):
            applicable.append('predicted')
        direct_facet = None
        if match and o.get('restriction_basis') == 'experimental':
            applicable.append('experimental')
            direct_facet = 'experimental'
        if match and o.get('monoallelic_verified') is True:
            applicable.append('monoallelic')
            direct_facet = 'monoallelic'
        eligible = o.get('eligible') is True
        facet = max(applicable, key=lambda name: (weights[name], name))
        contribution = q * t * weights[facet] if eligible else 0.
        presence = max(presence, contribution)
        if o.get('unit_id'):
            kind = o.get('unit_kind')
            if kind not in ('native_donor', 'engineered_origin'):
                raise ValueError('A resolved MS unit requires its explicit lineage kind')
            old = units.get(o['unit_id'])
            if old and old['kind'] != kind:
                raise ValueError('One resolved MS unit has conflicting material contexts')
            units[o['unit_id']] = dict(kind=kind, support=max(
                contribution, old['support'] if old else 0.))
        for ref in _texts(o.get('references')):
            group = references.root('reference:' + ref)
            ref_support[group] = max(ref_support.get(group, 0.), contribution)
        direct |= bool(eligible and direct_facet and q * t * weights[direct_facet] > 0)
        details.append(dict(observation_id=o['observation_id'], facet=facet,
                            applicable_facets=applicable, contribution=contribution,
                            eligible=eligible, exclusion_reasons=o.get('exclusion_reasons', [])))
    resolved = sum(u['support'] for u in units.values())
    bibliographic = parameters.reference_discount * sum(ref_support.values())
    # These are evidence features. The named Topiary criteria below compute
    # M/S/G and the prediction transformation from these three inputs.
    return dict(resolved=resolved, presence=presence, reference=bibliographic,
                direct=direct, units=units, reference_groups={
                    r: references.root(r) for r in sorted(references.parent)
                    if r.startswith('reference:')},
                reference_support=ref_support, observations=details)


def ms_criteria(parameters, base_expression, *, log_base_expression=None):
    """Topiary named criteria for the single MS-dependent power transform."""
    cap = parameters.max_relaxation
    g = 'criterion("ms_relaxation_raw")'
    if cap is not None:
        g = f'minimum({g}, {cap!r})'
    adjusted = ('criterion("ms_base") ** (1 / criterion("ms_relaxation"))'
                if log_base_expression is None else
                f'(({log_base_expression}) / criterion("ms_relaxation")).exp()')
    return [
        dict(name='ms_base', expression=base_expression, role='score'),
        dict(name='ms_mass', expression='maximum(ms_resolved, maximum(ms_presence, ms_reference))', role='score'),
        dict(name='ms_strength', expression='(1 + criterion("ms_mass")).log() / 0.6931471805599453', role='score'),
        dict(name='ms_relaxation_raw', expression=f'1 + {parameters.gain!r} * criterion("ms_strength")', role='score'),
        dict(name='ms_relaxation', expression=g, role='score'),
        dict(name='ms_adjusted', expression=adjusted, role='score'),
    ]


def validate_config(config):
    """Validate acquisition options separately from numerical policy knobs."""
    options = dict(config or {})
    allowed = {'parameters', 'lookup', 'package', 'base_score', 'soft_threshold',
               'log_base_score', 'ms_only_score'}
    if set(options) - allowed:
        raise ValueError(f'Unknown evidence.ms settings: {sorted(set(options) - allowed)}')
    MSParameters.from_dict(options.get('parameters'))
    if options.get('lookup', 'auto') not in ('auto', 'off', 'required'):
        raise ValueError('evidence.ms.lookup must be auto, off or required')
    threshold = options.get('soft_threshold')
    if threshold is not None and (not math.isfinite(threshold) or not 0 <= threshold <= 1):
        raise ValueError('evidence.ms.soft_threshold must be in [0, 1]')
    if options.get('ms_only_score') is not None and (
            not isinstance(options['ms_only_score'], (int, float))
            or not math.isfinite(options['ms_only_score'])
            or not 0 <= options['ms_only_score'] <= 1):
        raise ValueError('evidence.ms.ms_only_score must be an explicit finite score in [0, 1]')
    for name in ('base_score', 'log_base_score'):
        if options.get(name) is not None:
            from topiary.ranking import parse
            parse(options[name])
    return options


def score_frame(frame, cfg, *, group_keys, alleles=None, kind_support=None,
                default_methods=None, default_versions=None):
    """Run the named MS policy, retaining baseline, adjustment and admission.

    Raw normalized observations travel in ``hitlist_ms_records_json``. Replay
    recomputes facets from that captured evidence, never from installed data.
    A custom base requires an explicit bounded soft threshold; arbitrary user
    filters remain hard requirements. All arithmetic uses the Topiary DSL.
    """
    import numpy as np
    import pandas as pd
    from topiary import resolve_selection_policy, TopiaryResult, evaluate_selection_policy
    from topiary.ranking import EvalContext, parse
    from mhcgnomes import parse as parse_mhc

    config = validate_config(cfg.ms_evidence)
    parameters = MSParameters.from_dict(config.get('parameters'))
    context_options = dict(group_keys=group_keys, alleles=alleles,
                           kind_support=kind_support, default_methods=default_methods,
                           default_versions=default_versions)
    context = EvalContext(frame, **context_options)
    if config.get('base_score') is not None:
        if config.get('soft_threshold') is None:
            raise ValueError('A custom MS base_score requires an explicit soft_threshold')
        base_expression = config['base_score']
        cutoff = config['soft_threshold']
        operator = '>='
    else:
        if cfg.score_expr is not None or cfg.scoring_mode != 'affinity':
            raise ValueError('MS scaling of a custom score/rank requires base_score and soft_threshold')
        midpoint, width = cfg.logistic_epitope_score_midpoint, cfg.logistic_epitope_score_width
        base_expression = f'affinity.value.logistic_normalized({midpoint!r}, {width!r})'
        # Apply the SAME producer transform to the raw cutoff, preserving the
        # strict affinity cutoff and independent inclusive minimum-score gate.
        scalar = f'({cfg.binding_affinity_cutoff!r}).logistic_normalized({midpoint!r}, {width!r})'
        cutoff = float(parse(scalar).eval(context).iloc[0])
        operator = '>'
    base = parse(base_expression).eval(context).reindex(context.group_index)
    if config.get('log_base_score') is not None:
        log_base = parse(config['log_base_score']).eval(context).reindex(context.group_index)
        if (log_base.isna().ne(base.isna()).any() or log_base.gt(0).any()
                or not np.isclose(np.exp(log_base), base, rtol=1e-8, atol=0, equal_nan=True).all()):
            raise ValueError('log_base_score must be the natural logarithm of base_score, including missing values')
    present = base.dropna()
    if not np.isfinite(present).all() or not present.between(0, 1).all():
        raise ValueError('MS base score must be finite and in [0, 1], or explicitly missing')
    if 'pMHC_affinity' in set(frame.get('kind', [])):
        qualification = parse(f'affinity.value <= {parameters.prediction_ic50_cutoff!r}').eval(context)
    else:
        qualification = pd.Series(False, index=context.group_index)

    # The lookup annotation is occurrence-invariant but tied to exact peptide
    # and MHC species. Group expansion for the patient genotype must preserve it.
    captured = {}
    if 'hitlist_ms_records_json' in frame:
        for peptide, text in frame[['peptide', 'hitlist_ms_records_json']].itertuples(index=False, name=None):
            if not isinstance(text, str):
                continue
            if peptide in captured and captured[peptide] != text:
                raise ValueError('Conflicting captured MS observations for one peptide')
            captured[peptide] = text
    rows, evidence = [], {}
    peptide_index, allele_index = group_keys.index('peptide'), group_keys.index('allele')
    for key in context.group_index:
        peptide, allele = key[peptide_index], key[allele_index]
        if not allele:
            continue
        molecule = parse_mhc(allele, infer_class2_pairing=False)
        observed = json.loads(captured.get(peptide, '[]'))
        # A viral source organism does not change the restricting molecule's species.
        observations = [o for o in observed if o.get('mhc_species') == molecule.species.name
                        and o.get('mhc_class') in (None, '', 'I' if molecule.is_class1 else 'II')]
        qualifies = (molecule.is_class1 and
                     qualification.get(key) is not pd.NA and
                     pd.notna(qualification.get(key)) and bool(qualification.get(key)))
        features = support_features(observations, allele,
                                    prediction_qualifies=qualifies, parameters=parameters)
        evidence[key] = features
        rows.append(dict(zip(group_keys, key), ms_base_value=base.get(key),
                         ms_resolved=features['resolved'], ms_presence=features['presence'],
                         ms_reference=features['reference'], ms_direct=features['direct']))
    annotations = pd.DataFrame(rows)
    if annotations.empty:
        return pd.Series(dtype=float)
    # Expand only score-group annotations, retaining predictor rows and names.
    augmented = frame.merge(annotations, on=list(group_keys), how='outer', validate='many_to_one')
    augmented.attrs = dict(frame.attrs)
    criteria = ms_criteria(parameters, 'ms_base_value', log_base_expression=config.get('log_base_score'))
    def gate(value):
        return f'({value} {operator} {cutoff!r}) & ({value} >= {cfg.min_epitope_score!r})'
    baseline_gate = gate('criterion("ms_base")')
    adjusted_gate = gate('criterion("ms_adjusted")')
    direct = '(ms_direct > 0) & (criterion("ms_relaxation") > 1)'
    if not parameters.direct_admission:
        direct = '0 > 1'
    criteria += [dict(name='ms_baseline_gate', expression=baseline_gate, role='eligibility'),
                 dict(name='ms_adjusted_gate', expression=adjusted_gate, role='eligibility'),
                 dict(name='ms_direct_admission', expression=direct, role='eligibility')]
    numeric = 'criterion("ms_adjusted")'
    if config.get('ms_only_score') is not None:
        # Explicit, named missing-prediction policy; exact zero is never missing.
        criteria.append(dict(name='ms_only_score', expression=repr(config['ms_only_score']), role='score'))
        numeric = f'criterion("ms_adjusted").fillna({config["ms_only_score"]!r})'
    admitted = '(criterion("ms_adjusted_gate") | criterion("ms_direct_admission"))'
    if cfg.filter_expr is not None:
        admitted = f'({cfg.filter_expr}) & {admitted}'
    policy = resolve_selection_policy(dict(name='ms-evidence-v1', criteria=criteria,
        score_by=numeric, filter_by=admitted, default_methods=default_methods))
    from dataclasses import replace
    policy = replace(policy, default_versions=default_versions)
    evaluation = evaluate_selection_policy(TopiaryResult(augmented, form='long'), policy,
        group_keys=group_keys, alleles=alleles, kind_support=kind_support)
    selected = evaluation.selected
    if not selected.empty:
        finite = selected.score.dropna()
        if not np.isfinite(finite).all() or not finite.between(0, 1).all():
            raise ValueError('MS adjusted/MS-only scores must be finite and in [0, 1]')
    scores = (selected.set_index(list(group_keys)).score.astype(float)
              if not selected.empty else pd.Series(dtype=float))
    scores = scores.loc[np.isfinite(scores)]
    baseline_decisions = parse(f'({base_expression}) {operator} {cutoff!r}').eval(context)
    computed = {tuple(getattr(row, name) for name in group_keys):
                {c['criterion']: c['value'] for c in row.criteria}
                for row in evaluation.occurrences.itertuples()}
    records = {}
    for key, adjusted in scores.items():
        b = base.get(key)
        passed = bool(pd.notna(b) and b >= cfg.min_epitope_score
                      and pd.notna(baseline_decisions.get(key)) and baseline_decisions.get(key))
        b0 = float(b) if passed else 0.
        r = dict(evidence[key], base=float(b) if pd.notna(b) else None,
                 adjusted=float(adjusted), baseline=b0,
                 mass=computed[key]['ms_mass'], strength=computed[key]['ms_strength'],
                 relaxation=computed[key]['ms_relaxation'],
                 score_basis='prediction' if pd.notna(b) else 'ms_only')
        records.setdefault(key[:-1], {})[key[-1]] = r
    scores.attrs.update(policy_evaluation=evaluation, ms_scoring=records)
    return scores
