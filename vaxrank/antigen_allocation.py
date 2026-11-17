"""Allocate unique, already-admitted antigen sequences without changing scores."""

from collections import Counter
from dataclasses import asdict, dataclass, field

import msgspec

from .vaccine_antigen import ANTIGEN_KINDS

ANTIGEN_SELECTION_TYPES = ANTIGEN_KINDS | {'ncORF'}


def normalize_selection_type(value):
    return next((kind for kind in ANTIGEN_SELECTION_TYPES if kind.lower() == value.lower()), value)


def selection_type(source, source_evidence=None):
    if not isinstance(source, dict):
        from .coding_context import coding_contexts
        source = dict(kind=source.kind, coding_contexts=coding_contexts(source, source_evidence))
    if source['kind'] == 'mutation' and any(
            context.get('kind') == 'extra_orf' for context in source.get('coding_contexts', ())):
        return 'ncORF'
    return source['kind']


def select_ranked_sources(ranked, policy=None, limit=None, *, source_evidence=None, audit=None):
    """Apply category switches before the shared ranked-source budget.

    Positive type caps count unique assembled sequences later. This stage only
    removes disabled categories; the native dataset retains their evidence.
    """
    policy = AntigenAllocation.resolve(policy)
    selected, decisions = [], []
    for source, peptides in ranked:
        eligible = [p for p in peptides if policy.max_per_kind.get(
            selection_type(p.antigen, source_evidence)) != 0]
        reason = ('antigen_type_disabled' if not eligible else 'ranked_source_limit'
                  if limit is not None and len(selected) >= limit else 'selected')
        if reason == 'selected':
            selected.append((source, eligible))
        if audit is not None:
            from .input_coverage import report_source_key
            decisions.append(dict(source_key=report_source_key(source), reason=reason))
    if audit is not None:
        audit['ranked_source_selection'] = decisions
    return selected


@dataclass(frozen=True)
class AntigenAllocation:
    # Historical modes remain readable for exact native replay. New selection
    # is always ranked; allocate_antigens rejects the retired balanced mode.
    mode: str = 'ranked'
    min_per_kind: dict[str, int] = field(default_factory=dict)
    max_per_kind: dict[str, int] = field(default_factory=dict)

    def __post_init__(self):
        if self.mode not in ('balanced', 'ranked'):
            raise ValueError('Unknown recorded antigen allocation mode: %r' % self.mode)
        for limits in (self.min_per_kind, self.max_per_kind):
            for kind, count in limits.items():
                if kind not in ANTIGEN_SELECTION_TYPES:
                    raise ValueError('Unknown antigen type %r; choose from %s' % (kind, sorted(ANTIGEN_SELECTION_TYPES)))
                if type(count) is not int or count < 0:
                    raise ValueError('Antigen type limits must be nonnegative integers')
        for kind, minimum in self.min_per_kind.items():
            if minimum > self.max_per_kind.get(kind, minimum):
                raise ValueError('Antigen minimum exceeds maximum for %s' % kind)

    @classmethod
    def resolve(cls, value):
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if set(value) - set(cls.__dataclass_fields__):
            raise ValueError('Unknown antigen allocation settings: %s' % sorted(set(value) - set(cls.__dataclass_fields__)))
        return msgspec.convert(value, cls)


def allocate_antigens(antigens, policy, audit):
    """Return a deterministic packing order after sequence deduplication.

    One product is attributed to the kind of its first enabled ranked source. All other
    source kinds remain in provenance, but cannot earn multiple allocation units.
    Minimums are goals; packing may truncate them. No eligible type is fabricated.
    """
    policy = AntigenAllocation.resolve(policy)
    if policy.mode != 'ranked':
        raise ValueError(
            'Balanced antigen allocation has been removed. Exact saved products can still be replayed; '
            'to redesign, set antigen_allocation.mode to ranked in the configuration.')
    counts = Counter()
    pending = list(enumerate(antigens))
    ordered, decisions = [], []
    while pending:
        eligible = [(index, antigen) for index, antigen in pending
                    if counts[selection_type(antigen.sources[0])] < policy.max_per_kind.get(
                        selection_type(antigen.sources[0]), len(antigens))]
        if not eligible:
            break
        def key(item):
            index, antigen = item
            kind = selection_type(antigen.sources[0])
            deficit = counts[kind] < policy.min_per_kind.get(kind, 0)
            return (not deficit, counts[kind] if deficit else 0, index)
        index, antigen = min(eligible, key=key)
        kind = selection_type(antigen.sources[0])
        reason = ('representation_goal' if counts[kind] < policy.min_per_kind.get(kind, 0)
                  else 'ranked_order')
        counts[kind] += 1
        ordered.append(antigen)
        decisions.append(dict(name=antigen.name, kind=kind, original_order=index + 1,
                              allocation_order=len(ordered), reason=reason))
        pending.remove((index, antigen))
    decisions.extend(dict(name=a.name, kind=selection_type(a.sources[0]), original_order=i + 1,
                          allocation_order=None, reason='antigen_type_cap') for i, a in pending)
    audit['allocation'] = dict(policy=asdict(policy), unit='unique_antigen_sequence',
                               attribution='first_enabled_ranked_source', decisions=decisions)
    return ordered


def finalize_allocation(audit, products):
    """Record actual emitted counts, including goals unmet after length/cap packing."""
    allocation = audit.get('allocation')
    if allocation is None:
        return
    emitted = {name for product in products for name in product.antigen_names}
    counts = Counter(d['kind'] for d in allocation['decisions'] if d['name'] in emitted)
    available = Counter(d['kind'] for d in allocation['decisions'])
    allocation['emitted_counts'] = dict(sorted(counts.items()))
    allocation['unmet_goals'] = [dict(kind=kind, requested=goal, emitted=counts[kind],
        available=available[kind], reason=('insufficient_eligible_sequences' if available[kind] < goal
                                         else 'construct_capacity'))
        for kind, goal in sorted(allocation['policy']['min_per_kind'].items()) if counts[kind] < goal]
