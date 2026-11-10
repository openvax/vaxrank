"""Allocate unique, already-admitted antigen sequences without changing scores."""

from collections import Counter
from dataclasses import asdict, dataclass, field

import msgspec

from .vaccine_antigen import ANTIGEN_KINDS


@dataclass(frozen=True)
class AntigenAllocation:
    mode: str = 'balanced'
    min_per_kind: dict[str, int] = field(default_factory=dict)
    max_per_kind: dict[str, int] = field(default_factory=dict)

    def __post_init__(self):
        if self.mode not in ('balanced', 'ranked'):
            raise ValueError('antigen allocation mode must be balanced or ranked')
        for limits in (self.min_per_kind, self.max_per_kind):
            for kind, count in limits.items():
                if kind not in ANTIGEN_KINDS:
                    raise ValueError('Unknown antigen type %r; choose from %s' % (kind, sorted(ANTIGEN_KINDS)))
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

    One product is attributed to the kind of its first ranked source. All other
    source kinds remain in provenance, but cannot earn multiple allocation units.
    Minimums are goals; packing may truncate them. No eligible type is fabricated.
    """
    policy = AntigenAllocation.resolve(policy)
    counts = Counter()
    pending = list(enumerate(antigens))
    ordered, decisions = [], []
    while pending:
        eligible = [(index, antigen) for index, antigen in pending
                    if counts[antigen.sources[0]['kind']] < policy.max_per_kind.get(
                        antigen.sources[0]['kind'], len(antigens))]
        if not eligible:
            break
        def key(item):
            index, antigen = item
            kind = antigen.sources[0]['kind']
            deficit = counts[kind] < policy.min_per_kind.get(kind, 0)
            return (not deficit, counts[kind] if deficit or policy.mode == 'balanced' else 0, index)
        index, antigen = min(eligible, key=key)
        kind = antigen.sources[0]['kind']
        reason = ('representation_goal' if counts[kind] < policy.min_per_kind.get(kind, 0)
                  else 'balanced_types' if policy.mode == 'balanced' else 'ranked_order')
        counts[kind] += 1
        ordered.append(antigen)
        decisions.append(dict(name=antigen.name, kind=kind, original_order=index + 1,
                              allocation_order=len(ordered), reason=reason))
        pending.remove((index, antigen))
    decisions.extend(dict(name=a.name, kind=a.sources[0]['kind'], original_order=i + 1,
                          allocation_order=None, reason='antigen_type_cap') for i, a in pending)
    audit['allocation'] = dict(policy=asdict(policy), unit='unique_antigen_sequence',
                               attribution='first_ranked_source', decisions=decisions)
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
