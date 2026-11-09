"""Source-independent RNA measurements, with explicit protein-window scope."""

import json
from dataclasses import dataclass, replace

from isovar import validate_rna_feature_support
from serializable import DataclassSerializable


RNA_BINDINGS = frozenset(
    "rna_" + name
    for unit in ("reads", "fragments")
    for name in (unit, "library_" + unit, "effective_length_" + unit,
                 unit + "_per_million", unit + "_per_effective_kb",
                 unit + "_per_effective_kb_per_million"))


@dataclass(frozen=True)
class RNAFeatureEvidence(DataclassSerializable):
    """One Isovar feature measurement and its declared applicability.

    ``event`` attaches local event evidence to this source sequence and retains
    it through slicing; it does not assert observation of downstream coding
    bases. ``window`` describes exactly ``protein_interval`` (zero-based,
    half-open amino acids). Its score binds only for that exact window.
    Multiple libraries/features remain separate; no pooling or best-sample
    policy is inferred. Normalized values are computed by the producer before
    any explicitly chosen longitudinal comparison.
    """

    measurement: dict
    scope: str
    protein_interval: tuple[int, int] | None = None

    def __post_init__(self):
        object.__setattr__(self, "measurement", json.loads(json.dumps(validate_rna_feature_support(self.measurement), allow_nan=False)))
        if self.scope not in ("event", "window"):
            raise ValueError("RNA evidence scope must be event or window")
        if self.scope == "event":
            if self.protein_interval is not None:
                raise ValueError("Event RNA support does not assert a protein-window span")
        else:
            interval = tuple(self.protein_interval or ())
            if (len(interval) != 2 or any(type(x) is not int for x in interval)
                    or not 0 <= interval[0] < interval[1]):
                raise ValueError("Window RNA support requires an amino-acid interval")
            object.__setattr__(self, "protein_interval", interval)

    def __hash__(self):
        return hash((self.measurement["measurement_id"], self.scope, self.protein_interval))


def validate_rna_support(records, protein_length):
    records = tuple(records)
    identities = []
    for record in records:
        if not isinstance(record, RNAFeatureEvidence):
            raise TypeError("RNA support requires RNAFeatureEvidence records")
        validate_rna_feature_support(record.measurement)
        if record.protein_interval is not None and record.protein_interval[1] > protein_length:
            raise ValueError("RNA window evidence extends beyond its protein sequence")
        identities.append((record.measurement["measurement_id"], record.scope, record.protein_interval))
    if len(identities) != len(set(identities)):
        raise ValueError("Duplicate RNA feature measurement attachment")
    return records


def slice_rna_support(records, start, end):
    """Retain local event evidence and rebase contained, exactly named windows."""
    result = []
    for record in records:
        if record.scope == "event":
            result.append(record)
        elif start <= record.protein_interval[0] and record.protein_interval[1] <= end:
            result.append(replace(record, protein_interval=tuple(x - start for x in record.protein_interval)))
    return tuple(result)


def applicable_rna_support(antigen):
    return tuple(record for record in getattr(antigen, "rna_support", ())
                 if record.scope == "event" or record.protein_interval == (0, len(antigen.amino_acids)))


def rna_support_bindings(antigen):
    """Missing or multiply scoped evidence remains unavailable, never zero."""
    records = applicable_rna_support(antigen)
    values = records[0].measurement["values"] if len(records) == 1 else {}
    return {name: values.get(name.removeprefix("rna_")) for name in RNA_BINDINGS}
