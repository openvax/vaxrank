"""Shared Isovar API, source-file and YAML controls.

Reconstruction keys keep their producer names. Only documented, JSON-valued
parameters of the selected public Isovar API are accepted; resolved defaults
are captured with the evidence instead of being copied into Vaxrank.
"""

from inspect import Parameter, signature
from typing import Any

import msgspec


class IsovarAdmission(msgspec.Struct, frozen=True, kw_only=True, forbid_unknown_fields=True):
    """Evidence gates, independent of tumor specificity and epitope ranking.

    Partial and ambiguous proteins are retained even when held out. Fragment
    counts retain their producer scope; they are never interpreted as TPM,
    independent molecules or evidence of protein translation.
    """

    allow_partial: bool = False
    allow_ambiguous: bool = False
    min_fragments: int = 2

    def __post_init__(self):
        if type(self.allow_partial) is not bool or type(self.allow_ambiguous) is not bool:
            raise ValueError('Isovar admission switches must be booleans')
        if type(self.min_fragments) is not int or self.min_fragments < 1:
            raise ValueError('Isovar min_fragments must be a positive integer')


class IsovarSourceOptions(msgspec.Struct, frozen=True, kw_only=True, forbid_unknown_fields=True):
    """One source family's controls, with strict producer-specific validation.

    ``reconstruction`` supplies the selected reconstruction function's scalar
    controls. ``reads`` configures ReadCollector for variant/SV sources;
    ``protein`` configures ProteinSequenceCreator for variant ORFs. Inputs
    such as BAMs, variants and reference models belong in ``inputs`` instead.
    """

    admission: IsovarAdmission = IsovarAdmission()
    reconstruction: dict[str, Any] = msgspec.field(default_factory=dict)
    reads: dict[str, Any] = msgspec.field(default_factory=dict)
    protein: dict[str, Any] = msgspec.field(default_factory=dict)


FAMILIES = ('variant_orfs', 'sv', 'fusion', 'viral', 'isoforms', 'protein_hypotheses')


def source_options(value=None):
    """Use the same strict contract for Python, JSON/YAML and CLI overrides."""
    if isinstance(value, IsovarSourceOptions):
        # Struct constructors do not perform type validation.
        value = msgspec.to_builtins(value)
    return msgspec.convert({} if value is None else value, type=IsovarSourceOptions, strict=True)


def merge_options(*values):
    """Merge explicit mappings from source < YAML/--set < manifest entry."""
    merged = {}
    for value in values:
        if isinstance(value, IsovarSourceOptions):
            value = msgspec.to_builtins(value)
        if value is None:
            continue
        if not isinstance(value, dict):
            raise ValueError('Isovar options must be an object')
        for key, section in value.items():
            if not isinstance(section, dict):
                raise ValueError('Isovar option sections must be objects')
            merged.setdefault(key, {}).update(section)
    return source_options(merged)


def producer_options(function, supplied, *, exclude=(), defaults=None, nullable=()):
    """Validate JSON-valued controls and freeze the actual producer defaults.

    Required arguments and explicitly excluded data objects cannot be smuggled
    through options. Callback/profile objects remain available in Isovar's own
    Python API; Vaxrank source files never import or execute caller code.
    """
    resolved = {
        name: parameter.default
        for name, parameter in signature(function).parameters.items()
        if name not in exclude and name != 'self'
        and parameter.default is not Parameter.empty
        and isinstance(parameter.default, (type(None), str, int, float, bool, tuple, list))
    }
    unknown = set(supplied) - set(resolved)
    if unknown:
        raise ValueError(f'{function.__name__}: unknown or inapplicable controls: {sorted(unknown)}')
    resolved.update(defaults or {})
    for key, value in supplied.items():
        if value is None and key in nullable:
            resolved[key] = None
            continue
        default = resolved[key]
        if default is not None:
            expected = type(default)
            valid = (isinstance(value, (list, tuple)) if expected in (list, tuple)
                     else type(value) in (int, float) if expected is float
                     else type(value) is expected)
            if not valid:
                raise ValueError(f'{function.__name__}.{key}: expected {expected.__name__}')
        if value is None and default is not None:
            raise ValueError(f'{function.__name__}.{key}: null is not supported')
        resolved[key] = value
    # Normalize tuples and reject non-JSON values/nonfinite floats.
    import json
    return json.loads(json.dumps(resolved, allow_nan=False))
