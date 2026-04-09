from __future__ import annotations

import warnings
from typing import Any

import msgspec

from .schema import VaxrankConfigSchema

_LEGACY_KEY_MAP = {
    "epitope_config": "epitopes",
    "vaccine_config": "vaccine_peptides",
}


def _drop_none(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: _drop_none(nested_value)
            for key, nested_value in value.items()
            if nested_value is not None
        }
    if isinstance(value, list):
        return [_drop_none(item) for item in value]
    return value


def _schema_to_dict(schema: VaxrankConfigSchema) -> dict[str, Any]:
    return _drop_none(msgspec.to_builtins(schema))


def _read_config_file(config_path: str | None) -> dict[str, Any]:
    if not config_path:
        return {}

    with open(config_path) as f:
        content = f.read()

    if not content.strip():
        return {}

    raw = msgspec.yaml.decode(content, type=dict)
    if raw is None:
        return {}
    if not isinstance(raw, dict):
        raise TypeError("Expected top-level YAML object to be a mapping")
    return raw


def _split_override(override: str) -> tuple[str, str]:
    if "=" not in override:
        raise ValueError(
            f"Invalid config override '{override}'. Expected dotted.path=value."
        )
    path, value = override.split("=", 1)
    path = path.strip()
    if not path:
        raise ValueError(
            f"Invalid config override '{override}'. Override path cannot be empty."
        )
    return path, value


def _decode_override_value(value: str) -> Any:
    value = value.strip()
    if value == "":
        return ""
    try:
        return msgspec.yaml.decode(value, type=Any)
    except msgspec.DecodeError:
        return value


def _set_nested_value(target: dict[str, Any], path: str, value: Any) -> None:
    keys = path.split(".")
    current = target
    for key in keys[:-1]:
        existing = current.get(key)
        if existing is None:
            existing = {}
            current[key] = existing
        elif not isinstance(existing, dict):
            raise TypeError(
                f"Cannot assign config override '{path}': '{key}' is not a mapping."
            )
        current = existing
    current[keys[-1]] = value


def load_vaxrank_config(
    args=None,
    *,
    config_path: str | None = None,
    set_overrides: list[str] | None = None,
    expr_overrides: list[str] | None = None,
    ordered_overrides: list[tuple[str, str]] | None = None,
) -> dict[str, Any]:
    if args is not None:
        config_path = getattr(args, "config", config_path)
        set_overrides = getattr(args, "config_set_overrides", set_overrides)
        expr_overrides = getattr(args, "config_expr_overrides", expr_overrides)
        ordered_overrides = getattr(args, "config_overrides", ordered_overrides)

    raw = _read_config_file(config_path)
    if "schema_version" in raw:
        raise ValueError(
            "schema_version is not supported yet. Use the current unversioned config schema."
        )

    for legacy_key, new_key in _LEGACY_KEY_MAP.items():
        if legacy_key in raw:
            warnings.warn(
                f"Config key '{legacy_key}' is deprecated, use '{new_key}' instead. "
                f"Support for '{legacy_key}' will be removed in a future release.",
                DeprecationWarning,
                stacklevel=2,
            )
            if new_key in raw:
                raise ValueError(
                    f"Cannot use both '{legacy_key}' and '{new_key}' in the same config. "
                    f"Rename '{legacy_key}' to '{new_key}'."
                )
            raw[new_key] = raw.pop(legacy_key)

    if ordered_overrides is not None:
        for kind, override in ordered_overrides:
            path, value_text = _split_override(override)
            if kind == "set":
                value = _decode_override_value(value_text)
            elif kind == "expr":
                value = value_text.strip()
            else:
                raise ValueError(f"Unknown config override kind '{kind}'")
            _set_nested_value(raw, path, value)
    else:
        for override in set_overrides or []:
            path, value_text = _split_override(override)
            _set_nested_value(raw, path, _decode_override_value(value_text))

        for override in expr_overrides or []:
            path, value_text = _split_override(override)
            _set_nested_value(raw, path, value_text.strip())

    schema = msgspec.convert(raw, VaxrankConfigSchema)
    return _schema_to_dict(schema)


_MISSING = object()

# Declarative mapping: (dotted config path) -> (EpitopeConfig kwarg name)
_EPITOPE_CONFIG_MAPPING: list[tuple[str, str]] = [
    ("epitopes.min_score", "min_epitope_score"),
    ("epitopes.scoring_mode", "scoring_mode"),
    ("epitopes.logistic_midpoint", "logistic_epitope_score_midpoint"),
    ("epitopes.logistic_width", "logistic_epitope_score_width"),
    ("epitopes.affinity_cutoff", "binding_affinity_cutoff"),
    ("epitopes.percentile_rank_cutoff", "percentile_rank_cutoff"),
]

# Declarative mapping: (dotted config path) -> (VaccineConfig kwarg name)
_VACCINE_CONFIG_MAPPING: list[tuple[str, str]] = [
    ("vaccine_peptides.generation.padding_around_mutation", "padding_around_mutation"),
    ("vaccine_peptides.keep.per_mutation", "max_vaccine_peptides_per_variant"),
    ("vaccine_peptides.keep.score_fraction_of_best", "score_fraction_of_best"),
    ("vaccine_peptides.manufacturability.max_c_terminal_hydropathy", "max_c_terminal_hydropathy"),
    ("vaccine_peptides.manufacturability.min_kmer_hydropathy", "min_kmer_hydropathy"),
    ("vaccine_peptides.manufacturability.max_kmer_hydropathy_low_priority", "max_kmer_hydropathy_low_priority"),
    ("vaccine_peptides.manufacturability.max_kmer_hydropathy_high_priority", "max_kmer_hydropathy_high_priority"),
]


def _resolve_dotted(config: dict[str, Any], dotted_path: str) -> Any:
    """Resolve a dotted path like 'a.b.c' to a value in a nested dict."""
    keys = dotted_path.split(".")
    current = config
    for key in keys:
        if not isinstance(current, dict) or key not in current:
            return _MISSING
        current = current[key]
    return current


def _extract_via_mapping(
    config: dict[str, Any],
    mapping: list[tuple[str, str]],
) -> dict[str, Any]:
    """Extract kwargs from config using a declarative mapping."""
    kwargs: dict[str, Any] = {}
    for config_path, kwarg_name in mapping:
        value = _resolve_dotted(config, config_path)
        if value is not _MISSING:
            kwargs[kwarg_name] = value
    return kwargs


def extract_epitope_config_kwargs(config: dict[str, Any]) -> dict[str, Any]:
    return _extract_via_mapping(config, _EPITOPE_CONFIG_MAPPING)


def extract_vaccine_config_kwargs(config: dict[str, Any]) -> dict[str, Any]:
    kwargs = _extract_via_mapping(config, _VACCINE_CONFIG_MAPPING)

    # Handle lengths specially: list -> single value
    lengths = _resolve_dotted(config, "vaccine_peptides.generation.lengths")
    if lengths is not _MISSING:
        if isinstance(lengths, list):
            if len(lengths) != 1:
                raise ValueError(
                    "Current Vaxrank implementation supports a single vaccine peptide "
                    "length. Provide exactly one value in "
                    "vaccine_peptides.generation.lengths."
                )
            kwargs["vaccine_peptide_length"] = lengths[0]
        else:
            kwargs["vaccine_peptide_length"] = lengths

    # Handle num_mutant_epitopes_to_keep from two possible locations;
    # error if both are set.
    from_vaccine = _resolve_dotted(
        config, "vaccine_peptides.keep.max_epitopes_per_candidate"
    )
    from_epitopes = _resolve_dotted(
        config, "epitopes.top_epitopes_per_candidate"
    )
    if from_vaccine is not _MISSING and from_epitopes is not _MISSING:
        raise ValueError(
            "Cannot set both vaccine_peptides.keep.max_epitopes_per_candidate "
            "and epitopes.top_epitopes_per_candidate. Use one or the other."
        )
    if from_vaccine is not _MISSING:
        kwargs["num_mutant_epitopes_to_keep"] = from_vaccine
    elif from_epitopes is not _MISSING:
        kwargs["num_mutant_epitopes_to_keep"] = from_epitopes

    return kwargs
