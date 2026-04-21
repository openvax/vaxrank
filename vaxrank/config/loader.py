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


def _read_single_config_file(config_path: str) -> dict[str, Any]:
    with open(config_path) as f:
        content = f.read()

    if not content.strip():
        return {}

    raw = msgspec.yaml.decode(content, type=dict)
    if raw is None:
        return {}
    if not isinstance(raw, dict):
        raise TypeError(
            f"Expected top-level YAML object in {config_path} to be a mapping"
        )
    return raw


def _deep_merge(base: dict[str, Any], overlay: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge ``overlay`` into ``base``. Overlay scalars replace
    base scalars; nested mappings merge key-by-key."""
    for key, overlay_value in overlay.items():
        base_value = base.get(key)
        if isinstance(base_value, dict) and isinstance(overlay_value, dict):
            _deep_merge(base_value, overlay_value)
        else:
            base[key] = overlay_value
    return base


def _read_config_files(config_paths: list[str] | str | None) -> dict[str, Any]:
    if not config_paths:
        return {}

    if isinstance(config_paths, str):
        config_paths = [config_paths]

    merged: dict[str, Any] = {}
    for path in config_paths:
        if not path:
            continue
        _deep_merge(merged, _read_single_config_file(path))
    return merged


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


def _promote_nested_manufacturability(raw: dict[str, Any]) -> None:
    """Emit a deprecation warning and lift `vaccine_peptides.manufacturability`
    to top-level `manufacturability`. Errors if both are set."""
    vaccine_section = raw.get("vaccine_peptides")
    if not isinstance(vaccine_section, dict):
        return
    nested = vaccine_section.get("manufacturability")
    if nested is None:
        return
    if raw.get("manufacturability") is not None:
        raise ValueError(
            "Cannot set both top-level 'manufacturability' and "
            "'vaccine_peptides.manufacturability'. The nested form is "
            "deprecated — use the top-level section."
        )
    warnings.warn(
        "Config key 'vaccine_peptides.manufacturability' is deprecated; "
        "move it to the top-level 'manufacturability' section. Support for "
        "the nested form will be removed in a future release.",
        DeprecationWarning,
        stacklevel=3,
    )
    raw["manufacturability"] = vaccine_section.pop("manufacturability")


def load_vaxrank_config(
    args=None,
    *,
    config_path: str | list[str] | None = None,
    set_overrides: list[str] | None = None,
    expr_overrides: list[str] | None = None,
    ordered_overrides: list[tuple[str, str]] | None = None,
) -> dict[str, Any]:
    if args is not None:
        config_path = getattr(args, "config", config_path)
        set_overrides = getattr(args, "config_set_overrides", set_overrides)
        expr_overrides = getattr(args, "config_expr_overrides", expr_overrides)
        ordered_overrides = getattr(args, "config_overrides", ordered_overrides)

    raw = _read_config_files(config_path)
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

    _promote_nested_manufacturability(raw)

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
    ("epitopes.filter_expr", "filter_expr"),
    ("epitopes.score_expr", "score_expr"),
]

# Declarative mapping: (dotted config path) -> (VaccineConfig kwarg name)
_VACCINE_CONFIG_MAPPING: list[tuple[str, str]] = [
    ("vaccine_peptides.preferred_length", "preferred_peptide_length"),
    ("vaccine_peptides.min_length", "min_peptide_length"),
    ("vaccine_peptides.max_length", "max_peptide_length"),
    ("vaccine_peptides.padding_around_mutation", "padding_around_mutation"),
    ("vaccine_peptides.per_mutation", "max_vaccine_peptides_per_variant"),
    ("vaccine_peptides.score_fraction_of_best", "score_fraction_of_best"),
    ("vaccine_peptides.combined_score_mode", "combined_score_mode"),
    ("vaccine_peptides.ranking_rules", "ranking_rules"),
    ("manufacturability.max_c_terminal_hydropathy", "max_c_terminal_hydropathy"),
    ("manufacturability.min_kmer_hydropathy", "min_kmer_hydropathy"),
    ("manufacturability.max_kmer_hydropathy_low_priority", "max_kmer_hydropathy_low_priority"),
    ("manufacturability.max_kmer_hydropathy_high_priority", "max_kmer_hydropathy_high_priority"),
    ("manufacturability.rules", "manufacturability_rules"),
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

    # VaccineConfig declares *_rules fields as Optional[tuple[str, ...]];
    # YAML gives us lists, so coerce.
    for rules_key in ("manufacturability_rules", "ranking_rules"):
        if rules_key in kwargs and kwargs[rules_key] is not None:
            kwargs[rules_key] = tuple(kwargs[rules_key])

    # When preferred_peptide_length is set but min/max are not,
    # default them to match preferred so the range is consistent.
    if "preferred_peptide_length" in kwargs:
        pref = kwargs["preferred_peptide_length"]
        kwargs.setdefault("min_peptide_length", pref)
        kwargs.setdefault("max_peptide_length", pref)

    # Handle num_mutant_epitopes_to_keep from two possible locations;
    # error if both are set.
    from_vaccine = _resolve_dotted(
        config, "vaccine_peptides.max_epitopes_per_candidate"
    )
    from_epitopes = _resolve_dotted(
        config, "epitopes.top_epitopes_per_candidate"
    )
    if from_vaccine is not _MISSING and from_epitopes is not _MISSING:
        raise ValueError(
            "Cannot set both vaccine_peptides.max_epitopes_per_candidate "
            "and epitopes.top_epitopes_per_candidate. Use one or the other."
        )
    if from_vaccine is not _MISSING:
        kwargs["num_mutant_epitopes_to_keep"] = from_vaccine
    elif from_epitopes is not _MISSING:
        kwargs["num_mutant_epitopes_to_keep"] = from_epitopes

    return kwargs
