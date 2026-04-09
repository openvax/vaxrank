from __future__ import annotations

from typing import Any

import msgspec

from .schema import VaxrankConfigSchema


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


def extract_epitope_config_kwargs(config: dict[str, Any]) -> dict[str, Any]:
    kwargs = dict(config.get("epitope_config") or {})

    epitopes = config.get("epitopes") or {}
    filters = epitopes.get("filters") or {}
    scoring = epitopes.get("scoring") or {}
    derived_fields = scoring.get("derived_fields") or {}
    affinity_score = derived_fields.get("affinity_score") or {}
    percentile_score = derived_fields.get("percentile_score") or {}

    if "min_score" in filters:
        kwargs["min_epitope_score"] = filters["min_score"]
    if "mode" in scoring:
        kwargs["scoring_mode"] = scoring["mode"]
    if "midpoint" in affinity_score:
        kwargs["logistic_epitope_score_midpoint"] = affinity_score["midpoint"]
    if "width" in affinity_score:
        kwargs["logistic_epitope_score_width"] = affinity_score["width"]
    if "cutoff" in affinity_score:
        kwargs["binding_affinity_cutoff"] = affinity_score["cutoff"]
    if "worst" in percentile_score:
        kwargs["percentile_rank_cutoff"] = percentile_score["worst"]

    return kwargs


def extract_vaccine_config_kwargs(config: dict[str, Any]) -> dict[str, Any]:
    kwargs = dict(config.get("vaccine_config") or {})

    vaccine_peptides = config.get("vaccine_peptides") or {}
    generation = vaccine_peptides.get("generation") or {}
    keep = vaccine_peptides.get("keep") or {}
    epitope_keep = (config.get("epitopes") or {}).get("keep") or {}

    lengths = generation.get("lengths")
    if lengths is not None:
        if isinstance(lengths, list):
            if len(lengths) > 1:
                raise ValueError(
                    "Current Vaxrank implementation supports a single vaccine peptide "
                    "length. Use one value in vaccine_peptides.generation.lengths."
                )
            if lengths:
                kwargs["vaccine_peptide_length"] = lengths[0]
        else:
            kwargs["vaccine_peptide_length"] = lengths

    if "padding_around_mutation" in generation:
        kwargs["padding_around_mutation"] = generation["padding_around_mutation"]
    if "per_mutation" in keep:
        kwargs["max_vaccine_peptides_per_variant"] = keep["per_mutation"]
    if "score_fraction_of_best" in keep:
        kwargs["score_fraction_of_best"] = keep["score_fraction_of_best"]
    if "max_epitopes_per_candidate" in keep:
        kwargs["num_mutant_epitopes_to_keep"] = keep["max_epitopes_per_candidate"]
    if "top_n_per_candidate" in epitope_keep:
        kwargs["num_mutant_epitopes_to_keep"] = epitope_keep["top_n_per_candidate"]

    manufacturability = vaccine_peptides.get("manufacturability") or {}
    if "max_c_terminal_hydropathy" in manufacturability:
        kwargs["max_c_terminal_hydropathy"] = manufacturability["max_c_terminal_hydropathy"]
    if "min_kmer_hydropathy" in manufacturability:
        kwargs["min_kmer_hydropathy"] = manufacturability["min_kmer_hydropathy"]
    if "max_kmer_hydropathy_low_priority" in manufacturability:
        kwargs["max_kmer_hydropathy_low_priority"] = manufacturability["max_kmer_hydropathy_low_priority"]
    if "max_kmer_hydropathy_high_priority" in manufacturability:
        kwargs["max_kmer_hydropathy_high_priority"] = manufacturability["max_kmer_hydropathy_high_priority"]

    return kwargs
