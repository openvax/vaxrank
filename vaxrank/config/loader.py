from __future__ import annotations

import warnings
from typing import Any

import msgspec

from .schema import VaxrankConfigSchema
from ..modalities import known_modalities

# Top-level YAML keys that have been renamed across versions. Each
# entry deprecates the OLD key and routes its value to the NEW key,
# erroring if both are set in the same config.
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


# Antigen-design knobs that, post-2.19, live under ``vaccine_peptides:``
# but pre-2.19 lived at the ``vaccine_constructs:`` top level. Used by
# the back-compat migration AND by ``extract_construct_kwargs`` to
# thread these into the per-modality construct kwargs (so the
# construct assemblers receive them at the same place they always have).
_ANTIGEN_DESIGN_KNOBS = (
    "antigen_content",
    "epitopes_per_antigen",
    "candidates_per_slot",
    "min_antigen_length_aa",
    "max_antigen_length_aa",
)


def _migrate_legacy_layout(raw: dict[str, Any]) -> None:
    """Translate the pre-2.19 YAML layout into the post-2.19 shape
    in-place, with one ``DeprecationWarning`` per migrated section.
    Errors when both layouts are mixed for the same logical section.

    Migrations:

    * ``vaccine_peptides.manufacturability`` →
      ``peptide.manufacturability`` (was previously renamed to
      top-level ``manufacturability``; pre-2.19 it had also lived
      under vaccine_peptides).
    * top-level ``manufacturability`` → ``peptide.manufacturability``.
    * ``vaccine_constructs.<modality>`` → top-level ``<modality>:``.
    * ``vaccine_constructs.<antigen-design knob>`` →
      ``vaccine_peptides.<knob>`` (e.g. ``antigen_content``,
      ``epitopes_per_antigen``, ``candidates_per_slot``,
      ``min_antigen_length_aa``, ``max_antigen_length_aa``).
    """
    # 1. ``vaccine_peptides.manufacturability`` → top-level for the
    # next stage to migrate further.
    vp = raw.get("vaccine_peptides")
    if isinstance(vp, dict) and vp.get("manufacturability") is not None:
        if raw.get("manufacturability") is not None:
            raise ValueError(
                "Cannot set both 'vaccine_peptides.manufacturability' and "
                "top-level 'manufacturability'. Both are deprecated — pick "
                "one and prefer 'peptide.manufacturability'.")
        warnings.warn(
            "Config key 'vaccine_peptides.manufacturability' is deprecated; "
            "move it to 'peptide.manufacturability'. Auto-migrating.",
            DeprecationWarning,
            stacklevel=4,
        )
        raw["manufacturability"] = vp.pop("manufacturability")

    # 2. top-level ``manufacturability`` → ``peptide.manufacturability``.
    if raw.get("manufacturability") is not None:
        peptide_section = raw.setdefault("peptide", {})
        if not isinstance(peptide_section, dict):
            raise TypeError(
                "Config key 'peptide' must be a mapping; got %r"
                % type(peptide_section).__name__)
        if peptide_section.get("manufacturability") is not None:
            raise ValueError(
                "Cannot set both top-level 'manufacturability' and "
                "'peptide.manufacturability'. Top-level form is "
                "deprecated — keep just 'peptide.manufacturability'.")
        warnings.warn(
            "Config key 'manufacturability' (top-level) is deprecated; "
            "move it under 'peptide.manufacturability'. Manufacturability "
            "rules describe peptide-synthesis difficulty and don't apply "
            "to mRNA constructs. Auto-migrating.",
            DeprecationWarning,
            stacklevel=4,
        )
        peptide_section["manufacturability"] = raw.pop("manufacturability")

    # 3. ``vaccine_constructs:`` wrapper → top-level modality
    # sections + antigen-design knobs into ``vaccine_peptides:``.
    vc = raw.pop("vaccine_constructs", None)
    if isinstance(vc, dict) and vc:
        warnings.warn(
            "Config key 'vaccine_constructs' is deprecated. Move "
            "modality-specific entries to top-level 'peptide:' / "
            "'mrna:' sections; move antigen-design knobs "
            "(antigen_content, epitopes_per_antigen, "
            "candidates_per_slot, min/max_antigen_length_aa) under "
            "'vaccine_peptides:'. Auto-migrating.",
            DeprecationWarning,
            stacklevel=4,
        )
        modalities = set(known_modalities())
        for key, value in vc.items():
            if key in modalities:
                target = raw.setdefault(key, {})
                if not isinstance(target, dict):
                    raise TypeError(
                        "Config key %r must be a mapping; got %r"
                        % (key, type(target).__name__))
                if not isinstance(value, dict):
                    raise TypeError(
                        "vaccine_constructs.%s must be a mapping; got %r"
                        % (key, type(value).__name__))
                # Modality value at the new top level wins on
                # collision — the user's explicit new-layout entry
                # is more specific than the legacy wrapper.
                for k, v in value.items():
                    target.setdefault(k, v)
            elif key in _ANTIGEN_DESIGN_KNOBS:
                vp_section = raw.setdefault("vaccine_peptides", {})
                if not isinstance(vp_section, dict):
                    raise TypeError(
                        "Config key 'vaccine_peptides' must be a "
                        "mapping; got %r" % type(vp_section).__name__)
                vp_section.setdefault(key, value)
            else:
                # Unknown key under the legacy wrapper — let msgspec
                # surface it as a forbid_unknown_fields error rather
                # than guess where it goes.
                raise ValueError(
                    "Unknown legacy 'vaccine_constructs.%s'; remove "
                    "it or move under one of the known modality "
                    "sections (%s)." % (
                        key, ", ".join(sorted(modalities))))


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

    _migrate_legacy_layout(raw)

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

# Declarative mapping: (dotted config path) → (EpitopeConfig kwarg name)
_EPITOPE_CONFIG_MAPPING: list[tuple[str, str]] = [
    ("epitopes.min_score", "min_epitope_score"),
    ("epitopes.scoring_mode", "scoring_mode"),
    ("epitopes.logistic_midpoint", "logistic_epitope_score_midpoint"),
    ("epitopes.logistic_width", "logistic_epitope_score_width"),
    ("epitopes.affinity_cutoff", "binding_affinity_cutoff"),
    ("epitopes.percentile_rank_cutoff", "percentile_rank_cutoff"),
    ("epitopes.filter_expr", "filter_expr"),
    ("epitopes.score_expr", "score_expr"),
    ("epitopes.default_methods", "default_methods"),
]

# Declarative mapping: (dotted config path) → (VaccineConfig kwarg name)
#
# Post-2.20 manufacturability lives on its own struct
# (:class:`vaxrank.manufacturability_config.ManufacturabilityConfig`)
# threaded into the ranker as a separate parameter — see
# ``_MANUFACTURABILITY_CONFIG_MAPPING`` below. ``VaccineConfig``
# carries only modality-agnostic ranking knobs.
_VACCINE_CONFIG_MAPPING: list[tuple[str, str]] = [
    ("vaccine_peptides.preferred_length", "preferred_peptide_length"),
    ("vaccine_peptides.min_length", "min_peptide_length"),
    ("vaccine_peptides.max_length", "max_peptide_length"),
    ("vaccine_peptides.padding_around_mutation", "padding_around_mutation"),
    ("vaccine_peptides.per_mutation", "max_vaccine_peptides_per_variant"),
    ("vaccine_peptides.score_fraction_of_best", "score_fraction_of_best"),
    ("vaccine_peptides.combined_score_mode", "combined_score_mode"),
    ("vaccine_peptides.combined_score_expr", "combined_score_expr"),
    ("vaccine_peptides.ranking_rules", "ranking_rules"),
    ("vaccine_peptides.require_target_epitopes_in_variant",
     "require_target_epitopes_in_variant"),
]


# Declarative mapping: (dotted config path) → (ManufacturabilityConfig kwarg name)
_MANUFACTURABILITY_CONFIG_MAPPING: list[tuple[str, str]] = [
    ("peptide.manufacturability.max_c_terminal_hydropathy",
     "max_c_terminal_hydropathy"),
    ("peptide.manufacturability.min_kmer_hydropathy",
     "min_kmer_hydropathy"),
    ("peptide.manufacturability.max_kmer_hydropathy_low_priority",
     "max_kmer_hydropathy_low_priority"),
    ("peptide.manufacturability.max_kmer_hydropathy_high_priority",
     "max_kmer_hydropathy_high_priority"),
    ("peptide.manufacturability.rules", "rules"),
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


def extract_construct_kwargs(
    config: dict[str, Any],
    modality: str,
) -> dict[str, Any]:
    """Resolve construct-assembly kwargs for ``modality`` ('peptide',
    'mrna', or any other registered modality) by reading the
    top-level modality section, then layering on the antigen-design
    knobs from ``vaccine_peptides:`` (so each modality picks them up
    unless it overrides them in its own section).

    Empty dict when no inputs apply — same shape as before so code
    paths that don't use the per-modality section behave identically.
    """
    out: dict[str, Any] = {}
    # Antigen-design knobs from vaccine_peptides: thread into every
    # modality. Modality-specific overrides in ``<modality>:`` win.
    vp = _resolve_dotted(config, 'vaccine_peptides')
    if isinstance(vp, dict):
        for k in _ANTIGEN_DESIGN_KNOBS:
            v = vp.get(k)
            if v is not None:
                out[k] = v
    section = _resolve_dotted(config, modality)
    if isinstance(section, dict):
        for k, v in section.items():
            # ``manufacturability`` is a sub-section in the YAML
            # schema; it never goes into the construct-config kwargs
            # (those are flat). The per-modality manufacturability
            # config rides through ``_VACCINE_CONFIG_MAPPING`` for
            # now (will move onto PeptideConstructConfig once the
            # in-code split lands).
            if k == 'manufacturability':
                continue
            if v is None:
                continue
            out[k] = v
    return out


def extract_manufacturability_config_kwargs(
    config: dict[str, Any],
) -> dict[str, Any]:
    """Extract ManufacturabilityConfig kwargs from the YAML
    ``peptide.manufacturability`` sub-section. Returns an empty dict
    when the sub-section is absent — caller should use
    ``ManufacturabilityConfig()`` defaults in that case."""
    kwargs = _extract_via_mapping(
        config, _MANUFACTURABILITY_CONFIG_MAPPING)
    if "rules" in kwargs and kwargs["rules"] is not None:
        kwargs["rules"] = tuple(kwargs["rules"])
    return kwargs


def extract_vaccine_config_kwargs(config: dict[str, Any]) -> dict[str, Any]:
    kwargs = _extract_via_mapping(config, _VACCINE_CONFIG_MAPPING)

    # VaccineConfig declares ranking_rules as Optional[tuple[str, ...]];
    # YAML gives us a list, so coerce. Manufacturability rules live
    # on ManufacturabilityConfig now (post-2.20) and are coerced in
    # ``extract_manufacturability_config_kwargs``.
    if "ranking_rules" in kwargs and kwargs["ranking_rules"] is not None:
        kwargs["ranking_rules"] = tuple(kwargs["ranking_rules"])

    # When preferred_peptide_length is set but min/max are not,
    # default them to match preferred so the range is consistent.
    if "preferred_peptide_length" in kwargs:
        pref = kwargs["preferred_peptide_length"]
        kwargs.setdefault("min_peptide_length", pref)
        kwargs.setdefault("max_peptide_length", pref)

    # Handle num_target_epitopes_to_keep from two possible locations;
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
        kwargs["num_target_epitopes_to_keep"] = from_vaccine
    elif from_epitopes is not _MISSING:
        kwargs["num_target_epitopes_to_keep"] = from_epitopes

    return kwargs
