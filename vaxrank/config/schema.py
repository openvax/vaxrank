from typing import Any

import msgspec

from .defaults import DEFAULT_SCHEMA_VERSION


class VaxrankConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    schema_version: int = DEFAULT_SCHEMA_VERSION
    inputs: dict[str, Any] = msgspec.field(default_factory=dict)
    epitope_config: dict[str, Any] = msgspec.field(default_factory=dict)
    vaccine_config: dict[str, Any] = msgspec.field(default_factory=dict)
    self_proteome: dict[str, Any] = msgspec.field(default_factory=dict)
    epitopes: dict[str, Any] = msgspec.field(default_factory=dict)
    mutations: dict[str, Any] = msgspec.field(default_factory=dict)
    vaccine_peptides: dict[str, Any] = msgspec.field(default_factory=dict)
    compat: dict[str, Any] = msgspec.field(default_factory=dict)
