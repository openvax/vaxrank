"""Resolve RNA reconstruction settings without changing peptide/DNA padding."""

from isovar.default_parameters import (
    MIN_PROTEIN_SEQUENCE_SUPPORT_FRACTION,
    MIN_VARIANT_SEQUENCE_COVERAGE,
    PROTEIN_SEQUENCE_PREFERENCE,
)


def resolve_isovar_args(args, vaccine_config, merged_config):
    """Apply explicit CLI > YAML > defaults before creating Isovar's creator.

    An explicit protein length takes precedence over legacy padding. With
    neither supplied, leave the length unset so Isovar derives its adaptive
    request from the resolved peptide size. Keep the effective settings on
    the namespace for saved run provenance; no RNA policy is applied to
    external predictions, cached reports, or the DNA-only fallback.
    """
    config = merged_config.get("isovar", {})
    defaults = {
        "protein_sequence_length": None,
        "protein_context_peptide_length": vaccine_config.preferred_peptide_length,
        "protein_sequence_preference": PROTEIN_SEQUENCE_PREFERENCE,
        "min_protein_sequence_support_fraction": MIN_PROTEIN_SEQUENCE_SUPPORT_FRACTION,
        "min_variant_sequence_coverage": MIN_VARIANT_SEQUENCE_COVERAGE,
    }
    for name, default in defaults.items():
        value = getattr(args, name, None)
        if value is None:
            value = config.get(name)
        if value is None:
            value = default
        setattr(args, name, value)

    padding = getattr(args, "padding_around_mutation", None)
    if padding is None:
        padding = merged_config.get("vaccine_peptides", {}).get("padding_around_mutation")
    if args.protein_sequence_length is None and padding is not None:
        args.protein_sequence_length = vaccine_config.preferred_peptide_length + 2 * padding
