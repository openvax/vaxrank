# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Single source of truth for the set of vaccine modalities vaxrank
knows how to design (today: ``peptide`` and ``mrna``; future:
``dna``, …).

Every place that needs the list of modalities — CLI ``--vaccine-type``
choices, the YAML config schema, the dispatch table inside
``_emit_outputs``, the per-modality config loader — reads from here.
Adding a new modality is one entry: register the construct-config
class. The schema, CLI, and dispatch all pick up the new modality
automatically.
"""

from __future__ import annotations

# Filled lazily so the ``construct_config`` classes (peptide /
# mrna), which pull in vaccine_library / mhctools, aren't imported
# at every ``from vaxrank.modalities import known_modalities``.
_REGISTRY = None


def _build_registry():
    from .peptide import PeptideConstructConfig
    from .mrna import RNAConstructConfig

    return {
        'peptide': {
            'construct_config': PeptideConstructConfig,
        },
        'mrna': {
            'construct_config': RNAConstructConfig,
        },
    }


def _registry():
    global _REGISTRY
    if _REGISTRY is None:
        _REGISTRY = _build_registry()
    return _REGISTRY


def known_modalities():
    """Names of registered modalities, in registration order."""
    return tuple(_registry())


def construct_config_class(modality):
    """The ``*ConstructConfig`` class for the named modality."""
    return _registry()[modality]['construct_config']
