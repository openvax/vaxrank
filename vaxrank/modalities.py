# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Single source of truth for the set of vaccine modalities vaxrank
knows how to design (today: ``peptide`` and ``mrna``; future:
``dna``, …).
"""

from __future__ import annotations

_MODALITIES = ('peptide', 'mrna')


def known_modalities():
    """Names of registered modalities, in registration order."""
    return _MODALITIES
