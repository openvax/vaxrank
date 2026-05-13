# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Pure-Python reference implementation of the legacy per-prediction
epitope score, used only by tests as an independent oracle to verify
the DSL-evaluated score in :mod:`vaxrank.epitope_dsl`. Not imported
by any production code — production goes through the DSL exclusively.
"""

import numpy as np

DEFAULT_LOGISTIC_EPITOPE_SCORE_MIDPOINT = 350.0
DEFAULT_LOGISTIC_EPITOPE_SCORE_WIDTH = 150.0
DEFAULT_BINDING_AFFINITY_CUTOFF = 5000.0
DEFAULT_PERCENTILE_RANK_CUTOFF = 10.0


def legacy_score_one(ic50, percentile_rank, *,
                     midpoint=DEFAULT_LOGISTIC_EPITOPE_SCORE_MIDPOINT,
                     width=DEFAULT_LOGISTIC_EPITOPE_SCORE_WIDTH,
                     ic50_cutoff=DEFAULT_BINDING_AFFINITY_CUTOFF,
                     scoring_mode="affinity",
                     percentile_rank_cutoff=DEFAULT_PERCENTILE_RANK_CUTOFF):
    """Reference scorer. Mirrors what the DSL's default ``score_expr``
    produces row-by-row; parity is asserted in
    ``tests/test_epitope_dsl.py::test_default_score_matches_legacy_*``.
    """
    if scoring_mode == "percentile_rank":
        if percentile_rank is None:
            return 0.0
        rank = float(percentile_rank)
        if rank >= percentile_rank_cutoff:
            return 0.0
        return max(0.0, 1.0 - rank / percentile_rank_cutoff)

    if ic50 >= ic50_cutoff:
        return 0.0
    rescaled = (float(ic50) - midpoint) / width
    logistic = 1.0 / (1.0 + np.exp(rescaled))
    normalizer = 1.0 / (1.0 + np.exp(-midpoint / width))
    return logistic / normalizer
