# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np
from serializable import DataclassSerializable

from .config.defaults import (
    DEFAULT_BINDING_AFFINITY_CUTOFF,
    DEFAULT_LOGISTIC_EPITOPE_SCORE_MIDPOINT,
    DEFAULT_LOGISTIC_EPITOPE_SCORE_WIDTH,
    DEFAULT_PERCENTILE_RANK_CUTOFF,
)
from .manufacturability import (
    ManufacturabilityScores,
    compute_manufacturability_tuple,
)
from .ranking import compute_ranking_tuple
from .vaccine_config import COMBINED_SCORE_MODES, DEFAULT_COMBINED_SCORE_MODE


def _legacy_score_one(ic50, percentile_rank, *,
                      midpoint=DEFAULT_LOGISTIC_EPITOPE_SCORE_MIDPOINT,
                      width=DEFAULT_LOGISTIC_EPITOPE_SCORE_WIDTH,
                      ic50_cutoff=DEFAULT_BINDING_AFFINITY_CUTOFF,
                      scoring_mode="affinity",
                      percentile_rank_cutoff=DEFAULT_PERCENTILE_RANK_CUTOFF):
    """Per-prediction logistic score used by ``VaccinePeptide`` to sum
    scores across the leaf ``mhctools.Prediction`` records inside
    each ``CandidateEpitope``.

    Kept as a free function so the math stays in one place. The
    topiary DSL's default ``score_expr`` produces byte-identical
    output (pinned by ``test_default_score_matches_legacy_*``); a
    follow-up will route ``VaccinePeptide``'s score aggregation
    through ``score_predictions`` directly and this helper goes away.
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


@dataclass
class VaccinePeptide(DataclassSerializable):
    """
    VaccinePeptide combines the sequence information of MutantProteinFragment
    with MHC binding predictions for subsequences of the protein fragment.

    The resulting lists of mutant and wildtype ``CandidateEpitope`` objects
    are sorted by best mutant affinity.

    Parameters
    ----------
    mutant_protein_fragment : MutantProteinFragment

    epitopes : list of CandidateEpitope

    num_mutant_epitopes_to_keep : int or None
        If None or 0 then keep all mutant epitopes.

    epitope_score_params : dict or None
        Parameters for the per-prediction score (midpoint, width,
        ic50_cutoff, scoring_mode, percentile_rank_cutoff). Matches
        the topiary DSL default-score knobs byte-for-byte.

    sort_epitopes_by : str
        Field of the mutant's best affinity ``mhctools.Prediction``
        used for sorting epitopes overlapping mutation in ascending
        order. Can be either 'ic50' (Prediction.value) or
        'percentile_rank'.

    manufacturability_thresholds : dict or None
        Hydropathy thresholds for peptide synthesis difficulty scoring.
        Keys match the parameters of
        ``peptide_synthesis_difficulty_score_tuple``. When None the
        compiled-in defaults from ``config.defaults`` are used.

    manufacturability_rules : sequence of str or None
        Ordered list of manufacturability rule names to drive the
        lexicographic sort tuple. Names must appear in
        ``MANUFACTURABILITY_RULE_REGISTRY``. When None the legacy
        10-rule default order is used (byte-identical with prior
        releases).

    combined_score_mode : str or None
        How to fold expression into the final ``combined_score``. One
        of ``sqrt_reads_times_epitope`` (default, legacy),
        ``reads_times_epitope``, or ``epitope_only``. Does not affect
        the ``expression_score`` property, which always reports
        ``sqrt(n_alt_reads)``.

    ranking_rules : sequence of str or None
        Ordered list of ranking rule names driving the lexicographic
        sort tuple on ``lexicographic_sort_key``. Names must appear in
        ``RANKING_RULE_REGISTRY``; the special sentinel
        ``"manufacturability"`` expands to the peptide's manufacturability
        sort tuple. When None the 7-entry default order is used
        (byte-identical with prior releases).
    """

    mutant_protein_fragment: Any
    epitopes: list = field(default_factory=list)
    num_mutant_epitopes_to_keep: Optional[int] = None
    epitope_score_params: Optional[dict] = None
    sort_epitopes_by: str = "ic50"
    manufacturability_thresholds: Optional[dict] = None
    manufacturability_rules: Optional[tuple] = None
    combined_score_mode: Optional[str] = None
    # When set, supersedes ``combined_score_mode``. A string DSL
    # expression evaluated per-VP at scoring time — see
    # :mod:`vaxrank.combined_score_dsl` for the grammar and
    # bindings. ``None`` keeps the enum-mode behavior.
    combined_score_expr: Optional[str] = None
    ranking_rules: Optional[tuple] = None

    # Derived attributes computed in __post_init__ — not part of the
    # serialized form. `init=False` keeps them out of the generated
    # __init__ signature; default=None gives them a sane pre-computed
    # value for any pathway that inspects fields before __post_init__.
    mutant_epitopes: list = field(default_factory=list, init=False, repr=False)
    wildtype_epitopes: list = field(default_factory=list, init=False, repr=False)
    wildtype_epitope_score: float = field(default=0.0, init=False, repr=False)
    mutant_epitope_score: float = field(default=0.0, init=False, repr=False)
    manufacturability_scores: Any = field(default=None, init=False, repr=False)

    def __post_init__(self):
        # Normalize the optional collection/tuple fields.
        self.epitope_score_params = self.epitope_score_params or {}
        self.manufacturability_thresholds = self.manufacturability_thresholds or {}
        if self.manufacturability_rules is not None:
            self.manufacturability_rules = tuple(self.manufacturability_rules)
        if self.ranking_rules is not None:
            self.ranking_rules = tuple(self.ranking_rules)

        # Validate combined_score_mode; resolve None to the legacy default.
        resolved_mode = self.combined_score_mode or DEFAULT_COMBINED_SCORE_MODE
        if resolved_mode not in COMBINED_SCORE_MODES:
            raise ValueError(
                f"combined_score_mode must be one of {COMBINED_SCORE_MODES}, "
                f"got '{resolved_mode}'"
            )
        self.combined_score_mode = resolved_mode

        # Pre-parse + validate the optional DSL expression once at
        # construction time so syntax errors surface early (not at
        # ranking time, after we've already done all the loading
        # work). The parsed AST is stashed on the instance for
        # ``combined_score`` to evaluate per call.
        if self.combined_score_expr is not None:
            from .combined_score_dsl import parse_combined_score_expr
            self._combined_score_expr_ast = parse_combined_score_expr(
                self.combined_score_expr)
        else:
            self._combined_score_expr_ast = None

        # Sort key over CandidateEpitope: pick the strongest pMHC_affinity
        # record across all predictors / alleles inside the CandidateEpitope.
        # Legacy semantic: sort by the SINGLE best (peptide, allele,
        # predictor) row's ic50 or percentile_rank — which here is
        # the min over the flat leaf set. ``best_affinity()`` would
        # raise on multi-predictor Epitopes (cross-predictor ranking
        # ambiguous), so we iterate ``predictions_flat()`` + filter
        # by kind and take ``min``. Epitopes with no affinity record
        # sort to the end via +inf.
        sort_by = self.sort_epitopes_by

        def _epitope_sort_key(e):
            affinity_leaves = [
                p for p in e.predictions_flat()
                if p.kind == 'pMHC_affinity']
            if not affinity_leaves:
                return float("inf")
            if sort_by == "percentile_rank":
                ranks = [
                    p.percentile_rank for p in affinity_leaves
                    if p.percentile_rank is not None
                    and not (isinstance(p.percentile_rank, float)
                             and np.isnan(p.percentile_rank))]
                return min(ranks) if ranks else float("inf")
            return min(p.value for p in affinity_leaves)

        # only keep the top k epitopes
        self.mutant_epitopes = sorted(
            [
                e for e in self.epitopes
                if e.overlaps_mutation and not e.occurs_in_reference
            ],
            key=_epitope_sort_key,
        )
        if self.num_mutant_epitopes_to_keep:
            self.mutant_epitopes = (
                self.mutant_epitopes[: self.num_mutant_epitopes_to_keep]
            )

        self.wildtype_epitopes = sorted(
            [
                e for e in self.epitopes
                if not e.overlaps_mutation or e.occurs_in_reference
            ],
            key=_epitope_sort_key,
        )

        # Score each CandidateEpitope by summing the per-prediction logistic
        # score across all of its mutant affinity records. Iterates
        # ``predictions_flat()`` + kind-filters rather than
        # ``predictions_for(kind=...)`` so multi-predictor data
        # doesn't trip the predictor-disambiguation check.
        params = self.epitope_score_params

        def _epitope_total_score(e):
            return sum(
                _legacy_score_one(
                    ic50=p.value,
                    percentile_rank=p.percentile_rank,
                    **params)
                for p in e.predictions_flat()
                if p.kind == 'pMHC_affinity')

        self.wildtype_epitope_score = sum(
            _epitope_total_score(e) for e in self.wildtype_epitopes)
        self.mutant_epitope_score = sum(
            _epitope_total_score(e) for e in self.mutant_epitopes)

        self.manufacturability_scores = ManufacturabilityScores.from_amino_acids(
            self.mutant_protein_fragment.amino_acids
        )

    def peptide_synthesis_difficulty_score_tuple(
            self,
            max_c_terminal_hydropathy=1.5,
            min_kmer_hydropathy=0.0,
            max_kmer_hydropathy_low_priority=1.5,
            max_kmer_hydropathy_high_priority=2.5,
            rules=None):
        """
        Generates a tuple of scores used for lexicographic sorting of vaccine
        peptides. Each entry corresponds to one rule in ``rules`` (defaults
        to ``DEFAULT_MANUFACTURABILITY_RULES`` — the legacy 10-rule order).

        The most important default criterion is to minimize the number of
        cysteines in the sequence (to prevent disulfide bonds), followed by
        keeping C-terminal 7-mer GRAVY below 1.5 and no 7-mer window above
        2.5 (Kyte & Doolittle 1982). Remaining tie-breakers target
        problematic terminal residues, Asp-Pro bonds, and lower-priority
        hydropathy bounds.

        Users can reorder, disable, or subset these by passing a ``rules``
        list (names from ``MANUFACTURABILITY_RULE_REGISTRY``) — typically
        set once on the config and threaded through.

        (Sort criteria determined through conversations with manufacturer)
        """
        thresholds = {
            "max_c_terminal_hydropathy": max_c_terminal_hydropathy,
            "min_kmer_hydropathy": min_kmer_hydropathy,
            "max_kmer_hydropathy_low_priority": max_kmer_hydropathy_low_priority,
            "max_kmer_hydropathy_high_priority": max_kmer_hydropathy_high_priority,
        }
        return compute_manufacturability_tuple(
            self.manufacturability_scores,
            thresholds,
            rules=rules,
        )

    def lexicographic_sort_key(self):
        """
        Build the tuple used to sort vaccine peptides lexicographically.

        Rules come from ``self.ranking_rules`` if supplied, otherwise
        ``DEFAULT_RANKING_RULES`` (see ``vaxrank.ranking``). The special
        ``"manufacturability"`` sentinel expands in place using this
        peptide's own manufacturability rules + thresholds. The default
        order is byte-identical to what this method produced before the
        rules were extracted — the parity is pinned by
        ``tests/test_ranking.py::test_default_rules_match_legacy_tuple``.
        """
        return compute_ranking_tuple(self, rules=self.ranking_rules)

    def contains_mutant_epitopes(self):
        return len(self.mutant_epitopes) > 0

    @property
    def expression_score(self):
        # Honest expression metric regardless of combined_score_mode — this
        # is what reports display as "Expression score". The mode changes
        # how expression folds into `combined_score`, not the metric itself.
        return np.sqrt(self.mutant_protein_fragment.n_alt_reads)

    @property
    def combined_score(self):
        # When the DSL expression is set, it supersedes the enum
        # mode entirely — same precedence relationship as
        # ``epitopes.score_expr`` overriding the scalar logistic
        # knobs in :mod:`vaxrank.epitope_dsl`.
        if self._combined_score_expr_ast is not None:
            from .combined_score_dsl import evaluate_combined_score
            return evaluate_combined_score(
                self._combined_score_expr_ast, self)
        epitope = self.mutant_epitope_score
        if self.combined_score_mode == "reads_times_epitope":
            return float(self.mutant_protein_fragment.n_alt_reads) * epitope
        if self.combined_score_mode == "epitope_only":
            return epitope
        # Default: sqrt_reads_times_epitope (legacy)
        return self.expression_score * epitope

    def to_dict(self):
        # The persisted form combines the filtered mutant + wildtype lists
        # back into a single `epitopes` list. Also trims fields that match
        # their defaults so the JSON stays small and older readers don't
        # see keys they don't understand.
        epitopes = self.mutant_epitopes + self.wildtype_epitopes
        d = {
            "mutant_protein_fragment": self.mutant_protein_fragment,
            "epitopes": epitopes,
            "num_mutant_epitopes_to_keep": self.num_mutant_epitopes_to_keep,
            "epitope_score_params": self.epitope_score_params,
            "sort_epitopes_by": self.sort_epitopes_by,
        }
        if self.manufacturability_thresholds:
            d["manufacturability_thresholds"] = self.manufacturability_thresholds
        if self.manufacturability_rules is not None:
            d["manufacturability_rules"] = list(self.manufacturability_rules)
        if self.combined_score_mode != DEFAULT_COMBINED_SCORE_MODE:
            d["combined_score_mode"] = self.combined_score_mode
        if self.ranking_rules is not None:
            d["ranking_rules"] = list(self.ranking_rules)
        return d
