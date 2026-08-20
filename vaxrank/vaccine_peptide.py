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

from .combined_score_dsl import (
    evaluate_combined_score,
    parse_combined_score_expr,
)
from .manufacturability import (
    ManufacturabilityScores,
    compute_manufacturability_tuple,
)
from .ranking import compute_ranking_tuple
from .vaccine_config import DEFAULT_COMBINED_SCORE_EXPR


@dataclass
class VaccinePeptide(DataclassSerializable):
    """
    VaccinePeptide combines the sequence information of MutantProteinFragment
    with MHC binding predictions for subsequences of the protein fragment.

    The resulting lists of target and self ``CandidateEpitope`` objects
    are sorted by best (target) affinity.

    Parameters
    ----------
    mutant_protein_fragment : MutantProteinFragment

    epitopes : list of CandidateEpitope

    num_target_epitopes_to_keep : int or None
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

    combined_score_expr : str or None
        DSL expression evaluated against this VaccinePeptide to
        produce ``combined_score``. See
        :mod:`vaxrank.combined_score_dsl` for grammar and bindings
        (``target_epitope_score``, ``n_alt_reads``,
        ``expression_score``, ...). ``None`` resolves to
        :data:`vaxrank.vaccine_config.DEFAULT_COMBINED_SCORE_EXPR`.
        There is no separate mode enum — the expression IS the
        mechanism.

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
    num_target_epitopes_to_keep: Optional[int] = None
    sort_epitopes_by: str = "ic50"
    manufacturability_thresholds: Optional[dict] = None
    manufacturability_rules: Optional[tuple] = None
    # DSL expression that produces this peptide's ``combined_score``.
    # ``None`` resolves to :data:`DEFAULT_COMBINED_SCORE_EXPR`. There
    # is no fallback mode enum — this string is the single mechanism.
    combined_score_expr: Optional[str] = None
    ranking_rules: Optional[tuple] = None

    # Derived attributes computed in __post_init__ — not part of the
    # serialized form. `init=False` keeps them out of the generated
    # __init__ signature; default=None gives them a sane pre-computed
    # value for any pathway that inspects fields before __post_init__.
    target_epitopes: list = field(default_factory=list, init=False, repr=False)
    self_epitopes: list = field(default_factory=list, init=False, repr=False)
    self_epitope_score: float = field(default=0.0, init=False, repr=False)
    target_epitope_score: float = field(default=0.0, init=False, repr=False)
    manufacturability_scores: Any = field(default=None, init=False, repr=False)

    def __post_init__(self):
        # Normalize the optional collection/tuple fields.
        self.manufacturability_thresholds = self.manufacturability_thresholds or {}
        if self.manufacturability_rules is not None:
            self.manufacturability_rules = tuple(self.manufacturability_rules)
        if self.ranking_rules is not None:
            self.ranking_rules = tuple(self.ranking_rules)

        # Resolve the score expression: user-supplied or the canonical
        # default. Pre-parse once so syntax errors surface at config /
        # construction time, not after a full ranking run. The AST is
        # the single artifact ``combined_score`` evaluates against —
        # there is no fallback Python branch for the default formula.
        if self.combined_score_expr is None:
            self.combined_score_expr = DEFAULT_COMBINED_SCORE_EXPR
        self._combined_score_expr_ast = parse_combined_score_expr(
            self.combined_score_expr)

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

        def _usable_prediction_values(predictions, attr):
            values = []
            for p in predictions:
                value = getattr(p, attr)
                if value is None:
                    continue
                try:
                    if np.isnan(value):
                        continue
                except TypeError:
                    pass
                values.append(value)
            return values

        def _epitope_sort_key(e):
            affinity_leaves = [
                p for p in e.predictions_flat()
                if p.kind == 'pMHC_affinity']
            if not affinity_leaves:
                return float("inf")
            if sort_by == "percentile_rank":
                ranks = _usable_prediction_values(
                    affinity_leaves, "percentile_rank")
                return min(ranks) if ranks else float("inf")
            ic50s = _usable_prediction_values(affinity_leaves, "value")
            return min(ic50s) if ic50s else float("inf")

        # Universal target/self split — source-agnostic. An epitope
        # whose exact sequence appears in the patient's reference
        # proteome (``occurs_in_reference``) is unsafe by default
        # (would cross-react with self on normal tissue). For
        # mutation-derived candidates this matches the historical
        # behavior: peptides identical to WT land in self_epitopes
        # because they're in the reference. Viral peptides absent
        # from the reference land in target_epitopes. CTAs land in
        # self_epitopes today. Antigen-aware consumers introduced by
        # issue #303 use ``occurs_in_non_CTA_reference`` for CTA sources.
        self.target_epitopes = sorted(
            [e for e in self.epitopes if not e.occurs_in_reference],
            key=_epitope_sort_key,
        )
        if self.num_target_epitopes_to_keep:
            self.target_epitopes = (
                self.target_epitopes[: self.num_target_epitopes_to_keep]
            )

        self.self_epitopes = sorted(
            [e for e in self.epitopes if e.occurs_in_reference],
            key=_epitope_sort_key,
        )

        # The per-epitope total score lives on ``CandidateEpitope`` as
        # ``epitope_score`` (sum of its ``per_allele_scores`` populated
        # at predict time by the configured ``EpitopeConfig.score_expr``
        # — see :mod:`vaxrank.epitope_dsl`). VaccinePeptide is downstream
        # of scoring; it never recomputes from ic50 / percentile_rank.
        # This is what makes the DSL the single source of truth: change
        # the per-prediction formula in one place and the change flows
        # all the way through to ``combined_score``.
        self.self_epitope_score = sum(
            e.epitope_score for e in self.self_epitopes)
        self.target_epitope_score = sum(
            e.epitope_score for e in self.target_epitopes)

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

    def contains_target_epitopes(self):
        return len(self.target_epitopes) > 0

    @property
    def expression_score(self):
        # Honest expression metric. ``sqrt(n_alt_reads)`` is a binding
        # the default combined-score expression uses; if the user
        # writes a different ``combined_score_expr`` that doesn't
        # reference it, this property is still queryable for reports.
        return np.sqrt(self.mutant_protein_fragment.n_alt_reads)

    @property
    def combined_score(self):
        # Single code path: evaluate the (pre-parsed) DSL expression.
        # When the user didn't override, ``combined_score_expr`` was
        # resolved to ``DEFAULT_COMBINED_SCORE_EXPR`` in
        # ``__post_init__`` — the default formula is *evaluated*
        # through the DSL, not branched around.
        return evaluate_combined_score(
            self._combined_score_expr_ast, self)

    def to_dict(self):
        # The persisted form combines the filtered target + self lists
        # back into a single `epitopes` list. Also trims fields that match
        # their defaults so the JSON stays small and older readers don't
        # see keys they don't understand.
        epitopes = self.target_epitopes + self.self_epitopes
        d = {
            "mutant_protein_fragment": self.mutant_protein_fragment,
            "epitopes": epitopes,
            "num_target_epitopes_to_keep": self.num_target_epitopes_to_keep,
            "sort_epitopes_by": self.sort_epitopes_by,
        }
        if self.manufacturability_thresholds:
            d["manufacturability_thresholds"] = self.manufacturability_thresholds
        if self.manufacturability_rules is not None:
            d["manufacturability_rules"] = list(self.manufacturability_rules)
        if self.combined_score_expr != DEFAULT_COMBINED_SCORE_EXPR:
            d["combined_score_expr"] = self.combined_score_expr
        if self.ranking_rules is not None:
            d["ranking_rules"] = list(self.ranking_rules)
        return d
