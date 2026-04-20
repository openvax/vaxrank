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


from operator import attrgetter

import numpy as np
from serializable import Serializable

from .manufacturability import (
    ManufacturabilityScores,
    compute_manufacturability_tuple,
)
from .vaccine_config import COMBINED_SCORE_MODES, DEFAULT_COMBINED_SCORE_MODE


class VaccinePeptide(Serializable):
    """
    VaccinePeptide combines the sequence information of MutantProteinFragment
    with MHC binding predictions for subsequences of the protein fragment.

    The resulting lists of mutant and wildtype epitope predictions
    are sorted by affinity.
    """

    def __init__(
            self,
            mutant_protein_fragment,
            epitope_predictions,
            num_mutant_epitopes_to_keep=None,
            epitope_score_params=None,
            sort_predictions_by='ic50',
            manufacturability_thresholds=None,
            manufacturability_rules=None,
            combined_score_mode=None):
        """
        Parameters
        ----------
        mutant_protein_fragment : MutantProteinFragment

        epitope_predictions : list of EpitopePrediction

        num_mutant_epitopes_to_keep : int or None
            If None or 0 then keep all mutant epitopes.

        epitope_score_params : dict or None
            Parameters passed to EpitopePrediction.logistic_epitope_score.

        sort_predictions_by : str
            Field of EpitopePrediction used for sorting epitope predictions
            overlapping mutation in ascending order. Can be either 'ic50'
            or 'percentile_rank'.

        manufacturability_thresholds : dict or None
            Hydropathy thresholds for peptide synthesis difficulty scoring.
            Keys match the parameters of
            ``peptide_synthesis_difficulty_score_tuple``.  When None the
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
        """
        self.mutant_protein_fragment = mutant_protein_fragment
        self.epitope_predictions = epitope_predictions
        self.num_mutant_epitopes_to_keep = num_mutant_epitopes_to_keep
        self.epitope_score_params = epitope_score_params or {}
        self.sort_predictions_by = sort_predictions_by
        self.manufacturability_thresholds = manufacturability_thresholds or {}
        self.manufacturability_rules = (
            tuple(manufacturability_rules)
            if manufacturability_rules is not None
            else None
        )
        resolved_mode = combined_score_mode or DEFAULT_COMBINED_SCORE_MODE
        if resolved_mode not in COMBINED_SCORE_MODES:
            raise ValueError(
                f"combined_score_mode must be one of {COMBINED_SCORE_MODES}, "
                f"got '{resolved_mode}'"
            )
        self.combined_score_mode = resolved_mode

        sort_key = attrgetter(sort_predictions_by)

        # only keep the top k epitopes
        self.mutant_epitope_predictions = sorted([
            p for p in epitope_predictions
            if p.overlaps_mutation and not p.occurs_in_reference
        ], key=sort_key)
        if num_mutant_epitopes_to_keep:
            self.mutant_epitope_predictions = \
                self.mutant_epitope_predictions[:num_mutant_epitopes_to_keep]

        self.wildtype_epitope_predictions = sorted([
            p for p in epitope_predictions
            if not p.overlaps_mutation or p.occurs_in_reference
        ], key=sort_key)

        def epitope_score(prediction):
            return prediction.logistic_epitope_score(**self.epitope_score_params)

        self.wildtype_epitope_score = sum(
            epitope_score(p)
            for p in self.wildtype_epitope_predictions)
        # only keep the top k epitopes for the purposes of the score
        self.mutant_epitope_score = sum(
            epitope_score(p)
            for p in self.mutant_epitope_predictions)

        self.manufacturability_scores = \
            ManufacturabilityScores.from_amino_acids(
                self.mutant_protein_fragment.amino_acids)

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
        Create tuple of scores so that candidates get sorted lexicographically
        by multiple criteria. Make sure to make the wildtype epitope
        score positive (since we want fewer wildtype epitopes) but the others
        negative (since we want more of them).
        """
        # since we're sorting in decreasing order, numbers which we want
        # to be larger must have their signs flipped
        essential_score_tuple = (
            # Sum of normalized MHC binding affinities of subsequences
            # round to 5 digits to avoid floating point errors from
            # serving as tie-breakers
            -round(self.mutant_epitope_score, 6),

            # Number of reads supporting the variant
            -self.mutant_protein_fragment.n_alt_reads
        )
        manufacturability_score_tuple = self.peptide_synthesis_difficulty_score_tuple(
            rules=self.manufacturability_rules,
            **self.manufacturability_thresholds)
        extra_score_tuple = (
            # Number of reads supporting the particular protein sequence
            # sequence we're using for this vaccine peptide. Currently
            # all vaccine peptides are drawn from the same larger sequence
            # so this score shouldn't change.
            -self.mutant_protein_fragment.n_alt_reads_supporting_protein_sequence,

            # Minimize the sum of non-mutant MHC binding scores,
            # round to prevent floating point errors from serving as
            # tie-breakers
            round(self.wildtype_epitope_score, 6),

            # All else being equal, we prefer to maximize the number of
            # mutant amino acids
            -self.mutant_protein_fragment.n_mutant_amino_acids,

            # If nothing else can serve as a tie break then try to center
            # the mutation in the vaccine peptide.
            -self.mutant_protein_fragment.mutation_distance_from_edge
        )
        return (
            essential_score_tuple +
            manufacturability_score_tuple +
            extra_score_tuple
        )

    def contains_mutant_epitopes(self):
        return len(self.mutant_epitope_predictions) > 0

    @property
    def expression_score(self):
        # Honest expression metric regardless of combined_score_mode — this
        # is what reports display as "Expression score". The mode changes
        # how expression folds into `combined_score`, not the metric itself.
        return np.sqrt(self.mutant_protein_fragment.n_alt_reads)

    @property
    def combined_score(self):
        epitope = self.mutant_epitope_score
        if self.combined_score_mode == "reads_times_epitope":
            return float(self.mutant_protein_fragment.n_alt_reads) * epitope
        if self.combined_score_mode == "epitope_only":
            return epitope
        # Default: sqrt_reads_times_epitope (legacy)
        return self.expression_score * epitope

    def to_dict(self):
        epitope_predictions = self.mutant_epitope_predictions + self.wildtype_epitope_predictions
        d = {
            "mutant_protein_fragment": self.mutant_protein_fragment,
            "epitope_predictions": epitope_predictions,
            "num_mutant_epitopes_to_keep": self.num_mutant_epitopes_to_keep,
            "epitope_score_params": self.epitope_score_params,
            "sort_predictions_by": self.sort_predictions_by,
        }
        if self.manufacturability_thresholds:
            d["manufacturability_thresholds"] = self.manufacturability_thresholds
        if self.manufacturability_rules is not None:
            d["manufacturability_rules"] = list(self.manufacturability_rules)
        if self.combined_score_mode != DEFAULT_COMBINED_SCORE_MODE:
            d["combined_score_mode"] = self.combined_score_mode
        return d

    @classmethod
    def from_dict(cls, d):
        d = d.copy()
        if "sort_predictions_by" not in d:
            d["sort_predictions_by"] = "ic50"
        if "epitope_score_params" not in d:
            d["epitope_score_params"] = None
        if "manufacturability_thresholds" not in d:
            d["manufacturability_thresholds"] = None
        if "manufacturability_rules" not in d:
            d["manufacturability_rules"] = None
        if "combined_score_mode" not in d:
            d["combined_score_mode"] = None
        return cls(**d)
