# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from dataclasses import dataclass
from typing import ClassVar, Optional

import numpy as np
from serializable import DataclassSerializable

from .config.defaults import (
    DEFAULT_BINDING_AFFINITY_CUTOFF,
    DEFAULT_LOGISTIC_EPITOPE_SCORE_MIDPOINT,
    DEFAULT_LOGISTIC_EPITOPE_SCORE_WIDTH,
    DEFAULT_PERCENTILE_RANK_CUTOFF,
)


@dataclass
class EpitopePrediction(DataclassSerializable):
    # TODO:
    #   - rename to CandidateEpitope
    #   - add groups of predictions:
    #      * mhc_affinity_ic50
    #      * mhc_affinity_score
    #      * mhc_affinity_percentile_rank
    #      * mhc_stability_hours
    #      * mhc_stability_score
    #      * mhc_stability_percentile_rank
    #      * antigen_processing_score
    #      * antigen_processing_percentile_rank
    #      * pmhc_presentation_score
    #      * pmhc_presentation_percentile_rank
    #      * immunogenicity_score
    #      * immunogenicity_percentile_rank
    #  - also change wt_peptide_sequence to:
    #      * closest_reference_peptide_sequence_in_any_protein
    #      * closest_reference_peptide_sequence_in_same_protein
    #  - and calculate:
    #      * edit_distance_to_closest_reference_peptide_in_any_protein
    #      * edit_distance_to_closest_reference_peptide_in_same_protein
    #      * edit_distance_to_closest_reference_peptide_in_any_protein_PMBEC
    #      * edit_distance_to_closest_reference_peptide_in_same_protein_PMBEC
    #  - how to avoid duplicating every prediction for MT vs. WT-any vs. WT-same?
    allele: str
    peptide_sequence: str
    wt_peptide_sequence: str
    ic50: float
    wt_ic50: float
    percentile_rank: float
    prediction_method_name: str
    overlaps_mutation: bool
    source_sequence: str
    offset: int
    occurs_in_reference: bool
    # Empty when the loader doesn't know the predictor version; set by
    # load_lens from column prefixes like "mhcflurry_2.1.1.aff".
    predictor_version: str = ""

    # `length` used to be a constructor arg; since 1.1.0 it is derived from
    # `peptide_sequence`. Drop it from any old JSON blobs we happen to load.
    _SERIALIZABLE_KEYWORD_ALIASES: ClassVar[dict[str, Optional[str]]] = {
        "length": None,
    }

    @property
    def length(self) -> int:
        return len(self.peptide_sequence)

    def logistic_epitope_score(
            self,
            midpoint=DEFAULT_LOGISTIC_EPITOPE_SCORE_MIDPOINT,
            width=DEFAULT_LOGISTIC_EPITOPE_SCORE_WIDTH,
            ic50_cutoff=DEFAULT_BINDING_AFFINITY_CUTOFF,
            scoring_mode="affinity",
            percentile_rank_cutoff=DEFAULT_PERCENTILE_RANK_CUTOFF):
        """
        Map from binding predictions to score where 1.0 = strong binder,
        0.0 = weak binder.

        Parameters
        ----------
        scoring_mode : str
            "affinity" (default): logistic transform of IC50
            "percentile_rank": inverted percentile rank (lower rank = higher score)

        .. deprecated:: 2.1.0
            ``predict_epitopes`` now filters and scores via the topiary 5.0.0
            DSL (see :func:`vaxrank.epitope_dsl.build_score_node`). This
            method remains for :class:`~vaxrank.vaccine_peptide.VaccinePeptide`
            and report/IO paths that still score one prediction at a time;
            it will be migrated and removed in a follow-up release. Configure
            ``filter_expr`` / ``score_expr`` on ``EpitopeConfig`` for custom
            DSL scoring.
        """
        if scoring_mode == "percentile_rank":
            if self.percentile_rank is None:
                return 0.0
            # Percentile rank: 0 = best, 100 = worst
            # Convert to score: 0..1 where 1 = best
            # Epitopes with rank > 10% are considered weak
            rank = float(self.percentile_rank)
            if rank >= percentile_rank_cutoff:
                return 0.0
            return max(0.0, 1.0 - rank / percentile_rank_cutoff)

        # Default: affinity-based logistic scoring
        if self.ic50 >= ic50_cutoff:
            return 0.0

        rescaled = (float(self.ic50) - midpoint) / width
        logistic = 1.0 / (1.0 + np.exp(rescaled))
        normalizer = 1.0 / (1.0 + np.exp(-midpoint / width))
        return logistic / normalizer

    def slice_source_sequence(self, start_offset, end_offset):
        """

        Parameters
        ----------
        start_offset : int

        end_offset : int

        Return EpitopePrediction object with source sequence and offset
        adjusted. If this slicing would shorten the mutant peptide, then
        return None.
        """
        if self.offset < start_offset:
            # this peptide starts before the requested slice through the
            # source sequence
            return None

        if self.offset + self.length > end_offset:
            # this peptide goes beyond the end of the requested slice
            # through the source sequence
            return None

        return EpitopePrediction(
            allele=self.allele,
            peptide_sequence=self.peptide_sequence,
            wt_peptide_sequence=self.wt_peptide_sequence,
            ic50=self.ic50,
            wt_ic50=self.wt_ic50,
            percentile_rank=self.percentile_rank,
            prediction_method_name=self.prediction_method_name,
            overlaps_mutation=self.overlaps_mutation,
            source_sequence=self.source_sequence[start_offset:end_offset],
            offset=self.offset - start_offset,
            occurs_in_reference=self.occurs_in_reference)
