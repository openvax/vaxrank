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

import numpy as np
from serializable import Serializable


class EpitopePrediction(Serializable):
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
    def __init__(
            self,
            allele,
            peptide_sequence,
            wt_peptide_sequence,
            ic50,
            wt_ic50,
            percentile_rank,
            prediction_method_name,
            overlaps_mutation,
            source_sequence,
            offset,
            occurs_in_reference):
        self.allele = allele
        self.peptide_sequence = peptide_sequence
        self.wt_peptide_sequence = wt_peptide_sequence
        self.length = len(peptide_sequence)
        self.ic50 = ic50
        self.wt_ic50 = wt_ic50
        self.percentile_rank = percentile_rank
        self.prediction_method_name = prediction_method_name
        self.overlaps_mutation = overlaps_mutation
        self.source_sequence = source_sequence
        self.offset = offset
        self.occurs_in_reference = occurs_in_reference

    @classmethod
    def from_dict(cls, d):
        """
        Deserialize EpitopePrediction from a dictionary of keywords.
        """
        d = d.copy()
        if "length" in d:
            # length argument removed in version 1.1.0
            del d["length"]
        return cls(**d)

    def logistic_epitope_score(
            self,
            midpoint=350.0,
            width=150.0,
            ic50_cutoff=5000.0,
            scoring_mode="affinity"):
        """
        Map from binding predictions to score where 1.0 = strong binder,
        0.0 = weak binder.

        Parameters
        ----------
        scoring_mode : str
            "affinity" (default): logistic transform of IC50
            "percentile_rank": inverted percentile rank (lower rank = higher score)
        """
        if scoring_mode == "percentile_rank":
            if self.percentile_rank is None:
                return 0.0
            # Percentile rank: 0 = best, 100 = worst
            # Convert to score: 0..1 where 1 = best
            # Epitopes with rank > 10% are considered weak
            rank = float(self.percentile_rank)
            if rank >= 10.0:
                return 0.0
            return max(0.0, 1.0 - rank / 10.0)

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
