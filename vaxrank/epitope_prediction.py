# Copyright (c) 2016. Mount Sinai School of Medicine
#
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

from __future__ import absolute_import, print_function, division
from collections import namedtuple

import numpy as np

EpitopePredictionBase = namedtuple(
    "EpitopePrediction", [
        "allele",
        "peptide_sequence",
        "length",
        "ic50",
        "percentile_rank",
        "prediction_method_name",
        "overlaps_mutation",
        "source_sequence",
        "offset",
    ])

class EpitopePrediction(EpitopePredictionBase):

    @classmethod
    def from_mhctools_binding_prediction(
            cls,
            binding_prediction,
            overlaps_mutation):
        return cls(
            allele=binding_prediction.allele,
            peptide_sequence=binding_prediction.peptide,
            length=len(binding_prediction.peptide),
            ic50=binding_prediction.value,
            percentile_rank=binding_prediction.percentile_rank,
            prediction_method_name=binding_prediction.prediction_method_name,
            overlaps_mutation=overlaps_mutation,
            source_sequence=binding_prediction.source_sequence,
            offset=binding_prediction.offset)

    def logistic_score(
            self,
            midpoint=350.0,
            width=150.0,
            ic50_cutoff=2000.0):
        """
        Map from IC50 values to score where 1.0 = strong binder, 0.0 = weak binder
        Default midpoint and width for logistic determined by max likelihood fit
        for data from Alessandro Sette's 1994 paper:

           "The relationship between class I binding affinity
            and immunogenicity of potential cytotoxic T cell epitopes.

        TODO: Use a large dataset to find MHC binding range predicted to #
        correlate with immunogenicity
        """
        if self.ic50 >= ic50_cutoff:
            return 0.0

        rescaled = (float(self.ic50) - midpoint) / width
        # simplification of 1.0 - logistic(x) = logistic(-x)
        logistic = 1.0 / (1.0 + np.exp(rescaled))

        # since we're scoring IC50 values, let's normalize the output
        # so IC50 near 0.0 always returns a score of 1.0
        normalizer = 1.0 / (1.0 + np.exp(-midpoint / width))

        return logistic / normalizer

def predict_epitopes(mhc_predictor, protein_fragment):
    results = []
    mhctools_binding_predictions = mhc_predictor.predict(
        {"protein_fragment": protein_fragment.amino_acids})
    # convert from mhctools.BindingPrediction objects to EpitopePrediction
    # which differs primarily by also having a boolean field
    # 'overlaps_mutation' that indicates whether the epitope overlaps
    # mutant amino acids or both sides of a deletion
    for binding_prediction in mhctools_binding_predictions:
        peptide_start_offset = binding_prediction.offset
        peptide_end_offset = binding_prediction.offset + binding_prediction.length

        overlaps_mutation = protein_fragment.interval_overlaps_mutation(
            start_offset=peptide_start_offset,
            end_offset=peptide_end_offset)
        results.append(EpitopePrediction.from_mhctools_binding_prediction(
            binding_prediction,
            overlaps_mutation=overlaps_mutation))
    return results

def slice_epitope_predictions(
        epitope_predictions,
        start_offset,
        end_offset):
    """
    Return subsert of EpitopePrediction objects which overlap the given interval
    and slice through their source sequences and adjust their offset.
    """
    return [
        EpitopePrediction(
            allele=p.allele,
            peptide_sequence=p.peptide_sequence,
            length=p.length,
            ic50=p.ic50,
            percentile_rank=p.percentile_rank,
            prediction_method_name=p.prediction_method_name,
            overlaps_mutation=p.overlaps_mutation,
            source_sequence=p.source_sequence[start_offset:end_offset],
            offset=p.offset - start_offset)
        for p in epitope_predictions
        if p.offset >= start_offset and p.offset + p.length <= end_offset
    ]
