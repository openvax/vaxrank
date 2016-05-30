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

import numpy as np

class EpitopeScorer(object):
    def __init__(
            self,
            ic50_cutoff=500.0,
            transformation_function=None,
            allele_weights=None):
        """
        Parameters
        ----------
        ic50_cutoff : float
            Value above which we ignore epitopes

        transformation_function : fn : float -> float (optional)
            Given the predicted binding quantity, transform it into a value
            between 0.0 and 1.0

        allele_weights : dict (optional)
            Associates each MHC allele name with an importance weight
        """
        self.ic50_cutoff = ic50_cutoff
        self.transformation_function = transformation_function
        self.allele_weights = allele_weights

    def allele_weight(self, allele_name):
        if self.allele_weights:
            assert allele_name in self.allele_weights, \
                "Allele not in weight dictionary: %s" % (allele_name,)
            return self.allele_weights[allele_name]
        else:
            return 1.0

    def binding_value_score(self, ic50, allele_name=None):
        """
        Transform a binding prediction for a particular allele into a
        rescaled score (non-binders get a score of 0.0, strongest
        predicted binders should have a score of 1.0)

        Parameters
        ----------
        ic50 : float

        allele_name : string (optional)
        """
        if self.transformation_function:
            score = self.transformation_function(ic50)
        else:
            score = 1.0 if ic50 < self.ic50_cutoff else 0.0
        allele_weight = self.allele_weight(allele_name) if allele_name else 1.0
        return allele_weight * score

    def binding_record_score(self, record):
        return self.binding_value_score(
            ic50=record["IC50"],
            allele_name=record["allele"])

    def sum_binding_record_scores(self, epitope_binding_prediction_records):
        """

        Parameters
        ----------

        epitope_binding_prediction_records : dictionary of binding predictions
            Keys are alleles, values are binding prediction records with
            fields such as 'MHC_IC50' and 'MHC_Percentile_Rank'

        Returns sum of normalized binding prediction scores.
        """
        if isinstance(epitope_binding_prediction_records, dict):
            epitope_binding_prediction_records = list(
                epitope_binding_prediction_records.values())

        return sum(
            self.binding_record_score(record)
            for record
            in epitope_binding_prediction_records
        )

    def epitope_score(self, epitope):
        binding_records = epitope['MHC_Allele_Scores']
        return self.sum_binding_record_scores(binding_records)


class DecreasingLogisticFunction(object):
    def __init__(self, midpoint, width):
        assert width > 0
        self.midpoint = float(midpoint)
        self.width = float(width)

    def __call__(self, value):
        rescaled = (float(value) - self.midpoint) / self.width
        # simplification of 1.0 - logistic(x) = logistic(-x)
        logistic = 1.0 / (1.0 + np.exp(rescaled))

        # since we're scoring IC50 values, let's normalize the output
        # so IC50 near 0.0 always returns a score of 1.0
        normalizer = 1.0 / (1.0 + np.exp(-self.midpoint / self.width))

        return logistic / normalizer

# add up all the epitopes with IC50 <= 500nM
simple_class1_ic50_epitope_scorer = EpitopeScorer(cutoff=500.0)

# default midpoint and width for logistic determined by max likelihood fit
# for data from Alessandro Sette's 1994 paper:
#
#   "The relationship between class I binding affinity
#    and immunogenicity of potential cytotoxic T cell epitopes.
#
# TODO: Use a large dataset to find MHC binding range predicted to #
# correlate with immunogenicity
logistic_fn = DecreasingLogisticFunction(midpoint=350.0, width=150.0)

logistic_class1_ic50_epitope_scorer = EpitopeScorer(
    cutoff=2000.0,
    transformation_function=logistic_fn,
)
