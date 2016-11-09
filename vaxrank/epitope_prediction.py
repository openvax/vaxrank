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
from collections import namedtuple, defaultdict, OrderedDict
import logging
import os

from appdirs import user_cache_dir
import numpy as np
import shellinford
import six


logger = logging.getLogger(__name__)

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
        "occurs_in_reference",
    ])

class EpitopePrediction(EpitopePredictionBase):

    @classmethod
    def from_mhctools_binding_prediction(
            cls,
            binding_prediction,
            overlaps_mutation,
            occurs_in_reference):
        return cls(
            allele=binding_prediction.allele,
            peptide_sequence=binding_prediction.peptide,
            length=len(binding_prediction.peptide),
            ic50=binding_prediction.value,
            percentile_rank=binding_prediction.percentile_rank,
            prediction_method_name=binding_prediction.prediction_method_name,
            overlaps_mutation=overlaps_mutation,
            source_sequence=binding_prediction.source_sequence,
            offset=binding_prediction.offset,
            occurs_in_reference=occurs_in_reference)

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

def fm_index_path(genome):
    """
    Returns a path for cached reference peptides, for the given genome.
    """
    cache_dir = user_cache_dir('vaxrank')
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)

    return os.path.join(cache_dir, '%s_%d_%d.fm' % (
        genome.species.latin_name, genome.release, 2 if six.PY2 else 3))

def index_contains_kmer(fm, kmer):
    """
    Checks FM index for kmer, returns true if found.
    """
    found = False
    for _ in fm.search(kmer):
        found = True
        break
    return found

def load_reference_peptides_index(genome, force_reload=False):
    """
    Loads the FM index containing reference peptides.

    Parameters
    ----------
    genome : pyensembl.EnsemblRelease
        Input genome to load for reference peptides

    force_reload : bool, optional
        If true, will recompute index for this genome even if it already exists.

    Returns
    -------
    fm : shellinford.FMIndex
        Index populated with reference peptides from the genome
    """
    path = fm_index_path(genome)
    if force_reload or not os.path.exists(path):
        logger.info("Loading FM index to %s", path)
        fm = shellinford.FMIndex()
        fm.build(
            (t.protein_sequence for t in genome.transcripts() if t.is_protein_coding),
            path)
        logger.info("Done loading FM index")
        return fm
    return shellinford.FMIndex(filename=path)

def predict_epitopes(mhc_predictor, protein_fragment, min_epitope_score=0, genome=None):
    """
    Returns an OrderedDict of EpitopePrediction objects, keyed by a (peptide sequence, allele)
    tuple, that have a normalized score greater than min_epitope_score.

    Uses the input genome to evaluate whether the epitope occurs in reference.
    """
    results = OrderedDict()
    fm = load_reference_peptides_index(genome)

    mhctools_binding_predictions = mhc_predictor.predict(
        {"protein_fragment": protein_fragment.amino_acids})
    # convert from mhctools.BindingPrediction objects to EpitopePrediction
    # which differs primarily by also having a boolean field
    # 'overlaps_mutation' that indicates whether the epitope overlaps
    # mutant amino acids or both sides of a deletion
    num_total = 0
    num_occurs_in_reference = 0
    for binding_prediction in mhctools_binding_predictions:
        num_total += 1
        peptide_start_offset = binding_prediction.offset
        peptide_end_offset = binding_prediction.offset + binding_prediction.length

        peptide = binding_prediction.peptide
        occurs_in_reference = index_contains_kmer(fm, peptide)
        if occurs_in_reference:
            logger.debug('Peptide %s occurs in reference', peptide)
            num_occurs_in_reference += 1
        overlaps_mutation = protein_fragment.interval_overlaps_mutation(
            start_offset=peptide_start_offset,
            end_offset=peptide_end_offset)
        epitope_prediction = EpitopePrediction.from_mhctools_binding_prediction(
            binding_prediction,
            overlaps_mutation=overlaps_mutation,
            occurs_in_reference=occurs_in_reference)
        if epitope_prediction.logistic_score() >= min_epitope_score:
            key = (epitope_prediction.peptide_sequence, epitope_prediction.allele)
            results[key] = epitope_prediction

    logger.info('%d out of %d peptides occur in reference', num_occurs_in_reference, num_total)
    return results

def slice_epitope_predictions(
        epitope_predictions,
        start_offset,
        end_offset):
    """
    Return subset of EpitopePrediction objects which overlap the given interval
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
            offset=p.offset - start_offset,
            occurs_in_reference=p.occurs_in_reference)
        for p in epitope_predictions
        if p.offset >= start_offset and p.offset + p.length <= end_offset
    ]
