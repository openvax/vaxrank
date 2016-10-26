# Copyright (c) 2016. Mount Sinai School of Medicine
#
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

from __future__ import absolute_import, print_function, division

from argparse import ArgumentParser
from nose.tools import eq_, ok_

from mhctools import RandomBindingPredictor
from pyensembl import EnsemblRelease
from vaxrank.epitope_prediction import predict_epitopes
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.vaccine_peptide import VaccinePeptide


def test_reference_peptide_logic():
    protein_fragment = MutantProteinFragment(
    	variant=None,  # irrelevant to test
    	gene_name=u'Wdr13',
    	amino_acids='KLQGHSAPVLDVIVNCDESLLASSD',
    	mutant_amino_acid_start_offset=12,
    	mutant_amino_acid_end_offset=13,
    	n_overlapping_reads=71,
    	n_alt_reads=25,
    	n_ref_reads=46,
    	n_alt_reads_supporting_protein_sequence=2)

    epitope_predictions = predict_epitopes(
        mhc_predictor=RandomBindingPredictor(["H-2-Kb", "H-2-Db"]),
        protein_fragment=protein_fragment,
        genome=EnsemblRelease(species="mouse"))

    # occurs in protein ENSMUSP00000033506
    prediction_occurs_in_reference = epitope_predictions['NCDESLLAS']
    prediction_does_not_occur_in_reference = epitope_predictions['LDVIVNCDE']
    ok_(prediction_occurs_in_reference.occurs_in_reference)
    ok_(not prediction_does_not_occur_in_reference.occurs_in_reference)

    # construct a simple vaccine peptide having these two predictions, which makes it easy to check
    # for mutant/WT scores from single contributors
    vaccine_peptide = VaccinePeptide(
    	protein_fragment,
    	[prediction_occurs_in_reference, prediction_does_not_occur_in_reference])

    eq_(prediction_occurs_in_reference.logistic_score(),
    	vaccine_peptide.wildtype_epitope_score)
    eq_(prediction_does_not_occur_in_reference.logistic_score(),
    	vaccine_peptide.mutant_epitope_score)
