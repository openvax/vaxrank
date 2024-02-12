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

from nose.tools import eq_, ok_

from mhctools import RandomBindingPredictor
from pyensembl import genome_for_reference_name
from varcode import Variant
from vaxrank.epitope_prediction import predict_epitopes, EpitopePrediction
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.vaccine_peptide import VaccinePeptide

mouse_genome = genome_for_reference_name("GRCm38")

def test_reference_peptide_logic():

    wdr13_transcript = mouse_genome.transcripts_by_name("Wdr13-201")[0]

    protein_fragment = MutantProteinFragment(
        variant=Variant('X', '8125624', 'C', 'A'),
        gene_name='Wdr13',
        amino_acids='KLQGHSAPVLDVIVNCDESLLASSD',
        mutant_amino_acid_start_offset=12,
        mutant_amino_acid_end_offset=13,
        n_overlapping_reads=71,
        n_alt_reads=25,
        n_ref_reads=46,
        n_alt_reads_supporting_protein_sequence=2,
        supporting_reference_transcripts=[wdr13_transcript])

    epitope_predictions = predict_epitopes(
        mhc_predictor=RandomBindingPredictor(["H-2-Kb"]),
        protein_fragment=protein_fragment,
        genome=mouse_genome)

    # occurs in protein ENSMUSP00000033506
    prediction_occurs_in_reference = epitope_predictions[('NCDESLLAS', 'H-2-Kb')]
    prediction_does_not_occur_in_reference = epitope_predictions[('LDVIVNCDE', 'H-2-Kb')]
    ok_(prediction_occurs_in_reference.occurs_in_reference)
    ok_(not prediction_does_not_occur_in_reference.occurs_in_reference)

    # construct a simple vaccine peptide having these two predictions, which makes it easy to check
    # for mutant/WT scores from single contributors
    vaccine_peptide = VaccinePeptide(
        protein_fragment,
        [prediction_occurs_in_reference, prediction_does_not_occur_in_reference])

    eq_(prediction_occurs_in_reference.logistic_epitope_score(),
        vaccine_peptide.wildtype_epitope_score)
    eq_(prediction_does_not_occur_in_reference.logistic_epitope_score(),
        vaccine_peptide.mutant_epitope_score)

def test_mhc_predictor_error():
    wdr13_transcript = mouse_genome.transcripts_by_name("Wdr13-201")[0]

    protein_fragment = MutantProteinFragment(
        variant=Variant('X', '8125624', 'C', 'A'),
        gene_name='Wdr13',
        amino_acids='KLQGHSAPVLDVIVNCDESLLASSD',
        mutant_amino_acid_start_offset=12,
        mutant_amino_acid_end_offset=13,
        n_overlapping_reads=71,
        n_alt_reads=25,
        n_ref_reads=46,
        n_alt_reads_supporting_protein_sequence=2,
        supporting_reference_transcripts=[wdr13_transcript])

    # throws an error for each prediction, make sure vaxrank doesn't fall down
    class FakeMHCPredictor:
        def predict_subsequences(self, x):
            raise ValueError('I throw an error in your general direction')

    epitope_predictions = predict_epitopes(
        mhc_predictor=FakeMHCPredictor(),
        protein_fragment=protein_fragment,
        genome=mouse_genome)

    eq_(0, len(epitope_predictions))

def test_EpitopePrediction_json_serialization():
    e = EpitopePrediction(
        allele="HLA-A*02:01",
        peptide_sequence="SIINFEQL",
        ic50=2.0,
        wt_peptide_sequence="SIINFEKL",
        wt_ic50=2000.0,
        percentile_rank=0.3,
        prediction_method_name="ImaginationMHCpan",
        overlaps_mutation=True,
        source_sequence="SSIINFEQL",
        offset=1,
        occurs_in_reference=False)
    json = e.to_json()
    e2 = EpitopePrediction.from_json(json)
    eq_(e, e2)
