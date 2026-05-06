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

from mhctools import RandomBindingPredictor
from pyensembl import genome_for_reference_name
from varcode import Variant
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_logic import predict_epitopes
from vaxrank.epitope_prediction import EpitopePrediction
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.vaccine_peptide import VaccinePeptide

from .common import eq_, ok_

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

    # Use min_epitope_score=0 to ensure no epitopes are filtered by score
    # (RandomBindingPredictor generates random IC50 values that can cause filtering)
    epitope_config = EpitopeConfig(min_epitope_score=0)
    epitope_predictions = predict_epitopes(
        mhc_predictor=RandomBindingPredictor(["H-2-Kb"]),
        protein_fragment=protein_fragment,
        epitope_config=epitope_config,
        genome=mouse_genome)

    # ``predict_epitopes`` returns a flat list (post-2.24, #261) so
    # multi-predictor runs preserve every prediction. Single-predictor
    # tests like this one look up by (peptide, allele) themselves.
    by_key = {
        (p.peptide_sequence, p.allele): p for p in epitope_predictions}
    # occurs in protein ENSMUSP00000033506
    prediction_occurs_in_reference = by_key[('NCDESLLAS', 'H-2-Kb')]
    prediction_does_not_occur_in_reference = by_key[('LDVIVNCDE', 'H-2-Kb')]
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

def test_predict_epitopes_returns_one_row_per_predictor_for_multimodel():
    """Multi-model TopiaryPredictor (post-2.24, #261) emits one
    EpitopePrediction per (peptide, allele, predictor) instead of
    overwriting via dict-keyed-on-(peptide,allele). Verifies the
    flat-list return shape preserves every predictor's view."""
    from unittest.mock import patch
    import pandas as pd
    from topiary import TopiaryPredictor

    wdr13_transcript = mouse_genome.transcripts_by_name("Wdr13-201")[0]
    protein_fragment = MutantProteinFragment(
        variant=Variant('X', '8125624', 'C', 'A'),
        gene_name='Wdr13',
        amino_acids='KLQGHSAPVLDVIVNCDESLLASSD',
        mutant_amino_acid_start_offset=12,
        mutant_amino_acid_end_offset=13,
        n_overlapping_reads=71, n_alt_reads=25, n_ref_reads=46,
        n_alt_reads_supporting_protein_sequence=2,
        supporting_reference_transcripts=[wdr13_transcript])

    # Stub a Topiary DF with two predictors emitting the same
    # (peptide, allele) pair — the situation that pre-2.24 collapsed
    # to one row in the dict-keyed return value.
    rows = []
    for predictor in ('mhcflurry', 'netmhcpan'):
        for offset, peptide in [(0, 'KLQGHSAP'), (1, 'LQGHSAPV')]:
            rows.append({
                'peptide': peptide, 'peptide_length': 8,
                'peptide_offset': offset, 'allele': 'H-2-Kb',
                'affinity': 100.0, 'percentile_rank': 0.5,
                'value': 100.0, 'prediction_method_name': predictor,
                'source_sequence_name': 'Wdr13',
                'kind': 'pMHC_affinity',
                'predictor_version': '',
                'score': 0.5,
            })
    fake_df = pd.DataFrame(rows)

    class _StubTopiary(TopiaryPredictor):
        def __init__(self):
            pass  # don't call super().__init__ — we don't need real models

        def predict_from_named_sequences(self, _named_sequences):
            return fake_df

    # Multi-predictor mode: tell the topiary DSL which model to
    # default ``affinity`` to (or use bracket syntax in score_expr).
    # This is the real-world contract for multi-predictor scoring.
    epitope_config = EpitopeConfig(
        min_epitope_score=0,
        default_methods={'pMHC_affinity': 'mhcflurry'})
    with patch('vaxrank.epitope_logic.ReferenceProteome'):
        preds = predict_epitopes(
            mhc_predictor=_StubTopiary(),
            protein_fragment=protein_fragment,
            epitope_config=epitope_config,
            genome=mouse_genome)
    # 2 predictors × 2 peptides × 1 allele = 4 predictions; pre-fix
    # this would collapse to 2 (one per (peptide, allele) key).
    assert len(preds) == 4
    methods = sorted({p.prediction_method_name for p in preds})
    assert methods == ['mhcflurry', 'netmhcpan']
    # Same (peptide, allele) appears twice, once per predictor.
    by_pep_allele = {}
    for p in preds:
        by_pep_allele.setdefault(
            (p.peptide_sequence, p.allele), []).append(p.prediction_method_name)
    for key, names in by_pep_allele.items():
        assert sorted(names) == ['mhcflurry', 'netmhcpan'], (
            "Expected both predictors per (peptide, allele); got %r for %r"
            % (names, key))


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
        default_peptide_lengths = [9]
        def predict_subsequences(self, x):
            raise ValueError('I throw an error in your general direction')
        def predict_proteins_dataframe(self, x):
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
