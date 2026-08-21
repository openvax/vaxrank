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
from varcode import Variant
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.epitope_logic import predict_epitopes
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.vaccine_peptide import VaccinePeptide
from tests._legacy_score_reference import legacy_score_one as _legacy_score_one

from .common import eq_, ok_


def _by_pep_allele(epitopes):
    """Map ``(peptide, allele) -> (CandidateEpitope, mutant Prediction leaf)``
    for single-predictor test inputs."""
    out = {}
    for e in epitopes:
        for p in e.predictions_flat():
            if p.kind != 'pMHC_affinity':
                continue
            out[(e.sequence, p.allele)] = (e, p)
    return out


def test_reference_peptide_logic(mouse_genome):

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

    epitope_config = EpitopeConfig(min_epitope_score=0)
    epitopes = predict_epitopes(
        mhc_predictor=RandomBindingPredictor(["H-2-Kb"]),
        protein_fragment=protein_fragment,
        epitope_config=epitope_config,
        genome=mouse_genome)

    # Single-predictor test: one Prediction per (peptide, allele).
    by_key = _by_pep_allele(epitopes)
    # occurs in protein ENSMUSP00000033506
    e_in_ref, p_in_ref = by_key[('NCDESLLAS', 'H-2-Kb')]
    e_not_in_ref, p_not_in_ref = by_key[('LDVIVNCDE', 'H-2-Kb')]
    ok_(e_in_ref.occurs_in_reference)
    ok_(not e_not_in_ref.occurs_in_reference)
    assert e_in_ref.self_reference_match.occurs
    assert not e_not_in_ref.self_reference_match.occurs
    assert e_in_ref.self_reference_match.source_provenance_complete
    assert e_not_in_ref.self_reference_match.source_provenance_complete
    assert any(
        source.protein_id == "ENSMUSP00000033506"
        for source in e_in_ref.self_reference_match.sources
    )
    assert e_not_in_ref.self_reference_match.sources == ()
    assert e_in_ref.overlaps_targetable == e_in_ref.overlaps_mutation
    assert e_not_in_ref.overlaps_targetable == e_not_in_ref.overlaps_mutation

    # Build a VaccinePeptide from these two Epitopes; the mutant /
    # wildtype split comes out of __post_init__.
    vaccine_peptide = VaccinePeptide(
        protein_fragment, [e_in_ref, e_not_in_ref])
    assert vaccine_peptide.antigen.kind == "mutation"

    eq_(_legacy_score_one(p_in_ref.value, p_in_ref.percentile_rank),
        vaccine_peptide.self_epitope_score)
    eq_(_legacy_score_one(p_not_in_ref.value, p_not_in_ref.percentile_rank),
        vaccine_peptide.target_epitope_score)


def test_predict_epitopes_records_distinct_non_cta_reference_flag(mouse_genome):
    """The CTA-aware flag must come from its own reference index."""
    from unittest.mock import patch

    wdr13_transcript = mouse_genome.transcripts_by_name("Wdr13-201")[0]
    protein_fragment = MutantProteinFragment(
        variant=Variant("X", "8125624", "C", "A"),
        gene_name="Wdr13",
        amino_acids="KLQGHSAPVLDVIVNCDESLLASSD",
        mutant_amino_acid_start_offset=12,
        mutant_amino_acid_end_offset=13,
        n_overlapping_reads=71,
        n_alt_reads=25,
        n_ref_reads=46,
        n_alt_reads_supporting_protein_sequence=2,
        supporting_reference_transcripts=[wdr13_transcript],
    )

    with patch("vaxrank.epitope_logic.ReferenceProteome") as reference_cls:
        reference_cls.return_value.contains.return_value = True
        reference_cls.from_genome.return_value.contains.return_value = False
        epitopes = predict_epitopes(
            mhc_predictor=RandomBindingPredictor(["H-2-Kb"]),
            protein_fragment=protein_fragment,
            epitope_config=EpitopeConfig(min_epitope_score=0),
            genome=mouse_genome,
        )

    assert epitopes
    assert all(epitope.occurs_in_reference for epitope in epitopes)
    assert all(not epitope.occurs_in_non_CTA_reference for epitope in epitopes)
    reference_cls.from_genome.assert_called_once_with(
        mouse_genome, exclude_cta_genes=True)

def test_predict_epitopes_returns_one_row_per_predictor_for_multimodel(mouse_genome):
    """Multi-model TopiaryPredictor (post-2.24, #261) keeps each
    predictor's view of every (peptide, allele) pair. In the new
    CandidateEpitope shape, one CandidateEpitope per (peptide, source, offset) carries
    multiple per-(predictor, allele) leaf Predictions — verifies that
    no leaf record is lost in the grouping."""
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

    epitope_config = EpitopeConfig(
        min_epitope_score=0,
        default_methods={'pMHC_affinity': 'mhcflurry'})
    with patch('vaxrank.epitope_logic.ReferenceProteome'):
        epitopes = predict_epitopes(
            mhc_predictor=_StubTopiary(),
            protein_fragment=protein_fragment,
            epitope_config=epitope_config,
            genome=mouse_genome)
    # 2 unique peptide positions; each carries 2-predictor × 1-allele
    # leaf records (4 leaves total). Pre-fix dict-keyed return would
    # have lost one predictor per (peptide, allele).
    assert len(epitopes) == 2
    leaves = [p for e in epitopes for p in e.predictions_flat()]
    assert len(leaves) == 4
    methods = sorted({p.predictor_name for p in leaves})
    assert methods == ['mhcflurry', 'netmhcpan']
    for e in epitopes:
        per_method = {p.predictor_name for p in e.predictions_flat()}
        assert per_method == {'mhcflurry', 'netmhcpan'}


def test_wt_prediction_uses_named_peptides_keyed_by_mutant(mouse_genome):
    """WT comparator path: ``predict_from_named_peptides`` is called with
    each entry named after its mutant peptide, and the resulting row's
    own ``peptide`` column populates the CandidateEpitope's WT — not a
    side dict. Locks in the one-hop lookup invariant."""
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

    # 9-mer at offset 8 covers [8, 17) — overlaps mutation at [12, 13).
    mutant_peptide = protein_fragment.amino_acids[8:17]
    mutant_df = pd.DataFrame([{
        'peptide': mutant_peptide, 'peptide_length': 9,
        'peptide_offset': 8, 'allele': 'H-2-Kb',
        'affinity': 100.0, 'percentile_rank': 0.5,
        'value': 100.0, 'prediction_method_name': 'mhcflurry',
        'source_sequence_name': 'Wdr13', 'kind': 'pMHC_affinity',
        'predictor_version': '', 'score': 0.5,
    }])

    # Distinct WT string so we can verify the consumer reads it from
    # the row, not from a side dict.
    stub_wt_sequence = 'WTSEQUENC'
    captured_calls = []

    def _stub_predict_wt(name_to_peptide):
        captured_calls.append(dict(name_to_peptide))
        return pd.DataFrame([{
            'source_sequence_name': name, 'peptide': stub_wt_sequence,
            'peptide_length': 9, 'peptide_offset': 0,
            'allele': 'H-2-Kb', 'affinity': 1000.0,
            'percentile_rank': 5.0, 'value': 1000.0,
            'prediction_method_name': 'mhcflurry',
            'predictor_version': '', 'score': 0.1,
            'kind': 'pMHC_affinity',
        } for name in name_to_peptide])

    class _StubTopiary(TopiaryPredictor):
        def __init__(self):
            pass

        def predict_from_named_sequences(self, _named_sequences):
            return mutant_df

        def predict_from_named_peptides(self, name_to_peptide):
            return _stub_predict_wt(name_to_peptide)

    cfg = EpitopeConfig(
        min_epitope_score=0,
        default_methods={'pMHC_affinity': 'mhcflurry'})
    with patch('vaxrank.epitope_logic.ReferenceProteome'):
        epitopes = predict_epitopes(
            mhc_predictor=_StubTopiary(),
            protein_fragment=protein_fragment,
            epitope_config=cfg,
            genome=mouse_genome)

    # WT call happened exactly once, named by the mutant peptide.
    assert len(captured_calls) == 1
    assert list(captured_calls[0].keys()) == [mutant_peptide]

    # The CandidateEpitope's WT peptide came from the row's
    # ``peptide`` column, not from any side dict.
    assert len(epitopes) == 1
    wt = epitopes[0].wt
    assert wt is not None
    assert wt.sequence == stub_wt_sequence
    wt_leaves = wt.predictions_flat()
    assert len(wt_leaves) == 1
    assert wt_leaves[0].value == 1000.0


def test_mhc_predictor_error(mouse_genome):
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
