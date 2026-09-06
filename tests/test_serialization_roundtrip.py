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

"""Unit-level JSON / pickle / dict roundtrip coverage for the classes
migrated from ``serializable.Serializable`` to
``serializable.DataclassSerializable``. The existing integration tests
only exercise these code paths when an MHC predictor is available."""

import pickle

import pytest
from mhctools.pred import Prediction
from varcode import Variant

from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.candidate_epitope import COMPARATOR_WT, CandidateEpitope, Peptide
from vaxrank.vaccine_antigen import (
    SelfReferenceMatch,
    SelfReferenceSource,
    VaccineAntigen,
)
from vaxrank.vaccine_peptide import VaccinePeptide


def _sample_variant():
    # varcode.Variant with ensembl=None constructs without touching the
    # Ensembl DB — enough for serialization roundtrips.
    return Variant(contig="1", start=1000, ref="A", alt="G", ensembl=None)


def _sample_fragment(amino_acids="MKVTHFYL" * 3, n_alt_reads=50):
    return MutantProteinFragment(
        variant=_sample_variant(),
        gene_name="BRCA1",
        amino_acids=amino_acids,
        mutant_amino_acid_start_offset=4,
        mutant_amino_acid_end_offset=5,
        supporting_reference_transcripts=[],
        n_overlapping_reads=100,
        n_alt_reads=n_alt_reads,
        n_ref_reads=50,
        n_alt_reads_supporting_protein_sequence=n_alt_reads,
    )


def _sample_epitope(
    peptide_sequence="SIINFEQL",
    ic50=2.0,
    overlaps_mutation=True,
    occurs_in_reference=False,
    overlaps_targetable=None,
    self_reference_match=None,
):
    mutant_pred = Prediction(
        kind='pMHC_affinity', predictor_name='ImaginationMHCpan',
        predictor_version='', allele='HLA-A*02:01',
        peptide=peptide_sequence, value=ic50, score=0.0,
        percentile_rank=0.3)
    wt_pred = Prediction(
        kind='pMHC_affinity', predictor_name='ImaginationMHCpan',
        predictor_version='', allele='HLA-A*02:01',
        peptide='SIINFEKL', value=2000.0, score=0.0,
        percentile_rank=None)
    return CandidateEpitope.from_peptide(
        Peptide(
            sequence=peptide_sequence,
            source_sequence='SSIINFEQL', offset=1,
            predictions=(mutant_pred,)),
        comparators={COMPARATOR_WT: Peptide(
            sequence='SIINFEKL',
            source_sequence='SSIINFEQL', offset=1,
            predictions=(wt_pred,))},
        overlaps_mutation=overlaps_mutation,
        overlaps_targetable=overlaps_targetable,
        occurs_in_reference=occurs_in_reference,
        self_reference_match=self_reference_match)


class _FakeVariant:
    """Minimal Variant stand-in — enough for VaccinePeptide's sort key to
    function without depending on varcode."""

    def __init__(self, gene_name="BRAF"):
        self.gene_name = gene_name


class _FakeFragment:
    """Duck-typed MutantProteinFragment for VaccinePeptide JSON tests that
    don't exercise MutantProteinFragment's own serialization path (which
    requires varcode.Variant)."""

    def __init__(self, amino_acids="A" * 25, n_alt_reads=9):
        self.amino_acids = amino_acids
        self.mutant_amino_acid_start_offset = 12
        self.mutant_amino_acid_end_offset = 13
        self.n_alt_reads = n_alt_reads
        self.n_alt_reads_supporting_protein_sequence = n_alt_reads
        self.n_mutant_amino_acids = 1
        self.mutation_distance_from_edge = 5


# ---- CandidateEpitope -----------------------------------------------------------


def test_epitope_json_roundtrip():
    """CandidateEpitope and its nested Peptide graph use one dataclass schema."""
    self_reference_match = SelfReferenceMatch(
        peptide="SIINFEQL",
        occurs=True,
        antigen_kind="mutation",
        excluded_gene_ids=("ENSG00000141510",),
        sources=(SelfReferenceSource(
            gene_id="ENSG00000146648",
            transcript_id="ENST00000275493",
            protein_id="ENSP00000275493",
            gene_name="EGFR",
            species="human",
        ),),
        source_provenance_complete=True,
        genome_release="GRCh38",
    )
    epitope = _sample_epitope(self_reference_match=self_reference_match)
    restored = CandidateEpitope.from_json(epitope.to_json())
    assert restored == epitope


def test_epitope_pickle_roundtrip():
    """``CandidateEpitope`` and its nested ``Peptide`` survive pickle."""
    e = _sample_epitope()
    assert pickle.loads(pickle.dumps(e)) == e


def test_epitope_length_is_derived():
    """``length`` is a convenience property over the mutant context's
    peptide sequence."""
    e = _sample_epitope(peptide_sequence="SIINFEQLL")
    assert e.length == 9
    assert e.sequence == "SIINFEQLL"


# ---- MutantProteinFragment --------------------------------------------


def test_mutant_protein_fragment_json_roundtrip():
    f = _sample_fragment()
    restored = MutantProteinFragment.from_json(f.to_json())
    assert restored == f


def test_mutant_protein_fragment_pickle_roundtrip():
    f = _sample_fragment()
    assert pickle.loads(pickle.dumps(f)) == f


def test_mutant_protein_fragment_computed_properties_after_roundtrip():
    """@property methods (n_mutant_amino_acids, mutation_distance_from_edge,
    n_other_reads) must still compute correctly after a JSON roundtrip —
    guards against anyone accidentally making them stored fields."""
    f = _sample_fragment()
    restored = MutantProteinFragment.from_json(f.to_json())
    assert restored.n_mutant_amino_acids == f.n_mutant_amino_acids
    assert restored.mutation_distance_from_edge == f.mutation_distance_from_edge
    assert restored.n_other_reads == f.n_other_reads
    assert len(restored) == len(f)


# ---- VaccinePeptide ----------------------------------------------------


def test_vaccine_peptide_post_init_derives_epitope_lists():
    """Filtering and sorting moved from __init__ to __post_init__ but
    must produce the same mutant / wildtype split."""
    mutant = _sample_epitope(peptide_sequence="MUTANT", ic50=5.0)
    wildtype = _sample_epitope(
        peptide_sequence="REFERNC",
        ic50=50.0,
        overlaps_mutation=False,
        occurs_in_reference=True,
    )
    vp = VaccinePeptide(
        mutant_protein_fragment=_FakeFragment(),
        epitopes=[mutant, wildtype],
    )
    assert vp.target_epitopes == [mutant]
    assert vp.self_epitopes == [wildtype]


def test_vaccine_peptide_post_init_validates_combined_score_expr():
    """__post_init__ pre-parses the expression so syntax errors surface at
    construction time rather than during ranking. 3.1 single-mechanism
    contract — the expression IS the validated input, no enum fallback."""
    with pytest.raises(ValueError, match="combined_score_expr"):
        VaccinePeptide(
            mutant_protein_fragment=_FakeFragment(),
            epitopes=[],
            combined_score_expr="not !! valid",
        )


def test_vaccine_peptide_to_dict_omits_derived_fields():
    """The custom to_dict only emits constructor args, not the
    init=False dataclass fields (target_epitopes, etc.)."""
    vp = VaccinePeptide(
        mutant_protein_fragment=_FakeFragment(),
        epitopes=[_sample_epitope()],
    )
    d = vp.to_dict()
    for derived in (
        "target_epitopes",
        "self_epitopes",
        "non_target_epitopes",
        "target_epitope_score",
        "self_epitope_score",
        "manufacturability_scores",
    ):
        assert derived not in d


def test_vaccine_peptide_tuple_coercion_on_rules():
    """List of rule names coerced to tuple even after migration."""
    vp = VaccinePeptide(
        mutant_protein_fragment=_FakeFragment(),
        epitopes=[_sample_epitope()],
        manufacturability_rules=["cysteine_count", "cterm_hydropathy"],
    )
    assert vp.manufacturability_rules == ("cysteine_count", "cterm_hydropathy")
    assert isinstance(vp.manufacturability_rules, tuple)


def test_vaccine_peptide_json_roundtrip_with_real_fragment():
    """Full to_json / from_json with a real MutantProteinFragment and
    flat records. This is the path the CLI writes out — we want
    __post_init__ to re-derive the mutant/wildtype split after load and
    scores to re-compute consistently."""
    fragment = _sample_fragment()
    antigen = VaccineAntigen.from_mutant_protein_fragment(fragment)
    mutant = _sample_epitope(
        peptide_sequence="SIINFEKL",
        ic50=500.0,
        overlaps_targetable=True,
        self_reference_match=antigen.self_reference_match("SIINFEKL", False),
    )
    wildtype = _sample_epitope(
        peptide_sequence="TESTWILD",
        ic50=50.0,
        overlaps_mutation=False,
        occurs_in_reference=True,
    )
    vp = VaccinePeptide(
        mutant_protein_fragment=fragment,
        epitopes=[mutant, wildtype],
        manufacturability_rules=("cysteine_count", "cterm_hydropathy"),
        combined_score_expr="target_epitope_score",
    )
    restored = VaccinePeptide.from_json(vp.to_json())

    # Wire-level invariants
    assert restored.mutant_protein_fragment == vp.mutant_protein_fragment
    assert restored.antigen == vp.antigen
    assert restored.target_epitopes[0].self_reference_match == (
        vp.target_epitopes[0].self_reference_match
    )
    assert restored.manufacturability_rules == vp.manufacturability_rules
    assert restored.combined_score_expr == vp.combined_score_expr

    # Derived state must match — __post_init__ re-ran on load
    assert restored.target_epitopes == vp.target_epitopes
    assert restored.self_epitopes == vp.self_epitopes
    assert restored.target_epitope_score == pytest.approx(vp.target_epitope_score)
    assert restored.combined_score == pytest.approx(vp.combined_score)


def test_vaccine_peptide_pickle_roundtrip_with_real_fragment():
    vp = VaccinePeptide(
        mutant_protein_fragment=_sample_fragment(),
        epitopes=[_sample_epitope()],
    )
    restored = pickle.loads(pickle.dumps(vp))
    assert restored.mutant_protein_fragment == vp.mutant_protein_fragment
    assert restored.combined_score == pytest.approx(vp.combined_score)
