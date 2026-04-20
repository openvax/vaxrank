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
from varcode import Variant

from vaxrank.epitope_prediction import EpitopePrediction
from vaxrank.mutant_protein_fragment import MutantProteinFragment
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
):
    return EpitopePrediction(
        allele="HLA-A*02:01",
        peptide_sequence=peptide_sequence,
        wt_peptide_sequence="SIINFEKL",
        ic50=ic50,
        wt_ic50=2000.0,
        percentile_rank=0.3,
        prediction_method_name="ImaginationMHCpan",
        overlaps_mutation=overlaps_mutation,
        source_sequence="SSIINFEQL",
        offset=1,
        occurs_in_reference=occurs_in_reference,
    )


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
        self.n_alt_reads = n_alt_reads
        self.n_alt_reads_supporting_protein_sequence = n_alt_reads
        self.n_mutant_amino_acids = 1
        self.mutation_distance_from_edge = 5


# ---- EpitopePrediction -------------------------------------------------


def test_epitope_prediction_json_roundtrip():
    e = _sample_epitope()
    restored = EpitopePrediction.from_json(e.to_json())
    assert restored == e


def test_epitope_prediction_pickle_roundtrip():
    e = _sample_epitope()
    assert pickle.loads(pickle.dumps(e)) == e


def test_epitope_prediction_length_is_derived():
    """`length` moved from init-arg to @property in the migration; it
    still matches len(peptide_sequence) and doesn't appear in to_dict."""
    e = _sample_epitope(peptide_sequence="SIINFEQLL")
    assert e.length == 9
    assert "length" not in e.to_dict()


def test_epitope_prediction_from_dict_drops_legacy_length_field():
    """Old JSON blobs that include the retired `length` init arg still
    load, via _SERIALIZABLE_KEYWORD_ALIASES."""
    legacy_dict = {
        "allele": "HLA-A*02:01",
        "peptide_sequence": "SIINFEQL",
        "wt_peptide_sequence": "SIINFEKL",
        "ic50": 2.0,
        "wt_ic50": 2000.0,
        "percentile_rank": 0.3,
        "prediction_method_name": "ImaginationMHCpan",
        "overlaps_mutation": True,
        "source_sequence": "SSIINFEQL",
        "offset": 1,
        "occurs_in_reference": False,
        "length": 8,  # legacy field, should be silently dropped
    }
    e = EpitopePrediction.from_dict(legacy_dict)
    assert e.peptide_sequence == "SIINFEQL"
    assert e.length == 8  # derived from peptide_sequence


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
    """Filtering and sorting moved from __init__ to __post_init__ but must
    produce the same mutant / wildtype split."""
    mutant = _sample_epitope(peptide_sequence="MUTANT", ic50=5.0)
    wildtype = _sample_epitope(
        peptide_sequence="REFERNC",
        ic50=50.0,
        overlaps_mutation=False,
        occurs_in_reference=True,
    )
    vp = VaccinePeptide(
        mutant_protein_fragment=_FakeFragment(),
        epitope_predictions=[mutant, wildtype],
    )
    assert vp.mutant_epitope_predictions == [mutant]
    assert vp.wildtype_epitope_predictions == [wildtype]


def test_vaccine_peptide_post_init_validates_combined_score_mode():
    """__post_init__ now raises on unknown modes — the contract added in
    2.5.0 carries through the dataclass migration."""
    with pytest.raises(ValueError, match="combined_score_mode"):
        VaccinePeptide(
            mutant_protein_fragment=_FakeFragment(),
            epitope_predictions=[],
            combined_score_mode="garbage",
        )


def test_vaccine_peptide_to_dict_omits_derived_fields():
    """The custom to_dict only emits constructor args, not the
    init=False dataclass fields (mutant_epitope_predictions, etc.)."""
    vp = VaccinePeptide(
        mutant_protein_fragment=_FakeFragment(),
        epitope_predictions=[_sample_epitope()],
    )
    d = vp.to_dict()
    for derived in (
        "mutant_epitope_predictions",
        "wildtype_epitope_predictions",
        "mutant_epitope_score",
        "wildtype_epitope_score",
        "manufacturability_scores",
    ):
        assert derived not in d


def test_vaccine_peptide_tuple_coercion_on_rules():
    """List of rule names coerced to tuple even after migration."""
    vp = VaccinePeptide(
        mutant_protein_fragment=_FakeFragment(),
        epitope_predictions=[_sample_epitope()],
        manufacturability_rules=["cysteine_count", "cterm_hydropathy"],
    )
    assert vp.manufacturability_rules == ("cysteine_count", "cterm_hydropathy")
    assert isinstance(vp.manufacturability_rules, tuple)


def test_vaccine_peptide_json_roundtrip_with_real_fragment():
    """Full to_json / from_json with a real MutantProteinFragment and
    EpitopePredictions. This is the path the CLI writes out — we want
    __post_init__ to re-derive the mutant/wildtype split after load and
    scores to re-compute consistently."""
    mutant = _sample_epitope(peptide_sequence="SIINFEKL", ic50=500.0)
    wildtype = _sample_epitope(
        peptide_sequence="TESTWILD",
        ic50=50.0,
        overlaps_mutation=False,
        occurs_in_reference=True,
    )
    vp = VaccinePeptide(
        mutant_protein_fragment=_sample_fragment(),
        epitope_predictions=[mutant, wildtype],
        manufacturability_rules=("cysteine_count", "cterm_hydropathy"),
        combined_score_mode="epitope_only",
    )
    restored = VaccinePeptide.from_json(vp.to_json())

    # Wire-level invariants
    assert restored.mutant_protein_fragment == vp.mutant_protein_fragment
    assert restored.manufacturability_rules == vp.manufacturability_rules
    assert restored.combined_score_mode == vp.combined_score_mode

    # Derived state must match — __post_init__ re-ran on load
    assert restored.mutant_epitope_predictions == vp.mutant_epitope_predictions
    assert restored.wildtype_epitope_predictions == vp.wildtype_epitope_predictions
    assert restored.mutant_epitope_score == pytest.approx(vp.mutant_epitope_score)
    assert restored.combined_score == pytest.approx(vp.combined_score)


def test_vaccine_peptide_pickle_roundtrip_with_real_fragment():
    vp = VaccinePeptide(
        mutant_protein_fragment=_sample_fragment(),
        epitope_predictions=[_sample_epitope()],
    )
    restored = pickle.loads(pickle.dumps(vp))
    assert restored.mutant_protein_fragment == vp.mutant_protein_fragment
    assert restored.combined_score == pytest.approx(vp.combined_score)
