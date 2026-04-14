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

"""Regression tests for VaxrankResults.variant_counts() — pins the
nested-subset invariant described by the patient-info report labels
(total ⊇ coding ⊇ coding+RNA ⊇ coding+RNA+ligand)."""

from types import SimpleNamespace

from vaxrank.vaxrank_results import VaxrankResults


class _FakeVariant:
    """Identity-hashable stand-in for varcode.Variant (SimpleNamespace
    is unhashable, which breaks the `variant in dict` lookup in
    variant_properties())."""
    def __init__(self):
        self.contig = "1"
        self.start = 1
        self.ref = "A"
        self.alt = "G"


def _isovar_result(*, coding, rna, gene_name="GENE1"):
    variant = _FakeVariant()
    return variant, SimpleNamespace(
        variant=variant,
        top_gene_name=gene_name,
        predicted_effect_modifies_protein_sequence=coding,
        has_mutant_protein_sequence_from_rna=rna,
    )


def _results(isovar_results, variants_with_peptides=()):
    return VaxrankResults(
        isovar_results=isovar_results,
        variant_to_vaccine_peptides_dict={v: [] for v in variants_with_peptides},
        ranked_vaccine_peptides=[],
    )


def test_variant_counts_nested_subsets():
    # Four variants covering every coding/rna/peptide combination so the
    # AND-chain matters for each step of the nesting.
    v1_variant, v1 = _isovar_result(coding=True, rna=True)   # coding + RNA + peptide
    v2_variant, v2 = _isovar_result(coding=True, rna=True)   # coding + RNA, no peptide
    v3_variant, v3 = _isovar_result(coding=True, rna=False)  # coding only
    v4_variant, v4 = _isovar_result(coding=False, rna=True)  # RNA but not coding

    results = _results([v1, v2, v3, v4], variants_with_peptides=[v1_variant])
    counts = results.variant_counts()

    assert counts["num_total_variants"] == 4
    assert counts["num_coding_effect_variants"] == 3
    assert counts["num_variants_with_rna_support"] == 2
    assert counts["num_variants_with_vaccine_peptides"] == 1


def test_variant_counts_non_coding_rna_does_not_inflate():
    # Regression for #205: a non-coding variant with RNA support must
    # not be counted in num_variants_with_rna_support.
    _, v_coding = _isovar_result(coding=True, rna=False)
    _, v_noncoding_rna = _isovar_result(coding=False, rna=True)

    counts = _results([v_coding, v_noncoding_rna]).variant_counts()

    assert counts["num_coding_effect_variants"] == 1
    assert counts["num_variants_with_rna_support"] == 0
    assert counts["num_variants_with_rna_support"] <= counts["num_coding_effect_variants"]


def test_variant_counts_empty():
    counts = _results([]).variant_counts()
    assert counts == {
        "num_total_variants": 0,
        "num_coding_effect_variants": 0,
        "num_variants_with_rna_support": 0,
        "num_variants_with_vaccine_peptides": 0,
    }
