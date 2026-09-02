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

from dataclasses import fields

from mhctools import RandomBindingPredictor
from topiary import FRAGMENTS

from vaxrank.cli import make_vaxrank_arg_parser, run_vaxrank_from_parsed_args
from vaxrank.mutant_protein_fragment import (
    RNA_EVIDENCE_METHOD_ALIGNMENT,
    SEQUENCE_SOURCE_ISOVAR_ASSEMBLY,
    MutantProteinFragment,
)

from .common import eq_, almost_eq_
from .testing_helpers import data_path

random_binding_predictor = RandomBindingPredictor(["H-2-Kb", "H-2-Db"])


def test_generate_subsequences_preserves_evidence_and_provenance():
    from varcode import Variant

    fragment = MutantProteinFragment(
        variant=Variant("1", 1000, "A", "G", ensembl=None),
        gene_name="GENE1",
        amino_acids="ABCDEFGHIJKLMNOPQRSTUVWXYZ",
        mutant_amino_acid_start_offset=12,
        mutant_amino_acid_end_offset=13,
        supporting_reference_transcripts=["ENST000001"],
        n_overlapping_reads=101,
        n_alt_reads=31,
        n_ref_reads=61,
        n_alt_reads_supporting_protein_sequence=29,
        placeholder_alleles=True,
        rna_evidence_method="rna_alignment",
        rna_evidence_subject="fragments",
        n_overlapping_fragments=55,
        n_alt_fragments=21,
        n_ref_fragments=30,
        n_alt_fragments_supporting_protein_sequence=20,
        n_dna_overlapping=80,
        n_dna_alt=32,
        n_dna_ref=48,
        dna_vaf=0.4,
        dna_evidence_method="dna_depth_x_vaf",
        dna_evidence_subject="reads",
        sequence_source="isovar_assembly")

    offset, subsequence = list(fragment.generate_subsequences(15))[5]

    assert offset == 5
    assert subsequence.amino_acids == fragment.amino_acids[5:20]
    assert subsequence.mutant_amino_acid_start_offset == 7
    assert subsequence.mutant_amino_acid_end_offset == 8

    changed_fields = {
        "amino_acids",
        "mutant_amino_acid_start_offset",
        "mutant_amino_acid_end_offset",
    }
    for field in fields(MutantProteinFragment):
        if field.name not in changed_fields:
            assert getattr(subsequence, field.name) == getattr(fragment, field.name)


def check_mutant_amino_acids(variant, mutant_protein_fragment):
    predicted_effect = mutant_protein_fragment.predicted_effect()
    expected_amino_acids = predicted_effect.aa_alt
    vaxrank_mutant_amino_acids = mutant_protein_fragment.amino_acids[
        mutant_protein_fragment.mutant_amino_acid_start_offset:
        mutant_protein_fragment.mutant_amino_acid_end_offset]

    eq_(expected_amino_acids, vaxrank_mutant_amino_acids,
        "Expected amino acids '%s' for %s but got '%s' from vaxrank in '%s' %d:%d" % (
            expected_amino_acids,
            predicted_effect,
            vaxrank_mutant_amino_acids,
            mutant_protein_fragment.amino_acids,
            mutant_protein_fragment.mutant_amino_acid_start_offset,
            mutant_protein_fragment.mutant_amino_acid_end_offset))
    assert all(
        t.gene.name in variant.gene_names
        for t in
        mutant_protein_fragment.supporting_reference_transcripts), \
        "Wrong gene names for %s" % (mutant_protein_fragment.supporting_reference_transcripts,)


def check_native_rna_evidence(
        fragment,
        *,
        reads,
        fragments,
        protein_support):
    """Check native Isovar evidence after vaccine-peptide slicing."""
    assert (
        fragment.n_overlapping_reads,
        fragment.n_alt_reads,
        fragment.n_ref_reads,
    ) == reads
    assert (
        fragment.n_overlapping_fragments,
        fragment.n_alt_fragments,
        fragment.n_ref_fragments,
    ) == fragments
    assert (
        fragment.n_alt_reads_supporting_protein_sequence,
        fragment.n_alt_fragments_supporting_protein_sequence,
    ) == protein_support
    assert (
        fragment.n_rna_overlapping,
        fragment.n_rna_alt,
        fragment.n_rna_ref,
        fragment.n_rna_supporting_protein_sequence,
    ) == (
        fragments[0],
        fragments[1],
        fragments[2],
        protein_support[1],
    )
    assert fragment.rna_evidence_method == RNA_EVIDENCE_METHOD_ALIGNMENT
    assert fragment.rna_evidence_subject == FRAGMENTS
    assert fragment.sequence_source == SEQUENCE_SOURCE_ISOVAR_ASSEMBLY


def test_mutant_amino_acids_in_mm10_chrX_8125624_refC_altA_pS460I():
    # there are two co-occurring variants in the RNAseq data but since
    # they don't happen in the same codon then we're considering the Varcode
    # annotation to be correct
    # TODO:
    #  deal with phasing of variants explicitly so that both
    #  variant positions are considered mutated
    arg_parser = make_vaxrank_arg_parser()
    args = arg_parser.parse_args([
        "--vcf", data_path("b16.f10/b16.f10.Wdr13.vcf"),
        "--bam", data_path("b16.f10/b16.combined.sorted.bam"),
        "--vaccine-peptide-length", "15",
        "--padding-around-mutation", "5",
        "--mhc-predictor", "random",
        "--mhc-alleles", "HLA-A*02:01",
    ])
    results = run_vaxrank_from_parsed_args(args)
    ranked_list = results.ranked_vaccine_peptides

    for variant, vaccine_peptides in ranked_list:
        eq_(
            1,
            len(vaccine_peptides),
            "Expected 1 vaccine peptide for variant '%s' but got %d" % (
                variant, len(vaccine_peptides)))
        vaccine_peptide = vaccine_peptides[0]
        mutant_protein_fragment = vaccine_peptide.mutant_protein_fragment
        check_mutant_amino_acids(variant, mutant_protein_fragment)
        check_native_rna_evidence(
            mutant_protein_fragment,
            reads=(73, 25, 48),
            fragments=(65, 21, 44),
            protein_support=(25, 21))

def test_mutant_amino_acids_in_mm10_chr9_82927102_refGT_altTG_pT441H():
    # In the Isovar repository this test is weird because the VCF only
    # mentions the G>T variant but doesn't include the subsequent nucleotide
    # change T>G. To avoid having to think about phasing of variants I changed
    # the VCF in vaxrank to contain a GT>TG variant.
    arg_parser = make_vaxrank_arg_parser()
    args = arg_parser.parse_args([
        "--vcf", data_path("b16.f10/b16.f10.Phip.vcf"),
        "--bam", data_path("b16.f10/b16.combined.sorted.bam"),
        "--vaccine-peptide-length", "15",
        "--padding-around-mutation", "5",
        "--mhc-predictor", "random",
        "--mhc-alleles", "HLA-A*02:01",
    ])
    results = run_vaxrank_from_parsed_args(args)
    ranked_list = results.ranked_vaccine_peptides

    for variant, vaccine_peptides in ranked_list:
        vaccine_peptide = vaccine_peptides[0]
        mutant_protein_fragment = vaccine_peptide.mutant_protein_fragment
        check_mutant_amino_acids(
            variant,
            mutant_protein_fragment)
        check_native_rna_evidence(
            mutant_protein_fragment,
            reads=(56, 17, 39),
            fragments=(52, 15, 37),
            protein_support=(17, 15))

def test_keep_top_k_epitopes():
    arg_parser = make_vaxrank_arg_parser()
    keep_k_epitopes = 3
    args = arg_parser.parse_args([
        "--vcf", data_path("b16.f10/b16.f10.Phip.vcf"),
        "--bam", data_path("b16.f10/b16.combined.sorted.bam"),
        "--vaccine-peptide-length", "15",
        "--padding-around-mutation", "5",
        "--num-epitopes-per-vaccine-peptide", str(keep_k_epitopes),
        "--mhc-predictor", "netmhc",
        "--mhc-alleles", "HLA-A*02:01",
    ])
    results = run_vaxrank_from_parsed_args(args)

    ranked_list = results.ranked_vaccine_peptides

    from tests._legacy_score_reference import legacy_score_one as _legacy_score_one
    for variant, vaccine_peptides in ranked_list:
        vaccine_peptide = vaccine_peptides[0]
        eq_(keep_k_epitopes, len(vaccine_peptide.target_epitopes))
        # Recompute the expected score and verify the top-k slice from
        # ``ranked_vaccine_peptides()`` reached the per-CandidateEpitope sum.
        # Each CandidateEpitope contributes the sum of its per-allele leaf
        # ``pMHC_affinity`` records' logistic score.
        target_epitope_score = sum(
            _legacy_score_one(leaf.value, leaf.percentile_rank)
            for e in vaccine_peptide.target_epitopes
            for leaf in e.predictions_flat()
            if leaf.kind == 'pMHC_affinity')
        almost_eq_(target_epitope_score, vaccine_peptide.target_epitope_score)

def test_rna_vaf_property():
    from varcode import Variant
    fragment = MutantProteinFragment(
        variant=Variant('X', '8125624', 'C', 'A'),
        gene_name='Wdr13',
        amino_acids='KLQGHSAPVLDVIVNCDESLLASSD',
        mutant_amino_acid_start_offset=12,
        mutant_amino_acid_end_offset=13,
        n_overlapping_reads=80,
        n_alt_reads=25,
        n_ref_reads=75,
        n_alt_reads_supporting_protein_sequence=2,
        supporting_reference_transcripts=[])
    eq_(0.25, fragment.rna_vaf)

    zero_fragment = MutantProteinFragment(
        variant=Variant('X', '8125624', 'C', 'A'),
        gene_name='Wdr13',
        amino_acids='KLQGHSAPVLDVIVNCDESLLASSD',
        mutant_amino_acid_start_offset=12,
        mutant_amino_acid_end_offset=13,
        n_overlapping_reads=0,
        n_alt_reads=0,
        n_ref_reads=0,
        n_alt_reads_supporting_protein_sequence=0,
        supporting_reference_transcripts=[])
    assert zero_fragment.rna_vaf is None


def test_mutant_protein_fragment_serialization():
    arg_parser = make_vaxrank_arg_parser()
    keep_k_epitopes = 3
    args = arg_parser.parse_args([
        "--vcf", data_path("b16.f10/b16.f10.Phip.vcf"),
        "--bam", data_path("b16.f10/b16.combined.sorted.bam"),
        "--vaccine-peptide-length", "15",
        "--padding-around-mutation", "5",
        "--num-epitopes-per-vaccine-peptide", str(keep_k_epitopes),
        "--mhc-predictor", "netmhc",
        "--mhc-alleles", "HLA-A*02:01",
    ])
    results = run_vaxrank_from_parsed_args(args)

    ranked_list = results.ranked_vaccine_peptides

    for _, vaccine_peptides in ranked_list:
        mutant_protein_fragment = vaccine_peptides[0].mutant_protein_fragment
        json_str = mutant_protein_fragment.to_json()
        deserialized = MutantProteinFragment.from_json(json_str)
        eq_(mutant_protein_fragment, deserialized)
