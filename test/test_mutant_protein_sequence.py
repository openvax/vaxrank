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

from nose.tools import eq_, assert_almost_equal
from vaxrank.cli import make_vaxrank_arg_parser, run_vaxrank_from_parsed_args
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from mhctools import RandomBindingPredictor

from .testing_helpers import data_path

random_binding_predictor = RandomBindingPredictor(["H-2-Kb", "H-2-Db"])


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

    for variant, vaccine_peptides in ranked_list:
        vaccine_peptide = vaccine_peptides[0]
        eq_(keep_k_epitopes, len(vaccine_peptide.mutant_epitope_predictions))
        # recompute the expected score, make sure the top-k argument from ranked_vaccine_peptides()
        # propagated as expected
        mutant_epitope_score = sum(
            p.logistic_epitope_score() for p in vaccine_peptide.mutant_epitope_predictions)
        assert_almost_equal(mutant_epitope_score, vaccine_peptide.mutant_epitope_score)

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
