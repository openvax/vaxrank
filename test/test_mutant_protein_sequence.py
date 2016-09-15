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
import logging

from nose.tools import eq_
from vaxrank.core_logic import ranked_vaccine_peptides
from mhctools import RandomBindingPredictor
from isovar.cli.variant_sequences import make_variant_sequences_arg_parser
from isovar.cli.rna_reads import allele_reads_generator_from_args

from .common import data_path

def check_vaxrank_agrees_with_varcode(variant, mutant_protein_fragment):
    predicted_effect = variant.effects().top_priority_effect()

    varcode_mutant_amino_acids = predicted_effect.aa_alt

    vaxrank_mutant_amino_acids = mutant_protein_fragment.amino_acids[
        mutant_protein_fragment.mutant_amino_acid_start_offset:
        mutant_protein_fragment.mutant_amino_acid_end_offset]

    eq_(varcode_mutant_amino_acids, vaxrank_mutant_amino_acids,
        "Expected amino acids '%s' for %s but got '%s' from vaxrank in '%s' %d:%d" % (
            varcode_mutant_amino_acids,
            predicted_effect,
            vaxrank_mutant_amino_acids,
            mutant_protein_fragment.amino_acids,
            mutant_protein_fragment.mutant_amino_acid_start_offset,
            mutant_protein_fragment.mutant_amino_acid_end_offset))


def test_mutant_amino_acids_agree_with_varcode():
    logging.basicConfig(level=logging.WARN)

    bam_path = data_path("b16.f10/b16.combined.sorted.bam")
    vcf_path = data_path("b16.f10/b16.vcf")
    arg_parser = make_variant_sequences_arg_parser()
    args_list = [
        "--vcf", vcf_path,
        "--bam", bam_path
    ]
    args = arg_parser.parse_args(args_list)
    reads_generator = allele_reads_generator_from_args(args)
    ranked_list = ranked_vaccine_peptides(
        reads_generator=reads_generator,
        mhc_predictor=RandomBindingPredictor(["H-2-Kb", "H-2-Db"]),
        vaccine_peptide_length=15,
        padding_around_mutation=5,
        max_vaccine_peptides_per_variant=1,
        min_reads_supporting_cdna_sequence=1)

    for variant, vaccine_peptides in ranked_list:
        if len(vaccine_peptides) == 0:
            continue
        vaccine_peptide = vaccine_peptides[0]
        mutant_protein_fragment = vaccine_peptide.mutant_protein_fragment

        yield (
            check_vaxrank_agrees_with_varcode,
            variant,
            mutant_protein_fragment
        )
