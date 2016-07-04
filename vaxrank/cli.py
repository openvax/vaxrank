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
import sys
import logging

from isovar.args.variant_sequences import make_variant_sequences_arg_parser
from isovar.args.rna_reads import allele_reads_generator_from_args
from isovar.protein_sequences import (
    reads_generator_to_protein_sequences_generator,
)
from topiary.commandline_args.mhc import (
    add_mhc_args,
    mhc_alleles_from_args,
    mhc_binding_predictor_from_args,
)

from .old_vaccine_peptides import generate_candidate_vaccine_peptides
from .mutant_protein_fragment import MutantProteinFragment

# inherit all commandline options from Isovar
arg_parser = make_variant_sequences_arg_parser(
    prog="vaxrank",
    description=(
        "Select personalized vaccine peptides from cancer variants, "
        "expression data, and patient HLA type."),
)

add_mhc_args(arg_parser)

arg_parser.add_argument(
    "--output-csv",
    default="vaccine_peptides.csv",
    help="Name of CSV file which contains predicted sequences")

vaccine_peptide_group = arg_parser.add_argument_group("Vaccine peptide options")
vaccine_peptide_group.add_argument(
    "--vaccine-peptide-length",
    default=25,
    type=int,
    help="Number of amino acids in the vaccine peptides (default %(default)s)")

vaccine_peptide_group.add_argument(
    "--padding-around-mutation",
    default=0,
    type=int,
    help=(
        "Number of off-center windows around the mutation to consider "
        "as vaccine peptides (default %(default)s)"
    ))

def main(args_list=None):
    """
    Script to generate vaccine peptide predictions from somatic cancer variants,
    patient HLA type, and tumor RNA-seq data.

    Example usage:
        vaxrank
            --vcf somatic.vcf \
            --bam rnaseq.bam \
            --vaccine-peptide-length 25 \
            --output-csv vaccine-peptides.csv
    """
    if args_list is None:
        args_list = sys.argv[1:]

    logging.basicConfig(level=logging.DEBUG)
    args = arg_parser.parse_args(args_list)
    print(args)

    mhc_alleles = mhc_alleles_from_args(args)
    print("MHC alleles: %s" % (mhc_alleles,))
    mhc_predictor = mhc_binding_predictor_from_args(args)

    # generator that for each variant gathers all RNA reads, both those
    # supporting the variant and reference alleles
    reads_generator = allele_reads_generator_from_args(args)

    # total number of amino acids is the vaccine peptide length plus the
    # number of off-center windows around the mutation
    protein_fragment_sequence_length = (
        args.vaccine_peptide_length + 2 * args.padding_around_mutation)

    protein_sequences_generator = reads_generator_to_protein_sequences_generator(
        reads_generator,
        transcript_id_whitelist=None,
        protein_sequence_length=protein_fragment_sequence_length,
        min_reads_supporting_cdna_sequence=args.min_reads_supporting_variant_sequence,
        max_protein_sequences_per_variant=1)

    variant_to_mutant_protein_fragments = {}

    for variant, isovar_protein_sequences in protein_sequences_generator:
        isovar_protein_sequences = list(isovar_protein_sequences)
        if len(isovar_protein_sequences) == 0:
            logging.info("No protein sequences for %s" % (variant,))
            continue
        isovar_protein_sequence = isovar_protein_sequences[0]
        mutant_protein_fragment = MutantProteinFragment.from_isovar_protein_sequence(
            isovar_protein_sequence)
        variant_to_mutant_protein_fragments[variant] = mutant_protein_fragment
        epitope_predictions = mhc_predictor.predict(variant_to_amino_acid_sequences_dict)

    for variant, protein_sequence in variant_to_protein_sequence_objects_dict.items():
        generate_candidate_vaccine_peptides(
            protein_sequence.amino_acid,
            epitopes=epitope_predictions,
            mutation_start=protein_sequence.variant_aa_interval_start,
            mutation_end=protein_sequence.variant_aa_interval_end,
            epitope_scorer=None,
            result_length=25,
            padding=5)

    """
    df_vaccine_peptides = vaccine_peptides_dataframe_from_args(args)
    print(df_vaccine_peptides)
    if args.output_csv:
        df_vaccine_peptides.to_csv(args.output_csv)
    """
