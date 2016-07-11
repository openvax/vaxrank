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

import roman


def ascii_report_from_ranked_vaccine_peptides(
        ranked_variants_with_vaccine_peptides,
        mhc_alleles,
        variants,
        bam_path):
    # create dictionary mapping variants to coding effects
    variants_to_top_coding_effect_dict = {
        variant: effect_collection.top_priority_effect()
        for (variant, effect_collection)
        in variants.effects().drop_silent_and_noncoding().groupby_variant().items()
    }
    lines = [
        "VCF (somatic variants) path: %s" % variants.path,
        "BAM (RNAseq reads) path: %s" % bam_path,
        "MHC alleles: %s" % (" ".join(mhc_alleles)),
        "Total number of somatic variants: %d" % (len(variants),),
        "Somatic variants with predicted coding effects: %d" % (
            len(variants_to_top_coding_effect_dict),),
        "---",
    ]

    for i, (variant, vaccine_peptides) in enumerate(
            ranked_variants_with_vaccine_peptides):
        variant_short_description = variant.short_description
        if len(vaccine_peptides) == 0:
            logging.info(
                "Skipping %s, no vaccine peptides" % variant_short_description)
            continue

        gene_name = vaccine_peptides[0].mutant_protein_fragment.gene_name
        lines.append("\n%d) %s (%s)" % (
            i + 1,
            variant_short_description,
            gene_name,))
        lines.append(
            "\tTop score: %0.2f" % (vaccine_peptides[0].combined_score))
        lines.append(
            "\tPredicted effect: %s" % (
                variants_to_top_coding_effect_dict.get(variant)))
        lines.append(
            "\tReads supporting variant allele: %d" % (
                vaccine_peptides[0].mutant_protein_fragment.n_alt_reads))
        lines.append(
            "\tReads supporting reference allele: %d" % (
                vaccine_peptides[0].mutant_protein_fragment.n_ref_reads))
        lines.append(
            "\tReads supporting other alleles: %d" % (
                vaccine_peptides[0].mutant_protein_fragment.n_other_reads))

        lines.append("\tVaccine Peptides:")
        for j, vaccine_peptide in enumerate(vaccine_peptides):
            mutant_protein_fragment = vaccine_peptide.mutant_protein_fragment
            amino_acids = mutant_protein_fragment.amino_acids
            mutation_start = mutant_protein_fragment.mutant_amino_acid_start_offset
            mutation_end = mutant_protein_fragment.mutant_amino_acid_end_offset
            aa_before_mutation = amino_acids[:mutation_start]
            aa_mutant = amino_acids[mutation_start:mutation_end]
            aa_after_mutation = amino_acids[mutation_end:]
            lines.append("\t\t%s. %s_%s_%s (score = %0.2f)" % (
                roman.toRoman(j + 1).lower(),
                aa_before_mutation,
                aa_mutant,
                aa_after_mutation,
                vaccine_peptide.combined_score))
            lines.append(
                "\t\t   - Length: %d" % len(amino_acids))
            lines.append(
                "\t\t   - Expression score: %0.2f" % (vaccine_peptide.expression_score))
            lines.append(
                "\t\t   - Mutant epitope score: %0.2f" % (
                    vaccine_peptide.mutant_epitope_score))
            lines.append(
                "\t\t   - Wildtype epitope score: %0.2f" % (
                    vaccine_peptide.wildtype_epitope_score))

            lines.append(
                "\t\t   - Reads fully spanning cDNA sequence(s): %d" % (
                    mutant_protein_fragment.n_alt_reads_supporting_protein_sequence))
            lines.append(
                "\t\t   - Mutant amino acids: %d" % (
                    mutant_protein_fragment.n_mutant_amino_acids))
            lines.append(
                "\t\t   - Mutation distance from edge: %d" % (
                    mutant_protein_fragment.mutation_distance_from_edge))
            lines.append("\t\t   - Predicted mutant epitopes:")
            for epitope_prediction in vaccine_peptide.epitope_predictions:
                if epitope_prediction.overlaps_mutation and epitope_prediction.ic50 <= 2000:
                    lines.append("\t\t\t * %s %0.2f (%s)" % (
                        epitope_prediction.peptide_sequence,
                        epitope_prediction.ic50,
                        epitope_prediction.allele))
    return "\n".join(lines)
