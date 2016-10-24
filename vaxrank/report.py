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
import os
import tempfile

import jinja2
import pdfkit
import roman


logger = logging.getLogger(__name__)


JINJA_ENVIRONMENT = jinja2.Environment(
    loader=jinja2.FileSystemLoader(os.path.dirname(__file__)),
    extensions=['jinja2.ext.autoescape'],
    autoescape=False,
    trim_blocks=True,
    lstrip_blocks=True,
)


def compute_template_data(
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

    template_data = {
        'vcf_paths': "; ".join(variants.sources),
        'bam_path': bam_path,
        'mhc_alleles': " ".join(mhc_alleles),
        'num_total_somatic_variants': len(variants),
        'num_somatic_variants_with_predicted_coding_effects':
            len(variants_to_top_coding_effect_dict),
    }

    # list ranked variants with their peptides
    variants = []
    for i, (variant, vaccine_peptides) in enumerate(
            ranked_variants_with_vaccine_peptides):
        variant_short_description = variant.short_description
        if len(vaccine_peptides) == 0:
            logger.info("Skipping %s, no vaccine peptides", variant_short_description)
            continue
        gene_name = vaccine_peptides[0].mutant_protein_fragment.gene_name
        peptides = []
        variant_dict = {
            'num': i + 1,
            'short_description': variant_short_description,
            'gene_name': gene_name,
            'top_score': vaccine_peptides[0].combined_score,
            'predicted_effect': variants_to_top_coding_effect_dict.get(variant),
            'reads_supporting_variant_allele':
                vaccine_peptides[0].mutant_protein_fragment.n_alt_reads,
            'reads_supporting_reference_allele':
                vaccine_peptides[0].mutant_protein_fragment.n_ref_reads,
            'reads_supporting_other_alleles':
                vaccine_peptides[0].mutant_protein_fragment.n_other_reads,
            'peptides': peptides,
        }

        # compile peptide info
        for j, vaccine_peptide in enumerate(vaccine_peptides):
            epitopes = []
            # compile epitope info
            for epitope_prediction in vaccine_peptide.epitope_predictions:
                if epitope_prediction.overlaps_mutation:
                    score = epitope_prediction.logistic_score()
                    epitope_dict = {
                        'sequence': epitope_prediction.peptide_sequence,
                        'ic50': epitope_prediction.ic50,
                        'normalized_binding_score': score,
                        'allele': epitope_prediction.allele,
                    }
                    epitopes.append(epitope_dict)

            mutant_protein_fragment = vaccine_peptide.mutant_protein_fragment
            amino_acids = mutant_protein_fragment.amino_acids
            mutation_start = mutant_protein_fragment.mutant_amino_acid_start_offset
            mutation_end = mutant_protein_fragment.mutant_amino_acid_end_offset
            aa_before_mutation = amino_acids[:mutation_start]
            aa_mutant = amino_acids[mutation_start:mutation_end]
            aa_after_mutation = amino_acids[mutation_end:]
            manufacturability_scores = vaccine_peptide.manufacturability_scores
            peptide_dict = {
                'num': roman.toRoman(j + 1).lower(),
                'aa_before_mutation': aa_before_mutation,
                'aa_mutant': aa_mutant,
                'aa_after_mutation': aa_after_mutation,
                'score': vaccine_peptide.combined_score,
                'length': len(amino_acids),
                'expression_score': vaccine_peptide.expression_score,
                'mutant_epitope_score': vaccine_peptide.mutant_epitope_score,
                'wildtype_epitope_score': vaccine_peptide.wildtype_epitope_score,
                'num_alt_reads_supporting_protein_sequence':
                    mutant_protein_fragment.n_alt_reads_supporting_protein_sequence,
                'num_mutant_amino_acids':
                    mutant_protein_fragment.n_mutant_amino_acids,
                'mutation_distance_from_edge':
                    mutant_protein_fragment.mutation_distance_from_edge,
                'epitopes': epitopes,
                'difficult_n_terminal_residue':
                    int(manufacturability_scores.difficult_n_terminal_residue),
                'c_terminal_cysteine':
                    int(manufacturability_scores.c_terminal_cysteine),
                'c_terminal_proline':
                    int(manufacturability_scores.c_terminal_proline),
                'n_terminal_asparagine':
                    int(manufacturability_scores.n_terminal_asparagine),
                'asparagine_proline_bond_count':
                    manufacturability_scores.asparagine_proline_bond_count,
                'cysteine_count':
                    manufacturability_scores.cysteine_count,
                'cterm_7mer_gravy_score':
                    manufacturability_scores.cterm_7mer_gravy_score,
                'max_7mer_gravy_score':
                    manufacturability_scores.max_7mer_gravy_score
            }

            peptides.append(peptide_dict)
        variants.append(variant_dict)

    template_data['variants'] = variants
    return template_data


def _make_report(
        template_data,
        file_handle,
        template_path):
    template = JINJA_ENVIRONMENT.get_template(template_path)
    report = template.render(template_data)
    file_handle.write(report)

def make_ascii_report(
        template_data,
        ascii_report_path):
    with open(ascii_report_path, "w") as f:
        _make_report(template_data, f, 'templates/template.txt')
    logger.info('Wrote ASCII report to %s', ascii_report_path)


def make_html_report(
        template_data,
        html_report_path):
    with open(html_report_path, "w") as f:
        _make_report(template_data, f, 'templates/template.html')
    logger.info('Wrote HTML report to %s', html_report_path)

def make_pdf_report(
        template_data,
        pdf_report_path):
    with tempfile.NamedTemporaryFile(mode='w', suffix='.html') as f:
        _make_report(template_data, f, 'templates/template.html')
        f.flush()

        options = {
            'zoom': 0.6,
            'margin-top': '20mm'
        }
        pdfkit.from_file(f.name, pdf_report_path, options=options)
    logger.info('Wrote PDF report to %s', pdf_report_path)
