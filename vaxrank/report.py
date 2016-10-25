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

from __future__ import absolute_import, division
from collections import OrderedDict
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
        bam_path,
        output_values):
    # create dictionary mapping variants to coding effects
    variants_to_top_coding_effect_dict = {
        variant: effect_collection.top_priority_effect()
        for (variant, effect_collection)
        in variants.effects().drop_silent_and_noncoding().groupby_variant().items()
    }

    ############################### Patient info ###############################

    patient_info = OrderedDict()
    patient_info['Patient ID'] = output_values['patient_id']
    patient_info['VCF (somatic variants) path(s)'] = "; ".join(variants.sources)
    patient_info['BAM (RNAseq reads) path'] = bam_path
    patient_info['MHC alleles'] = " ".join(mhc_alleles)
    patient_info['Total number of somatic variants'] = len(variants)
    patient_info['Somatic variants with predicted coding effects'] = len(
        variants_to_top_coding_effect_dict)

    # list ranked variants with their peptides
    variants = []
    for i, (variant, vaccine_peptides) in enumerate(
            ranked_variants_with_vaccine_peptides):
        variant_short_description = variant.short_description
        if len(vaccine_peptides) == 0:
            logger.info("Skipping %s, no vaccine peptides", variant_short_description)
            continue

        ############################### Variant info ###############################

        variant_data = OrderedDict()
        mutant_protein_fragment = vaccine_peptides[0].mutant_protein_fragment
        top_score = round(vaccine_peptides[0].combined_score, 4)
        variant_data['Gene name'] = mutant_protein_fragment.gene_name
        variant_data['Top score'] = top_score
        variant_data['Reads supporting variant allele'] = mutant_protein_fragment.n_alt_reads
        variant_data['Reads supporting reference allele'] = mutant_protein_fragment.n_ref_reads
        variant_data['Reads supporting other alleles'] = mutant_protein_fragment.n_other_reads

        ############################### Predicted effect info ###############################

        predicted_effect = variants_to_top_coding_effect_dict.get(variant)
        effect_data = OrderedDict()
        effect_data['Effect type'] = predicted_effect.__class__.__name__
        effect_data['Transcript name'] = predicted_effect.transcript_name
        effect_data['Transcript ID'] = predicted_effect.transcript_id
        effect_data['Effect description'] = predicted_effect.short_description

        peptides = []
        for j, vaccine_peptide in enumerate(vaccine_peptides):
            mutant_protein_fragment = vaccine_peptide.mutant_protein_fragment
            amino_acids = mutant_protein_fragment.amino_acids
            mutation_start = mutant_protein_fragment.mutant_amino_acid_start_offset
            mutation_end = mutant_protein_fragment.mutant_amino_acid_end_offset
            aa_before_mutation = amino_acids[:mutation_start]
            aa_mutant = amino_acids[mutation_start:mutation_end]
            aa_after_mutation = amino_acids[mutation_end:]

            header_display_data = {
                'num': roman.toRoman(j + 1).lower(),
                'aa_before_mutation': aa_before_mutation,
                'aa_mutant': aa_mutant,
                'aa_after_mutation': aa_after_mutation,
            }

            ############################### Peptide info ###############################

            peptide_data = OrderedDict()
            peptide_data['Transcript name'] = predicted_effect.transcript_name
            peptide_data['Length'] = len(amino_acids)
            peptide_data['Expression score'] = round(vaccine_peptide.expression_score, 4)
            peptide_data['Mutant epitope score'] = round(vaccine_peptide.mutant_epitope_score, 4)
            peptide_data['Combined score'] = round(vaccine_peptide.combined_score, 4)
            peptide_data['Reads fully spanning cDNA sequence(s)'] = \
                mutant_protein_fragment.n_alt_reads_supporting_protein_sequence
            peptide_data['Mutant amino acids'] = mutant_protein_fragment.n_mutant_amino_acids
            peptide_data['Mutation distance from edge'] = \
                mutant_protein_fragment.mutation_distance_from_edge

            ############################### Manufacturability info ###############################
            
            manufacturability_data = OrderedDict()
            scores = vaccine_peptide.manufacturability_scores
            manufacturability_data['C-terminal 7mer GRAVY score'] = \
                    round(scores.cterm_7mer_gravy_score, 4)
            manufacturability_data['Max 7mer GRAVY score'] = round(scores.max_7mer_gravy_score, 4)
            manufacturability_data['N-terminal Glutamine, Glutamic Acid, or Cysteine'] = int(
                scores.difficult_n_terminal_residue)
            manufacturability_data['C-terminal Cysteine'] = int(
                scores.c_terminal_cysteine)
            manufacturability_data['C-terminal Proline'] = int(scores.c_terminal_proline)
            manufacturability_data['Total number of Cysteine residues'] = scores.cysteine_count
            manufacturability_data['N-terminal Asparagine'] = int(scores.n_terminal_asparagine)
            manufacturability_data['Number of Asparagine-Proline bonds'] = \
                scores.asparagine_proline_bond_count
            
            ############################### Epitopes info ###############################

            epitopes = []
            for epitope_prediction in vaccine_peptide.epitope_predictions:
                if epitope_prediction.overlaps_mutation:
                    epitope_data = OrderedDict()
                    epitope_data['Sequence'] = epitope_prediction.peptide_sequence
                    epitope_data['IC50'] = epitope_prediction.ic50
                    epitope_data['Normalized binding score'] = round(
                        epitope_prediction.logistic_score(), 4)
                    epitope_data['Allele'] = epitope_prediction.allele
                    epitopes.append(epitope_data)

            peptide_dict = {
                'header_display_data': header_display_data,
                'peptide_data': peptide_data,
                'manufacturability_data': manufacturability_data,
                'epitopes': epitopes,
            }
            peptides.append(peptide_dict)

        variant_dict = {
            'num': i + 1,
            'short_description': variant_short_description,
            'variant_data': variant_data,
            'effect_data': effect_data,
            'peptides': peptides,
        }
        variants.append(variant_dict)

    template_data = {
        'patient_info': patient_info,
        'variants': variants,
    }
    template_data.update(output_values)
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
