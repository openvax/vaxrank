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
from collections import namedtuple, OrderedDict
from copy import copy
import logging
import os
import tempfile

import jinja2
import pandas as pd
import pdfkit
import roman
from varcode.effects import top_priority_effect

from .manufacturability import ManufacturabilityScores


logger = logging.getLogger(__name__)


JINJA_ENVIRONMENT = jinja2.Environment(
    loader=jinja2.FileSystemLoader(os.path.dirname(__file__)),
    extensions=['jinja2.ext.autoescape'],
    autoescape=False,
    trim_blocks=True,
    lstrip_blocks=True,
)


PatientInfo = namedtuple(
    "PatientInfo", (
        "patient_id",
        "vcf_paths",
        "bam_path",
        "mhc_alleles",
        "num_somatic_variants",
        "num_coding_effect_variants",
))


class TemplateDataCreator(object):
    def __init__(
            self,
            ranked_variants_with_vaccine_peptides,
            patient_info,
            final_review,
            reviewers,
            args_for_report,
            input_json_file):
        """
        Construct a TemplateDataCreator object, from the output of the vaxrank pipeline.
        """
        self.ranked_variants_with_vaccine_peptides = ranked_variants_with_vaccine_peptides
        self.patient_info = patient_info

        # filter output-related command-line args: we want to display everything else
        args_to_display_in_report = {
            k: v for k, v in args_for_report.items() if not k.startswith("output")
        }

        self.template_data = {
            'args': sorted(args_to_display_in_report.items()),
            'reviewers': reviewers.split(','),
            'final_review': final_review,
            'input_json_file': input_json_file,
        }

    def _patient_info(self):
        """
        Returns an OrderedDict with patient info.
        """
        patient_info = OrderedDict([
            ('Patient ID', self.patient_info.patient_id),
            ('VCF (somatic variants) path(s)', '; '.join(self.patient_info.vcf_paths)),
            ('BAM (RNAseq reads) path', self.patient_info.bam_path),
            ('MHC alleles', ' '.join(self.patient_info.mhc_alleles)),
            ('Total number of somatic variants', self.patient_info.num_somatic_variants),
            ('Somatic variants with predicted coding effects',
                self.patient_info.num_coding_effect_variants),
        ])
        return patient_info

    def _variant_data(self, top_vaccine_peptide):
        """
        Returns an OrderedDict with info used to populate variant info section.
        """
        variant_data = OrderedDict()
        mutant_protein_fragment = top_vaccine_peptide.mutant_protein_fragment
        top_score = round(top_vaccine_peptide.combined_score, 4)
        variant_data = OrderedDict([
            ('Gene name', mutant_protein_fragment.gene_name),
            ('Top score', top_score),
            ('Reads supporting variant allele', mutant_protein_fragment.n_alt_reads),
            ('Reads supporting reference allele', mutant_protein_fragment.n_ref_reads),
            ('Reads supporting other alleles', mutant_protein_fragment.n_other_reads),
        ])
        return variant_data

    def _effect_data(self, predicted_effect):
        """
        Returns an OrderedDict with info about the given predicted effect.
        """
        effect_data = OrderedDict([
            ('Effect type', predicted_effect.__class__.__name__),
            ('Transcript name',predicted_effect.transcript_name),
            ('Transcript ID', predicted_effect.transcript_id),
            ('Effect description', predicted_effect.short_description),
        ])
        return effect_data

    def _peptide_header_display_data(self, vaccine_peptide, rank):
        """
        Returns a dictionary with info used to populate the header section of a peptide table.

        Parameters
        ----------
        vaccine_peptide : VaccinePeptide
          The given peptide to convert to display form

        rank : int
          Rank of vaccine peptide in list
        """
        mutant_protein_fragment = vaccine_peptide.mutant_protein_fragment
        amino_acids = mutant_protein_fragment.amino_acids
        mutation_start = mutant_protein_fragment.mutant_amino_acid_start_offset
        mutation_end = mutant_protein_fragment.mutant_amino_acid_end_offset
        aa_before_mutation = amino_acids[:mutation_start]
        aa_mutant = amino_acids[mutation_start:mutation_end]
        aa_after_mutation = amino_acids[mutation_end:]

        header_display_data = {
            'num': roman.toRoman(rank + 1).lower(),
            'aa_before_mutation': aa_before_mutation,
            'aa_mutant': aa_mutant,
            'aa_after_mutation': aa_after_mutation,
        }
        return header_display_data

    def _peptide_data(self, vaccine_peptide, transcript_name):
        """
        Returns a dictionary with info used to populate peptide table contents.

        Parameters
        ----------
        vaccine_peptide : VaccinePeptide
          The given peptide to convert to display form

        transcript_name : str
            RNA transcript name (should match that displayed in effect section)
        """
        mutant_protein_fragment = vaccine_peptide.mutant_protein_fragment
        amino_acids = mutant_protein_fragment.amino_acids
        peptide_data = OrderedDict([
            ('Transcript name', transcript_name),
            ('Length', len(amino_acids)),
            ('Expression score', round(vaccine_peptide.expression_score, 4)),
            ('Mutant epitope score', round(vaccine_peptide.mutant_epitope_score, 4)),
            ('Combined score', round(vaccine_peptide.combined_score, 4)),
            ('Reads fully spanning cDNA sequence(s)',
                mutant_protein_fragment.n_alt_reads_supporting_protein_sequence),
            ('Mutant amino acids', mutant_protein_fragment.n_mutant_amino_acids),
            ('Mutation distance from edge',
                mutant_protein_fragment.mutation_distance_from_edge),
        ])
        return peptide_data

    def _manufacturability_data(self, vaccine_peptide):
        """
        Returns an OrderedDict with manufacturability data for the given peptide.
        """
        scores = vaccine_peptide.manufacturability_scores
        manufacturability_data = OrderedDict([
            ('C-terminal 7mer GRAVY score', round(scores.cterm_7mer_gravy_score, 4)),
            ('Max 7mer GRAVY score', round(scores.max_7mer_gravy_score, 4)),
            ('N-terminal Glutamine, Glutamic Acid, or Cysteine',
                int(scores.difficult_n_terminal_residue)),
            ('C-terminal Cysteine', int(scores.c_terminal_cysteine)),
            ('C-terminal Proline', int(scores.c_terminal_proline)),
            ('Total number of Cysteine residues', scores.cysteine_count),
            ('N-terminal Asparagine', int(scores.n_terminal_asparagine)),
            ('Number of Asparagine-Proline bonds', scores.asparagine_proline_bond_count),
        ])
        return manufacturability_data

    def _epitope_data(self, epitope_prediction):
        """
        Returns an OrderedDict with epitope data from the given prediction.
        """
        epitope_data = OrderedDict([
            ('Sequence', epitope_prediction.peptide_sequence),
            ('IC50', epitope_prediction.ic50),
            ('Normalized binding score', round(epitope_prediction.logistic_score(), 4)),
            ('Allele', epitope_prediction.allele),
        ])
        return epitope_data

    def compute_template_data(self):
        patient_info = self._patient_info()

        # list ranked variants with their peptides
        variants = []
        for i, (variant, vaccine_peptides) in enumerate(
                self.ranked_variants_with_vaccine_peptides):
            variant_short_description = variant.short_description
            if len(vaccine_peptides) == 0:
                logger.info("Skipping gene(s) %s, variant %s: no vaccine peptides",
                    variant.gene_names, variant_short_description)
                continue

            variant_data = self._variant_data(vaccine_peptides[0])
            # TODO(julia): is this right?
            predicted_effect = top_priority_effect([
                variant.effect_on_transcript(t) for t in variant.transcripts])
            effect_data = self._effect_data(predicted_effect)

            peptides = []
            for j, vaccine_peptide in enumerate(vaccine_peptides):
                header_display_data = self._peptide_header_display_data(vaccine_peptide, j)
                peptide_data = self._peptide_data(vaccine_peptide, predicted_effect.transcript_name)
                manufacturability_data = self._manufacturability_data(vaccine_peptide)

                epitopes = []
                for epitope_prediction in vaccine_peptide.epitope_predictions:
                    if epitope_prediction.overlaps_mutation:
                        epitope_data = self._epitope_data(epitope_prediction)
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

        self.template_data.update({
            'patient_info': patient_info,
            'variants': variants,
        })
        return self.template_data


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

def new_columns():
    columns = OrderedDict([
        ("amino_acids", []),
        ("chr", []),
        ("pos", []),
        ("ref", []),
        ("alt", []),
        ("variant_rank", []),
        ("peptide_rank", []),
        ("mutation_start", []),
        ("mutation_end", []),
        ("combined_score", []),
        ("mutant_epitope_score", []),
    ])
    for field in ManufacturabilityScores._fields:
        columns[field] = []
    return columns

def _sanitize(val):
    """
    Converts values into display-friendly
    """
    if type(val) == bool:
        val = int(val)
    elif type(val) == float:
        val = round(val, 4)
    return val

def make_csv_report(
        ranked_variants_with_vaccine_peptides,
        report_dir_path,
        combined_report_path=None):
    if report_dir_path and not os.path.exists(report_dir_path):
        os.makedirs(report_dir_path)

    frames = []
    for i, (variant, vaccine_peptides) in enumerate(ranked_variants_with_vaccine_peptides):
        if not vaccine_peptides:
            continue
        filename = '%d_%s_chr%s_%d_%s_%s.csv' % (
            i + 1, vaccine_peptides[0].mutant_protein_fragment.gene_name,
            variant.contig, variant.start, variant.ref, variant.alt)
        path = os.path.join(report_dir_path, filename)
        columns = new_columns()
        for j, vaccine_peptide in enumerate(vaccine_peptides):
            columns["chr"].append(variant.contig)
            columns["pos"].append(variant.original_start)
            columns["ref"].append(variant.original_ref)
            columns["alt"].append(variant.original_alt)
            columns["variant_rank"].append(i + 1)
            columns["peptide_rank"].append(j + 1)
            columns["amino_acids"].append(vaccine_peptide.mutant_protein_fragment.amino_acids)
            columns["mutation_start"].append(
                vaccine_peptide.mutant_protein_fragment.mutant_amino_acid_start_offset)
            columns["mutation_end"].append(
                vaccine_peptide.mutant_protein_fragment.mutant_amino_acid_end_offset)
            columns["combined_score"].append(round(vaccine_peptide.combined_score, 4))
            columns["mutant_epitope_score"].append(round(vaccine_peptide.mutant_epitope_score, 4))
            for field in ManufacturabilityScores._fields:
                columns[field].append(
                    _sanitize(getattr(vaccine_peptide.manufacturability_scores, field)))
        df = pd.DataFrame(columns, columns=columns.keys())
        frames.append(df)
        if report_dir_path:
            df.to_csv(path, index=False)
            logger.info('Wrote CSV to %s', path)

    if combined_report_path:
        all_dfs = pd.concat(frames)
        # move rank columns to the front of the lines, for easy visual grouping
        colnames = all_dfs.columns.tolist()
        colnames.insert(0, colnames.pop(colnames.index('peptide_rank')))
        colnames.insert(0, colnames.pop(colnames.index('variant_rank')))
        all_dfs = all_dfs.reindex(columns=colnames)

        all_dfs.to_csv(combined_report_path, index=False)
        logger.info('Wrote combined CSV to %s', combined_report_path)

