# Copyright (c) 2016-2018. Mount Sinai School of Medicine
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
from importlib import import_module
import logging
import os
import sys
import tempfile

from astropy.io import ascii as asc
import jinja2
import pandas as pd
import pdfkit
import requests
import roman
from varcode import load_vcf_fast

from .manufacturability import ManufacturabilityScores

logger = logging.getLogger(__name__)


JINJA_ENVIRONMENT = jinja2.Environment(
    loader=jinja2.FileSystemLoader(os.path.dirname(__file__)),
    extensions=['jinja2.ext.autoescape'],
    autoescape=False,
    trim_blocks=True,
    lstrip_blocks=True,
)


PatientInfo = namedtuple("PatientInfo", (
    "patient_id",
    "vcf_paths",
    "bam_path",
    "mhc_alleles",
    "num_somatic_variants",
    "num_coding_effect_variants",
    "num_variants_with_rna_support",
    "num_variants_with_vaccine_peptides",
))


class TemplateDataCreator(object):
    def __init__(
            self,
            ranked_variants_with_vaccine_peptides,
            patient_info,
            final_review,
            reviewers,
            args_for_report,
            input_json_file,
            cosmic_vcf_filename=None):
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
            'reviewers': reviewers.split(',') if reviewers else [],
            'final_review': final_review,
            'input_json_file': input_json_file,
            # these report sections are optional
            'include_manufacturability': args_for_report['manufacturability'],
            'include_wt_epitopes': args_for_report['wt_epitopes'],
        }

        # map from peptide objects to their COSMIC IDs if they exist
        if cosmic_vcf_filename:
            logger.info('Loading COSMIC data...')
            self.cosmic_variant_collection = load_vcf_fast(
                cosmic_vcf_filename, allow_extended_nucleotides=True, include_info=False)
            logger.info('COSMIC data loaded.')
        else:
            self.cosmic_variant_collection = None

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
            ('Somatic variants with predicted coding effects and RNA support',
                self.patient_info.num_variants_with_rna_support),
            ('Somatic variants with predicted coding effects, RNA support and predicted MHC '
                'ligands',
                self.patient_info.num_variants_with_vaccine_peptides),
        ])
        return patient_info

    def _variant_data(self, variant, top_vaccine_peptide):
        """
        Returns an OrderedDict with info used to populate variant info section.
        """
        variant_data = OrderedDict()
        mutant_protein_fragment = top_vaccine_peptide.mutant_protein_fragment
        top_score = _sanitize(top_vaccine_peptide.combined_score)
        igv_locus = "chr%s:%d" % (variant.contig, variant.start)
        variant_data = OrderedDict([
            ('IGV locus', igv_locus),
            ('Gene name', mutant_protein_fragment.gene_name),
            ('Top score', top_score),
            ('RNA reads supporting variant allele', mutant_protein_fragment.n_alt_reads),
            ('RNA reads supporting reference allele', mutant_protein_fragment.n_ref_reads),
            ('RNA reads supporting other alleles', mutant_protein_fragment.n_other_reads),
        ])
        return variant_data

    def _effect_data(self, predicted_effect):
        """
        Returns an OrderedDict with info about the given predicted effect.
        """
        effect_data = OrderedDict([
            ('Effect type', predicted_effect.__class__.__name__),
            ('Transcript name', predicted_effect.transcript_name),
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
            ('Expression score', _sanitize(vaccine_peptide.expression_score)),
            ('Mutant epitope score', _sanitize(vaccine_peptide.mutant_epitope_score)),
            ('Combined score', _sanitize(vaccine_peptide.combined_score)),
            ('Max coding sequence coverage',
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
            ('C-terminal 7mer GRAVY score', _sanitize(scores.cterm_7mer_gravy_score)),
            ('Max 7mer GRAVY score', _sanitize(scores.max_7mer_gravy_score)),
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
        # if the WT peptide is too short, it's possible that we're missing a prediction for it
        if epitope_prediction.wt_ic50 is not None:
            wt_ic50_str = '%.2f nM' % epitope_prediction.wt_ic50
        else:
            wt_ic50_str = 'No prediction'
        epitope_data = OrderedDict([
            ('Sequence', epitope_prediction.peptide_sequence),
            ('IC50', '%.2f nM' % epitope_prediction.ic50),
            ('Score', _sanitize(epitope_prediction.logistic_epitope_score())),
            ('Allele', epitope_prediction.allele.replace('HLA-', '')),
            ('WT sequence', epitope_prediction.wt_peptide_sequence),
            ('WT IC50', wt_ic50_str),
        ])
        return epitope_data

    def _query_cosmic(self, variant):
        if not self.cosmic_variant_collection:
            return None
        if variant in self.cosmic_variant_collection.metadata:
            # IDs in the DB are of the form 'COSM725245'
            cosmic_id = self.cosmic_variant_collection.metadata[variant]['id']
            link_for_report = \
                "http://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=%s" % cosmic_id[4:]
            logger.info("Link for report: %s", link_for_report)
            return link_for_report

        logger.info("Variant not in COSMIC")
        return None

    def _query_wustl(self, predicted_effect, gene_name):
        """
        Returns a link to the WUSTL page for this variant, if present.
        """
        amino_acids = predicted_effect.short_description
        api_url = "http://www.docm.info/api/v1/variants.json?amino_acids=%s&genes=%s" % (
            amino_acids, gene_name.upper())
        logger.info("WUSTL link: %s", api_url)

        try:
            contents = requests.get(api_url).json()
            if len(contents) > 0:
                hgvs = contents[0]['hgvs']
                link_for_report = "http://docm.genome.wustl.edu/variants/%s" % hgvs
                logger.info("Link for report: %s", link_for_report)
                return link_for_report
        except requests.exceptions.ConnectionError as e:
            logger.warn('ConnectionError reaching WUSTL: %s', e)
            return None

        logger.info('Variant not found in WUSTL')
        return None

    def _databases(self, variant, predicted_effect, gene_name):
        databases = {}
        wustl_link = self._query_wustl(predicted_effect, gene_name)
        if wustl_link:
            databases['WUSTL'] = wustl_link

        cosmic_link = self._query_cosmic(variant)
        if cosmic_link:
            databases['COSMIC'] = cosmic_link

        return databases

    def compute_template_data(self):
        patient_info = self._patient_info()

        # list ranked variants with their peptides
        variants = []
        num = 0
        for (variant, vaccine_peptides) in self.ranked_variants_with_vaccine_peptides:
            variant_short_description = variant.short_description
            if len(vaccine_peptides) == 0:
                logger.info(
                    "Skipping gene(s) %s, variant %s: no vaccine peptides",
                    variant.gene_names, variant_short_description)
                continue
            num += 1

            top_peptide = vaccine_peptides[0]
            variant_data = self._variant_data(variant, top_peptide)
            predicted_effect = top_peptide.mutant_protein_fragment.predicted_effect()
            effect_data = self._effect_data(predicted_effect)

            databases = self._databases(
                variant, predicted_effect, top_peptide.mutant_protein_fragment.gene_name)

            peptides = []
            for j, vaccine_peptide in enumerate(vaccine_peptides):
                if not vaccine_peptide.contains_mutant_epitopes():
                    logger.info('No epitopes for peptide: %s', vaccine_peptide)
                    continue

                header_display_data = self._peptide_header_display_data(vaccine_peptide, j)
                peptide_data = self._peptide_data(vaccine_peptide, predicted_effect.transcript_name)
                manufacturability_data = self._manufacturability_data(vaccine_peptide)

                epitopes = []
                wt_epitopes = []
                for mutant_epitope_prediction in vaccine_peptide.mutant_epitope_predictions:
                    epitopes.append(self._epitope_data(mutant_epitope_prediction))

                for wt_epitope_prediction in vaccine_peptide.wildtype_epitope_predictions:
                    epitope_data = self._epitope_data(wt_epitope_prediction)
                    key_list = ['Allele', 'IC50', 'Sequence']
                    wt_epitopes.append({key: epitope_data[key] for key in key_list})

                # hack: make a nicely-formatted fixed width table for epitopes, used in ASCII report
                with tempfile.TemporaryFile(mode='r+') as temp:
                    asc.write(epitopes, temp, format='fixed_width_two_line', delimiter_pad=' ')
                    temp.seek(0)
                    ascii_epitopes = temp.read()

                ascii_wt_epitopes = None
                if len(wt_epitopes) > 0:
                    with tempfile.TemporaryFile(mode='r+') as temp:
                        asc.write(
                            wt_epitopes, temp, format='fixed_width_two_line', delimiter_pad=' ')
                        temp.seek(0)
                        ascii_wt_epitopes = temp.read()
                peptide_dict = {
                    'header_display_data': header_display_data,
                    'peptide_data': peptide_data,
                    'manufacturability_data': manufacturability_data,
                    'epitopes': epitopes,
                    'ascii_epitopes': ascii_epitopes,
                    'wt_epitopes': wt_epitopes,
                    'ascii_wt_epitopes': ascii_wt_epitopes,
                }
                peptides.append(peptide_dict)

            # if there are no peptides for this variant, exclude from report
            if len(peptides) == 0:
                logger.info('No peptides for variant: %s', variant)
                continue

            variant_dict = {
                'num': num,
                'short_description': variant_short_description,
                'variant_data': variant_data,
                'effect_data': effect_data,
                'peptides': peptides,
                'databases': databases,
            }
            variants.append(variant_dict)

        # add package metadata to the report
        package_versions = {}
        for name in ['vaxrank', 'isovar', 'mhctools', 'varcode', 'pyensembl']:
            module = import_module(name)
            version = getattr(module, '__version__')
            package_versions[name] = version

        self.template_data.update({
            'patient_info': patient_info,
            'variants': variants,
            'package_versions': package_versions,
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
            'zoom': 0.55,
            'margin-top': '20mm'
        }

        if sys.platform in ('linux', 'linux2'):
            # pdfkit uses wkhtmltopdf, which doesn't work on headless servers;
            # recommended workaround is to use xvfb, as documented here:
            # https://github.com/wkhtmltopdf/wkhtmltopdf/issues/2037#issuecomment-62019521
            from xvfbwrapper import Xvfb
            logger.info('Running pdfkit inside xvfb wrapper')
            with Xvfb():
                pdfkit.from_file(f.name, pdf_report_path, options=options)

        else:
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
        ("expression_score", []),
        ("mutant_epitope_score", []),
    ])
    for field in ManufacturabilityScores._fields:
        columns[field] = []
    return columns

def _str_sig_figs(input, n_sig_figs):
    return '{:g}'.format(float('{:.{p}g}'.format(input, p=n_sig_figs)))

def _sanitize(val):
    """
    Converts values into display-friendly
    """
    if type(val) == bool:
        return int(val)
    else:
        return _str_sig_figs(val, 5)

def resize_columns(worksheet, amino_acids_col, pos_col):
    """
    Resizes amino acid and mutant position columns in the Excel sheet so that they don't
    have to be expanded.
    """
    worksheet.set_column('%s:%s' % (amino_acids_col, amino_acids_col), 40)
    worksheet.set_column('%s:%s' % (pos_col, pos_col), 12)

def make_minimal_neoepitope_report(
        ranked_variants_with_vaccine_peptides,
        num_epitopes_per_peptide=None,
        excel_report_path=None):
    """
    Creates a simple Excel spreadsheet containing one neoepitope per row

    Parameters
    ----------
    ranked_variants_with_vaccine_peptides :
      Ranked list of (variant, list of its vaccine peptides)

    num_epitopes_per_peptide : int
      The number of epitopes to include for each vaccine peptide; these are sorted before cutoff.
      If None, all epitopes will be included in the output

    excel_report_path : str
      Path to which to write the output Excel file
    """
    rows = []
    # each row in the spreadsheet is a neoepitope
    for (variant, vaccine_peptides) in ranked_variants_with_vaccine_peptides:
        for vaccine_peptide in vaccine_peptides:
            # only include mutant epitopes
            for epitope_prediction in vaccine_peptide.mutant_epitope_predictions:
                row = OrderedDict([
                    ('Allele', epitope_prediction.allele),
                    ('Mutant peptide sequence', epitope_prediction.peptide_sequence),
                    ('Score', vaccine_peptide.mutant_epitope_score),
                    ('Predicted mutant pMHC affinity', '%.2f nM' % epitope_prediction.ic50),
                    ('Variant allele RNA read count',
                        vaccine_peptide.mutant_protein_fragment.n_alt_reads),
                    ('Wildtype sequence', epitope_prediction.wt_peptide_sequence),
                    ('Predicted wildtype pMHC affinity',
                        '%.2f nM' % epitope_prediction.wt_ic50),
                    ('Gene name', vaccine_peptide.mutant_protein_fragment.gene_name),
                    ('Genomic variant', variant.short_description),
                ])
                rows.append(row)

    if len(rows) > 0:
        df = pd.DataFrame.from_dict(rows)
        writer = pd.ExcelWriter(excel_report_path, engine='xlsxwriter')
        df.to_excel(writer, sheet_name='Neoepitopes', index=False)

        # resize columns to be not crappy
        worksheet = writer.sheets['Neoepitopes']
        worksheet.set_column('%s:%s' % ('B', 'B'), 23)
        worksheet.set_column('%s:%s' % ('D', 'D'), 27)
        worksheet.set_column('%s:%s' % ('E', 'E'), 26)
        worksheet.set_column('%s:%s' % ('F', 'F'), 17)
        worksheet.set_column('%s:%s' % ('G', 'G'), 30)
        worksheet.set_column('%s:%s' % ('H', 'H'), 9)
        worksheet.set_column('%s:%s' % ('I', 'I'), 18)
        writer.save()
        logger.info('Wrote XLSX neoepitope report file to %s', excel_report_path)


def make_csv_report(
        ranked_variants_with_vaccine_peptides,
        excel_report_path=None,
        csv_report_path=None):
    """
    Writes out CSV/XLSX reports as needed.
    """
    # make a bunch of pd frames, save them in an OrderedDict with keys being descriptive Excel
    # sheet names (will be used later for making the Excel report if needed)
    frames = OrderedDict()
    for i, (variant, vaccine_peptides) in enumerate(ranked_variants_with_vaccine_peptides):
        any_vaccine_peptides = False
        if not vaccine_peptides:
            continue

        sheet_name = '%d_%s_chr%s_%d_%s_%s' % (
            i + 1, vaccine_peptides[0].mutant_protein_fragment.gene_name,
            variant.contig, variant.start, variant.ref, variant.alt)
        columns = new_columns()
        for j, vaccine_peptide in enumerate(vaccine_peptides):

            # if there are no predicted epitopes, exclude this peptide from the report
            if not vaccine_peptide.contains_mutant_epitopes():
                logger.info('No epitopes for peptide: %s', vaccine_peptide)
                continue

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
            columns["combined_score"].append(_sanitize(vaccine_peptide.combined_score))
            columns["expression_score"].append(_sanitize(vaccine_peptide.expression_score))
            columns["mutant_epitope_score"].append(_sanitize(vaccine_peptide.mutant_epitope_score))
            for field in ManufacturabilityScores._fields:
                columns[field].append(
                    _sanitize(getattr(vaccine_peptide.manufacturability_scores, field)))
            any_vaccine_peptides = True

        if not any_vaccine_peptides:
            continue

        df = pd.DataFrame(columns, columns=columns.keys())
        frames[sheet_name] = df

    if not frames:
        logger.info('No data for CSV or XLSX report')
        return

    all_dfs = pd.concat(frames.values())
    if csv_report_path:
        all_dfs.to_csv(csv_report_path, index=False)
        logger.info('Wrote CSV report file to %s', csv_report_path)

    if excel_report_path:
        writer = pd.ExcelWriter(excel_report_path, engine='xlsxwriter')

        # copy the variant rank column to position 0, make first sheet called "All"
        all_dfs[''] = all_dfs['variant_rank']
        colnames = all_dfs.columns.tolist()
        colnames.insert(0, colnames.pop(colnames.index('')))
        all_dfs = all_dfs.reindex(columns=colnames)
        all_dfs.to_excel(writer, sheet_name='All', index=False)
        resize_columns(writer.sheets['All'], 'B', 'D')

        # add one sheet per variant
        for sheet_name, df in frames.items():
            # trim sheet names to 31 characters due to limit in Excel
            # should still be unique since they start with the variant
            # index
            shortened_sheet_name = sheet_name[:31]
            df.to_excel(writer, sheet_name=shortened_sheet_name, index=False)
            resize_columns(writer.sheets[shortened_sheet_name], 'A', 'C')

        writer.save()
        logger.info('Wrote manufacturer XLSX file to %s', excel_report_path)
