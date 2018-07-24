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
from collections import defaultdict, namedtuple, OrderedDict
import logging

from numpy import isclose
import pandas as pd

from isovar.protein_sequences import (
    reads_generator_to_protein_sequences_generator,
)

from .mutant_protein_fragment import MutantProteinFragment
from .epitope_prediction import predict_epitopes, slice_epitope_predictions
from .vaccine_peptide import VaccinePeptide


logger = logging.getLogger(__name__)

class VaxrankCoreLogic:
    def __init__(
            self,
            variants,
            reads_generator,
            mhc_predictor,
            vaccine_peptide_length,
            padding_around_mutation,
            max_vaccine_peptides_per_variant,
            min_alt_rna_reads,
            min_variant_sequence_coverage,
            variant_sequence_assembly,
            num_mutant_epitopes_to_keep=10000,
            min_epitope_score=0.0,
            interferon_gamma_response_csv=None,
            class1_mhc_presentation_pathway_csv=None,
            cancer_driver_genes_csv=None,
            cancer_driver_variants_csv=None):
        """
        Construct a VaxrankCoreLogic object.

        Parameters
        ----------
        variants : VariantCollection
            Variant objects to evaluate for vaccine inclusion

        reads_generator : generator
            Yields sequence of varcode.Variant objects paired with sequences of AlleleRead objects
            that support that variant.

        mhc_predictor : mhctools.BasePredictor
            Object with predict_peptides method, used for making pMHC binding predictions

        vaccine_peptide_length : int
            Length of vaccine SLP to construct

        padding_around_mutation : int
            Number of off-center windows around the mutation to consider as vaccine peptides.

        max_vaccine_peptides_per_variant : int
            Number of vaccine peptides to generate for each mutation.

        min_alt_rna_reads : int
            Drop variant sequences at loci with fewer than this number of reads supporting the alt
            allele.

        min_variant_sequence_coverage : int
            Trim variant sequences to positions supported by at least this number of RNA reads.

        variant_sequence_assembly : int
            If True, then assemble variant cDNA sequences based on overlap of RNA reads. If False,
            then variant cDNA sequences must be fully spanned and contained within RNA reads.

        num_mutant_epitopes_to_keep : int, optional
            Number of top-ranking epitopes for each vaccine peptide to include in computation.

        min_epitope_score : float, optional
            Ignore peptides with binding predictions whose normalized score is less than this.

        interferon_gamma_response_csv : str, optional
            Local path to interferon-gamma response CSV file.

        class1_mhc_presentation_pathway_csv : str, optional
            Local path to MHC class I presentation pathway CSV file.

        cancer_driver_genes_csv : str, optional
            Local path to cancer driver genes CSV file.

        cancer_driver_variants_csv : str, optional
            Local path to cancer driver variants CSV file.

        """
        self.variants = variants
        self.reads_generator = reads_generator
        self.mhc_predictor = mhc_predictor
        self.vaccine_peptide_length = vaccine_peptide_length
        self.padding_around_mutation = padding_around_mutation
        self.max_vaccine_peptides_per_variant = max_vaccine_peptides_per_variant
        self.min_alt_rna_reads = min_alt_rna_reads
        self.min_variant_sequence_coverage = min_variant_sequence_coverage
        self.variant_sequence_assembly = variant_sequence_assembly
        self.num_mutant_epitopes_to_keep = num_mutant_epitopes_to_keep
        self.min_epitope_score = min_epitope_score

        # load in gene dataframes we want to search
        self.interferon_gamma_response = None
        self.class1_mhc_presentation_pathway = None
        self.cancer_driver_genes = None
        self.cancer_driver_variants = None

        if interferon_gamma_response_csv:
            self.interferon_gamma_response = pd.read_csv(interferon_gamma_response_csv)
        if class1_mhc_presentation_pathway_csv:
            self.class1_mhc_presentation_pathway = pd.read_csv(class1_mhc_presentation_pathway_csv)
        if cancer_driver_genes_csv:
            self.cancer_driver_genes = pd.read_csv(cancer_driver_genes_csv)
        if cancer_driver_variants_csv:
            self.cancer_driver_variants = pd.read_csv(cancer_driver_variants_csv)

    def vaccine_peptides_for_variant(
            self,
            variant,
            isovar_protein_sequence):
        """
        Returns sorted list of VaccinePeptide objects.
        """
        protein_fragment = MutantProteinFragment.from_isovar_protein_sequence(
            variant=variant,
            protein_sequence=isovar_protein_sequence)

        logger.info(
            "Mutant protein fragment for %s: %s",
            variant,
            protein_fragment)

        epitope_predictions = predict_epitopes(
            mhc_predictor=self.mhc_predictor,
            protein_fragment=protein_fragment,
            min_epitope_score=self.min_epitope_score,
            genome=variant.ensembl).values()

        candidate_vaccine_peptides = []

        for offset, candidate_fragment in protein_fragment.top_k_subsequences(
                subsequence_length=self.vaccine_peptide_length,
                k=2 * self.padding_around_mutation + 1):
            subsequence_epitope_predictions = slice_epitope_predictions(
                epitope_predictions,
                start_offset=offset,
                end_offset=offset + len(candidate_fragment))
            # filter out peptides that have no epitopes
            if not subsequence_epitope_predictions:
                logger.info(
                    "No epitope predictions for mutant protein fragment %s",
                    candidate_fragment)
                continue

            assert all(
                p.source_sequence == candidate_fragment.amino_acids
                for p in subsequence_epitope_predictions)

            candidate_vaccine_peptide = VaccinePeptide(
                mutant_protein_fragment=candidate_fragment,
                epitope_predictions=subsequence_epitope_predictions,
                num_mutant_epitopes_to_keep=self.num_mutant_epitopes_to_keep)

            logger.debug(
                "%s, combined score: %0.4f",
                candidate_vaccine_peptide,
                candidate_vaccine_peptide.combined_score)
            candidate_vaccine_peptides.append(candidate_vaccine_peptide)

        n_total_candidates = len(candidate_vaccine_peptides)
        if n_total_candidates == 0:
            logger.info("No candidate peptides for variant %s", variant.short_description)
            return []

        max_score = max(vp.combined_score for vp in candidate_vaccine_peptides)
        if isclose(max_score, 0.0):
            filtered_candidate_vaccine_peptides = candidate_vaccine_peptides
        else:
            # only keep candidate vaccines that are within 1% of the maximum
            # combined score
            filtered_candidate_vaccine_peptides = [
                vp
                for vp in candidate_vaccine_peptides
                if vp.combined_score / max_score > 0.99
            ]

        logger.info(
            "Keeping %d/%d vaccine peptides for %s",
            len(filtered_candidate_vaccine_peptides),
            n_total_candidates,
            variant)
        filtered_candidate_vaccine_peptides.sort(key=VaccinePeptide.lexicographic_sort_key)

        if len(filtered_candidate_vaccine_peptides) > 0:
            logger.debug("Top vaccine peptides for %s:", variant)
            for i, vaccine_peptide in enumerate(filtered_candidate_vaccine_peptides):
                logger.debug(
                    "%d) %s (combined score = %0.4f)",
                    i + 1,
                    vaccine_peptide,
                    vaccine_peptide.combined_score)
        return filtered_candidate_vaccine_peptides[:self.max_vaccine_peptides_per_variant]

    def generate_vaccine_peptides(self):
        """
        Returns a tuple of two values:
        - dictionary mapping each variant to list of VaccinePeptide objects
        - dictionary containing some variant counts for report display
        """

        # total number of amino acids is the vaccine peptide length plus the
        # number of off-center windows around the mutation
        protein_fragment_sequence_length = (
            self.vaccine_peptide_length + 2 * self.padding_around_mutation)

        protein_sequences_generator = reads_generator_to_protein_sequences_generator(
            self.reads_generator,
            transcript_id_whitelist=None,
            protein_sequence_length=protein_fragment_sequence_length,
            min_alt_rna_reads=self.min_alt_rna_reads,
            min_variant_sequence_coverage=self.min_variant_sequence_coverage,
            variant_sequence_assembly=self.variant_sequence_assembly,
            max_protein_sequences_per_variant=1)

        variant_metadata_dicts = []

        result_dict = {}
        counts_dict = defaultdict(int)
        for variant, isovar_protein_sequences in protein_sequences_generator:
            # import pdb; pdb.set_trace()
            # figure out which source(s) this variant came from
            # is this variant coding?
            # does it have RNA support?
            # do we predict it to result in a strong MHC binder/does it have mutant epitopes?
            # is the variant's gene in a driver list, or is the variant a known driver some way?

            # TODO: check for variants that are in a coding region, silent, but near splice sites
            # (ExonicSpliceSite)

            # some default values
            variant_metadata_dict = OrderedDict([
                ('chr', variant.contig),
                ('start', variant.start),
                ('ref', variant.ref),
                ('alt', variant.alt),
                ('is_coding_nonsynonymous', False),
                ('rna_support', False),
                ('mhc_binder', False),
                ('interferon_gamma_response', False),
                ('class1_mhc_presentation_pathway', False),
                ('cancer_driver_gene', False),
                ('cancer_driver_variant', False),
            ])
            variant_metadata_dicts.append(variant_metadata_dict)

            # check this variant against some databases
            effect = variant.effects().top_priority_effect().short_description
            gene_ids = variant.gene_ids

            # check each gene ID against databases
            def is_present(df, value, column):
                return len(df.loc[df[column] == value]) > 0

            for gene_id in gene_ids:
                if self.interferon_gamma_response is not None and \
                        is_present(self.interferon_gamma_response, gene_id, 'Ensembl Gene ID'):
                    variant_metadata_dict['interferon_gamma_response'] = True
                if self.class1_mhc_presentation_pathway is not None and \
                        is_present(self.class1_mhc_presentation_pathway, gene_id, 'Ensembl Gene ID'):
                    variant_metadata_dict['class1_mhc_presentation_pathway'] = True
                if self.cancer_driver_genes is not None and \
                        is_present(self.cancer_driver_genes, gene_id, 'Ensembl Gene ID'):
                    variant_metadata_dict['cancer_driver_gene'] = True
                if self.cancer_driver_variants is not None and \
                        len(self.cancer_driver_variants.loc[
                            (self.cancer_driver_variants['Ensembl Gene ID'] == gene_id) & 
                            (self.cancer_driver_variants['Mutation'] == effect)]) > 0:
                    variant_metadata_dict['cancer_driver_variant'] = True

            counts_dict['num_total_variants'] += 1
            if len(variant.effects().drop_silent_and_noncoding()) > 0:
                counts_dict['num_coding_effect_variants'] += 1
                variant_metadata_dict['is_coding_nonsynonymous'] = True
            isovar_protein_sequences = list(isovar_protein_sequences)
            if len(isovar_protein_sequences) == 0:
                # this means the variant RNA support is below threshold
                logger.info("No protein sequences for %s", variant)
                continue

            variant_metadata_dict['rna_support'] = True
            counts_dict['num_variants_with_rna_support'] += 1

            # use the first protein sequence - why?
            isovar_protein_sequence = isovar_protein_sequences[0]
            vaccine_peptides = self.vaccine_peptides_for_variant(
                variant=variant,
                isovar_protein_sequence=isovar_protein_sequence)

            # do any of this variant's vaccine peptides contain mutant epitopes?
            any_mutant_epitopes = False
            for vaccine_peptide in vaccine_peptides:
                if vaccine_peptide.contains_mutant_epitopes():
                    any_mutant_epitopes = True
                    break
            if any_mutant_epitopes:
                variant_metadata_dict['mhc_binder'] = True
                counts_dict['num_variants_with_vaccine_peptides'] += 1
            result_dict[variant] = vaccine_peptides

        for key, value in counts_dict.items():
            logger.info('%s: %d', key, value)

        df = pd.DataFrame(variant_metadata_dicts)
        df.to_csv('all_variants.csv', index=False)

        return result_dict, counts_dict

    def ranked_vaccine_peptides(
            self,
            num_mutant_epitopes_to_keep=10000,
            min_epitope_score=0):
        """
        This function is an entrypoint into the core vaccine peptide ranking logic.
        Returns a tuple of two values:
        - sorted list whose first element is a Variant and whose second
        element is a list of VaccinePeptide objects
        - dictionary containing "funnel" variant counts for report display. Keys are e.g.
        "num_variants_with_rna_support", values are integer counts.
        """
        variants_to_vaccine_peptides_dict, variant_counts_dict = self.generate_vaccine_peptides()
        result_list = list(variants_to_vaccine_peptides_dict.items())
        # TODO: move this sort key into its own function, also make less nuts
        result_list.sort(
            key=lambda x: x[1][0].combined_score if len(x[1]) > 0 else 0.0,
            reverse=True)
        return result_list, variant_counts_dict
