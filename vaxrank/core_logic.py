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
from .gene_pathway_check import GenePathwayCheck


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
            gene_pathway_check=None):
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

        gene_pathway_check : GenePathwayCheck object, optional
            If provided, will check against known pathways/gene sets/variant sets and include the
            info in the all-variants output file.

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
        self.gene_pathway_check = gene_pathway_check

        # dictionary which will contain some overall variant counts for report display
        self.counts_dict = {
            'num_total_variants': len(variants),
            'num_coding_effect_variants': 0,
            'num_variants_with_rna_support': 0,
            'num_variants_with_vaccine_peptides': 0,
        }

        # dictionary: varcode.Variant -> dictionary of properties we want to analyze later, e.g.
        # whether this variant is part of a pathway of interest, is a strong MHC binder, etc.
        self.variant_properties_dict = {}
        for variant in self.variants:
            self.variant_properties_dict[variant] = OrderedDict((
                ('contig', variant.contig),
                ('start', variant.start),
                ('ref', variant.ref),
                ('alt', variant.alt),
                ('is_coding_nonsynonymous', False),
                ('rna_support', False),
                ('mhc_binder', False),
            ))

        # dictionary: varcode.Variant -> list of VaccinePeptide objects
        self.variant_peptides_dict = {}

        self.variants_processed = False

    def vaccine_peptides_for_variant(
            self,
            variant,
            isovar_protein_sequence):
        """
        Returns sorted list of VaccinePeptide objects. If there are no suitable vaccine peptides
        (no strong MHC binder subsequences), returns an empty list.

        At this point, we know the variant has RNA support, as per the isovar_protein_sequence.
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

    def process_variants(self):
        """
        Constructs a dictionary mapping each vaccine variant to list of VaccinePeptide objects, for
        variants where we have candidate vaccine peptides. Also populates some other dictionaries.
        """

        # total number of amino acids is the vaccine peptide length plus the
        # number of off-center windows around the mutation
        protein_fragment_sequence_length = (
            self.vaccine_peptide_length + 2 * self.padding_around_mutation)

        """
        These sequences are only the ones that overlap the variant and support the mutation.
        Right now, this generator yields:
        - (variant, mutant protein sequences) if there's enough alt RNA support
        - (variant, None) if the variant is silent or there are ref reads overlapping the variant
        locus but inadequate alt RNA support.
        - does not return the variant if there's no RNA support for ref or alt - we may miss some
        coding variants this way unless we check for them explicitly

        Future intended behavior: returns all passing variants, with a protein sequences generator
        that is non empty only if there are enough alt RNA reads supporting the variant
        """
        protein_sequences_generator = reads_generator_to_protein_sequences_generator(
            self.reads_generator,
            transcript_id_whitelist=None,
            protein_sequence_length=protein_fragment_sequence_length,
            min_alt_rna_reads=self.min_alt_rna_reads,
            min_variant_sequence_coverage=self.min_variant_sequence_coverage,
            variant_sequence_assembly=self.variant_sequence_assembly,
            max_protein_sequences_per_variant=1)

        for variant, isovar_protein_sequences in protein_sequences_generator:
            # add pathway info if available
            if self.gene_pathway_check is not None:
                pathway_dict = self.gene_pathway_check.make_variant_dict(variant)
                self.variant_properties_dict[variant].update(pathway_dict)

            if len(variant.effects().drop_silent_and_noncoding()) > 0:
                self.counts_dict['num_coding_effect_variants'] += 1
                self.variant_properties_dict[variant]['is_coding_nonsynonymous'] = True
            isovar_protein_sequences = list(isovar_protein_sequences)
            if len(isovar_protein_sequences) == 0:
                # this means the variant RNA support is below threshold
                logger.info("No protein sequences for %s", variant)
                continue

            self.variant_properties_dict[variant]['rna_support'] = True
            self.counts_dict['num_variants_with_rna_support'] += 1

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
                # TODO: need to know if there are strong binders in the predicted mutant fragment
                # that's just based on the predicted effect, rather than being assembled from RNA.
                # For now, the 'mhc_binder' property will depend on whether there's RNA support
                self.variant_properties_dict[variant]['mhc_binder'] = True
                self.counts_dict['num_variants_with_vaccine_peptides'] += 1
                self.variant_peptides_dict[variant] = vaccine_peptides

        self.variants_processed = True

    def ranked_vaccine_peptides(
            self,
            num_mutant_epitopes_to_keep=10000,
            min_epitope_score=0):
        """
        This function returns a sorted list whose first element is a Variant and whose second
        element is a list of VaccinePeptide objects.
        """
        if not self.variants_processed:
            raise ValueError("Must call process_variants() first")
            
        result_list = list(self.variant_peptides_dict.items())
        # TODO: move this sort key into its own function, also make less nuts
        result_list.sort(
            key=lambda x: x[1][0].combined_score if len(x[1]) > 0 else 0.0,
            reverse=True)
        return result_list

    def variant_properties(self):
        """
        This function returns a list of 
        """
        if not self.variants_processed:
            raise ValueError("Must call process_variants() first")
        # return a list of dictionaries representing all variants
        return self.variant_properties_dict.values()

    def variant_counts(self):
        """
        This function returns a dictionary containing "funnel" variant counts for report display.
        Keys are e.g. "num_variants_with_rna_support", values are integer counts.
        """
        if not self.variants_processed:
            raise ValueError("Must call process_variants() first")
        for key, value in self.counts_dict.items():
            logger.info('%s: %d', key, value)
        return self.counts_dict
