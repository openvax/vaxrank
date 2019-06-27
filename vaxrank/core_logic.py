# Copyright (c) 2016-2019. Mount Sinai School of Medicine
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
from collections import OrderedDict
import logging

from numpy import isclose

from isovar import run_isovar


from .mutant_protein_fragment import MutantProteinFragment
from .epitope_prediction import predict_epitopes, slice_epitope_predictions
from .vaccine_peptide import VaccinePeptide

logger = logging.getLogger(__name__)


class VaxrankCoreLogic(object):
    def __init__(
            self,
            read_collector,
            protein_sequence_creator,
            filter_thresholds,
            mhc_predictor,
            vaccine_peptide_length,
            padding_around_mutation,
            max_vaccine_peptides_per_variant,
            num_mutant_epitopes_to_keep=10000,
            min_epitope_score=0.0,
            gene_pathway_check=None):
        """
        Parameters
        ----------
        read_collector : isovar.ReadCollector
            Isovar object used to extract reads from the tumor RNA, holds
            options related to read processing.

        protein_sequence_creator : isovar.ProteinSequenceCreator
            Isovar object used to assemble mutant protein sequences from RNA
            reads, holds options related to assembly.

        filter_thresholds : dict
            Dictionary mapping names of IsovarResult properties (with either
            "min_" or "max" beforehand) to threshold cutoffs,
            e.g. {"max_num_ref_reads": 10}

        mhc_predictor : mhctools.BasePredictor
            Object with predict_peptides method, used for making pMHC binding
            predictions

        vaccine_peptide_length : int
            Length of vaccine SLP to construct

        padding_around_mutation : int
            Number of off-center windows around the mutation to consider as vaccine
            peptides.

        num_mutant_epitopes_to_keep : int, optional
            Number of top-ranking epitopes for each vaccine peptide to include in
            computation.

        min_epitope_score : float, optional
            Ignore peptides with binding predictions whose normalized score is less
            than this.

        gene_pathway_check : GenePathwayCheck, optional
            If provided, will check against known pathways/gene sets/variant sets and
            include the info in the all-variants output file.
        """
        self.read_collector = read_collector
        self.protein_sequence_creator = protein_sequence_creator
        self.filter_thresholds = filter_thresholds
        self.mhc_predictor = mhc_predictor
        self.vaccine_peptide_length = vaccine_peptide_length
        self.padding_around_mutation = padding_around_mutation
        self.max_vaccine_peptides_per_variant = max_vaccine_peptides_per_variant
        self.num_mutant_epitopes_to_keep = num_mutant_epitopes_to_keep
        self.min_epitope_score = min_epitope_score
        self.gene_pathway_check = gene_pathway_check

        # will be a dictionary: varcode.Variant -> IsovarResult
        self._variant_to_isovar_result_dict = None

        # will be a dictionary: varcode.Variant -> list(VaccinePeptide)
        self._variant_to_vaccine_peptides_dict = None

    @property
    def variant_to_isovar_result_dict(self):
        if self._variant_to_isovar_result_dict is None:
            raise ValueError("You must call VaxrankCoreLogic.run_isovar")
        return self._variant_to_isovar_result_dict

    @property
    def variant_to_vaccine_peptides_dict(self):
        if self._variant_to_vaccine_peptides_dict is None:
            raise ValueError("You must call VaxrankCoreLogic.select_vaccine_peptides")
        return self._variant_to_vaccine_peptides_dict

    def vaccine_peptides_for_variant(self, isovar_result):
        """
        Parameters
        ----------
        isovar_result : isovar.IsovarResult

        Returns sorted list of VaccinePeptide objects. If there are no suitable
        vaccine peptides returns an empty list.
        """
        variant = isovar_result.variant
        protein_fragment = MutantProteinFragment.from_isovar_result(isovar_result)

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

            assert all([
                p.source_sequence == candidate_fragment.amino_acids
                for p in subsequence_epitope_predictions])

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
        n_filtered = len(filtered_candidate_vaccine_peptides)
        logger.info(
            "Keeping %d/%d vaccine peptides for %s",
            n_filtered,
            n_total_candidates,
            variant)

        if n_filtered == 0:
            return []

        filtered_candidate_vaccine_peptides.sort(key=VaccinePeptide.lexicographic_sort_key)

        logger.debug("Top vaccine peptides for %s:", variant)
        for i, vaccine_peptide in enumerate(filtered_candidate_vaccine_peptides):
            logger.debug(
                "%d) %s (combined score = %0.4f)",
                i + 1,
                vaccine_peptide,
                vaccine_peptide.combined_score)

        return filtered_candidate_vaccine_peptides[:self.max_vaccine_peptides_per_variant]

    def run_isovar(self, variants, alignment_file):
        """
        This function populates a dictionary of Variant objects to a single
        IsovarResult (called `_variant_to_isovar_result_dict`), which will be later
        used to construct VaccinePeptides and compute statistics about the
        variants. If this function has been previously called, the result will
        be cached.

        Parameters
        ----------
        variants : varcode.VariantCollection or str

        alignment_file : pysam.AlignmentFile or str
        """
        if self._isovar_protein_sequence_dict is None:
            # total number of amino acids is the vaccine peptide length plus the
            # number of off-center windows around the mutation
            protein_fragment_sequence_length = (
                self.vaccine_peptide_length + 2 * self.padding_around_mutation)
            """
            These sequences are only the ones that overlap the variant and support the mutation.
            Right now, this generator yields:
            - (variant, mutant protein sequences) if there's enough alt RNA support
            - (variant, None) if the variant is silent or there are ref reads overlapping the
            variant locus but inadequate alt RNA support.
            - does not return the variant if there's no RNA support for ref or alt - we may miss
            some coding variants this way unless we check for them explicitly

            Future intended behavior: returns all passing variants, with a protein sequences
            generator that is non empty if there are enough alt RNA reads supporting the variant
            """

            protein_sequences_generator = run_isovar(
                variants,
                alignment_file,
                read_collector=self.read_collector,
                protein_sequence_creator=self.protein_sequence_creator,
                filter_thresholds=self.filter_thresholds)

                self.reads_generator,
                transcript_id_whitelist=None,
                protein_sequence_length=protein_fragment_sequence_length,
                min_alt_rna_reads=self.min_alt_rna_reads,
                min_variant_sequence_coverage=self.min_variant_sequence_coverage,
                variant_sequence_assembly=self.variant_sequence_assembly,
                max_protein_sequences_per_variant=1)

            self._isovar_protein_sequence_dict = {}
            for variant, isovar_protein_sequences in protein_sequences_generator:
                if len(isovar_protein_sequences) == 0:
                    # variant RNA support is below threshold
                    logger.info("No protein sequences for %s", variant)
                    continue

                # use the first protein sequence - why?
                self._isovar_protein_sequence_dict[variant] = isovar_protein_sequences[0]

        return self._isovar_protein_sequence_dict

    @property
    def vaccine_peptides(self, isovar_results):
        """
        Determine vaccine peptides for each variant in a list of IsovarResult
        objects that have non-coding effects and sufficient RNA to determine
        a mutant sequence.

        Parameters
        ----------
        isovar_results : list of IsovarResult

        Returns
        -------
        Dictionary of varcode.Variant objects to a list of VaccinePeptides.
        """
        vaccine_peptides_dict = OrderedDict()
        for isovar_result in isovar_results:
            vaccine_peptides = self.vaccine_peptides_for_single_variant(
                isovar_result,
                require_mutant_mhc_ligands=True)

            if len(vaccine_peptides) == 0:
                continue

            vaccine_peptides_dict[isovar_result.variant] = vaccine_peptides

        return vaccine_peptides_dict

    def ranked_vaccine_peptides(
            self,
            num_mutant_epitopes_to_keep=10000,
            min_epitope_score=0):
        """
        This function returns a sorted list whose first element is a Variant and whose second
        element is a list of VaccinePeptide objects.
        """
        variant_peptides_dict = self.vaccine_peptides
        result_list = list(variant_peptides_dict.items())
        # TODO: move this sort key into its own function, also make less nuts
        result_list.sort(
            key=lambda x: x[1][0].combined_score if len(x[1]) > 0 else 0.0,
            reverse=True)
        return result_list

    def variant_properties(self, isovar_results):
        """
        Parameters
        ----------
        isovar_results : list of isovar.IsovarResult

        Returns
        -------
        Dictionary from varcode.Variant to dictionary of properties we want to
        analyze later, e.g. whether this variant is part of a pathway of interest,
        is a strong MHC binder, etc.
        """
        variant_properties_dict = OrderedDict()
        for variant in self.variants:
            gene_name = ''
            if variant.gene_names:
                gene_name = variant.effects().top_priority_effect().gene_name
            variant_dict = OrderedDict((
                ('contig', variant.contig),
                ('start', variant.start),
                ('ref', variant.ref),
                ('alt', variant.alt),
                ('is_coding_nonsynonymous', False),
                ('rna_support', False),
                ('mhc_binder', False),
                ('gene_name', gene_name),
            ))
            if self.gene_pathway_check is not None:
                pathway_dict = self.gene_pathway_check.make_variant_dict(variant)
                variant_dict.update(pathway_dict)
            if len(variant.effects().drop_silent_and_noncoding()) > 0:
                variant_dict['is_coding_nonsynonymous'] = True
            if variant in self.isovar_protein_sequence_dict:
                variant_dict['rna_support'] = True
            # TODO: compute MHC binder status for variants that don't have RNA support
            if variant in self.vaccine_peptides:
                variant_dict['mhc_binder'] = True
            variant_properties_dict[variant] = variant_dict
        return list(variant_properties_dict.values())

    def variant_counts(self, isovar_results):
        """
        Parameters
        ----------
        isovar_results : list of isovar.IsovarResult

        Returns
        -------
        Dictionary from keys such as 'num_total_variants' to int
        """
        # dictionary which will contain some overall variant counts for report display
        counts_dict = {
            'num_total_variants': len(isovar_results),
            'num_coding_effect_variants': 0,
            'num_variants_with_rna_support': 0,
            'num_variants_with_vaccine_peptides': 0,
        }
        for isovar_result in isovar_results:
            isovar_result.predicted_effect
            if len(variant.effects().drop_silent_and_noncoding()) > 0:
                counts_dict['num_coding_effect_variants'] += 1
            if variant in self.isovar_protein_sequence_dict:
                counts_dict['num_variants_with_rna_support'] += 1
            if variant in self.vaccine_peptides:
                counts_dict['num_variants_with_vaccine_peptides'] += 1
        return counts_dict
