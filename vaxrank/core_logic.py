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

from isovar.protein_sequences import (
    reads_generator_to_protein_sequences_generator,
)



logger = logging.getLogger(__name__)


class VaxrankCoreLogic(object):
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
        Parameters
        ----------
        variants : VariantCollection
            Variant objects to evaluate for vaccine inclusion

        reads_generator : generator
            Yields sequence of varcode.Variant objects paired with sequences of
            AlleleRead objects that support that variant.

        mhc_predictor : mhctools.BasePredictor
            Object with predict_peptides method, used for making pMHC binding
            predictions

        vaccine_peptide_length : int
            Length of vaccine SLP to construct

        padding_around_mutation : int
            Number of off-center windows around the mutation to consider as vaccine
            peptides.

        max_vaccine_peptides_per_variant : int
            Number of vaccine peptides to generate for each mutation.

        min_alt_rna_reads : int
            Drop variant sequences at loci with fewer than this number of reads
            supporting the alt allele.

        min_variant_sequence_coverage : int
            Trim variant sequences to positions supported by at least this number
            of RNA reads.

        variant_sequence_assembly : int
            If True, then assemble variant cDNA sequences based on overlap of RNA
            reads. If False, then variant cDNA sequences must be fully spanned and
            contained within RNA reads.

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

        # will be a dictionary: varcode.Variant -> isovar protein sequence object
        self._isovar_protein_sequence_dict = None

        # will be a dictionary: varcode.Variant -> list(VaccinePeptide)
        self._vaccine_peptides_dict = None


    @property
    def isovar_protein_sequence_dict(self):
        """
        This function computes a dictionary of Variant objects to a
        single isovar protein sequence that will be used to try to construct
        VaccinePeptides. If this function has been previously
        called, the result will be cached.
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
            protein_sequences_generator = reads_generator_to_protein_sequences_generator(
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
