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

def create_variant_to_protein_sequence_dict(
        vaccine_peptide_length=25,):
    """
    This function computes a dictionary of Variant objects to a
    single isovar protein sequence that will be used to try to construct
    VaccinePeptides.

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

    """


    variant_to_protein_sequence_dict = {}
    # total number of amino acids is the vaccine peptide length plus the
    # number of off-center windows around the mutation
    protein_fragment_sequence_length = (
            vaccine_peptide_length + 2 * padding_around_mutation)

    # These sequences are only the ones that overlap the variant and support
    # the mutation.
    # Right now, this generator yields:
    # -  (variant, mutant protein sequences) if there's enough alt RNA support
    #
    # -  (variant, None) if the variant is silent or there are ref reads overlapping the
    #    variant locus but inadequate alt RNA support.
    # -  does not return the variant if there's no RNA support for ref or alt,
    #
    #    we may miss some coding variants this way unless we check for them
    #    explicitly
    #
    # Future intended behavior: returns all passing variants, with a protein
    # sequences generator that is non empty if there are enough alt RNA reads
    # supporting the variant
    protein_sequences_generator = reads_generator_to_protein_sequences_generator(
        reads_generator,
        transcript_id_whitelist=None,
        protein_sequence_length=protein_fragment_sequence_length,
        min_alt_rna_reads=min_alt_rna_reads,
        min_variant_sequence_coverage=min_variant_sequence_coverage,
        variant_sequence_assembly=variant_sequence_assembly,
        max_protein_sequences_per_variant=1)

    for variant, isovar_protein_sequences in protein_sequences_generator:
        if len(isovar_protein_sequences) == 0:
            # variant RNA support is below threshold
            logger.info("No protein sequences for %s", variant)
            continue

        # use the first protein sequence - why?
        variant_to_protein_sequence_dict[variant] = \
            isovar_protein_sequences[0]
    return variant_to_protein_sequence_dict
