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

from collections import namedtuple

# using a namedtuple base class for the immutable fields of a MutantProteinFragment
# since it makes it clearer what the essential information is and provides
# useful comparison/hashing methods

MutantProteinFragment_Fields = namedtuple("MutantProteinFragment", (
    # varcode.Variant
    "variant",
    # gene and transcript(s) which were used to translate one or more
    # variant cDNA sequences into the following amino acids
    "gene_name",

    ###
    # Translated protein sequence, aggregated from possibly multiple
    # synonymous coding sequences
    ###

    "amino_acids",
    # offsets of amino acids which differ due to the mutation
    "mutant_amino_acid_start_offset",
    "mutant_amino_acid_end_offset",

    ###
    # RNA evidence
    ###

    # number of reads overlapping the variant locus
    "n_overlapping_reads",
    # number of reads supporting the variant
    "n_alt_reads",
    # number of reads supporting the reference allele
    "n_ref_reads",

    # number of RNA reads fully spanning the cDNA sequence(s) from which we
    # translated this amino acid sequence.
    "n_alt_reads_supporting_protein_sequence",
))

class MutantProteinFragment(MutantProteinFragment_Fields):
    @classmethod
    def from_isovar_protein_sequence(cls, protein_sequence):
        return cls(
            variant=protein_sequence.variant,
            gene_name=";".join(protein_sequence.gene),
            amino_acids=protein_sequence.amino_acids,
            mutant_amino_acid_start_offset=protein_sequence.variant_aa_interval_start,
            mutant_amino_acid_end_offset=protein_sequence.variant_aa_interval_end,
            n_overlapping_reads=len(protein_sequence.overlapping_reads),
            n_alt_reads=len(protein_sequence.alt_reads),
            n_ref_reads=len(protein_sequence.ref_reads),
            n_alt_reads_supporting_protein_sequence=len(
                protein_sequence.alt_reads_supporting_protein_sequence))

    @property
    def n_mutant_amino_acids(self):
        return self.mutant_amino_acid_end_offset - self.mutant_amino_acid_start_offset

    @property
    def is_deletion(self):
        return self.n_mutant_amino_acids == 0 and self.variant.is_deletion

    @property
    def n_other_reads(self):
        """
        Number of reads supporting alleles which are neither ref nor alt
        """
        return self.n_overlapping_reads - (self.n_ref_reads + self.n_alt_reads)
