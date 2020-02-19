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

from varcode.effects import top_priority_effect
from serializable import Serializable

logger = logging.getLogger(__name__)


class MutantProteinFragment(Serializable):
    def __init__(
            self,
            variant,
            gene_name,
            amino_acids,
            mutant_amino_acid_start_offset,
            mutant_amino_acid_end_offset,
            supporting_reference_transcripts,
            n_overlapping_reads,
            n_alt_reads,
            n_ref_reads,
            n_alt_reads_supporting_protein_sequence):
        """
        Parameters
        ----------
        variant : varcode.Variant
            Somatic mutation.

        gene_name : str
            Gene from which we used a transcript to translate this mutation.

        amino_acids : str
            Translated protein sequence, aggregated from possibly multiple
            synonymous coding sequences.

        mutant_amino_acid_start_offset : int
            Starting offset of amino acids which differ due to the mutation

        mutant_amino_acid_end_offset : int
            End offset of amino acids which differ due to the mutation

        supporting_reference_transcripts : list of pyensembl.Transcript
            PyEnsembl Transcript objects for reference transcripts which
            were used to establish the reading frame of coding sequence(s)
            detected from RNA.

        n_overlapping_reads : int
            Number of reads overlapping the variant locus.

        n_alt_reads  : int
            Number of reads supporting the variant.

        n_ref_reads : int
            Number of reads supporting the reference allele.

        n_alt_reads_supporting_protein_sequence : int
            Number of RNA reads fully spanning the cDNA sequence(s) from which
            we translated this amino acid sequence.
        """
        self.variant = variant
        self.gene_name = gene_name
        self.amino_acids = amino_acids
        self.mutant_amino_acid_start_offset = mutant_amino_acid_start_offset
        self.mutant_amino_acid_end_offset = mutant_amino_acid_end_offset
        self.supporting_reference_transcripts = \
            supporting_reference_transcripts
        self.n_overlapping_reads = n_overlapping_reads
        self.n_alt_reads = n_alt_reads
        self.n_ref_reads = n_ref_reads
        self.n_alt_reads_supporting_protein_sequence = \
            n_alt_reads_supporting_protein_sequence

    @classmethod
    def from_isovar_result(cls, isovar_result):
        """
        Create a MutantProteinFragment from an isovar.IsovarResult object

        Parameters
        ----------
        isovar_result : isovar.IsovarResult

        Returns
        -------
        MutantProteinFragment
        """
        protein_sequence = isovar_result.top_protein_sequence
        if protein_sequence is None:
            return None
        return cls(
            variant=isovar_result.variant,
            gene_name=protein_sequence.gene_name,
            amino_acids=protein_sequence.amino_acids,
            mutant_amino_acid_start_offset=protein_sequence.mutation_start_idx,
            mutant_amino_acid_end_offset=protein_sequence.mutation_end_idx,

            # TODO: distinguish reads and fragments in Vaxrank?
            n_overlapping_reads=isovar_result.num_total_fragments,
            n_alt_reads=isovar_result.num_alt_fragments,
            n_ref_reads=isovar_result.num_ref_fragments,
            n_alt_reads_supporting_protein_sequence=protein_sequence.num_supporting_fragments,
            supporting_reference_transcripts=protein_sequence.transcripts)

    def __len__(self):
        return len(self.amino_acids)

    @property
    def n_mutant_amino_acids(self):
        return (
            self.mutant_amino_acid_end_offset - self.mutant_amino_acid_start_offset)

    @property
    def mutation_distance_from_edge(self):
        distance_from_left = self.mutant_amino_acid_start_offset
        distance_from_right = len(self) - self.mutant_amino_acid_end_offset
        return min(distance_from_left, distance_from_right)

    @property
    def is_deletion(self):
        return self.n_mutant_amino_acids == 0 and self.variant.is_deletion

    @property
    def n_other_reads(self):
        """
        Number of reads supporting alleles which are neither ref nor alt
        """
        return self.n_overlapping_reads - (self.n_ref_reads + self.n_alt_reads)

    def interval_overlaps_mutation(self, start_offset, end_offset):
        """
        Does the given start_offset:end_offset interval overlap the mutated
        region of this MutantProteinFragment? Interval offsets are expected
        to be base-0 half-open (start is inclusive, end is exclusive).
        """
        return (
            start_offset < self.mutant_amino_acid_end_offset and
            end_offset > self.mutant_amino_acid_start_offset)

    def generate_subsequences(self, subsequence_length):
        """
        Yields (int, MutantProteinFragment) pairs, where the integer
        indicates the offset into the amino acid sequences.
        """
        n_total_amino_acids = len(self.amino_acids)
        if n_total_amino_acids <= subsequence_length:
            yield (0, self)
        else:
            for subsequence_start_offset in range(
                    0,
                    n_total_amino_acids - subsequence_length + 1):
                subsequence_end_offset = subsequence_start_offset + subsequence_length
                amino_acids = self.amino_acids[
                    subsequence_start_offset:subsequence_end_offset]
                mutant_amino_acid_start_offset = max(
                    0,
                    self.mutant_amino_acid_start_offset - subsequence_start_offset)
                mutant_amino_acid_end_offset = min(
                    len(amino_acids),
                    max(
                        0,
                        self.mutant_amino_acid_end_offset - subsequence_start_offset))
                n_supporting_reads = self.n_alt_reads_supporting_protein_sequence
                subsequence_mutant_protein_fragment = MutantProteinFragment(
                    variant=self.variant,
                    gene_name=self.gene_name,
                    amino_acids=amino_acids,
                    mutant_amino_acid_start_offset=mutant_amino_acid_start_offset,
                    mutant_amino_acid_end_offset=mutant_amino_acid_end_offset,
                    n_overlapping_reads=self.n_overlapping_reads,
                    n_ref_reads=self.n_ref_reads,
                    n_alt_reads=self.n_alt_reads,
                    n_alt_reads_supporting_protein_sequence=n_supporting_reads,
                    supporting_reference_transcripts=self.supporting_reference_transcripts)
                yield subsequence_start_offset, subsequence_mutant_protein_fragment

    def sorted_subsequences(
            self,
            subsequence_length,
            limit=None,
            sort_key=lambda x: (
                -x[1].mutation_distance_from_edge,
                -x[1].n_mutant_amino_acids)):
        """
        Returns subsequences, paired with their offset from the start of the
        protein fragment. The default sort criterion is maximizing the
        mutation distance from the edge of the sequence and secondarily
        maximizing the number of mutant amino acids.
        """
        subsequences = list(self.generate_subsequences(subsequence_length))
        subsequences.sort(key=sort_key)
        if limit:
            subsequences = subsequences[:limit]
        return subsequences

    def predicted_effect(self):
        effects = [
            self.variant.effect_on_transcript(t) for t in
            self.supporting_reference_transcripts
        ]
        predicted_effect = top_priority_effect(effects)
        return predicted_effect

    def global_start_pos(self):
        # position of mutation start relative to the full amino acid sequence
        global_mutation_start_pos = self.predicted_effect().aa_mutation_start_offset
        if global_mutation_start_pos is None:
            logger.error(
                'Could not find mutation start pos for variant %s',
                self.variant)
            return -1

        # get the global position of the mutant protein fragment: shift left by the amount of
        # the relative mutant start position
        return (
            global_mutation_start_pos - self.mutant_amino_acid_start_offset
        )
