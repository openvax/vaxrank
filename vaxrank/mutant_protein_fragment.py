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


import logging
from dataclasses import dataclass
from typing import Any

from serializable import DataclassSerializable
from varcode.effects import top_priority_effect

logger = logging.getLogger(__name__)


def _find_mutation_region(ref_protein, mut_protein):
    """
    Find the mutated region by diffing reference and mutant protein sequences.
    Returns (start, end) offsets in the mutant protein (0-based, half-open).
    """
    min_len = min(len(ref_protein), len(mut_protein))

    # Find first difference
    start = min_len
    for i in range(min_len):
        if ref_protein[i] != mut_protein[i]:
            start = i
            break

    # Find where sequences re-align from the right
    ref_i = len(ref_protein)
    mut_i = len(mut_protein)
    while ref_i > start and mut_i > start:
        if ref_protein[ref_i - 1] == mut_protein[mut_i - 1]:
            ref_i -= 1
            mut_i -= 1
        else:
            break

    # Proteins differ only in length (e.g. stop-loss extension)
    if start == min_len and len(ref_protein) != len(mut_protein):
        return min_len, len(mut_protein)

    return start, mut_i


@dataclass
class MutantProteinFragment(DataclassSerializable):
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

    variant: Any
    gene_name: str
    amino_acids: str
    mutant_amino_acid_start_offset: int
    mutant_amino_acid_end_offset: int
    supporting_reference_transcripts: list
    n_overlapping_reads: int
    n_alt_reads: int
    n_ref_reads: int
    n_alt_reads_supporting_protein_sequence: int

    # True when ``variant.ref`` / ``variant.alt`` are placeholder
    # nucleotides synthesized by the loader because the upstream
    # source (e.g. a 2-part LENS ``chr:pos`` row) didn't supply real
    # alleles. Downstream code that interprets ref/alt as biology —
    # varcode effect annotation, isovar reannotation, anything that
    # touches the variant's transcript context — MUST detect this
    # flag and skip or refuse rather than silently produce nonsense.
    # Construct assembly only uses (chr, pos) and the AA fragment, so
    # is unaffected. Defaults False since the VCF/BAM pipeline always
    # supplies real alleles.
    placeholder_alleles: bool = False

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

    @classmethod
    def from_variant_dna(cls, variant, protein_sequence_length):
        """
        Create a MutantProteinFragment from DNA annotation alone (no RNA),
        using varcode's MutantTranscript.  Used as an opt-in fallback when
        isovar has no RNA support for a variant.

        Transcript selection: among protein-modifying effects, pick the one
        with the longest mutant protein, breaking ties by lex-sorted
        transcript ID.  (Effects are already priority-sorted by varcode.)

        Parameters
        ----------
        variant : varcode.Variant
        protein_sequence_length : int
            Desired length of the protein fragment window around the mutation.

        Returns
        -------
        MutantProteinFragment or None
        """
        from varcode.mutant_transcript import apply_variant_to_transcript

        effects = variant.effects()
        coding_effects = [
            e for e in effects
            if hasattr(e, 'transcript') and e.modifies_protein_sequence
        ]
        if not coding_effects:
            logger.debug(
                "No protein-modifying effects for %s, skipping DNA fallback",
                variant.short_description)
            return None

        # Among same-type effects (same priority tier), prefer longest protein
        # then lex-sorted transcript ID.  variant.effects() is already sorted
        # by priority, so the first element's type defines the top tier.
        best_type = type(coding_effects[0])
        same_tier = [e for e in coding_effects if type(e) is best_type]
        same_tier.sort(key=lambda e: (
            -len(e.transcript.protein_sequence or ''),
            e.transcript.id,
        ))
        best_effect = same_tier[0]
        transcript = best_effect.transcript

        mt = apply_variant_to_transcript(variant, transcript)
        if mt is None or mt.mutant_protein_sequence is None:
            logger.debug(
                "MutantTranscript construction failed for %s on %s",
                variant.short_description, transcript.id)
            return None

        mut_protein = mt.mutant_protein_sequence
        ref_protein = transcript.protein_sequence
        if not ref_protein:
            return None

        # Locate the mutation in the protein
        mut_start, mut_end = _find_mutation_region(ref_protein, mut_protein)
        if mut_start >= len(mut_protein):
            return None

        # Extract a window of protein_sequence_length centered on the mutation
        mutation_midpoint = (mut_start + mut_end) // 2
        half_window = protein_sequence_length // 2
        window_start = max(0, mutation_midpoint - half_window)
        window_end = min(len(mut_protein), window_start + protein_sequence_length)
        window_start = max(0, window_end - protein_sequence_length)

        amino_acids = mut_protein[window_start:window_end]
        local_mut_start = max(0, mut_start - window_start)
        local_mut_end = min(len(amino_acids), mut_end - window_start)

        return cls(
            variant=variant,
            gene_name=transcript.gene.name,
            amino_acids=amino_acids,
            mutant_amino_acid_start_offset=local_mut_start,
            mutant_amino_acid_end_offset=local_mut_end,
            n_overlapping_reads=0,
            n_alt_reads=0,
            n_ref_reads=0,
            n_alt_reads_supporting_protein_sequence=0,
            supporting_reference_transcripts=[transcript],
        )

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

    @property
    def rna_vaf(self):
        """
        RNA variant allele frequency.

        Computed as ``n_alt_reads / (n_alt_reads + n_ref_reads)`` —
        ``n_other_reads`` (reads matching neither ref nor alt) is
        excluded from the denominator so the value reports an
        alt-vs-ref ratio rather than alt-out-of-total-overlapping.
        Matches the convention used by
        ``vaxrank_results.variant_properties`` for JSON/CSV output.
        Differs from ``isovar.IsovarResult.fraction_alt_fragments``,
        which divides by the total fragment count including "other".

        Returns None when no ref+alt reads were observed.
        """
        denom = self.n_alt_reads + self.n_ref_reads
        if denom == 0:
            return None
        return self.n_alt_reads / denom

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
        """Top-priority varcode effect across the supporting transcripts.

        Returns ``None`` when no Transcript objects are available —
        common on external-input paths (LENS / pVACseq) where a
        transcript ID didn't resolve in the configured pyensembl
        release, or for ERV / non-genic antigens that carry no
        transcript context at all. Also returns ``None`` when varcode
        rejects every transcript (e.g. ``ReferenceMismatchError`` —
        the variant was called against a different genome build than
        the configured pyensembl release; common when the user hasn't
        set ``--ensembl-release`` to match LENS / pVACseq's
        reference). Callers (template-report renderers) tolerate
        ``None`` rather than crash.

        Only catches the specific exceptions varcode raises for
        ref/alt mismatches (``ReferenceMismatchError``) and missing /
        invalid transcript data (``ValueError``, ``KeyError``); other
        exceptions propagate so genuine bugs aren't swallowed.
        """
        if not self.supporting_reference_transcripts:
            return None
        # Imported here to keep the module-load cost down and avoid a
        # circular-ish import (varcode is heavy).
        from varcode.errors import ReferenceMismatchError
        effects = []
        for t in self.supporting_reference_transcripts:
            try:
                effects.append(self.variant.effect_on_transcript(t))
            except (ReferenceMismatchError, ValueError, KeyError) as e:
                logger.debug(
                    "varcode.effect_on_transcript failed for %s on %s: "
                    "%s — skipping that transcript.",
                    self.variant, t, e)
        if not effects:
            return None
        return top_priority_effect(effects)

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
