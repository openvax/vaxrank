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
from topiary import FRAGMENTS

from .varcode_effects import (
    OUTCOME_SELECTION_HIGHEST_PRIORITY,
    OUTCOME_SELECTION_MOST_LIKELY,
    OUTCOME_SELECTION_MULTI_OUTCOME,
    select_varcode_effect_outcome,
)

logger = logging.getLogger(__name__)


def find_mutation_region(reference_protein, mutant_protein):
    """Return the differing region between reference and mutant proteins.

    The returned offsets use zero-based, half-open coordinates in the mutant
    protein. Insertions span their inserted residues; deletions have an empty
    span at the deletion point.
    """
    min_len = min(len(reference_protein), len(mutant_protein))

    # Find first difference
    start = min_len
    for i in range(min_len):
        if reference_protein[i] != mutant_protein[i]:
            start = i
            break

    # Find where sequences re-align from the right
    ref_i = len(reference_protein)
    mut_i = len(mutant_protein)
    while ref_i > start and mut_i > start:
        if reference_protein[ref_i - 1] == mutant_protein[mut_i - 1]:
            ref_i -= 1
            mut_i -= 1
        else:
            break

    # Proteins differ only in length (e.g. stop-loss extension)
    if start == min_len and len(reference_protein) != len(mutant_protein):
        return min_len, len(mutant_protein)

    return start, mut_i


# Read-count vocabulary, all of it topiary's. Two orthogonal axes:
#
#   rna_evidence_method    how the number was obtained -- rna_reads,
#                        rna_depth_x_vaf, cds_overlap_reads, ...
#   rna_evidence_subject   what it counted -- reads or fragments
#
# vaxrank asked for an ``rna_fragments`` *method* and that was the wrong
# shape (openvax/topiary#239): a subject is not a derivation, and a term
# alongside the methods would have left n_alt_reads ambiguous on the
# isovar path while looking resolved. ``rna_reads`` deliberately fixes no
# subject -- it says a count came from an alignment, not whether the
# aligner counted reads or fragments -- so the producer states the
# subject separately.
#
# Read from topiary rather than restated, so a rename upstream fails at
# import instead of stamping a term nothing recognizes.
def _sequence_source(name):
    """Return a topiary sequence-source term, failing loudly if it is gone.

    Reading the name from topiary rather than restating it, so a rename
    upstream surfaces here as an import-time error instead of silently
    stamping fragments with a term nothing recognizes.
    """
    from topiary import SEQUENCE_SOURCES

    if name not in SEQUENCE_SOURCES:
        raise ValueError(
            "topiary.SEQUENCE_SOURCES no longer defines %r, which vaxrank "
            "stamps on the fragments it builds" % name)
    return name


SEQUENCE_SOURCE_ISOVAR_ASSEMBLY = _sequence_source("isovar_assembly")
SEQUENCE_SOURCE_VARCODE_TRANSLATION = _sequence_source("varcode_translation")


def _rna_evidence_method(name):
    """Return a topiary RNA-evidence method term, failing loudly if it is gone.

    The guard earned itself: ``rna_reads`` was renamed to ``rna_alignment``
    in topiary 5.45.0, and this raised at import across the whole test suite
    rather than stamping fragments with a term nothing recognizes.
    """
    from topiary import READ_COUNT_METHODS

    if name not in READ_COUNT_METHODS:
        raise ValueError(
            "topiary.READ_COUNT_METHODS no longer defines %r, which vaxrank "
            "stamps on the fragments it builds" % name)
    return name


# Counted from an alignment. Says how, not what — the subject is separate,
# because an aligner can count either reads or fragments.
RNA_EVIDENCE_METHOD_ALIGNMENT = _rna_evidence_method("rna_alignment")


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

    # How the read counts above were obtained, in topiary's vocabulary:
    # "rna_reads" (counted from an alignment), "rna_depth_x_vaf" (depth
    # times VAF, rounded — an estimate), "cds_overlap_reads" (counted, but
    # of CDS overlap rather than variant support), "tpm_x_dna_vaf" (an
    # expression proxy). Blank when the source did not say.
    #
    # The number alone cannot distinguish a measurement from an estimate,
    # and the combined-score DSL weights them identically —
    # sqrt(n_alt_reads) does not know that 51 was arithmetic rather than
    # observation. Carrying the label is what lets a reader tell.
    rna_evidence_method: str = ""

    # What those counts counted: topiary's FRAGMENTS or READS, blank when
    # the source did not say. Separate from the method because how a number
    # was obtained and what it counted are independent — an alignment-derived
    # count may be either. Five fragments and five reads are different bars,
    # and converting between them needs library information no source
    # carries, so this is stated rather than normalized.
    rna_evidence_subject: str = ""

    # The same counts in fragments, where the source has them. isovar
    # reports every count in both units; vaxrank stored only the fragment
    # numbers, in the fields named for reads, and then described that as a
    # subject the caller had to ask about. The reads were never
    # unavailable — carrying both is what makes the subject a fact about
    # the source rather than a workaround for a lie (openvax/topiary#239).
    n_overlapping_fragments: Any = None
    n_alt_fragments: Any = None
    n_ref_fragments: Any = None
    n_alt_fragments_supporting_protein_sequence: Any = None

    # How the protein sequence was obtained: "lens_pep_context",
    # "pvacseq_epitope", "varcode_translation", "isovar_assembly",
    # "caller_supplied". Distinct from the antigen's biological kind —
    # this is method, not biology. "Assembled from RNA" versus "translated
    # from reference" is the first question when auditing a ranking, and
    # was previously answerable only for a whole file.
    sequence_source: str = ""

    @staticmethod
    def slp_window_around_mutation(
            protein_aa, mut_start, mut_end, target_length):
        """Center a ``target_length``-aa SLP window on the mutation span.

        Returns ``(windowed_aa, new_mut_start, new_mut_end)``. When
        ``protein_aa`` is already at-or-below ``target_length``,
        returns the inputs unchanged. Mutation spans wider than the
        target are anchored on the mutation start and the trailing
        flank is truncated.

        The pipeline path (Isovar / varcode) gets correctly-sized
        fragments by construction — they're built from a 25mer-wide
        RNA-derived translation. External-input paths inherit
        whatever protein-context length the upstream tool chose
        (LENS's ``pep_context`` is sometimes a 100+aa protein
        prefix), so they call this helper to land on a fragment of
        the same canonical size. Both paths produce a
        ``MutantProteinFragment`` of the same shape afterwards.
        """
        if not protein_aa or target_length <= 0:
            return protein_aa, mut_start, mut_end
        if len(protein_aa) <= target_length:
            return protein_aa, mut_start, mut_end
        mut_len = max(1, mut_end - mut_start)
        if mut_len >= target_length:
            win_start = max(0, mut_start)
            win_end = min(len(protein_aa), win_start + target_length)
        else:
            mut_mid = (mut_start + mut_end) // 2
            half = target_length // 2
            win_start = max(0, mut_mid - half)
            win_end = win_start + target_length
            if win_end > len(protein_aa):
                win_end = len(protein_aa)
                win_start = max(0, win_end - target_length)
        windowed = protein_aa[win_start:win_end]
        new_start = max(0, mut_start - win_start)
        new_end = min(len(windowed), mut_end - win_start)
        return windowed, new_start, new_end

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

            # isovar reports every count in both units, so carry both.
            # vaxrank used to store the fragment numbers in the read fields
            # and call the discrepancy a subject the caller had to ask
            # about — but the reads were there the whole time.
            n_overlapping_reads=isovar_result.num_total_reads,
            n_alt_reads=isovar_result.num_alt_reads,
            n_ref_reads=isovar_result.num_ref_reads,
            n_alt_reads_supporting_protein_sequence=protein_sequence.num_supporting_reads,
            n_overlapping_fragments=isovar_result.num_total_fragments,
            n_alt_fragments=isovar_result.num_alt_fragments,
            n_ref_fragments=isovar_result.num_ref_fragments,
            n_alt_fragments_supporting_protein_sequence=protein_sequence.num_supporting_fragments,
            # Counted from an alignment; the subject says which unit the
            # n_rna_* accessors report, which is fragments when both exist.
            rna_evidence_method=RNA_EVIDENCE_METHOD_ALIGNMENT,
            rna_evidence_subject=FRAGMENTS,
            sequence_source=SEQUENCE_SOURCE_ISOVAR_ASSEMBLY,
            supporting_reference_transcripts=protein_sequence.transcripts)

    @classmethod
    def from_variant_dna(cls, variant, protein_sequence_length):
        """
        Create a MutantProteinFragment from DNA annotation alone (no RNA),
        using varcode's MutantTranscript.  Used as an opt-in fallback when
        isovar has no RNA support for a variant.

        Transcript selection: ask varcode 5 for splice outcome sets, collapse
        each to its highest-priority concrete outcome, then choose the most
        protein-disruptive effect. Within that tier, pick the longest mutant
        protein, breaking ties by lex-sorted transcript ID.

        Parameters
        ----------
        variant : varcode.Variant
        protein_sequence_length : int
            Desired length of the protein fragment window around the mutation.

        Returns
        -------
        MutantProteinFragment or None
        """
        from varcode.effects.effect_ordering import effect_priority
        from varcode.mutant_transcript import apply_variant_to_transcript

        effects = variant.effects(splice_outcomes=True)
        coding_effects = [
            select_varcode_effect_outcome(e, OUTCOME_SELECTION_HIGHEST_PRIORITY)
            for e in effects
        ]
        coding_effects = [
            e for e in coding_effects
            if (
                e is not None and
                hasattr(e, 'transcript') and
                e.modifies_protein_sequence)
        ]
        if not coding_effects:
            logger.debug(
                "No protein-modifying effects for %s, skipping DNA fallback",
                variant.short_description)
            return None

        # Prefer the most protein-disruptive concrete outcome, then choose
        # the longest protein / lex-sorted transcript within that tier.
        best_priority = max(effect_priority(e) for e in coding_effects)
        same_tier = [
            e for e in coding_effects
            if effect_priority(e) == best_priority
        ]
        same_tier.sort(key=lambda e: (
            -len(
                getattr(e, 'mutant_protein_sequence', None) or
                e.transcript.protein_sequence or ''),
            e.transcript.id,
        ))
        best_effect = same_tier[0]
        transcript = best_effect.transcript

        mut_protein = getattr(best_effect, 'mutant_protein_sequence', None)
        if mut_protein is None:
            mt = getattr(best_effect, 'mutant_transcript', None)
            if mt is not None:
                mut_protein = mt.mutant_protein_sequence
        if mut_protein is None:
            mt = apply_variant_to_transcript(variant, transcript)
            if mt is not None:
                mut_protein = mt.mutant_protein_sequence
        if mut_protein is None:
            logger.debug(
                "MutantTranscript construction failed for %s on %s",
                variant.short_description, transcript.id)
            return None

        ref_protein = transcript.protein_sequence
        if not ref_protein:
            return None

        # Locate the mutation in the protein
        mut_start, mut_end = find_mutation_region(ref_protein, mut_protein)
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
            # No RNA was consulted, so there is nothing to label: zero
            # here means "no reads counted", and rna_evidence_method stays
            # blank rather than claiming a derivation produced the zeros.
            n_overlapping_reads=0,
            n_alt_reads=0,
            n_ref_reads=0,
            n_alt_reads_supporting_protein_sequence=0,
            sequence_source=SEQUENCE_SOURCE_VARCODE_TRANSLATION,
            supporting_reference_transcripts=[transcript],
        )

    # ── RNA evidence, in whichever unit the source actually has ──────────
    #
    # Fragments where a source reports them, reads otherwise. Fragments win
    # because two mates of one fragment are one piece of evidence, so a
    # read count overstates paired-end support by roughly two — the wrong
    # thing to feed a diminishing-returns weight like sqrt().
    #
    # These replace a ``count_in(name, subject)`` accessor, which existed
    # only because fragments were stored under the reads name and something
    # had to explain that reads were unavailable. They were not
    # unavailable. Naming both is what isovar always did.

    def _rna(self, fragments_value, reads_value):
        return fragments_value if fragments_value is not None else reads_value

    @property
    def n_rna_alt(self):
        """Alt-supporting RNA evidence, fragments preferred."""
        return self._rna(self.n_alt_fragments, self.n_alt_reads)

    @property
    def n_rna_ref(self):
        """Reference-supporting RNA evidence, fragments preferred."""
        return self._rna(self.n_ref_fragments, self.n_ref_reads)

    @property
    def n_rna_overlapping(self):
        """RNA evidence covering the position, fragments preferred."""
        return self._rna(self.n_overlapping_fragments,
                         self.n_overlapping_reads)

    @property
    def n_rna_supporting_protein_sequence(self):
        """RNA evidence spanning the assembled sequence, fragments preferred."""
        return self._rna(self.n_alt_fragments_supporting_protein_sequence,
                         self.n_alt_reads_supporting_protein_sequence)

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

    def predicted_effect(
            self,
            outcome_selection=OUTCOME_SELECTION_HIGHEST_PRIORITY):
        """Top-priority varcode effect across the supporting transcripts.

        Varcode 5 can represent splice-disrupting variants as multi-outcome
        sets. By default vaxrank collapses those sets to the highest-priority
        concrete outcome so peptide mechanics keep working with a single
        protein effect. Use ``outcome_selection="most_likely"`` for producer
        order, or ``outcome_selection="multi_outcome"`` for report metadata.

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
        if outcome_selection not in {
                OUTCOME_SELECTION_HIGHEST_PRIORITY,
                OUTCOME_SELECTION_MOST_LIKELY,
                OUTCOME_SELECTION_MULTI_OUTCOME}:
            select_varcode_effect_outcome(None, outcome_selection)
        if not self.supporting_reference_transcripts:
            return None
        # Imported here to keep the module-load cost down and avoid a
        # circular-ish import (varcode is heavy).
        from varcode.errors import ReferenceMismatchError
        from varcode.splice_outcomes import enumerate_splice_outcomes
        effects = []
        for t in self.supporting_reference_transcripts:
            try:
                effect = self.variant.effect_on_transcript(t)
                effect = enumerate_splice_outcomes(effect)
                effect = select_varcode_effect_outcome(
                    effect,
                    outcome_selection=outcome_selection)
                if effect is not None:
                    effects.append(effect)
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
