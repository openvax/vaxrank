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

"""Downstream contracts required from Vaxrank's Isovar dependency."""

from isovar.allele_read import AlleleRead
from isovar.protein_sequence_creator import ProteinSequenceCreator
from isovar.reference_context import ReferenceContext
from isovar.variant_sequence_creator import VariantSequenceCreator
from varcode import Variant


def test_contained_rna_assembly_reaches_transcript_matching():
    """A longer incompatible assembly must not erase a compatible core."""
    variant = Variant("1", 100, "G", "C", "GRCh38")
    prefixes = ("AAA", "GGGAAA", "CCCGGGAAA", "TTTCCCGGGAAA")
    reads = [
        AlleleRead(
            prefix=prefix,
            allele=variant.alt,
            suffix="A" * 8,
            name=f"{prefix}-{i}",
        )
        for prefix in prefixes
        for i in range(2)
    ]

    sequence_creator = VariantSequenceCreator(
        min_variant_sequence_coverage=2,
        preferred_sequence_length=len(prefixes[-1]) * 2 + 1,
        variant_sequence_assembly=True,
        min_assembly_overlap_size=1,
    )
    candidates = sequence_creator.reads_to_variant_sequences(
        variant=variant,
        reads=reads,
    )
    reference_context = ReferenceContext(
        strand="+",
        sequence_before_variant_locus="A" * len(prefixes[-1]),
        sequence_at_variant_locus=variant.ref,
        sequence_after_variant_locus="A" * 8,
        offset_to_first_complete_codon=0,
        contains_start_codon=False,
        overlaps_start_codon=False,
        contains_five_prime_utr=False,
        amino_acids_before_variant="",
        variant=variant,
        transcripts=(),
    )
    protein_creator = ProteinSequenceCreator(
        protein_sequence_length=8,
        min_variant_sequence_coverage=2,
        min_transcript_prefix_length=3,
        max_transcript_mismatches=2,
    )

    translations = protein_creator.all_pairs_translations(
        variant_sequences=candidates,
        reference_contexts=[reference_context],
    )

    assert {candidate.prefix for candidate in candidates} == set(prefixes)
    # Isovar 1.7.9 also lets longer candidates trim to this same compatible
    # core. The downstream contract is that the original core reaches
    # translation exactly once; equivalent translations may accompany it and
    # are deduplicated later when Isovar groups protein sequences.
    assert {translation.amino_acids for translation in translations} == {"KQKK"}
    core_translation, = [
        translation
        for translation in translations
        if translation.untrimmed_variant_sequence.prefix == "AAA"
    ]
    assert core_translation.contains_mutation
