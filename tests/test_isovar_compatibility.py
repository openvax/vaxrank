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

from types import SimpleNamespace

from isovar.allele_read import AlleleRead
from isovar.protein_sequence_creator import ProteinSequenceCreator
from isovar.protein_sequence_helpers import group_equivalent_translations
from isovar.reference_context import ReferenceContext
from isovar.variant_sequence_creator import VariantSequenceCreator
from varcode import Variant

from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.vaccine_antigen import VaccineAntigen
from vaxrank.vaccine_config import VaccineConfig


def _isovar_result(variant, protein_sequence):
    """Minimal Isovar result carrying the fields Vaxrank consumes."""
    n_reads = protein_sequence.num_supporting_reads
    n_fragments = protein_sequence.num_supporting_fragments
    return SimpleNamespace(
        variant=variant,
        top_protein_sequence=protein_sequence,
        num_total_reads=n_reads,
        num_alt_reads=n_reads,
        num_ref_reads=0,
        num_total_fragments=n_fragments,
        num_alt_fragments=n_fragments,
        num_ref_fragments=0,
    )


def _reference_context(variant, prefix, suffix):
    return ReferenceContext(
        strand="+",
        sequence_before_variant_locus=prefix,
        sequence_at_variant_locus=variant.ref,
        sequence_after_variant_locus=suffix,
        offset_to_first_complete_codon=0,
        contains_start_codon=False,
        overlaps_start_codon=False,
        contains_five_prime_utr=False,
        amino_acids_before_variant="",
        variant=variant,
        transcripts=(),
    )


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
    reference_context = _reference_context(
        variant=variant,
        prefix="A" * len(prefixes[-1]),
        suffix="A" * 8,
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

    protein_sequence, = group_equivalent_translations(translations)
    fragment = MutantProteinFragment.from_isovar_result(
        _isovar_result(variant, protein_sequence)
    )
    assert fragment.n_alt_reads_supporting_protein_sequence == 8
    assert fragment.n_alt_fragments_supporting_protein_sequence == 8


def test_multibase_substitution_marks_every_targetable_amino_acid():
    """An alternate interval crossing a codon boundary changes both codons."""
    variant = Variant("1", 100, "GAA", "TCC", "GRCh38")
    prefix = "AAA" * 4 + "AT"
    suffix = "AGGG"
    reads = [
        AlleleRead(prefix=prefix, allele=variant.alt, suffix=suffix, name=str(i))
        for i in range(2)
    ]
    sequences = VariantSequenceCreator().reads_to_variant_sequences(variant, reads)
    translations = ProteinSequenceCreator().all_pairs_translations(
        sequences,
        [_reference_context(variant, prefix, suffix)],
    )
    protein_sequence, = group_equivalent_translations(translations)

    fragment = MutantProteinFragment.from_isovar_result(
        _isovar_result(variant, protein_sequence)
    )
    antigen = VaccineAntigen.from_mutant_protein_fragment(fragment)

    assert fragment.amino_acids == "KKKKIPG"
    assert (
        fragment.mutant_amino_acid_start_offset,
        fragment.mutant_amino_acid_end_offset,
    ) == (4, 6)
    assert antigen.interval_is_targetable(5, 6)


def test_vaxrank_length_keeps_context_for_long_alternate(monkeypatch):
    """A long alternate must not consume context required to establish its ORF."""
    variant = Variant("1", 100, "", "A" * 90, "GRCh38")
    prefix = "ACG" * 4
    suffix = "G" * 30
    reads = [
        AlleleRead(prefix=prefix, allele=variant.alt, suffix=suffix, name=str(i))
        for i in range(2)
    ]
    reference_context = _reference_context(variant, prefix, suffix)
    monkeypatch.setattr(
        "isovar.protein_sequence_creator.reference_contexts_for_variant",
        lambda *args, **kwargs: [reference_context],
    )
    vaccine_config = VaccineConfig()
    # Match the translation length the CLI derives from Vaxrank's defaults.
    protein_sequence_length = (
        vaccine_config.preferred_peptide_length
        + 2 * vaccine_config.padding_around_mutation
    )
    creator = ProteinSequenceCreator(protein_sequence_length=protein_sequence_length)

    translations = creator.translate_variant_reads(variant, reads)
    protein_sequence, = group_equivalent_translations(translations)
    fragment = MutantProteinFragment.from_isovar_result(
        _isovar_result(variant, protein_sequence)
    )

    assert len(fragment.amino_acids) == 35
    assert fragment.mutant_amino_acid_start_offset == 3
    assert fragment.mutant_amino_acid_end_offset == 33
