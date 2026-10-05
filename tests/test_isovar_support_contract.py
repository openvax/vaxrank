"""Original BAM evidence keeps assembly counts separate from catalog scores."""

import json

import pysam
import pytest
from isovar import (
    __version__ as isovar_version, IsovarResult, ProteinSequenceCreator,
    ReadCollector, export_protein_hypotheses,
)
from isovar.dna import reverse_complement_dna
from pyensembl import Genome
from varcode import Variant

from vaxrank.epitope_dataset import EpitopeDataset
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.native_serialization import from_native_json, to_native_json
from .test_core_logic_config import _make_epitope


CDNA = "ATG" + "GCT" * 4 + "AAA" + "GGC" * 4 + "TAA"
MUTANT = CDNA[:15] + "C" + CDNA[16:]


@pytest.fixture(params=["+", "-"])
def reference(tmp_path, request):
    strand = request.param
    attrs = ('gene_id "g"; gene_name "G"; transcript_id "t"; transcript_name "T"; '
             'exon_number "1"; exon_id "e"; protein_id "p"; '
             'gene_biotype "protein_coding"; transcript_biotype "protein_coding";')
    features = [("gene", 1001, 1033), ("transcript", 1001, 1033), ("exon", 1001, 1033),
                ("CDS", 1001, 1030), ("start_codon", 1001, 1003), ("stop_codon", 1031, 1033)]
    if strand == "-":
        features = [(name, 2034 - end, 2034 - start) for name, start, end in features]
    gtf = tmp_path / "ref.gtf"
    gtf.write_text("".join(f"1\ttest\t{name}\t{start}\t{end}\t.\t{strand}\t0\t{attrs}\n"
                           for name, start, end in features))
    cdna, protein = tmp_path / "cdna.fa", tmp_path / "protein.fa"
    cdna.write_text(">t\n" + CDNA + "\n")
    protein.write_text(">p transcript:t\nMAAAAKGGGG\n")
    genome = Genome(reference_name="support-test-" + strand, annotation_name="test",
                    annotation_version=1, gtf_path_or_url=str(gtf),
                    transcript_fasta_paths_or_urls=[str(cdna)],
                    protein_fasta_paths_or_urls=[str(protein)],
                    cache_directory_path=str(tmp_path / "cache"))
    genome.index()
    return strand, Variant("1", 1016 if strand == "+" else 1018,
                           "A" if strand == "+" else "T",
                           "C" if strand == "+" else "G", ensembl=genome)


def reconstruct(tmp_path, reference, sequences):
    """Collect actual source-scoped reads through Isovar's public APIs."""
    strand, variant = reference
    header = pysam.AlignmentHeader.from_dict({
        "SQ": [{"SN": "1", "LN": 20000}], "RG": [{"ID": "rg", "SM": "sample"}]})
    path = tmp_path / "reads.bam"
    with pysam.AlignmentFile(path, "wb", header=header) as output:
        for index, (sequence, operations, offset) in enumerate(sequences):
            if strand == "-":
                sequence, operations = reverse_complement_dna(sequence), operations[::-1]
                offset = len(CDNA) - offset - sum(n for n, op in operations if op != "I")
            read = pysam.AlignedSegment(header)
            read.query_name = f"fragment-{index}"
            read.query_sequence = sequence
            read.flag = 16 if strand == "-" else 0
            read.reference_id, read.reference_start = 0, 1000 + offset
            read.mapping_quality = 60
            read.cigarstring = "".join(f"{n}{op}" for n, op in operations)
            read.set_tag("RG", "rg")
            output.write(read)  # Missing QUAL remains unassessed.
    pysam.sort("-o", str(path), str(path))
    pysam.index(str(path))
    with pysam.AlignmentFile(path) as bam:
        evidence = ReadCollector().read_evidence_for_variant(variant, bam)
        creator = ProteinSequenceCreator(protein_sequence_length=20,
            min_transcript_prefix_length=3, min_variant_sequence_coverage=2,
            max_protein_sequences_per_variant=0, protein_sequence_preference="context")
        proteins = creator.sorted_protein_sequences_for_variant(variant, evidence)
    return IsovarResult(variant, evidence, predicted_effect=None,
                        sorted_protein_sequences=proteins,
                        protein_sequence_settings=creator.settings())


@pytest.mark.parametrize("operation", ["I", "D"])
@pytest.mark.parametrize("position", [6, 24])
def test_flanking_indel_support_is_additive_without_relabeling_assembly_counts(
        tmp_path, reference, operation, position):
    sequence = (MUTANT[:position] + "A" + MUTANT[position:] if operation == "I"
                else MUTANT[:position] + MUTANT[position + 1:])
    operations = [(position, "M"), (1, operation),
                  (len(MUTANT) - position - (operation == "D"), "M")]
    result = reconstruct(tmp_path, reference, [
        (MUTANT, [(len(MUTANT), "M")], 0),
        (MUTANT, [(len(MUTANT), "M")], 0), (sequence, operations, 0)])
    fragment = MutantProteinFragment.from_isovar_result(result)
    raw = export_protein_hypotheses([result], sample_id="sample", source="synthetic-bam")
    scored = export_protein_hypotheses([result], sample_id="sample", source="synthetic-bam",
                                       partial_read_support=True)
    event, = scored["events"]
    row, = event["protein_hypotheses"]
    assert fragment.amino_acids == "MAAAAQGGGG"
    assert fragment.n_alt_reads == fragment.n_alt_fragments == 3
    assert fragment.n_rna_supporting_protein_sequence == 2
    assert row["partial_read_support"]["fractional_fragments"] == 3
    assert event["partial_read_support"]["policy"]["comparison"] == "anchored_edit_compatibility.v2"
    assert {key: value for key, value in row.items() if key != "partial_read_support"} == raw[
        "events"][0]["protein_hypotheses"][0]
    assert MutantProteinFragment.from_isovar_result(result) == fragment
    assert fragment.sequence_source_version == isovar_version
    restored = from_native_json(to_native_json(fragment), MutantProteinFragment)
    assert restored == fragment
    assert restored.n_rna_supporting_protein_sequence == 2


def test_ambiguous_partial_fragments_do_not_become_fractional_raw_read_counts(tmp_path, reference):
    alternative = MUTANT[:24] + "TTC" + MUTANT[27:]
    sequences = [(seq, [(len(seq), "M")], 0) for seq in (MUTANT, MUTANT, alternative, alternative)]
    sequences += [(MUTANT[9:22], [(13, "M")], 9)] * 3
    result = reconstruct(tmp_path, reference, sequences)
    fragment = MutantProteinFragment.from_isovar_result(result)
    export = export_protein_hypotheses([result], sample_id="sample", source="synthetic-bam",
        partial_read_support=True, support_max_edits=0)
    scores = export["events"][0]["partial_read_support"]
    # Discovery also retains the shorter shared local protein, rather than
    # pretending the two full proteins are the complete candidate catalog.
    assert len(scores["protein_scores"]) == 3
    assert scores["scored_fragments"] == 7
    assert sorted(r["fractional_fragments"] for r in scores["protein_scores"]) == [2., 2., 3.]
    assert any(r["protein_weight"] == pytest.approx(1 / 3) for r in scores["fragments"])
    assert fragment.n_alt_reads == fragment.n_alt_fragments == 7
    assert fragment.n_rna_supporting_protein_sequence == 5
    epitope = _make_epitope(fragment.amino_acids[1:10], ic50=100., wt_ic50=1000.,
                            source_sequence=fragment.amino_acids, offset=1)
    dataset = EpitopeDataset.from_predictions([epitope])
    dataset.mutation_fragments[epitope.prediction_id] = fragment
    path = tmp_path / "native.tsv"
    dataset.save(path)
    reloaded = EpitopeDataset.load(path).mutation_fragments[epitope.prediction_id]
    assert reloaded.n_rna_supporting_protein_sequence == 5
    assert reloaded.sequence_source_version == isovar_version
    # An older fragment has genuinely unknown producer provenance.
    legacy = json.loads(to_native_json(fragment))
    legacy.pop("sequence_source_version")
    assert from_native_json(json.dumps(legacy), MutantProteinFragment).sequence_source_version == ""
