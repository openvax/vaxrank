"""Independent Sid construct definitions and shared original-RNA inputs."""

from dataclasses import replace
from hashlib import sha256
import json
from pathlib import Path

from vaxrank.construct_sequence import ConstructEvidence, ConstructSequence, ConstructSequenceEdit
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.vaccine_antigen import VaccineAntigen

from .osteosarc_selection_helpers import reconstruct_selection


DATA = Path(__file__).parent / "data" / "osteosarc" / "construct_audit"
DOCUMENTED = json.loads((DATA / "documented.json").read_text())
RNA = json.loads((DATA / "rna" / "manifest.json").read_text())


def verify_rna_inputs():
    for filename, digest in RNA["files"].items():
        if sha256((DATA / "rna" / filename).read_bytes()).hexdigest() != digest:
            raise ValueError("Sid RNA fixture checksum mismatch: " + filename)


def reconstruct_source(inputs, case, length=25):
    return reconstruct_selection(inputs, DOCUMENTED["variant_id"], length,
                                 bam_path=DATA / "rna" / case["bam"])[0]


def antigen_from_source(result, source_id, *, source_scope="selected_fixture"):
    fragment = MutantProteinFragment.from_isovar_result(result)
    antigen = VaccineAntigen.from_mutant_protein_fragment(fragment)
    # Human provenance is explicit in the published GRCh38 reference, not a
    # guess from gene naming or peptide identity. Keep RNA context unchanged.
    return replace(antigen, species="Homo sapiens",
                   source_identifier=antigen.source_identifier + ";rna_source=" + source_id
                   + ";rna_scope=" + source_scope)


def construct_from_result(result, record, source_id, *, source_scope="selected_fixture"):
    """Map a documented native window; never align or trim away a disagreement."""
    fragment = MutantProteinFragment.from_isovar_result(result)
    native = record["native_sequence"]
    if fragment.amino_acids.count(native) != 1:
        raise ValueError("Documented native window is absent or non-unique in RNA context")
    start = fragment.amino_acids.index(native)
    end = start + len(native)
    if [fragment.mutant_amino_acid_start_offset - start,
            fragment.mutant_amino_acid_end_offset - start] != record["native_mutation_interval"]:
        raise ValueError("RNA target coordinates disagree with the documented target")
    evidence = ConstructEvidence(
        source=DOCUMENTED["source"], evidence_level="documented",
        description="Published DYNC1H1 vaccine amino-acid sequence; chemistry unresolved",
        provider=record["provider"], vaccine_version=record["vaccine_version"] or "")
    rationale = record.get("addition_rationale")
    rationale = ConstructEvidence(**rationale, provider=record["provider"]) if rationale else None
    edits = tuple(ConstructSequenceEdit(offset, offset, addition, evidence, rationale)
                  for offset, addition in ((start, record["n_terminal_addition"]),
                                           (end, record["c_terminal_addition"])) if addition)
    antigen = antigen_from_source(result, source_id, source_scope=source_scope)
    return ConstructSequence(
        name=record["id"], sequence=record["sequence"], modality=record["modality"],
        evidence=evidence, native_antigen=antigen, native_start=start, native_end=end,
        sequence_edits=edits)
