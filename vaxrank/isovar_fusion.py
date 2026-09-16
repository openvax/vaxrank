"""Adapt trusted Isovar fusion reconstructions without choosing an isoform."""

from dataclasses import dataclass, replace
from hashlib import sha256
import json

from isovar.genetic_code import standard_genetic_code
from serializable import DataclassSerializable

from .vaccine_antigen import (
    ANTIGEN_KIND_FUSION, ATTESTATION_HELD_OUT, AminoAcidInterval,
    TargetableMask, TumorSpecificityAttestation, VaccineAntigen,
)


@dataclass(frozen=True)
class IsovarFusionAntigens(DataclassSerializable):
    """All coding hypotheses plus the original, possibly unresolved result."""

    reconstruction: dict
    antigens: tuple[VaccineAntigen, ...]

    @property
    def admitted_antigens(self) -> tuple[VaccineAntigen, ...]:
        return tuple(a for a in self.antigens if a.tumor_specificity.admits_construct)


def fusion_antigens_from_isovar(result, *, tumor_specificity=None, gene_name="", species=""):
    """Adapt a result from ``isovar.reconstruct_fusion`` (schema 1).

    Parameters
    ----------
    result : dict
        Trusted Isovar output, not a caller's unvalidated fusion declaration.
        The complete result is JSON-normalized into the return value and each antigen's
        source metadata, so provenance survives downstream serialization.
    tumor_specificity : TumorSpecificityAttestation, optional
        Independent evidence for tumor specificity. RNA reconstruction alone
        does not establish this. Ambiguous or partial-CDS hypotheses remain
        held out regardless of this attestation.
    gene_name, species : str
        Optional source labels; no identity is inferred from an event name.

    Returns
    -------
    IsovarFusionAntigens
        Every coding hypothesis, without ranking, and the original outcome.
        Unresolved results have no antigens, but keep their reasons/evidence.
    """
    payload = json.dumps(result, sort_keys=True)
    result = json.loads(payload)
    if result.get("schema_version") != 1:
        raise ValueError("Unsupported Isovar fusion result schema")
    status = result["status"]
    if status not in {"translated", "ambiguous", "unresolved_frame", "insufficient_support"}:
        raise ValueError("Unknown Isovar fusion reconstruction status")
    sequence = result["cdna_sequence"]
    if sha256(sequence.encode()).hexdigest() != result["sequence_sha256"]:
        raise ValueError("Isovar fusion sequence checksum mismatch")
    translations = result["translations"]
    if status == "translated" and (len(translations) != 1 or result["reasons"]):
        raise ValueError("Translated Isovar result must have one unambiguous coding hypothesis")
    if tumor_specificity is None:
        tumor_specificity = TumorSpecificityAttestation(
            status=ATTESTATION_HELD_OUT, evidence_kind="missing_tumor_specificity",
            evidence_source="Isovar RNA reconstruction", rationale_code="rna_is_not_tumor_specificity",
            requires_review=True,
        )
    antigens = []
    for index, protein in enumerate(translations):
        start = protein["translation_start"]
        complete = protein["complete_5prime"]
        if (type(start) is not int or not 0 <= start < len(sequence)
                or type(complete) is not bool or protein["cds_start"] != (start if complete else None)):
            raise ValueError("Inconsistent Isovar fusion CDS start")
        expected = standard_genetic_code.translate(sequence[start:], first_codon_is_start=complete)
        if tuple(expected) != (protein["amino_acids"], protein["ends_with_stop_codon"]):
            raise ValueError("Isovar fusion protein disagrees with its RNA sequence/frame")
        boundaries = protein["junction_in_translated_cds"]
        if (len(boundaries) != 2 or boundaries != [b - start for b in result["junction_interval"]]
                or not 0 < boundaries[0] <= boundaries[1] < 3 * len(protein["amino_acids"])):
            raise ValueError("Inconsistent Isovar fusion junction coordinates")
        # Separate boundaries avoid admitting a peptide wholly inside an insert.
        # A split codon itself contains bases from both sides of that boundary.
        intervals = sorted({(b // 3, (b + 2) // 3) for b in boundaries})
        attestation = tumor_specificity
        if status != "translated" or not complete:
            attestation = replace(attestation, status=ATTESTATION_HELD_OUT,
                rationale_code="isovar_" + (status if status != "translated" else "partial_cds"),
                requires_review=True, override_reason="")
        antigens.append(VaccineAntigen(
            kind=ANTIGEN_KIND_FUSION, amino_acids=protein["amino_acids"],
            targetable_mask=TargetableMask(tuple(AminoAcidInterval(a, b) for a, b in intervals)),
            tumor_specificity=attestation, gene_name=gene_name, species=species,
            transcript_ids=tuple(sorted(set(protein["donor_transcript_ids"])
                | set(result["compatible_transcripts"]["acceptor"]))),
            source_identifier=result["event_id"] + ":hypothesis:" + str(index + 1),
            source_metadata=(("isovar_fusion_result", payload),
                ("isovar_translation_index", str(index)),
                ("requested_tumor_specificity", tumor_specificity.to_json())),
        ))
    return IsovarFusionAntigens(reconstruction=result, antigens=tuple(antigens))
