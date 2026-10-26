"""Convert retained Vaxrank results to ORFdata without rerunning analysis."""

import hashlib
import json
import math
from collections.abc import Mapping
from dataclasses import replace

from isovar import from_protein_hypotheses as isovar_from_protein_hypotheses
from isovar import to_proteome as isovar_to_proteome
from orfdata import (
    ProducerInfo,
    ProducerLink,
    ProducerRecord,
    ProteinSequence,
    Proteome,
    ReferenceContext,
)
from serializable.helpers import to_serializable_repr

from .mutant_protein_fragment import (
    SEQUENCE_SOURCE_ISOVAR_ASSEMBLY,
    SEQUENCE_SOURCE_VARCODE_TRANSLATION,
    MutantProteinFragment,
)
from .vaxrank_results import VaxrankResults
from .version import __version__


SCHEMA = "vaxrank.results.v1"


def _json_value(value):
    """Return strict JSON data, retaining unknown non-finite values as null."""
    value = to_serializable_repr(value)

    def normalize(item):
        if isinstance(item, float) and not math.isfinite(item):
            return None
        if isinstance(item, dict):
            return {key: normalize(child) for key, child in item.items()}
        if isinstance(item, (list, tuple)):
            return [normalize(child) for child in item]
        return item

    value = normalize(value)
    json.dumps(value, allow_nan=False)
    return value


def _with_export_path(proteome):
    """Record that Vaxrank exported Isovar data without changing its identity."""
    proteins = tuple(
        replace(
            protein,
            producer=replace(
                protein.producer,
                exporter=None,
                exporters=tuple(sorted(set(protein.producer.exporters + ("vaxrank",)))),
            ),
        )
        if protein.producer is not None else protein
        for protein in proteome.proteins
    )
    records = tuple(
        replace(
            record,
            exporter=None,
            exporters=tuple(sorted(set(record.exporters + ("vaxrank",)))),
        )
        for record in proteome.records
    )
    return Proteome(
        proteome.reference,
        proteins,
        proteome.layers,
        proteome.excluded,
        proteome.genome,
        reference_policy=proteome.reference_policy,
        records=records,
        links=proteome.links,
    )


def _variant_fields(variant):
    fields = {}
    for name in ("reference_name", "contig", "start", "ref", "alt"):
        value = getattr(variant, name, None)
        if value is not None:
            fields[name] = value
    description = getattr(variant, "short_description", None)
    if description:
        fields["description"] = description
    return fields


def _same_variant(recorded, variant):
    if not isinstance(recorded, Mapping):
        return False
    current = _variant_fields(variant)
    compared = False
    for name in ("reference_name", "contig", "start", "ref", "alt"):
        if name in recorded and name in current:
            compared = True
            if recorded[name] != current[name]:
                return False
    return compared


def _isovar_protein(proteome, fragment):
    """Find the Isovar protein selected to make one retained Vaxrank fragment."""
    records = {record.key: record for record in proteome.records}
    proteins = proteome.by_accession()
    candidates = []
    event_links = {}
    protein_links = {}
    for link in proteome.links:
        if link.target_kind == "record" and link.relation == "has_hypothesis":
            event_links.setdefault(link.source, []).append(link.target)
        elif link.target_kind == "protein" and link.relation == "has_sequence":
            protein_links.setdefault(link.source, []).append(link.target)
    for event in proteome.records:
        if (event.namespace != "isovar" or event.kind != "event"
                or not _same_variant(event.payload.get("variant"), fragment.variant)):
            continue
        for hypothesis_key in event_links.get(event.key, ()):
            hypothesis = records[hypothesis_key]
            if hypothesis.payload.get("amino_acids") != fragment.amino_acids:
                continue
            for accession in protein_links.get(hypothesis_key, ()):
                candidates.append((hypothesis.payload.get("isovar_rank"), proteins[accession]))
    best = [protein for rank, protein in candidates if rank == 1]
    if len(best) == 1:
        return best[0]
    if len(candidates) == 1:
        return candidates[0][1]
    if candidates:
        raise ValueError("Retained Isovar fragment matches more than one source hypothesis")
    return None


def _source_reference(fragment, reference):
    declared = getattr(fragment.variant, "reference_name", None)
    if not declared:
        return None
    if declared != reference.genome_build:
        raise ValueError(
            f"Fragment reference {declared!r} conflicts with {reference.genome_build!r}"
        )
    return reference


def _transcript_ids(fragment):
    return tuple(sorted({
        transcript.id for transcript in fragment.supporting_reference_transcripts
        if getattr(transcript, "id", None)
    }))


def _gene_ids(fragment):
    identifiers = set()
    for transcript in fragment.supporting_reference_transcripts:
        gene_id = getattr(transcript, "gene_id", None)
        if gene_id is None:
            gene_id = getattr(getattr(transcript, "gene", None), "id", None)
        if gene_id:
            identifiers.add(gene_id)
    return tuple(sorted(identifiers))


def _fragment_source(fragment):
    source = fragment.sequence_source
    if source == SEQUENCE_SOURCE_ISOVAR_ASSEMBLY:
        return "reference_cds", "rna_assembly", "isovar"
    if source == SEQUENCE_SOURCE_VARCODE_TRANSLATION or fragment.dna_effect_selection is not None:
        return "reference_cds", "dna_prediction", "varcode"
    if source == "lens_pep_context":
        return "unknown", "supplied_fragment", "lens"
    if source == "pvacseq_epitope":
        return "unknown", "supplied_fragment", "pvacseq"
    return "unknown", "supplied_fragment", "unknown"


def _fragment_native_id(fragment, tool):
    variant = _variant_fields(fragment.variant)
    contig = str(variant.get("contig", "unknown"))
    start = str(variant.get("start", "unknown"))
    ref, alt = str(variant.get("ref", "?")), str(variant.get("alt", "?"))
    change = f"{ref}>{alt}" if len(ref) + len(alt) <= 20 else f"{len(ref)}bp>{len(alt)}bp"
    window_start = fragment.protein_start_offset
    window = (
        f"aa{window_start}-{window_start + len(fragment.amino_acids)}"
        if window_start is not None else f"aa?-len{len(fragment.amino_acids)}"
    )
    identity = {
        "variant": variant,
        "tool": tool,
        "sequence_source": fragment.sequence_source,
        "window": [window_start, len(fragment.amino_acids)],
        "mutation_interval": [
            fragment.mutant_amino_acid_start_offset,
            fragment.mutant_amino_acid_end_offset,
        ],
        "transcript_ids": _transcript_ids(fragment),
        "dna_effect_selection": fragment.dna_effect_selection,
    }
    if (fragment.dna_effect_selection is None
            and fragment.sequence_source != SEQUENCE_SOURCE_ISOVAR_ASSEMBLY):
        identity["sequence_sha256"] = hashlib.sha256(
            fragment.amino_acids.encode("ascii")
        ).hexdigest()
    digest = hashlib.sha256(
        json.dumps(identity, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()[:12]
    return f"{contig}:{start}:{change}:{window}:{digest}"


def _fragment_protein(fragment, sample_id, source, reference):
    origin, method, tool = _fragment_source(fragment)
    # Native Isovar/Varcode translations are tied to the caller-asserted
    # annotation. Imported fragments do not become reference-qualified merely
    # because Vaxrank used that reference while processing them.
    exact_reference = (
        _source_reference(fragment, reference)
        if tool in ("isovar", "varcode") else None
    )
    producer = ProducerInfo(
        origin=origin,
        method=method,
        tool=tool,
        native_id=_fragment_native_id(fragment, tool),
        tool_version=fragment.sequence_source_version or None,
        sample_id=sample_id,
        run_id=source,
        gene_ids=_gene_ids(fragment),
        gene_names=(fragment.gene_name,) if fragment.gene_name else (),
        transcript_ids=_transcript_ids(fragment),
        reference=exact_reference,
        start_complete=False if (fragment.protein_start_offset or 0) > 0 else None,
        stop_complete=None,
        exporters=("vaxrank",),
    )
    return ProteinSequence.from_producer(
        fragment.amino_acids,
        producer=producer,
        taxon_id=None if exact_reference is None else exact_reference.taxon_id,
    )


def _fragment_record(fragment, sample_id, source):
    _, _, tool = _fragment_source(fragment)
    start = fragment.protein_start_offset
    return ProducerRecord(
        namespace="vaxrank",
        schema_version="v1",
        kind="fragment",
        native_id=_fragment_native_id(fragment, tool),
        sample_id=sample_id,
        run_id=source,
        exporters=("vaxrank",),
        payload={
            "fragment": _json_value(fragment),
            "length_aa": len(fragment.amino_acids),
            "mutation_interval": [
                fragment.mutant_amino_acid_start_offset,
                fragment.mutant_amino_acid_end_offset,
            ],
            "protein_interval": (
                None if start is None else [start, start + len(fragment.amino_acids)]
            ),
            "sequence_source": fragment.sequence_source or None,
        },
    )


def to_proteome(results, *, sample_id, source, reference):
    """Convert retained Vaxrank results to an offline ORFdata proteome.

    This is an export operation: it does not predict MHC binding, rerank
    peptides, reconstruct full proteins from fragments, or change host layers.
    Isovar hypotheses keep their native source identity and gain ``vaxrank``
    only as an additional export path.
    """
    if not isinstance(results, VaxrankResults):
        raise TypeError("results must be VaxrankResults")
    if any(not isinstance(value, str) or not value.strip() for value in (sample_id, source)):
        raise ValueError("sample_id and source must be nonempty strings")
    if not isinstance(reference, ReferenceContext) or not reference.is_genomic:
        raise ValueError("reference must be a genomic ORFdata ReferenceContext")
    if any(not isinstance(fragment, MutantProteinFragment)
           for fragment in results.retained_protein_fragments):
        raise TypeError("retained_protein_fragments must contain MutantProteinFragment values")

    saved_isovar = results.isovar_protein_hypotheses
    if results.isovar_results and saved_isovar is not None:
        raise ValueError("Results cannot contain both live and saved Isovar hypotheses")
    if saved_isovar is not None:
        if not isinstance(saved_isovar, Mapping):
            raise ValueError("Saved Isovar hypotheses must be a versioned JSON object")
        if (saved_isovar.get("sample_id"), saved_isovar.get("source")) != (sample_id, source):
            raise ValueError("Saved Isovar sample_id and source disagree with the export request")
        proteome = _with_export_path(isovar_from_protein_hypotheses(
            saved_isovar, reference=reference))
    elif results.isovar_results:
        proteome = _with_export_path(isovar_to_proteome(
            results.isovar_results,
            sample_id=sample_id,
            source=source,
            reference=reference,
        ))
    else:
        proteome = Proteome(reference, ())

    export_record = ProducerRecord(
        namespace="vaxrank",
        schema_version="v1",
        kind="export",
        native_id="results",
        sample_id=sample_id,
        run_id=source,
        complete=True,
        exporters=("vaxrank",),
        payload={
            "schema": SCHEMA,
            "vaxrank_version": __version__,
            "retained_fragment_count": len(results.retained_protein_fragments),
            "ranked_variant_count": len(results.ranked_vaccine_peptides),
            "selection": {
                "variant_to_vaccine_peptides": _json_value(
                    results.variant_to_vaccine_peptides_dict
                ),
                "ranked_vaccine_peptides": _json_value(results.ranked_vaccine_peptides),
            },
            "interpretation": (
                "Recorded Vaxrank fragments and selection state; export does not rerun prediction"
            ),
        },
    )
    records = [export_record]
    links = [
        ProducerLink(export_record.key, record.key, "record", "uses_upstream_export")
        for record in proteome.records
        if record.namespace == "isovar" and record.kind == "export"
    ]
    proteins = []
    for fragment in results.retained_protein_fragments:
        fragment_record = _fragment_record(fragment, sample_id, source)
        records.append(fragment_record)
        links.append(ProducerLink(
            export_record.key, fragment_record.key, "record", "contains_fragment"
        ))
        protein = None
        if fragment.sequence_source == SEQUENCE_SOURCE_ISOVAR_ASSEMBLY:
            protein = _isovar_protein(proteome, fragment)
            if protein is None and (results.isovar_results or saved_isovar is not None):
                raise ValueError(
                    "Retained Isovar fragment does not match an exported Isovar hypothesis"
                )
        if protein is None:
            protein = _fragment_protein(fragment, sample_id, source, reference)
            proteins.append(protein)
        links.append(ProducerLink(
            fragment_record.key,
            protein.accession,
            "protein",
            "retains_fragment" if protein in proteome.proteins else "has_sequence",
            payload={
                "protein_start_offset": fragment.protein_start_offset,
                "mutation_interval": [
                    fragment.mutant_amino_acid_start_offset,
                    fragment.mutant_amino_acid_end_offset,
                ],
            },
        ))

    export_record = replace(
        export_record,
        payload={
            **export_record.to_dict()["payload"],
            "fragment_record_ids": [
                record.key for record in records if record.kind == "fragment"
            ],
        },
    )
    records[0] = export_record

    if proteins:
        compatibility = (
            "allow_unknown"
            if any(protein.taxon_id is None or protein.producer.reference is None
                   for protein in proteins)
            else None
        )
        proteome = proteome.with_proteins(proteins, compatibility=compatibility)
    return proteome.with_records(records, links=links)
