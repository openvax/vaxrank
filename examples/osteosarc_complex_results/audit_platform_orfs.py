"""Build a platform and assembly audit from the pinned Isovar osteosarc corpus.

This maintainer command consumes the checksummed full-product audit and the
bounded original-read BAM subsets in an Isovar checkout. It does not access
the network or copy large alignment files into the Vaxrank repository.
"""

import argparse
from collections import defaultdict
import gzip
from hashlib import sha256
import json
import logging
from pathlib import Path
import subprocess
import sys


PLATFORMS = ("ILMN", "ONT", "PacBio")
ATTRIBUTION_NOTES = {
    "NME1-chr17-51154443": (
        "The non-exact rows end one residue earlier than the isolated-edit "
        "window; other products recover the exact local sequence."),
    "NR2F2-chr15-96332299": (
        "The RNA window contains an unexplained three-nucleotide transcript "
        "deletion; no matched DNA variant has yet been linked to it."),
    "NTF3-chr12-5494381": (
        "RNA establishes the adjacent base change and compound AG>GT allele, "
        "changing Varcode's isolated Lys-to-Arg call to serine; germline "
        "versus somatic origin remains unresolved."),
    "ROBO2-chr3-77607853": (
        "One ILMN product contains an additional amino-acid substitution; "
        "other products recover the isolated-edit sequence."),
    "ZNF674-chrX-46501045": (
        "One ILMN product contains an additional amino-acid substitution; "
        "other products recover the isolated-edit sequence."),
    "ZNF764-chr16-30557995": (
        "Two ILMN products contain a divergent downstream protein context; "
        "other products recover the isolated-edit sequence."),
}


def digest(path):
    """Return the SHA-256 digest of a local input artifact."""
    checksum = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def platform_for_source(source):
    """Classify a recorded source without treating products as replicates."""
    text = " ".join(
        str(source.get(key, "")) for key in ("name", "key", "url"))
    text = text.lower()
    if "pacbio" in text:
        return "PacBio"
    if "ont" in text or "nanopore" in text:
        return "ONT"
    return "ILMN"


def _json_edit(edit):
    variant = edit.source_variant
    return {
        "transcript_id": edit.transcript_id,
        "cdna_start": edit.cdna_start,
        "cdna_end": edit.cdna_end,
        "alt_bases": edit.alt_bases,
        "source_variant": None if variant is None else {
            "contig": variant.contig,
            "start": variant.start,
            "ref": variant.ref,
            "alt": variant.alt,
        },
    }


def summarize_full_matrix(matrix):
    """Summarize the assembly-enabled 44-variant full-product audit."""
    sources = {source["source_id"]: source for source in matrix["sources"]}
    summaries = []
    discrepancies = defaultdict(lambda: {
        "platforms": set(), "rows": [], "actual_expected": set()})
    for platform in PLATFORMS:
        rows = [
            row for row in matrix["rows"]
            if platform_for_source(sources[row["source_id"]]) == platform]
        completed = [
            row for row in rows
            if row.get("defaults", {}).get("status") == "ok"]
        with_alt = {
            row["variant_id"] for row in completed
            if (row["defaults"].get("counts") or {}).get(
                "read_objects", {}).get("alt", 0) > 0}
        with_protein = {
            row["variant_id"] for row in completed
            if row["defaults"].get("proteins")}
        exact = {
            row["variant_id"] for row in completed
            if row["defaults"].get("outcome") == "expected_top_protein"}
        nonexact = {
            row["variant_id"] for row in completed
            if row["defaults"].get("outcome") == "different_top_protein"}
        summaries.append({
            "platform": platform,
            "catalogued_products": len({row["source_id"] for row in rows}),
            "products_with_completed_rows": len({
                row["source_id"] for row in completed}),
            "completed_source_variant_rows": len(completed),
            "variants_assessed": len({row["variant_id"] for row in completed}),
            "variants_with_alt_evidence": len(with_alt),
            "variants_with_validated_protein_window": len(with_protein),
            "variants_exact_in_at_least_one_product": len(exact),
            "variants_nonexact_in_at_least_one_product": len(nonexact),
            "variants_with_both_exact_and_nonexact_products": len(exact & nonexact),
        })
        for row in completed:
            if row["defaults"].get("outcome") != "different_top_protein":
                continue
            protein = row["defaults"]["proteins"][0]
            expected = {
                check["expected_amino_acids"]
                for check in protein["checks"]
                if check.get("expected_amino_acids") is not None}
            record = discrepancies[row["variant_id"]]
            record["platforms"].add(platform)
            record["rows"].append({
                "source_id": row["source_id"],
                "source_name": sources[row["source_id"]]["name"],
                "actual": protein["amino_acids"],
                "expected": sorted(expected),
                "supporting_fragments": protein["supporting_template_names"],
            })
            for expected_sequence in expected:
                record["actual_expected"].add(
                    (protein["amino_acids"], expected_sequence))
    discrepancy_records = []
    for variant_id, record in sorted(discrepancies.items()):
        context_only = all(
            expected.startswith(actual) or actual.startswith(expected)
            for actual, expected in record.pop("actual_expected"))
        discrepancy_records.append({
            "variant_id": variant_id,
            "platforms": sorted(record.pop("platforms")),
            "classification": (
                "incomplete_local_context" if context_only
                else "rna_protein_differs_from_isolated_edit"),
            "attribution_note": ATTRIBUTION_NOTES[variant_id],
            **record,
        })
    return summaries, discrepancy_records


def _load_isovar_helpers(isovar_repo):
    import isovar
    from isovar import run_isovar
    from isovar.read_collector import ReadCollector
    from isovar.protein_sequence_creator import ProteinSequenceCreator
    import pysam
    from varcode import Variant

    # Load the independent cohort oracle only after fixing the runtime Isovar
    # module to the installed release used by Vaxrank.
    sys.path.insert(0, str(isovar_repo))
    from tests.data.osteosarc.expansion.references import (
        apply_variant, load_reference, reference_genome)
    from tests.data.osteosarc.expansion.runner import protein_check

    return {
        "ProteinSequenceCreator": ProteinSequenceCreator,
        "ReadCollector": ReadCollector,
        "Variant": Variant,
        "apply_variant": apply_variant,
        "load_reference": load_reference,
        "protein_check": protein_check,
        "pysam": pysam,
        "reference_genome": reference_genome,
        "run_isovar": run_isovar,
        "runtime_module_sha256": digest(isovar.__file__),
        "runtime_version": isovar.__version__,
    }


def run_paired_corpus(corpus, sources, cache_directory, helpers):
    """Rerun identical original reads with Isovar assembly enabled/disabled."""
    manifests = {}
    genomes = {}
    models = {}
    for reference_name in sorted({case["reference"] for case in corpus["cases"]}):
        reference = corpus["root"] / "references" / reference_name
        manifests[reference_name], models[reference_name] = helpers[
            "load_reference"](reference)
        genomes[reference_name] = helpers["reference_genome"](
            reference, cache_directory / reference_name)

    cases = []
    for case in corpus["cases"]:
        record = case["variant"]
        reference_name = case["reference"]
        contig = record["chrom"].removeprefix("chr")
        variant = helpers["Variant"](
            "MT" if contig == "M" else contig,
            record["pos"], record["ref"], record["alt"],
            ensembl=genomes[reference_name])
        transcript_ids = manifests[reference_name]["variant_transcripts"][
            record["variant_id"]]
        expectations = {
            transcript_id: helpers["apply_variant"](
                record, models[reference_name][transcript_id])
            for transcript_id in transcript_ids}
        varcode_effects = []
        for effect in variant.effects():
            if effect.transcript_id not in transcript_ids:
                continue
            varcode_effects.append({
                "effect_class": type(effect).__name__,
                "short_description": effect.short_description,
                "transcript_id": effect.transcript_id,
                "transcript_name": effect.transcript_name,
                "has_exact_mutant_protein_sequence": bool(
                    getattr(effect, "mutant_protein_sequence", None)),
            })
        modes = {}
        for assembly in (True, False):
            creator = helpers["ProteinSequenceCreator"](
                variant_sequence_assembly=assembly)
            collector = helpers["ReadCollector"](
                merge_overlapping_fragments=True)
            with helpers["pysam"].AlignmentFile(
                    corpus["root"] / case["bam"]) as bam:
                result, = helpers["run_isovar"](
                    [variant], bam,
                    transcript_id_whitelist=set(transcript_ids),
                    read_collector=collector,
                    protein_sequence_creator=creator)
            top = (
                result.sorted_protein_sequences[0]
                if result.sorted_protein_sequences else None)
            if top is None:
                modes["assembly_on" if assembly else "assembly_off"] = {
                    "protein_status": "not_established",
                    "amino_acids": None,
                    "matches_isolated_edit": None,
                    "unexplained_transcript_edits": [],
                }
                continue
            checked = helpers["protein_check"](
                top, expectations, creator.protein_sequence_length,
                creator.protein_context_peptide_length)
            modes["assembly_on" if assembly else "assembly_off"] = {
                "protein_status": "validated_local_window",
                "amino_acids": top.amino_acids,
                "mutation_interval": [
                    top.mutation_start_idx, top.mutation_end_idx],
                "supporting_fragments": top.num_supporting_fragments,
                "matches_isolated_edit": checked["matches_expected"],
                "expected_amino_acids": sorted({
                    check["expected_amino_acids"]
                    for check in checked["checks"]}),
                "known_germline_transcript_edits": [
                    _json_edit(edit)
                    for edit in top.known_germline_transcript_edits],
                "known_somatic_transcript_edits": [
                    _json_edit(edit)
                    for edit in top.known_somatic_transcript_edits],
                "unexplained_transcript_edits": [
                    _json_edit(edit)
                    for edit in top.unexplained_transcript_edits],
            }
        cases.append({
            "case_id": case["case_id"],
            "variant_id": record["variant_id"],
            "gene": record["gene"],
            "source_id": case["source_id"],
            "platform": platform_for_source(sources[case["source_id"]]),
            "source_url": case["source_url"],
            "varcode_effects": varcode_effects,
            "modes": modes,
        })
    return cases


def summarize_paired_cases(cases):
    summaries = []
    for platform in PLATFORMS:
        platform_cases = [case for case in cases if case["platform"] == platform]
        for mode in ("assembly_on", "assembly_off"):
            established = [
                case for case in platform_cases
                if case["modes"][mode]["protein_status"] ==
                "validated_local_window"]
            exact = [
                case for case in established
                if case["modes"][mode]["matches_isolated_edit"]]
            summaries.append({
                "platform": platform,
                "mode": mode,
                "cases": len(platform_cases),
                "variants": len({case["variant_id"] for case in platform_cases}),
                "validated_local_protein_windows": len(established),
                "exact_isolated_edit_windows": len(exact),
                "nonexact_windows": len(established) - len(exact),
                "no_protein_window": len(platform_cases) - len(established),
            })
    return summaries


def summarize_assembly_effect(cases):
    """Compare paired protein availability and sequence for each platform."""
    summaries = []
    for platform in PLATFORMS:
        counts = {
            "same_sequence": 0,
            "different_sequence": 0,
            "assembly_only": 0,
            "no_assembly_only": 0,
            "neither": 0,
        }
        changed = []
        for case in (case for case in cases if case["platform"] == platform):
            assembled = case["modes"]["assembly_on"]["amino_acids"]
            unassembled = case["modes"]["assembly_off"]["amino_acids"]
            if assembled and unassembled:
                category = (
                    "same_sequence" if assembled == unassembled
                    else "different_sequence")
            elif assembled:
                category = "assembly_only"
            elif unassembled:
                category = "no_assembly_only"
            else:
                category = "neither"
            counts[category] += 1
            if category == "different_sequence":
                changed.append({
                    "variant_id": case["variant_id"],
                    "assembly_on_length": len(assembled),
                    "assembly_off_length": len(unassembled),
                    "assembly_on_sequence": assembled,
                    "assembly_off_sequence": unassembled,
                })
        summaries.append({"platform": platform, **counts,
                          "changed_cases": changed})
    return summaries


def summarize_varcode(cases):
    """Count unique loci where Varcode supplies an isolated-effect protein."""
    by_variant = defaultdict(list)
    for case in cases:
        by_variant[case["variant_id"]].extend(case["varcode_effects"])
    with_protein = {
        variant_id for variant_id, effects in by_variant.items()
        if any(effect["has_exact_mutant_protein_sequence"] for effect in effects)}
    return {
        "variants": len(by_variant),
        "variants_with_exact_isolated_effect_protein": len(with_protein),
        "variants_without_exact_isolated_effect_protein": sorted(
            set(by_variant) - with_protein),
        "interpretation": (
            "Varcode's exact sequence is the reference transcript plus the "
            "nominated edit; it is not a sample-specific RNA haplotype."),
    }


def write_markdown(payload, path):
    full = payload["full_matrix"]["platform_summary"]
    paired = payload["paired_corpus"]["platform_mode_summary"]
    lines = [
        "# Osteosarc RNA platform and ORF audit", "",
        "A protein result below is a validated local RNA-derived coding window, "
        "not a reconstructed full-length ORF. Products are not biological "
        "replicates, and counts are not platform-sensitivity estimates.", "",
        "## Full 44-variant audit (assembly enabled)", "",
        "| Platform | Products | Completed rows | ALT loci | Protein loci | "
        "Exact in ≥1 product | Non-exact in ≥1 product |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in full:
        lines.append(
            "| {platform} | {catalogued_products} | "
            "{completed_source_variant_rows} | {variants_with_alt_evidence} | "
            "{variants_with_validated_protein_window} | "
            "{variants_exact_in_at_least_one_product} | "
            "{variants_nonexact_in_at_least_one_product} |".format(**row))
    lines.extend([
        "", "Varcode produces an exact reference-plus-isolated-edit protein "
        "for all 44 loci. That sequence is an annotation result, not an "
        "observed sample-specific RNA haplotype.",
        "", "## Paired original-read corpus", "",
        "The same bounded BAM subsets were rerun with Isovar assembly on and "
        "off. This 49-case corpus is enriched for informative loci and is not "
        "an unbiased cohort sample.", "",
        "| Platform | Assembly | Cases | Protein windows | Exact | Non-exact | "
        "No window |", "| --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ])
    for row in paired:
        lines.append(
            "| {platform} | {mode} | {cases} | "
            "{validated_local_protein_windows} | "
            "{exact_isolated_edit_windows} | {nonexact_windows} | "
            "{no_protein_window} |".format(**row))
    lines.extend([
        "", "Assembly changed the returned local protein length in five of "
        "32 ILMN cases, but did not create or remove a protein window in this "
        "selected corpus. No paired ONT or PacBio case changed sequence; the "
        "PacBio denominator is only one case.",
        "", "## Non-exact RNA protein windows", "",
        "| Locus | Platform | Product rows | Interpretation |",
        "| --- | --- | ---: | --- |",
    ])
    for record in payload["full_matrix"]["nonexact_records"]:
        lines.append("| %s | %s | %d | %s |" % (
            record["variant_id"].split("-chr", 1)[0],
            " + ".join(record["platforms"]),
            len(record["rows"]),
            record["attribution_note"]))
    lines.extend([
        "", "## Attribution", "",
        "The current public Isovar path identifies the focal somatic edit and "
        "retains additional transcript edits as unexplained. It can phase "
        "multiple supplied somatic variants by shared fragment names, but this "
        "44-variant set contains no nearby nominated pair. A matched germline "
        "variant input is not wired into `run_isovar`; consequently no observed "
        "mismatch can yet be promoted to known germline by this audit.", "",
        "See the adjacent JSON for every source row, sequence, transcript edit, "
        "checksum, and denominator.", "",
    ])
    path.write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--isovar-repo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--markdown-output", type=Path)
    parser.add_argument("--cache-directory", required=True, type=Path)
    args = parser.parse_args()

    expansion = (
        args.isovar_repo / "tests/data/osteosarc/expansion")
    matrix_path = expansion / "audit/matrix.json.gz"
    corpus_root = expansion / "corpus"
    corpus_manifest_path = corpus_root / "manifest.json"
    with gzip.open(matrix_path, "rt") as handle:
        matrix = json.load(handle)
    corpus_manifest = json.loads(corpus_manifest_path.read_text())
    corpus_manifest["root"] = corpus_root
    sources = {source["source_id"]: source for source in matrix["sources"]}

    logging.disable(logging.CRITICAL)
    helpers = _load_isovar_helpers(args.isovar_repo)
    cases = run_paired_corpus(
        corpus_manifest, sources, args.cache_directory, helpers)
    full_summary, discrepancies = summarize_full_matrix(matrix)
    payload = {
        "schema_version": 1,
        "scope": {
            "full_matrix": (
                "44 vaccine variants across 164 catalogued RNA products; "
                "assembly enabled; products are not independent samples"),
            "paired_corpus": (
                "49 selected original-read source/variant cases; same BAM "
                "subsets rerun with assembly on and off; selection enriched"),
            "protein_definition": (
                "independently validated local RNA-derived coding window; "
                "not a full-length ORF"),
        },
        "provenance": {
            "isovar_repository": "https://github.com/openvax/isovar",
            "isovar_evidence_commit": subprocess.check_output(
                ["git", "-C", str(args.isovar_repo), "rev-parse", "HEAD"],
                text=True).strip(),
            "isovar_runtime_module_sha256": helpers[
                "runtime_module_sha256"],
            "isovar_runtime_version": helpers["runtime_version"],
            "matrix_sha256": digest(matrix_path),
            "corpus_manifest_sha256": digest(corpus_manifest_path),
            "matrix_schema_version": matrix["schema_version"],
        },
        "full_matrix": {
            "variant_count": len(matrix["variants"]),
            "catalogued_product_count": len(matrix["sources"]),
            "source_variant_row_count": len(matrix["rows"]),
            "platform_summary": full_summary,
            "nonexact_records": discrepancies,
        },
        "paired_corpus": {
            "case_count": len(cases),
            "platform_mode_summary": summarize_paired_cases(cases),
            "assembly_effect_summary": summarize_assembly_effect(cases),
            "cases": cases,
        },
        "varcode": summarize_varcode(cases),
        "attribution": {
            "known_germline_matches": 0,
            "reason_known_germline_is_zero": (
                "run_isovar has no matched-germline variant input; protein "
                "objects therefore retain non-focal differences as unexplained"),
            "somatic_phasing_capability": (
                "shared-fragment phase groups are computed for multiple "
                "somatic variants supplied in one run"),
            "somatic_phasing_in_44_variant_audit": (
                "not evaluable: no nominated pair is close enough to share a "
                "local RNA fragment"),
            "source_linked_compound_examples": [
                "NTF3 AG>GT changes the isolated A>G Lys-to-Arg prediction to serine",
                "MAP2 and CD109 require explicit compound-haplotype inputs in the extended Vaxrank examples",
            ],
        },
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    if args.markdown_output:
        write_markdown(payload, args.markdown_output)


if __name__ == "__main__":
    main()
