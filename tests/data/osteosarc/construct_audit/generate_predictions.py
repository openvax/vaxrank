"""Maintainer-only real-model audit; no test invokes models or downloads.

Run from the repository root with the published dependencies and local models:
python -m tests.data.osteosarc.construct_audit.generate_predictions --regional-cache ISOVAR_CACHE --output DIR
Requires an already installed full Ensembl human release 114, NetMHCpan 4.2c
and Pepsickle 0.1.3. Reproduction uses the recorded settings, not latest models.
"""

import argparse
from dataclasses import replace
import gzip
from hashlib import sha256
import importlib
import json
import logging
from pathlib import Path
import platform
import shutil
import tempfile

from mhctools import NetMHCpan42
from mhctools.peptidases import get_cleavage_model
from pyensembl import EnsemblRelease

from vaxrank import MHCRequest, audit_sequence_contexts, final_sequence_context, native_sequence_context
from vaxrank.native_serialization import to_native_json
from vaxrank.reference_proteome import (
    ensembl_dataset_cache_identity, oncoref_cta_source_gene_ids, self_reference_matches,
)

from tests.osteosarc_construct_helpers import (
    DATA, DOCUMENTED, RNA, antigen_from_source, construct_from_result,
    verify_rna_inputs,
)
from tests.osteosarc_selection_helpers import DATA as SELECTION_DATA, load_selection_inputs, reconstruct_selection
from .import_rna import digest


LENGTHS = (8, 9, 10, 11)
KINDS = ("pMHC_affinity", "pMHC_presentation")
ANCHOR = "1b66c15da594a3ef"  # Earliest bulk T0 source; all four windows present without overrides.
CHEMISTRY_BASIS = (
    "Conditional bare-peptide model input with free N/C termini, not evidence of "
    "the historical manufactured chemistry. No acetylation/amidation is documented."
)


def canonical_audits(audits):
    """Canonicalize observation order only; preserve every native model value.

    Backends may emit allele batches in hash-dependent order. Ordering by exact
    prediction identity makes offline artifacts reproducible across hash seeds.
    """
    return tuple(replace(audit, ligands=tuple(sorted(
        (replace(ligand, predictions=tuple(sorted(ligand.predictions, key=lambda p: p.identity)))
         for ligand in audit.ligands), key=lambda ligand: (ligand.start, len(ligand.peptide), ligand.peptide))))
        for audit in audits)


def reconstruct_contexts(inputs, regional_cache):
    """Keep every source identity, its true gates, and exact window limitations."""
    results, contexts, summaries = {}, [], []
    for case in RNA["cases"]:
        sid = case["source_id"]
        bam = regional_cache / "alignments" / sid / "regions-GRCh38.bam"
        actual_digest = digest(bam)
        receipts = [r for r in case["source"]["acquisition_receipts"]
                    if r["status"] == "ok" and r["bam_sha256"] == actual_digest]
        if len(receipts) != 1:
            raise ValueError("Full regional BAM does not match its released receipt: " + sid)
        result, _ = reconstruct_selection(inputs, DOCUMENTED["variant_id"], bam_path=bam)
        if result.top_protein_sequence is None:
            raise ValueError("Pilot source no longer reconstructs a protein: " + sid)
        results[sid] = result
        context = native_sequence_context(antigen_from_source(
            result, sid, source_scope="full_original_region"), name="rna-" + sid)
        contexts.append(context)
        windows = []
        for record in DOCUMENTED["records"]:
            native = record["native_sequence"]
            offsets = [i for i in range(len(context.sequence) - len(native) + 1)
                       if context.sequence[i:i + len(native)] == native]
            windows.append(dict(construct_id=record["id"], native_occurrence_offsets=offsets,
                                status="exact_unique_native_window" if len(offsets) == 1
                                else "absent_from_selected_rna_context" if not offsets else "ambiguous_native_window"))
        summaries.append(dict(
            source_id=sid, context_id=context.source_id, sequence=context.sequence,
            input_scope="full_original_region", regional_bam_sha256=actual_digest,
            modality=case["modality"], timepoint=case["timepoint"], attribution=case["attribution"],
            n_alt_reads=int(result.num_alt_reads), n_alt_fragments=int(result.num_alt_fragments),
            passes_all_filters=bool(result.passes_all_filters),
            filter_values={k: bool(v) for k, v in result.filter_values.items()},
            windows=windows, historical_selection_status="not_assessed_by_construct_audit"))
    for record in DOCUMENTED["records"]:
        construct = construct_from_result(results[ANCHOR], record, ANCHOR, source_scope="full_original_region")
        chemistry = (dict(n_term="free", c_term="free", chemistry_basis=CHEMISTRY_BASIS)
                     if record["modality"] == "peptide" else {})
        contexts.append(final_sequence_context(construct, **chemistry))
    return tuple(contexts), summaries


def model_metadata():
    executable_path = shutil.which("netMHCpan")
    if not executable_path:
        raise ValueError("NetMHCpan executable is not installed")
    executable = Path(executable_path).resolve()
    root = executable.parent
    binary = root / (platform.system() + "_" + platform.machine()) / "bin/netMHCpan-4.2"
    if not (root / "data").is_dir() or not binary.is_file():
        raise ValueError("Expected netmhc-bundle 4.2 layout for complete model provenance")
    assets = [binary, *(p for p in (root / "data").rglob("*") if p.is_file())]
    return dict(executable_name=executable.name, executable_sha256=digest(executable),
                model_assets_sha256={str(p.relative_to(root)): digest(p) for p in sorted(assets)})


def main(output, regional_cache):
    logging.disable(logging.WARNING)
    verify_rna_inputs()
    if output.exists():
        raise ValueError("Choose a new output directory; existing real caches are never overwritten")
    packages = {p: getattr(importlib.import_module(p), "__version__", None)
                for p in ("isovar", "mhctools", "topiary", "varcode", "pyensembl", "vaxrank", "oncoref")}
    if packages["isovar"] != "1.8.1":
        raise ValueError("This reference run requires published Isovar 1.8.1")
    with tempfile.TemporaryDirectory(prefix="sid-audit-reference-") as tmp:
        inputs = load_selection_inputs(Path(tmp))
        contexts, rna_results = reconstruct_contexts(inputs, regional_cache)
        metadata = model_metadata()
        requests = tuple(MHCRequest(k, "netMHCpan", "4.2c", a, size)
                         for k in KINDS for a in DOCUMENTED["prediction_alleles"] for size in LENGTHS)
        print("Running real models on", len(contexts), "complete source-aware contexts", flush=True)
        audits = audit_sequence_contexts(
            contexts,
            mhc_predictor=NetMHCpan42(alleles=DOCUMENTED["prediction_alleles"],
                                       default_peptide_lengths=LENGTHS, process_limit=1),
            mhc_requests=requests, pepsickle=True,
            peptidase_predictors=[get_cleavage_model("cpn-basic")])
        audits = canonical_audits(audits)
        for audit in audits:
            if (audit.mhc_status != "predictions_returned" or audit.unrequested_predictions
                    or any(missing for _, _, missing in audit.coverage)):
                raise ValueError("Incomplete real MHC output: " + audit.context.name + " " + audit.error_message)
            profile = audit.profiles[0]
            if profile.status != "observations_returned" or len(profile.sites) != len(audit.context.sequence) - 1:
                raise ValueError("Incomplete Pepsickle output: " + audit.context.name + " " + profile.error_message)
            print(audit.context.name, len(audit.context.sequence), "aa", len(audit.ligands), "ligands", flush=True)
        # Full human reference: never reuse the tiny RNA reconstruction subset
        # to establish absence of a self match. index() resolves existing files.
        genome = EnsemblRelease(114)
        genome.index()
        print("Matching full human Ensembl 114 proteins and every source gene", flush=True)
        peptides = sorted({ligand.peptide for audit in audits for ligand in audit.ligands})
        antigen = contexts[0].native_antigen
        if antigen.self_reference_excluded_gene_ids:
            raise ValueError("Mutation pilot must not exclude self-source genes")
        matches = self_reference_matches(peptides, antigen, genome)
        reference_paths = [genome.gtf_path, *(genome.protein_fasta_paths or [])]
        if any(p is None or not Path(p).is_file() for p in reference_paths) or len(reference_paths) < 2:
            raise ValueError("Full reference source paths are unavailable")
        reference = dict(
            release=114, assembly=genome.reference_name, species=genome.species.latin_name,
            dataset_identity=ensembl_dataset_cache_identity(genome),
            annotated_protein_ids=len(genome.protein_ids()),
            protein_coding_transcripts=len(genome.transcript_ids(biotype="protein_coding")),
            files={Path(p).name: digest(p) for p in reference_paths},
            source="https://ftp.ensembl.org/pub/release-114/",
            policy="All protein-coding transcripts with a protein sequence; no gene exclusions")
        if reference["dataset_identity"] is None or reference["protein_coding_transcripts"] < 80000:
            raise ValueError("Expected a fingerprinted full Ensembl 114 human reference, not a subset")
        cta_ids = sorted(oncoref_cta_source_gene_ids())
        cta = dict(
            package_version=packages["oncoref"], gene_ids=cta_ids,
            gene_ids_sha256=sha256(("\n".join(cta_ids) + "\n").encode()).hexdigest(),
            policy="oncoref.cta.cta_unfiltered_gene_ids; broad CTA candidate universe, not tumor-specific admission",
            interpretation="Non-CTA means outside this versioned candidate universe, not proof of normal-tissue expression or safety")
        output.mkdir(parents=True)
        payload = dict(audits=audits, self_matches=matches, rna_results=rna_results,
                       self_reference=reference, cta_catalog=cta)
        data = gzip.compress(to_native_json(payload).encode(), mtime=0)
        (output / "predictions.json.gz").write_bytes(data)
        metadata.update(
            package_versions=packages, output_sha256=sha256(data).hexdigest(),
            documented_sha256=digest(DATA / "documented.json"),
            rna_manifest_sha256=digest(DATA / "rna" / "manifest.json"),
            reconstruction_manifest_sha256=digest(SELECTION_DATA / "isovar" / "manifest.json"),
            generator_sha256=digest(Path(__file__)),
            peptide_lengths=list(LENGTHS), prediction_kinds=list(KINDS),
            alleles=DOCUMENTED["prediction_alleles"],
            reconstruction="Complete original regional BAMs, not selected fixtures. Vaxrank defaults: preferred vaccine length 25, adaptive RNA context, balanced, 0.85 compatible read-name support, independent two-read-object per-base floor; no gate overrides or DNA fallback",
            method="NetMHCpan 4.2c -BA affinity and EL, all observations; Pepsickle 0.1.3 default proteasome model; CPN initial-terminal qualitative rule",
            scope="Present-day model evidence, not historical selection agreement or clinical safety. Final peptide chemistry is explicitly conditional; complete historical mRNA products are unavailable.")
        (output / "manifest.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
        print("Saved", len(data), "bytes;", len(matches), "peptides,",
              sum(m.occurs for m in matches.values()), "exact self matches", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--regional-cache", type=Path, required=True)
    args = parser.parse_args()
    main(args.output, args.regional_cache)
