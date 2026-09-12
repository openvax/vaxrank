"""Maintainer-only real NetMHCpan cache generation, never run by CI.

Run from the repository root with published Isovar 1.8.1 and NetMHCpan 4.2:
    python -m tests.data.osteosarc.selection_validation.generate_predictions --output NEW_DIRECTORY
Does not generate or rewrite the independent documented sequence expectations.
"""

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import importlib
import json
import logging
from pathlib import Path
import platform
import shutil
import tempfile

import msgspec
from mhctools import NetMHCpan42
from topiary import CachedPredictor, TopiaryPredictor

from tests.osteosarc_selection_helpers import DATA, load_selection_inputs, reconstruct_selection
from vaxrank.mutant_protein_fragment import MutantProteinFragment
from vaxrank.epitope_config import EpitopeConfig
from vaxrank.vaccine_config import VaccineConfig


def main(output_directory):
    if output_directory.exists():
        raise ValueError("Choose a new output directory; existing real predictions are never overwritten")
    logging.disable(logging.INFO)
    documented = json.loads((DATA / "documented.json").read_text())
    contexts = {}
    requests = {(r["variant_id"], len(r["native_sequence"])) for r in documented["records"]}
    with tempfile.TemporaryDirectory(prefix="vaxrank-sid-prediction-") as directory:
        inputs = load_selection_inputs(Path(directory))
        requests.update((v, 25) for v in inputs[1])
        peptides = set()
        lengths = documented["prediction_peptide_lengths"]
        for variant_id, size in sorted(requests):
            result, _ = reconstruct_selection(inputs, variant_id, size)
            if result.top_protein_sequence is None:
                continue
            fragment = MutantProteinFragment.from_isovar_result(result)
            name = "%s-size%d" % (variant_id, size)
            contexts[name] = fragment.amino_acids
            # Include every candidate k-mer and every position-aligned WT
            # comparator, before filtering, not only the eventual winners.
            reference = fragment.predicted_effect().original_protein_sequence
            start = fragment.global_start_pos()
            for length in lengths:
                for offset in range(len(fragment) - length + 1):
                    peptides.add(fragment.amino_acids[offset:offset + length])
                    wt = reference[start + offset:start + offset + length]
                    if len(wt) == length:
                        peptides.add(wt)
        model = NetMHCpan42(alleles=documented["prediction_alleles"],
                           default_peptide_lengths=lengths, process_limit=1)
        predictor = TopiaryPredictor(models=[model])
        named = {"peptide-%04d" % i: p for i, p in enumerate(sorted(peptides))}
        frame = predictor.predict_from_named_peptides(named)
        cache = CachedPredictor.from_dataframe(frame)
        output_directory.mkdir(parents=True)
        output = output_directory / "netmhcpan42.tsv"
        cache.save(str(output))
        executable = Path(shutil.which("netMHCpan")).resolve()
        model_root = executable.parent
        binary = model_root / (platform.system() + "_" + platform.machine()) / "bin/netMHCpan-4.2"
        if not (model_root / "data").is_dir() or not binary.is_file():
            raise ValueError("Expected the netmhc-bundle 4.2 layout for model-asset provenance")
        asset_paths = [binary, *(p for p in (model_root / "data").rglob("*") if p.is_file())]
        assets = {str(p.relative_to(model_root)): sha256(p.read_bytes()).hexdigest()
                  for p in sorted(asset_paths)}
        metadata = dict(
            generated=datetime.now(timezone.utc).isoformat(timespec="seconds"),
            purpose="Present-day class-I comparison, not historical provider predictions",
            method="NetMHCpan 4.2c -BA binding-affinity and elution-score outputs, default model settings",
            predictor_name=cache.prediction_method_name, predictor_version=cache.predictor_version,
            package_versions={p: getattr(importlib.import_module(p), "__version__", None)
                              for p in ("isovar", "mhctools", "topiary", "varcode", "vaxrank")},
            executable_name=executable.name, executable_sha256=sha256(executable.read_bytes()).hexdigest(),
            model_assets_sha256=assets,
            generator_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
            documented_sha256=sha256((DATA / "documented.json").read_bytes()).hexdigest(),
            input_manifest_sha256=sha256((DATA / "isovar/manifest.json").read_bytes()).hexdigest(),
            output_sha256=sha256(output.read_bytes()).hexdigest(),
            contexts=contexts, requested_peptides=sorted(peptides),
            alleles=documented["prediction_alleles"], peptide_lengths=lengths,
            reconstruction="Vaxrank adaptive context for each native peptide size; balanced, support fraction 0.85, coverage floor 2, assembly enabled; no DNA fallback",
            vaccine_config=msgspec.to_builtins(VaccineConfig()),
            epitope_config=msgspec.to_builtins(EpitopeConfig()),
            self_reference_scope="131-transcript GRCh38 fixture only, not a whole-human safety audit")
        (output_directory / "predictions_manifest.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
        print("Cached", len(frame), "real prediction rows for", len(peptides), "peptides")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    main(parser.parse_args().output)
