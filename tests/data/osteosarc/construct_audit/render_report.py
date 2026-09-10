"""Render the independently documented Sid pilot from verified offline evidence.

python -m tests.data.osteosarc.construct_audit.render_report --output DIR
No prediction backend, genome download or RNA reprocessing is invoked here.
"""

import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path

import jinja2
from mhctools.pred import value_unit
import pandas as pd
from topiary.ranking import apply_filter, parse

from vaxrank.context_report import write_sequence_context_audits
from vaxrank.native_serialization import from_native_json
from tests.osteosarc_construct_helpers import DATA, DOCUMENTED, RNA, verify_rna_inputs
from tests.osteosarc_selection_helpers import DATA as SELECTION_DATA


STRONG_POLICY = "NetMHCpan 4.2c, separate affinity/EL percentile ranks < 0.5 (strict); no ranking changes"
POLICY_SOURCE = "https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/"


def load_evidence():
    verify_rna_inputs()
    directory = DATA / "predictions"
    metadata = json.loads((directory / "manifest.json").read_text())
    for path, key in ((directory / "predictions.json.gz", "output_sha256"),
                      (DATA / "documented.json", "documented_sha256"),
                      (DATA / "rna" / "manifest.json", "rna_manifest_sha256"),
                      (SELECTION_DATA / "isovar" / "manifest.json", "reconstruction_manifest_sha256")):
        if sha256(path.read_bytes()).hexdigest() != metadata[key]:
            raise ValueError("Sid audit input/output checksum mismatch: " + str(path))
    evidence = from_native_json(gzip.decompress((directory / "predictions.json.gz").read_bytes()).decode(), dict)
    inventory = json.loads(gzip.decompress((DATA / "rna" / "inventory.json.gz").read_bytes()))
    return evidence, metadata, inventory


def strong_predictions(audits):
    """Use the Topiary DSL, preserving kind/model/HLA/occurrence boundaries."""
    rows = [dict(
        source_sequence_name=a.context.source_id, peptide=ligand.peptide, peptide_offset=ligand.start,
        peptide_length=len(ligand.peptide), kind=p.kind, prediction_method_name=p.predictor_name,
        predictor_version=p.predictor_version, allele=p.allele, score=p.score,
        value=p.value, percentile_rank=p.percentile_rank)
        for a in audits for ligand in a.ligands for p in ligand.predictions]
    frame = pd.DataFrame(rows)
    keys = ["source_sequence_name", "peptide", "peptide_offset", "kind",
            "prediction_method_name", "predictor_version", "allele"]
    selected = []
    if frame.empty:
        return set()
    for kind, expression in (("pMHC_affinity", "affinity.percentile_rank < 0.5"),
                             ("pMHC_presentation", "presentation.percentile_rank < 0.5")):
        part = frame[frame["kind"] == kind]
        selected.extend(apply_filter(part, parse(expression), group_keys=keys)[keys]
                        .itertuples(index=False, name=None))
    return set(selected)


def comparison_rows(evidence):
    """Overlay independently intended ligands and retain every self-source gene."""
    strong = strong_predictions(evidence["audits"])
    cta_ids = set(evidence["cta_catalog"]["gene_ids"])
    rows = []
    for audit in evidence["audits"]:
        context = audit.context
        target = DOCUMENTED["intended_ligand"]
        offsets = [i for i in range(len(context.sequence) - len(target) + 1)
                   if context.sequence[i:i + len(target)] == target]
        if len(offsets) != 1:
            raise ValueError("Pilot intended ligand must map uniquely; do not guess a position")
        start, = offsets
        intervals = [(p, p.interval_evidence(start, start + len(target))) for p in audit.profiles]
        findings = []
        for ligand in audit.ligands:
            match = evidence["self_matches"][ligand.peptide]
            non_cta = tuple(s for s in match.sources if s.gene_id not in cta_ids)
            flags = tuple(p for p in ligand.predictions if (
                context.source_id, ligand.peptide, ligand.start, p.kind,
                p.predictor_name, p.predictor_version, p.allele) in strong)
            if non_cta and flags and (ligand.start, ligand.peptide) != (start, target):
                findings.append(dict(ligand=ligand, predictions=flags, self_match=match,
                                     non_cta_sources=non_cta,
                                     cta_sources=tuple(s for s in match.sources if s.gene_id in cta_ids),
                                     overlaps_intended=(ligand.start < start + len(target) and start < ligand.end)))
        rows.append(dict(audit=audit, target_start=start, target_end=start + len(target),
                         target_profiles=intervals, non_target_self_findings=findings))
    return rows


def render(output):
    evidence, metadata, inventory = load_evidence()
    rows = comparison_rows(evidence)
    environment = jinja2.Environment(
        loader=jinja2.FileSystemLoader(str(DATA)), autoescape=True,
        undefined=jinja2.StrictUndefined)
    html = environment.get_template("report.html.jinja").render(
        rows=rows, evidence=evidence, metadata=metadata, inventory=inventory,
        documented=DOCUMENTED, rna=RNA, policy=STRONG_POLICY, policy_source=POLICY_SOURCE,
        value_unit=value_unit,
        sources={s["source_id"]: s for s in inventory["sources"]})
    output.mkdir(parents=True, exist_ok=True)
    (output / "sid-construct-audit.html").write_text(html)
    write_sequence_context_audits(evidence["audits"], json_path=output / "full-contexts.json",
                                 html_path=output / "full-contexts.html")
    print(output / "sid-construct-audit.html")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    render(parser.parse_args().output)
