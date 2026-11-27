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

"""CLI handler for ``vaxrank mhc-status`` (openvax/vaxrank#936).

Builds a canonical, cross-species ``mhc_allele_status.json`` (plus optional
``tumor_cluster_ntpm.tsv`` and ``scrna_cluster_summary.json``) from any
combination of:
- Upstream companion tools (T1K, arcasHLA, HLApers, LOHHLA, DASH, LILAC, POLYSOLVER)
- Direct computational bulk RNA-seq personalized-allele discriminating k-mer ASE
- Single-cell RNA-seq tumor-cluster extraction, per-cluster gene nTPM, and
  per-cluster MHC allele UMI aggregation (scHLAcount / barcode-allele UMI tables)
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import sys

import argcomplete

from .help import HelpArgumentParser

logger = logging.getLogger(__name__)


def _build_parser(prog: str = "vaxrank mhc-status") -> argparse.ArgumentParser:
    parser = HelpArgumentParser(
        prog=prog,
        description=(
            "Build a cross-species mhc_allele_status.json (and optional scRNA "
            "tumor-cluster gene nTPM table) from upstream MHC tools, direct "
            "bulk RNA-seq k-mer ASE, and/or single-cell tumor-cluster analysis."
        ),
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Destination path for canonical mhc_allele_status.json.",
    )
    parser.add_argument(
        "--config",
        default=None,
        help="Optional per-trial YAML policy file (reads mhc_expression: and scrna_tumor:).",
    )
    parser.add_argument(
        "--tumor-purity",
        type=float,
        default=None,
        help="Estimated tumor purity in [0.0, 1.0] for bulk RNA ASE confounding guards.",
    )
    parser.add_argument(
        "--base-status",
        default=None,
        help="Optional existing mhc_allele_status JSON/TSV to merge underneath new assays.",
    )

    # Upstream companion tool ingestors
    for tool in ("t1k", "arcashla", "hlapers", "lohhla", "dash", "lilac", "polysolver"):
        parser.add_argument(
            f"--from-{tool}",
            dest=f"from_{tool}",
            default=None,
            metavar="PATH",
            help=f"Ingest MHC allele status / expression from {tool.upper()} output file.",
        )

    # Direct computational bulk RNA-seq k-mer ASE
    parser.add_argument(
        "--allele-fasta",
        default=None,
        metavar="FASTA",
        help="Personalized patient MHC allele nucleotide FASTA (cross-species: HLA/DLA/H2/SLA) "
             "for direct discriminating k-mer bulk RNA-seq ASE estimation.",
    )
    parser.add_argument(
        "--bulk-rna-reads",
        nargs="+",
        default=None,
        metavar="FASTQ",
        help="Bulk tumor RNA-seq FASTQ/FASTQ.gz/FASTA read file(s) paired with --allele-fasta.",
    )

    # Single-cell tumor-cluster extraction + per-cluster nTPM + MHC ASE
    parser.add_argument(
        "--scrna-counts",
        default=None,
        metavar="PATH",
        help="scRNA count matrix (10x MTX directory, TSV/CSV, or .h5ad).",
    )
    parser.add_argument(
        "--scrna-clusters",
        default=None,
        metavar="PATH",
        help="Per-barcode cluster assignment TSV/CSV (barcode, cluster).",
    )
    parser.add_argument(
        "--scrna-somatic-barcodes",
        default=None,
        metavar="PATH",
        help="Optional per-barcode somatic variant calls (e.g. VarTrix/cellsnp-lite TSV/CSV).",
    )
    parser.add_argument(
        "--scrna-mhc-umis",
        default=None,
        metavar="PATH",
        help="Optional per-barcode MHC allele UMIs (scHLAcount directory or TSV/CSV).",
    )
    parser.add_argument(
        "--scrna-chemistry",
        default=None,
        help="scRNA library chemistry (e.g. 10x_5p, 10x_3p, smartseq2).",
    )
    parser.add_argument(
        "--scrna-tumor-markers",
        default=None,
        help="Optional comma-separated tumor lineage marker gene symbols (e.g. EPCAM,KRT18).",
    )
    parser.add_argument(
        "--output-tumor-ntpm",
        default=None,
        metavar="PATH",
        help="Optional path to write per-gene tumor-cluster vs immune-control nTPM TSV.",
    )
    parser.add_argument(
        "--output-scrna-summary",
        default=None,
        metavar="PATH",
        help="Optional path to write scRNA cluster classification summary JSON.",
    )
    return parser


def run(argv: list[str] | None = None, prog: str = "vaxrank mhc-status") -> int:
    parser = _build_parser(prog=prog)
    argcomplete.autocomplete(parser)
    args = parser.parse_args(argv)

    from ..config.loader import load_vaxrank_config
    from ..mhc_allele_status import (
        MHCAlleleStatusBundle,
        load_mhc_allele_status,
        merge_mhc_allele_status_bundles,
    )
    from ..mhc_expression import (
        MHCExpressionPolicy,
        estimate_bulk_mhc_ase_from_reads,
        load_mhc_status_from_tool,
    )
    from ..scrna_tumor import (
        SCRNATumorClusterPolicy,
        analyze_scrna_tumor_clusters,
        load_cluster_assignments,
        load_scrna_counts,
        load_scrna_mhc_allele_umis,
        load_somatic_barcode_calls,
    )

    mhc_expr_cfg: dict = {}
    scrna_cfg: dict = {}
    if args.config:
        cfg = load_vaxrank_config(config_path=args.config)
        if isinstance(cfg.get("mhc_expression"), dict):
            mhc_expr_cfg.update(cfg["mhc_expression"])
        if isinstance(cfg.get("scrna_tumor"), dict):
            scrna_cfg.update(cfg["scrna_tumor"])
    if args.scrna_chemistry is not None:
        scrna_cfg["chemistry"] = args.scrna_chemistry
    if args.scrna_tumor_markers:
        scrna_cfg["tumor_markers"] = [
            m.strip() for m in args.scrna_tumor_markers.split(",") if m.strip()
        ]

    expr_policy = MHCExpressionPolicy.resolve(mhc_expr_cfg)
    scrna_policy = SCRNATumorClusterPolicy.resolve(scrna_cfg)

    merged: MHCAlleleStatusBundle | None = None
    if args.base_status:
        merged = load_mhc_allele_status(args.base_status)

    # 1. Upstream companion tools (DNA LOH/mutations + bulk RNA ASE)
    for tool in ("lohhla", "dash", "lilac", "polysolver", "t1k", "arcashla", "hlapers"):
        path = getattr(args, f"from_{tool}", None)
        if path:
            bundle = load_mhc_status_from_tool(
                path,
                tool,
                tumor_purity=args.tumor_purity,
                policy=expr_policy,
            )
            merged = merge_mhc_allele_status_bundles(merged, bundle)

    # 2. Direct computational bulk RNA-seq discriminating k-mer ASE
    if args.allele_fasta or args.bulk_rna_reads:
        if not (args.allele_fasta and args.bulk_rna_reads):
            parser.error("--allele-fasta and --bulk-rna-reads must be supplied together.")
        kmer_bundle = estimate_bulk_mhc_ase_from_reads(
            args.allele_fasta,
            args.bulk_rna_reads,
            tumor_purity=args.tumor_purity,
            policy=expr_policy,
        )
        merged = merge_mhc_allele_status_bundles(merged, kmer_bundle)

    # 3. Single-cell tumor-cluster extraction + per-cluster nTPM + MHC ASE
    if args.scrna_counts or args.scrna_clusters:
        if not (args.scrna_counts and args.scrna_clusters):
            parser.error("--scrna-counts and --scrna-clusters must be supplied together.")
        counts = load_scrna_counts(args.scrna_counts)
        clusters = load_cluster_assignments(args.scrna_clusters)
        somatic = (
            load_somatic_barcode_calls(args.scrna_somatic_barcodes)
            if args.scrna_somatic_barcodes
            else None
        )
        mhc_umis = (
            load_scrna_mhc_allele_umis(args.scrna_mhc_umis)
            if args.scrna_mhc_umis
            else None
        )
        scrna_res = analyze_scrna_tumor_clusters(
            counts,
            clusters,
            somatic_barcodes=somatic,
            mhc_allele_umis=mhc_umis,
            policy=scrna_policy,
        )
        merged = merge_mhc_allele_status_bundles(
            merged, scrna_res.to_mhc_allele_status_bundle()
        )
        if args.output_tumor_ntpm:
            scrna_res.write_tumor_cluster_ntpm_tsv(args.output_tumor_ntpm)
        if args.output_scrna_summary:
            scrna_res.write_summary_json(args.output_scrna_summary)

    if merged is None:
        parser.error(
            "Provide at least one MHC/scRNA evidence source "
            "(--from-t1k/--from-lohhla/..., --allele-fasta + --bulk-rna-reads, "
            "or --scrna-counts + --scrna-clusters)."
        )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as handle:
        json.dump(merged.to_dict(), handle, indent=2, allow_nan=False)

    print(
        f"Wrote MHC allele status bundle ({len(merged.alleles)} allele(s)"
        f"{', APM recorded' if merged.apm and merged.apm.any_recorded() else ''}) "
        f"to {out_path}",
        file=sys.stderr,
    )
    return 0
