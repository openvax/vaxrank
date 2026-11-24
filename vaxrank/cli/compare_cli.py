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

"""CLI handler for ``vaxrank compare <run_a> <run_b>`` (openvax/vaxrank#920)."""

from __future__ import annotations

import argparse
import json
import sys

from ..compare import (
    ComparePolicy,
    compare_runs,
    render_comparison_ascii,
    write_comparison_outputs,
)


def _parse_top_k(value: str) -> tuple[int, ...]:
    parts = [p.strip() for p in value.replace(" ", ",").split(",") if p.strip()]
    if not parts:
        raise argparse.ArgumentTypeError("--top-k must list at least one positive integer")
    try:
        nums = tuple(int(p) for p in parts)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"--top-k expects comma-separated integers; got {value!r}"
        ) from exc
    if any(n < 1 for n in nums):
        raise argparse.ArgumentTypeError("--top-k values must be >= 1")
    return nums


def _build_parser(prog: str = "vaxrank compare") -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Compare two vaxrank runs (output directories, ranked_candidates.csv, "
            "or manifest.json bundles) for rank preservation, Top-K target "
            "attribution, and peptide/mRNA construct diffs."
        ),
    )
    parser.add_argument("run_a", help="Baseline run directory, ranked_candidates.csv, or manifest.json")
    parser.add_argument("run_b", help="Comparison run directory, ranked_candidates.csv, or manifest.json")
    parser.add_argument(
        "--output-dir", "-o", default=None,
        help="Write comparison_summary.json, comparable_targets.csv, "
             "comparison_report.txt, comparison_report.html (and optional PDF) "
             "into this directory. Omit to print to stdout.",
    )
    parser.add_argument("--label-a", default=None, help="Display label for run A")
    parser.add_argument("--label-b", default=None, help="Display label for run B")
    parser.add_argument(
        "--config", default=None,
        help="Optional per-trial YAML policy file supplying compare: thresholds.",
    )
    parser.add_argument(
        "--top-k", type=_parse_top_k, default=None,
        help="Comma-separated Top-K cutoffs for overlap/Jaccard (default: 5,10,20).",
    )
    parser.add_argument(
        "--primary-top-k", type=int, default=None,
        help="Primary Top-K cutoff for gained/lost attribution tables (default: 10).",
    )
    parser.add_argument(
        "--pdf", action="store_true", default=False,
        help="Also render comparison_report.pdf via WeasyPrint (requires --output-dir).",
    )
    parser.add_argument(
        "--json", action="store_true", default=False,
        help="Print machine-readable comparison_summary JSON to stdout.",
    )
    return parser


def run(argv: list[str] | None = None, prog: str = "vaxrank compare") -> int:
    """Execute ``vaxrank compare <run_a> <run_b> [options]``."""
    parser = _build_parser(prog=prog)
    args = parser.parse_args(argv)

    policy_dict: dict = {}
    if args.config:
        from ..config import load_vaxrank_config
        cfg = load_vaxrank_config(config_path=args.config)
        if isinstance(cfg.get("compare"), dict):
            policy_dict.update(cfg["compare"])
    if args.top_k is not None:
        policy_dict["top_k_values"] = args.top_k
    if args.primary_top_k is not None:
        policy_dict["primary_top_k"] = args.primary_top_k
    if args.pdf:
        policy_dict["include_pdf"] = True

    policy = ComparePolicy.resolve(policy_dict or None)
    comparison = compare_runs(
        args.run_a,
        args.run_b,
        policy=policy,
        label_a=args.label_a,
        label_b=args.label_b,
    )

    if args.output_dir:
        written = write_comparison_outputs(comparison, args.output_dir, policy=policy)
        if args.json:
            sys.stdout.write(json.dumps(comparison, indent=2) + "\n")
        else:
            sys.stdout.write(render_comparison_ascii(comparison))
            sys.stdout.write("\nComparison outputs written:\n")
            for kind, path in written.items():
                sys.stdout.write(f"  {kind:<8s} {path}\n")
    else:
        if args.json:
            sys.stdout.write(json.dumps(comparison, indent=2) + "\n")
        else:
            sys.stdout.write(render_comparison_ascii(comparison))
    return 0
