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

import argparse

import msgspec

from ..config.loader import extract_epitope_config_kwargs, load_vaxrank_config
from ..epitope_config import EpitopeConfig


def add_epitope_prediction_args(arg_parser : argparse.ArgumentParser):
    epitope_prediction_args = arg_parser.add_argument_group("T-cell epitope prediction options")
    epitope_prediction_args.add_argument(
        "--min-epitope-score",
        type=float,
        help=(
            "Ignore predicted MHC ligands whose normalized binding score "
            "falls below this threshold."))
    epitope_prediction_args.add_argument(
        "--scoring-mode",
        choices=["affinity", "percentile_rank"],
        default=None,
        help=(
            "How to score epitopes. 'affinity' (default) uses IC50 binding "
            "affinity with logistic scoring. 'percentile_rank' uses the "
            "percentile rank directly (lower rank = better)."))

    # Per-Kind default methods. Apply whenever an input exposes multiple
    # models for a Kind (e.g. LENS with both mhcflurry and netmhcpan for
    # pMHC_affinity) — the unqualified ``affinity`` / ``stability`` /
    # ``presentation`` references in DSL formulas (and the synthesized
    # default score node) resolve to these. Bracketed references like
    # ``affinity['netmhcpan']`` still pick up their specific model. If
    # omitted and the data has multiple models, vaxrank auto-picks a
    # canonical method (mhcflurry > netmhcpan > alphabetical).
    default_method_args = arg_parser.add_argument_group(
        "Default predictor methods (for unqualified DSL references)")
    default_method_args.add_argument(
        "--default-affinity-predictor",
        default=None,
        metavar="METHOD",
        help="Method name to use for unqualified pMHC_affinity DSL "
             "references when multiple affinity models are present "
             "(e.g. 'mhcflurry').")
    default_method_args.add_argument(
        "--default-stability-predictor",
        default=None,
        metavar="METHOD",
        help="Method name to use for unqualified pMHC_stability DSL "
             "references (e.g. 'netmhcstabpan').")
    default_method_args.add_argument(
        "--default-presentation-predictor",
        default=None,
        metavar="METHOD",
        help="Method name to use for unqualified pMHC_presentation DSL "
             "references.")


def epitope_config_from_args(args : argparse.Namespace, merged_config=None) -> EpitopeConfig:
    """
    Extract config path and overrides from argument namespace.

    Parameters
    ----------
    args : argparse.Namespace
    merged_config : dict, optional
        Pre-loaded merged config dict.  When supplied the YAML file is
        not re-read, avoiding duplicate I/O when both epitope and
        vaccine configs are built from the same args.
    """
    if merged_config is None:
        merged_config = load_vaxrank_config(args)
    epitope_config_kwargs = extract_epitope_config_kwargs(merged_config)

    if args.min_epitope_score is not None:
        epitope_config_kwargs["min_epitope_score"] = args.min_epitope_score
    if getattr(args, 'scoring_mode', None) is not None:
        epitope_config_kwargs["scoring_mode"] = args.scoring_mode

    # Merge CLI --default-{affinity,stability,presentation}-predictor
    # into default_methods. CLI values override YAML.
    cli_defaults = {}
    for cli_attr, kind in (
        ("default_affinity_predictor", "pMHC_affinity"),
        ("default_stability_predictor", "pMHC_stability"),
        ("default_presentation_predictor", "pMHC_presentation"),
    ):
        val = getattr(args, cli_attr, None)
        if val is not None:
            cli_defaults[kind] = val
    if cli_defaults:
        merged = dict(epitope_config_kwargs.get("default_methods") or {})
        merged.update(cli_defaults)
        epitope_config_kwargs["default_methods"] = merged

    epitope_config = msgspec.convert(epitope_config_kwargs, EpitopeConfig)
    return epitope_config
