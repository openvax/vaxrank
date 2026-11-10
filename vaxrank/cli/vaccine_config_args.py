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

from ..config.loader import (
    extract_manufacturability_config_kwargs,
    extract_vaccine_config_kwargs,
    load_vaxrank_config,
)
from ..manufacturability_config import ManufacturabilityConfig
from ..vaccine_config import VaccineConfig
from ..vaccine_config import LENS_ANTIGEN_SOURCES
from ..vaccine_antigen import ANTIGEN_KINDS



def add_vaccine_peptide_args(arg_parser : argparse.ArgumentParser) -> None:
    default_vaccine_config = VaccineConfig()

    vaccine_peptide_group = arg_parser.add_argument_group("Vaccine peptide options")
    vaccine_peptide_group.add_argument(
        '--max-windows-per-antigen-source', type=int, default=None,
        help='Maximum distinct non-overlapping windows per antigen source (default: 3). '
             'Manufacturing alternatives do not consume additional windows.')
    allocation_group = arg_parser.add_argument_group(
        'Antigen type selection',
        'Counts are unique antigen sequences across the whole vaccine after deduplication: '
        'long windows or individual minimal epitopes. Minimums are representation goals; '
        'maximums are hard limits.')
    allocation_group.add_argument(
        '--antigen-allocation', choices=('balanced', 'ranked'), default=None,
        help='Selection order within the type limits: both modes first serve minimum goals; '
             'balanced then favors the least represented eligible type, while ranked keeps '
             'score/allele-coverage order (default: balanced). Scores are unchanged.')
    for kind in sorted(ANTIGEN_KINDS, key=str.lower):
        allocation_group.add_argument(
            '--min-%s-antigens' % kind.lower(), type=int, default=None, metavar='COUNT',
            help='Minimum requested %s antigen count (default: 0). A goal subject to eligibility '
                 'and vaccine capacity; unmet goals are reported.' % kind)
        allocation_group.add_argument(
            '--max-%s-antigens' % kind.lower(), type=int, default=None, metavar='COUNT',
            help='Maximum %s antigen count (default: no type-specific cap). '
                 'Set to 0 to exclude this type from the vaccine.' % kind)
    vaccine_peptide_group.add_argument(
        "--vaccine-peptide-length",
        default=None,
        type=int,
        help=(
            "Preferred number of amino acids in the vaccine peptides. "
            f"(default: {default_vaccine_config.preferred_peptide_length})"
        ))

    vaccine_peptide_group.add_argument(
        "--padding-around-mutation",
        default=None,
        type=int,
        help=(
            "Explicit legacy context request: vaccine peptide length + 2*padding. "
            "Without this option, RNA context adapts to peptide size via Isovar. "
            "--protein-sequence-length takes precedence for RNA. "
            f"DNA-only fallback padding defaults to {default_vaccine_config.padding_around_mutation}."
        ))

    vaccine_peptide_group.add_argument(
        "--max-vaccine-peptides-per-variant",
        "--max-vaccine-peptides-per-mutation",  # legacy alias
        dest="max_vaccine_peptides_per_variant",
        default=None,
        type=int,
        help=(
            "Legacy overlapping candidate limit (saved designs with max_windows_per_antigen_source=null). "
            "For new designs use --max-windows-per-antigen-source. "
            f"(default: {default_vaccine_config.max_vaccine_peptides_per_variant})"
        ))

    vaccine_peptide_group.add_argument(
        "--num-epitopes-per-vaccine-peptide",
        default=None,
        type=int,
        help=(
            "Maximum number of mutant epitopes to consider when scoring "
            "each vaccine peptide. Set to 0 to keep all epitopes. "
            f"(default: {default_vaccine_config.num_target_epitopes_to_keep})"
        ))

    vaccine_peptide_group.add_argument(
        "--include-antigen-source",
        dest="included_antigen_sources",
        action="append",
        choices=sorted(LENS_ANTIGEN_SOURCES),
        default=None,
        help=(
            "Add a LENS antigen source to construct ranking. Repeat for "
            "multiple sources. SNV, INDEL, and FUSION are included by "
            "default; SPLICE, CTA/SELF, and ERV require explicit inclusion."
        ),
    )
    vaccine_peptide_group.add_argument(
        "--exclude-antigen-source",
        dest="excluded_antigen_sources",
        action="append",
        choices=sorted(LENS_ANTIGEN_SOURCES),
        default=None,
        help=(
            "Remove a LENS antigen source from construct ranking. Repeat for "
            "multiple sources. Exclusion wins if a source is both included "
            "and excluded."
        ),
    )



    

def vaccine_config_from_args(args : argparse.Namespace, merged_config=None) -> VaccineConfig:
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
    vaccine_config_kwargs = extract_vaccine_config_kwargs(merged_config)
    if getattr(args, 'max_windows_per_antigen_source', None) is not None:
        vaccine_config_kwargs['max_windows_per_antigen_source'] = args.max_windows_per_antigen_source
    vaccine_config_kwargs['antigen_allocation'] = antigen_allocation_from_args(
        args, vaccine_config_kwargs.get('antigen_allocation'))

    if args.vaccine_peptide_length is not None:
        vaccine_config_kwargs["preferred_peptide_length"] = args.vaccine_peptide_length
        # When user specifies a single length via CLI, pin min/max to match
        vaccine_config_kwargs["min_peptide_length"] = args.vaccine_peptide_length
        vaccine_config_kwargs["max_peptide_length"] = args.vaccine_peptide_length
    if args.padding_around_mutation is not None:
        vaccine_config_kwargs["padding_around_mutation"] = args.padding_around_mutation
    max_vaccine_peptides_per_variant = getattr(
        args, "max_vaccine_peptides_per_variant", None)
    if max_vaccine_peptides_per_variant is None:
        max_vaccine_peptides_per_variant = getattr(
            args, "max_vaccine_peptides_per_mutation", None)
    if max_vaccine_peptides_per_variant is not None:
        vaccine_config_kwargs["max_vaccine_peptides_per_variant"] = (
            max_vaccine_peptides_per_variant
        )
    if args.num_epitopes_per_vaccine_peptide is not None:
        vaccine_config_kwargs["num_target_epitopes_to_keep"] = args.num_epitopes_per_vaccine_peptide
    included_antigen_sources = getattr(args, "included_antigen_sources", None)
    excluded_antigen_sources = getattr(args, "excluded_antigen_sources", None)
    if included_antigen_sources or excluded_antigen_sources:
        selected_sources = set(vaccine_config_kwargs.get(
            "included_antigen_sources",
            VaccineConfig().included_antigen_sources,
        ))
        selected_sources.update(included_antigen_sources or ())
        selected_sources.difference_update(excluded_antigen_sources or ())
        vaccine_config_kwargs["included_antigen_sources"] = tuple(
            source
            for source in sorted(LENS_ANTIGEN_SOURCES)
            if source in selected_sources
        )

    vaccine_config = msgspec.convert(vaccine_config_kwargs, VaccineConfig)
    return vaccine_config


def antigen_allocation_from_args(args, configured=None):
    """Merge shared allocation CLI overrides with the resolved YAML/replay policy."""
    from dataclasses import asdict
    from ..antigen_allocation import AntigenAllocation
    policy = asdict(AntigenAllocation.resolve(configured))
    if getattr(args, 'antigen_allocation', None) is not None:
        policy['mode'] = args.antigen_allocation
    for bound in ('min', 'max'):
        for kind in ANTIGEN_KINDS:
            count = getattr(args, '%s_%s_antigens' % (bound, kind.lower()), None)
            if count is not None:
                policy[bound + '_per_kind'][kind] = count
    return asdict(AntigenAllocation.resolve(policy))


def manufacturability_config_from_args(
        args: argparse.Namespace, merged_config=None) -> ManufacturabilityConfig:
    """Build a ``ManufacturabilityConfig`` from the YAML
    ``peptide.manufacturability`` sub-section. Always returns a real
    config (``ManufacturabilityConfig()`` defaults if the YAML
    sub-section is absent) — the ranker treats it as the source of
    peptide-synthesis difficulty thresholds + rules."""
    if merged_config is None:
        merged_config = load_vaxrank_config(args)
    kwargs = extract_manufacturability_config_kwargs(merged_config)
    return msgspec.convert(kwargs, ManufacturabilityConfig)
