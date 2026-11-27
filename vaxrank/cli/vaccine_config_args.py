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
from typing import TYPE_CHECKING

from ..config.defaults import (
    ANTIGEN_SELECTION_CHOICES,
    ANTIGEN_SELECTION_TYPES,
    DEFAULT_MAX_VACCINE_PEPTIDES_PER_VARIANT,
    DEFAULT_NUM_MUTANT_EPITOPES_TO_KEEP,
    DEFAULT_PADDING_AROUND_MUTATION,
    DEFAULT_PREFERRED_PEPTIDE_LENGTH,
    LENS_ANTIGEN_SOURCES,
    SHARED_SEQUENCE_RNA_POLICIES,
    SV_TYPES,
    VARIANT_TYPES,
    normalize_selection_type,
)

if TYPE_CHECKING:
    from ..manufacturability_config import ManufacturabilityConfig
    from ..vaccine_config import VaccineConfig



def add_vaccine_peptide_args(arg_parser : argparse.ArgumentParser) -> None:
    vaccine_peptide_group = arg_parser.add_argument_group("Vaccine peptide options")
    vaccine_peptide_group.add_argument(
        '--max-windows-per-antigen-source', type=int, default=None,
        help='Maximum distinct non-overlapping windows per antigen source (default: 3). '
             'Manufacturing alternatives do not consume additional windows.')
    vaccine_peptide_group.add_argument(
        '--shared-sequence-rna-policy',
        choices=sorted(SHARED_SEQUENCE_RNA_POLICIES),
        default=None,
        help=(
            'How candidate RNA support is evaluated when multiple origins share '
            'the exact same candidate amino-acid window: "source" (default, '
            'each origin uses its own RNA measurement) or "within_library_union" '
            '(compute exact within-library Isovar read/fragment set unions for '
            'rna_shared_* score bindings without pooling across libraries).'
        ),
    )
    allocation_group = arg_parser.add_argument_group(
        'Antigen type selection',
        'Counts are unique antigen sequences across the whole vaccine after deduplication: '
        'long windows or individual minimal epitopes. Minimums are representation goals; '
        'maximums are hard limits. Selection follows the existing score/allele-coverage '
        'order after serving explicit minimum goals, subject to maximum caps.')
    allocation_group.add_argument('--antigen-types', nargs='+', type=normalize_selection_type,
        choices=sorted(ANTIGEN_SELECTION_CHOICES), default=None,
        help='Antigen types eligible for vaccine products (default: all supplied types). '
             'Use variant for small and structural variants in any coding context. '
             'Historical mutation selects small variants outside recorded ncORFs. '
             'Source evidence remains available in reports.')
    allocation_group.add_argument('--variant-types', nargs='+', type=str.upper, default=None,
        choices=sorted(VARIANT_TYPES | {'SV:' + kind for kind in SV_TYPES}),
        help='Restrict genomic variant sources by recorded type (case insensitive). '
             'Examples: SNV INDEL; SV:DEL SV:INV. Other selected antigen categories are unaffected.')
    allocation_group.add_argument('--variant-contexts', nargs='+', type=str.lower, default=None,
        help='Restrict variant coding contexts: annotated, ncORF, an ORFdata orf_type such as '
             'uORF, or unknown. Multiple values are alternatives, case insensitive.')
    for kind in sorted(ANTIGEN_SELECTION_TYPES, key=str.lower):
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
            f"(default: {DEFAULT_PREFERRED_PEPTIDE_LENGTH})"
        ))

    vaccine_peptide_group.add_argument(
        "--padding-around-mutation",
        default=None,
        type=int,
        help=(
            "Explicit legacy context request: vaccine peptide length + 2*padding. "
            "Without this option, RNA context adapts to peptide size via Isovar. "
            "--protein-sequence-length takes precedence for RNA. "
            f"DNA-only fallback padding defaults to {DEFAULT_PADDING_AROUND_MUTATION}."
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
            f"(default: {DEFAULT_MAX_VACCINE_PEPTIDES_PER_VARIANT})"
        ))

    vaccine_peptide_group.add_argument(
        "--num-epitopes-per-vaccine-peptide",
        default=None,
        type=int,
        help=(
            "Maximum number of mutant epitopes to consider when scoring "
            "each vaccine peptide. Set to 0 to keep all epitopes. "
            f"(default: {DEFAULT_NUM_MUTANT_EPITOPES_TO_KEEP})"
        ))

    allocation_group.add_argument(
        "--include-antigen-type",
        dest="included_antigen_types",
        action="append",
        choices=sorted(ANTIGEN_SELECTION_CHOICES, key=str.lower), type=normalize_selection_type,
        default=None,
        help=(
            "Enable an eligible antigen type for vaccine selection (variants, "
            "fusions, CTA, LENS, assemblies). Repeat for several types. "
            "Lifts a zero cap set in YAML or by --max-TYPE-antigens 0; "
            "ncORF selects mutations in recorded noncanonical ORFs separately from ordinary mutations. "
            "Does not bypass source evidence or override LENS source labels "
            "(see --include-lens-source)."
        ))
    allocation_group.add_argument(
        "--exclude-antigen-type",
        dest="excluded_antigen_types",
        action="append",
        choices=sorted(ANTIGEN_SELECTION_CHOICES, key=str.lower), type=normalize_selection_type,
        default=None,
        help=(
            "Keep an antigen type out of the vaccine for any input. Same as "
            "--max-TYPE-antigens 0. Exclusion wins if a type is both "
            "included and excluded."
        ))
    lens_group = arg_parser.add_argument_group(
        "LENS source labels",
        "Apply only to --input-lens. LENS labels say where a peptide came "
        "from, not that it is a safe tumor-specific target, so SPLICE, "
        "CTA/SELF and ERV rows stay in the report but out of constructs "
        "until opted in.")
    lens_group.add_argument(
        "--include-lens-source", "--include-antigen-source",
        dest="included_lens_sources",
        action="append",
        choices=sorted(LENS_ANTIGEN_SOURCES), type=str.upper,
        default=None,
        help=(
            "Add a LENS source label to construct ranking. Repeat for "
            "multiple sources. SNV, INDEL, and FUSION are included by "
            "default; SPLICE, CTA/SELF, and ERV require explicit inclusion. "
            "Older spelling: --include-antigen-source."
        ),
    )
    lens_group.add_argument(
        "--exclude-lens-source", "--exclude-antigen-source",
        dest="excluded_lens_sources",
        action="append",
        choices=sorted(LENS_ANTIGEN_SOURCES), type=str.upper,
        default=None,
        help=(
            "Remove a LENS source label from construct ranking. Exclusion "
            "wins if a source is both included and excluded. Older "
            "spelling: --exclude-antigen-source."
        ),
    )




def vaccine_config_from_args(args: argparse.Namespace, merged_config=None) -> "VaccineConfig":
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
    import msgspec
    from ..config.loader import extract_vaccine_config_kwargs, load_vaxrank_config
    from ..vaccine_config import (
        PROXY_EXPRESSION_COMBINED_SCORE_EXPR,
        SOURCE_AGNOSTIC_COMBINED_SCORE_EXPR,
        VaccineConfig,
    )

    if merged_config is None:
        merged_config = load_vaxrank_config(args)
    vaccine_config_kwargs = extract_vaccine_config_kwargs(merged_config)
    if getattr(args, 'max_windows_per_antigen_source', None) is not None:
        vaccine_config_kwargs['max_windows_per_antigen_source'] = args.max_windows_per_antigen_source
    if getattr(args, 'shared_sequence_rna_policy', None) is not None:
        vaccine_config_kwargs['shared_sequence_rna_policy'] = args.shared_sequence_rna_policy
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
    included_lens_sources = getattr(args, "included_lens_sources", None)
    excluded_lens_sources = getattr(args, "excluded_lens_sources", None)
    if included_lens_sources or excluded_lens_sources:
        selected_sources = set(vaccine_config_kwargs.get(
            "included_lens_sources",
            VaccineConfig().included_lens_sources,
        ))
        selected_sources.update(included_lens_sources or ())
        selected_sources.difference_update(excluded_lens_sources or ())
        vaccine_config_kwargs["included_lens_sources"] = tuple(
            source
            for source in sorted(LENS_ANTIGEN_SOURCES)
            if source in selected_sources
        )

    if "combined_score_expr" not in vaccine_config_kwargs:
        if getattr(args, "proxy_expression_cancer_type", None):
            vaccine_config_kwargs["combined_score_expr"] = PROXY_EXPRESSION_COMBINED_SCORE_EXPR
        elif (
            getattr(args, "allow_dna_only_fallback", False)
            and not getattr(args, "bam", None)
            and any(getattr(args, name, None) for name in ("vcf", "maf", "variant", "json_variants"))
        ):
            vaccine_config_kwargs["combined_score_expr"] = SOURCE_AGNOSTIC_COMBINED_SCORE_EXPR

    vaccine_config = msgspec.convert(vaccine_config_kwargs, VaccineConfig)
    return vaccine_config


def antigen_allocation_from_args(args, configured=None):
    """Merge shared allocation CLI overrides with the resolved YAML/replay policy."""
    from ..antigen_allocation import AntigenAllocation, expand_selection_types
    policy = AntigenAllocation.resolve(configured).to_dict()
    selected = getattr(args, 'antigen_types', None)
    if selected is not None:
        selected = expand_selection_types(selected)
        for kind in ANTIGEN_SELECTION_TYPES:
            if kind not in selected:
                policy['max_per_kind'][kind] = 0
                policy['min_per_kind'].pop(kind, None)
            elif policy['max_per_kind'].get(kind) == 0:
                del policy['max_per_kind'][kind]
    for kind in expand_selection_types(getattr(args, 'included_antigen_types', None) or ()):
        if policy['max_per_kind'].get(kind) == 0:
            del policy['max_per_kind'][kind]
    for bound in ('min', 'max'):
        for kind in ANTIGEN_SELECTION_TYPES:
            count = getattr(args, '%s_%s_antigens' % (bound, kind.lower()), None)
            if count is not None:
                policy[bound + '_per_kind'][kind] = count
    for kind in expand_selection_types(getattr(args, 'excluded_antigen_types', None) or ()):
        policy['max_per_kind'][kind] = 0
    for name in ('variant_types', 'variant_contexts'):
        value = getattr(args, name, None)
        if value is not None:
            policy[name] = value
    return AntigenAllocation.resolve(policy).to_dict()


def manufacturability_config_from_args(
        args: argparse.Namespace, merged_config=None) -> "ManufacturabilityConfig":
    """Build a ``ManufacturabilityConfig`` from the YAML
    ``peptide.manufacturability`` sub-section. Always returns a real
    config (``ManufacturabilityConfig()`` defaults if the YAML
    sub-section is absent) — the ranker treats it as the source of
    peptide-synthesis difficulty thresholds + rules."""
    import msgspec
    from ..config.loader import extract_manufacturability_config_kwargs, load_vaxrank_config
    from ..manufacturability_config import ManufacturabilityConfig

    if merged_config is None:
        merged_config = load_vaxrank_config(args)
    kwargs = extract_manufacturability_config_kwargs(merged_config)
    return msgspec.convert(kwargs, ManufacturabilityConfig)
