"""msgspec schemas for the vaxrank YAML config.

Top-level shape (post-2.19):

    epitopes:               # filter / score MHC ligand predictions
    vaccine_peptides:       # rank + antigen-window selection
                            # (modality-agnostic; absorbs the
                            # antigen-design knobs that used to
                            # live at the ``vaccine_constructs:``
                            # top level — antigen_content,
                            # epitopes_per_antigen,
                            # candidates_per_slot,
                            # min/max_antigen_length_aa)
    peptide:                # peptide-vaccine assembly +
                            # ``manufacturability`` sub-section
    mrna:                   # mRNA-vaccine assembly
    # future: dna:, ...

The pre-2.19 layout (``vaccine_constructs:`` wrapper, top-level
``manufacturability:``, ``vaccine_peptides.manufacturability:``) is
still accepted by the loader for one deprecation cycle — see
:func:`vaxrank.config.loader._migrate_legacy_layout`. The schema
itself is the new shape.
"""

from typing import Optional

import msgspec


class EpitopesConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    min_score: Optional[float] = None
    scoring_mode: Optional[str] = None
    logistic_midpoint: Optional[float] = None
    logistic_width: Optional[float] = None
    affinity_cutoff: Optional[float] = None
    percentile_rank_cutoff: Optional[float] = None
    top_epitopes_per_candidate: Optional[int] = None
    filter_expr: Optional[str] = None
    score_expr: Optional[str] = None
    default_methods: Optional[dict[str, str]] = None


class ManufacturabilityConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    """Peptide-synthesis difficulty rules + thresholds.

    Lives under ``peptide.manufacturability:`` in the YAML (was
    top-level ``manufacturability:`` pre-2.19; the loader still
    accepts the old location with a deprecation warning).
    """
    max_c_terminal_hydropathy: Optional[float] = None
    min_kmer_hydropathy: Optional[float] = None
    max_kmer_hydropathy_low_priority: Optional[float] = None
    max_kmer_hydropathy_high_priority: Optional[float] = None
    rules: Optional[list[str]] = None


class VaccinePeptidesConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    """Modality-agnostic ranking + antigen-window selection.

    Absorbs the antigen-design knobs (``antigen_content``,
    ``epitopes_per_antigen``, ``candidates_per_slot``,
    ``min_antigen_length_aa`` / ``max_antigen_length_aa``) that
    lived at the ``vaccine_constructs:`` top level pre-2.19 — those
    are antigen *design*, consumed by every modality's writer, so
    they belong here next to the rest of the antigen-design knobs
    rather than bolted onto a "constructs" wrapper.

    Section name kept as ``vaccine_peptides`` for back-compat with
    pre-mRNA configs even though the contents are modality-agnostic.
    """
    preferred_length: Optional[int] = None
    min_length: Optional[int] = None
    max_length: Optional[int] = None
    padding_around_mutation: Optional[int] = None
    per_mutation: Optional[int] = None
    max_epitopes_per_candidate: Optional[int] = None
    score_fraction_of_best: Optional[float] = None
    combined_score_expr: Optional[str] = None
    ranking_rules: Optional[list[str]] = None
    require_target_epitopes_in_variant: Optional[bool] = None
    # Antigen-design knobs hoisted from ``vaccine_constructs:`` top
    # level (post-2.19).
    antigen_content: Optional[str] = None
    epitopes_per_antigen: Optional[int] = None
    candidates_per_slot: Optional[int] = None
    min_antigen_length_aa: Optional[int] = None
    max_antigen_length_aa: Optional[int] = None


class PeptideConstructConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    """Top-level ``peptide:`` section — peptide-vaccine assembly +
    peptide-only manufacturability.

    Pre-2.19 this lived under ``vaccine_constructs.peptide:``; the
    wrapper got dropped because the cross-modality knobs that used
    to justify it (``antigen_content`` etc.) actually belong with
    the antigen-design knobs in ``vaccine_peptides:``.
    """
    mode: Optional[str] = None
    linker: Optional[str] = None
    antigens_per_construct: Optional[int] = None
    max_constructs: Optional[int] = None
    n_terminal_acetyl: Optional[bool] = None
    c_terminal_amide: Optional[bool] = None
    scale_mg: Optional[float] = None
    purity_percent: Optional[float] = None
    counterion: Optional[str] = None
    manufacturability: Optional[ManufacturabilityConfigSchema] = None


class MrnaConstructConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    """Top-level ``mrna:`` section — mRNA-vaccine assembly.

    No manufacturability sub-section: solid-phase peptide synthesis
    rules don't apply to in-vivo-translated mRNA. If mRNA grows its
    own assembly-difficulty rules (codon-optimization friction, IVT
    yield, …), add them as a parallel sub-section under a different
    name.
    """
    linker: Optional[str] = None
    antigens_per_construct: Optional[int] = None
    max_constructs: Optional[int] = None
    max_length_nt: Optional[int] = None
    signal_peptide: Optional[str] = None
    include_mitd: Optional[bool] = None
    mitd: Optional[str] = None
    utr_5p: Optional[str] = None
    utr_3p: Optional[str] = None
    poly_a_length: Optional[int] = None
    poly_a_segmented: Optional[bool] = None
    poly_a_first_segment: Optional[int] = None
    poly_a_segment_linker: Optional[str] = None
    codon_species: Optional[str] = None
    codon_method: Optional[str] = None
    optimize_linkers: Optional[bool] = None
    junction_candidates: Optional[str] = None
    junction_rank_strong: Optional[float] = None
    junction_rank_mild: Optional[float] = None


class VaxrankConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    epitopes: Optional[EpitopesConfigSchema] = None
    vaccine_peptides: Optional[VaccinePeptidesConfigSchema] = None
    peptide: Optional[PeptideConstructConfigSchema] = None
    mrna: Optional[MrnaConstructConfigSchema] = None
