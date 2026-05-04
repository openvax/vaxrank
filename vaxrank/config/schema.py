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
    preferred_length: Optional[int] = None
    min_length: Optional[int] = None
    max_length: Optional[int] = None
    padding_around_mutation: Optional[int] = None
    per_mutation: Optional[int] = None
    max_epitopes_per_candidate: Optional[int] = None
    score_fraction_of_best: Optional[float] = None
    combined_score_mode: Optional[str] = None
    ranking_rules: Optional[list[str]] = None
    require_mutant_epitopes_in_variant: Optional[bool] = None
    manufacturability: Optional[ManufacturabilityConfigSchema] = None


class _SharedConstructConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    """Construct-assembly knobs shared between peptide and mRNA.

    Modality-specific subsections (``peptide``, ``mrna``) override
    these. Knobs whose defaults differ between modalities (linker,
    ``antigens_per_construct``, ``max_constructs``) live in the
    per-modality subsections, not here.
    """
    min_antigen_length_aa: Optional[int] = None
    max_antigen_length_aa: Optional[int] = None
    candidates_per_slot: Optional[int] = None
    antigen_content: Optional[str] = None
    epitopes_per_antigen: Optional[int] = None


class _PeptideConstructConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    """Peptide-specific construct-assembly knobs. Any knob also in
    ``shared`` overrides the shared value when set here."""
    # Shared knobs (override of vaccine_constructs.shared)
    min_antigen_length_aa: Optional[int] = None
    max_antigen_length_aa: Optional[int] = None
    candidates_per_slot: Optional[int] = None
    antigen_content: Optional[str] = None
    epitopes_per_antigen: Optional[int] = None
    # Peptide-only
    mode: Optional[str] = None
    linker: Optional[str] = None
    antigens_per_construct: Optional[int] = None
    max_constructs: Optional[int] = None
    n_terminal_acetyl: Optional[bool] = None
    c_terminal_amide: Optional[bool] = None
    scale_mg: Optional[float] = None
    purity_percent: Optional[float] = None
    counterion: Optional[str] = None


class _MrnaConstructConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    """mRNA-specific construct-assembly knobs. Any knob also in
    ``shared`` overrides the shared value when set here."""
    # Shared knobs (override of vaccine_constructs.shared)
    min_antigen_length_aa: Optional[int] = None
    max_antigen_length_aa: Optional[int] = None
    candidates_per_slot: Optional[int] = None
    antigen_content: Optional[str] = None
    epitopes_per_antigen: Optional[int] = None
    # mRNA-only
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


class VaccineConstructsConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    """Construct-assembly section of the YAML config.

    Layout:

        vaccine_constructs:
          shared:    # cross-modality defaults (length bounds, antigen_content, …)
          peptide:   # peptide-specific overrides (mode, linker, n-term acetyl, …)
          mrna:      # mRNA-specific overrides (signal_peptide, UTRs, polyA, …)

    Resolution order (highest precedence first):

        explicit CLI flag  >  per-modality config  >  shared config  >  built-in default

    A user can drive both modalities from one config file, and only
    needs to repeat a knob in a modality subsection when they want
    that modality to differ from the shared default.
    """
    shared: Optional[_SharedConstructConfigSchema] = None
    peptide: Optional[_PeptideConstructConfigSchema] = None
    mrna: Optional[_MrnaConstructConfigSchema] = None


class VaxrankConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    epitopes: Optional[EpitopesConfigSchema] = None
    vaccine_peptides: Optional[VaccinePeptidesConfigSchema] = None
    manufacturability: Optional[ManufacturabilityConfigSchema] = None
    vaccine_constructs: Optional[VaccineConstructsConfigSchema] = None
