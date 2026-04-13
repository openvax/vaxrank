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
    manufacturability: Optional[ManufacturabilityConfigSchema] = None


class VaxrankConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    epitopes: Optional[EpitopesConfigSchema] = None
    vaccine_peptides: Optional[VaccinePeptidesConfigSchema] = None
