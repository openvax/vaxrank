from typing import Optional, Union

import msgspec


class AffinityScoreConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    midpoint: Optional[float] = None
    width: Optional[float] = None
    cutoff: Optional[float] = None


class PercentileScoreConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    worst: Optional[float] = None


class DerivedFieldsConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    affinity_score: Optional[AffinityScoreConfigSchema] = None
    percentile_score: Optional[PercentileScoreConfigSchema] = None


class EpitopesFiltersConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    min_score: Optional[float] = None


class EpitopesScoringConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    mode: Optional[str] = None
    derived_fields: Optional[DerivedFieldsConfigSchema] = None


class EpitopesKeepConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    top_n_per_candidate: Optional[int] = None


class EpitopesConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    filters: Optional[EpitopesFiltersConfigSchema] = None
    scoring: Optional[EpitopesScoringConfigSchema] = None
    keep: Optional[EpitopesKeepConfigSchema] = None


class VaccinePeptideGenerationConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    lengths: Optional[Union[list[int], int]] = None
    padding_around_mutation: Optional[int] = None


class VaccinePeptideKeepConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    per_mutation: Optional[int] = None
    max_epitopes_per_candidate: Optional[int] = None
    score_fraction_of_best: Optional[float] = None


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
    generation: Optional[VaccinePeptideGenerationConfigSchema] = None
    keep: Optional[VaccinePeptideKeepConfigSchema] = None
    manufacturability: Optional[ManufacturabilityConfigSchema] = None


class VaxrankConfigSchema(
    msgspec.Struct,
    frozen=True,
    kw_only=True,
    forbid_unknown_fields=True,
):
    epitopes: Optional[EpitopesConfigSchema] = None
    vaccine_peptides: Optional[VaccinePeptidesConfigSchema] = None
