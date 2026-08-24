# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""Derive safety-critical normal-tissue proteins from oncoref HPA facts.

HPA v23 defines each cell-type ``Level`` from staining intensity and fraction
of stained cells, and assigns the gene-level IHC profile one of Enhanced,
Supported, Approved, or Uncertain reliability.  The default here—Medium/High
and any non-Uncertain reliability—is a conservative vaxrank safety policy over
those observations, not an HPA claim that lower expression is biologically
irrelevant or proof that a retained protein will be presented by MHC.

Primary definitions:
https://v23.proteinatlas.org/about/assays%2Bannotation
https://v23.proteinatlas.org/about/help
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field
from typing import Optional

import pandas as pd
from serializable import DataclassSerializable

from .identifiers import normalize_ensembl_gene_id

DEFAULT_HPA_LEVELS = ("High", "Medium")
DEFAULT_HPA_RELIABILITIES = ("Approved", "Enhanced", "Supported")
HPA_NORMAL_TISSUE_SOURCE = "hpa_normal_tissue"
PRAME_GENE_ID = "ENSG00000185686"


class TissueRiskDerivationError(RuntimeError):
    """Raised when HPA safety evidence is missing or inconsistent."""


def _text(value) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


@dataclass(frozen=True)
class SafetyTissueGroup(DataclassSerializable):
    """One oncoref-owned safety group resolved to exact HPA tissue names."""

    name: str
    tissues: tuple[str, ...]

    def __post_init__(self):
        if not self.name or not self.tissues:
            raise ValueError("Safety tissue group name and tissues are required")
        object.__setattr__(self, "tissues", tuple(sorted(set(self.tissues))))


@dataclass(frozen=True)
class TissueRiskPolicy(DataclassSerializable):
    """Resolved HPA selection rule; no tissue vocabulary is copied by vaxrank."""

    tissue_groups: tuple[SafetyTissueGroup, ...]
    allowed_levels: tuple[str, ...] = DEFAULT_HPA_LEVELS
    allowed_reliabilities: tuple[str, ...] = DEFAULT_HPA_RELIABILITIES

    def __post_init__(self):
        groups = tuple(sorted(self.tissue_groups, key=lambda group: group.name))
        names = [group.name for group in groups]
        if not groups or len(names) != len(set(names)):
            raise ValueError("Tissue-risk policy requires unique safety groups")
        levels = tuple(sorted(set(self.allowed_levels)))
        reliabilities = tuple(sorted(set(self.allowed_reliabilities)))
        if not levels or not reliabilities:
            raise ValueError("Tissue-risk levels and reliabilities cannot be empty")
        object.__setattr__(self, "tissue_groups", groups)
        object.__setattr__(self, "allowed_levels", levels)
        object.__setattr__(self, "allowed_reliabilities", reliabilities)


@dataclass(frozen=True)
class HPAReferenceProvenance(DataclassSerializable):
    """Version and content identity for the exact HPA IHC artifact used."""

    oncoref_version: str
    oncoref_data_version: str
    oncoref_source_matrix_version: str
    source_name: str
    hpa_version: str
    source_path: str
    source_bytes: int
    source_sha256: str
    verification_status: str

    def __post_init__(self):
        if len(self.source_sha256) != 64:
            raise ValueError("HPA provenance requires a SHA-256 digest")
        if self.source_bytes <= 0:
            raise ValueError("HPA provenance requires a non-empty source artifact")
        if self.verification_status not in {
            "verified", "checksum_mismatch", "manifest_unavailable"
        }:
            raise ValueError("Unknown HPA verification status")


@dataclass(frozen=True)
class TissueRiskEvidence(DataclassSerializable):
    """One qualifying HPA tissue/cell-type observation."""

    tissue_group: str
    hpa_tissue: str
    cell_type: str
    level: str
    reliability: str


@dataclass(frozen=True)
class TissueRiskProtein(DataclassSerializable):
    """A non-CTA source gene with qualifying normal-tissue IHC evidence."""

    gene_id: str
    gene_name: str
    evidence: tuple[TissueRiskEvidence, ...]

    def __post_init__(self):
        normalized = normalize_ensembl_gene_id(self.gene_id)
        if not normalized or not self.evidence:
            raise ValueError("Tissue-risk protein requires a gene ID and evidence")
        object.__setattr__(self, "gene_id", normalized)
        object.__setattr__(
            self,
            "evidence",
            tuple(sorted(set(self.evidence), key=lambda item: (
                item.tissue_group,
                item.hpa_tissue,
                item.cell_type,
                item.level,
                item.reliability,
            ))),
        )


@dataclass(frozen=True)
class HPAEvidenceCoverage(DataclassSerializable):
    """Counts that distinguish non-qualification from missing evidence."""

    total_rows: int
    selected_tissue_rows: int
    missing_gene_id_rows: int
    missing_level_rows: int
    missing_reliability_rows: int
    nonqualifying_level_rows: int
    nonqualifying_reliability_rows: int
    qualifying_rows: int
    cta_source_rows_excluded: int


@dataclass(frozen=True)
class TissueRiskAssessment(DataclassSerializable):
    """Immutable HPA tissue-risk gene set with complete policy provenance."""

    policy: TissueRiskPolicy
    proteins: tuple[TissueRiskProtein, ...]
    coverage: HPAEvidenceCoverage
    hpa_provenance: HPAReferenceProvenance
    cta_unfiltered_gene_count: int
    cta_unfiltered_gene_ids_sha256: str
    excluded_cta_source_gene_ids: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self):
        proteins = tuple(sorted(self.proteins, key=lambda protein: protein.gene_id))
        ids = [protein.gene_id for protein in proteins]
        if len(ids) != len(set(ids)):
            raise ValueError("Tissue-risk assessment has duplicate source genes")
        if len(self.cta_unfiltered_gene_ids_sha256) != 64:
            raise ValueError("CTA universe provenance requires a SHA-256 digest")
        if self.cta_unfiltered_gene_count <= 0:
            raise ValueError("CTA universe provenance cannot be empty")
        object.__setattr__(self, "proteins", proteins)
        object.__setattr__(
            self,
            "excluded_cta_source_gene_ids",
            tuple(sorted(set(self.excluded_cta_source_gene_ids))),
        )

    def to_report_dict(self) -> dict:
        return json.loads(json.dumps(asdict(self)))

    def evidence_rows(self) -> list[dict]:
        rows = []
        for protein in self.proteins:
            for evidence in protein.evidence:
                rows.append({
                    "gene_id": protein.gene_id,
                    "gene_name": protein.gene_name,
                    "tissue_group": evidence.tissue_group,
                    "hpa_tissue": evidence.hpa_tissue,
                    "cell_type": evidence.cell_type,
                    "level": evidence.level,
                    "reliability": evidence.reliability,
                    "hpa_version": self.hpa_provenance.hpa_version,
                    "hpa_source_sha256": self.hpa_provenance.source_sha256,
                })
        return rows


def resolve_tissue_risk_policy(
    tissue_group_names: Optional[tuple[str, ...]] = None,
    *,
    allowed_levels: tuple[str, ...] = DEFAULT_HPA_LEVELS,
    allowed_reliabilities: tuple[str, ...] = DEFAULT_HPA_RELIABILITIES,
) -> TissueRiskPolicy:
    """Resolve selected group names through oncoref's live vocabulary."""
    from oncoref.cta_tissues import SAFETY_TISSUE_GROUPS

    if tissue_group_names is None:
        tissue_group_names = tuple(sorted(SAFETY_TISSUE_GROUPS))
    unknown = set(tissue_group_names) - set(SAFETY_TISSUE_GROUPS)
    if unknown:
        raise ValueError(f"Unknown oncoref safety tissue groups: {sorted(unknown)}")
    return TissueRiskPolicy(
        tissue_groups=tuple(
            SafetyTissueGroup(
                name=name,
                tissues=tuple(SAFETY_TISSUE_GROUPS[name]),
            )
            for name in tissue_group_names
        ),
        allowed_levels=allowed_levels,
        allowed_reliabilities=allowed_reliabilities,
    )


def build_hpa_reference_provenance() -> HPAReferenceProvenance:
    """Hash and verify oncoref's pinned HPA normal-tissue source artifact."""
    import oncoref
    from oncoref import reference_data
    from oncoref.version import DATA_VERSION, SOURCE_MATRIX_VERSION

    hpa_version = reference_data.DEFAULT_HPA_VERSION
    path = reference_data.ensure(HPA_NORMAL_TISSUE_SOURCE, hpa_version)
    source_digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            source_digest.update(chunk)
    try:
        verified = reference_data.verify(HPA_NORMAL_TISSUE_SOURCE, hpa_version)
        verification_status = "verified" if verified else "checksum_mismatch"
    except reference_data.ReferenceDataError:
        verification_status = "manifest_unavailable"
    return HPAReferenceProvenance(
        oncoref_version=oncoref.__version__,
        oncoref_data_version=DATA_VERSION,
        oncoref_source_matrix_version=SOURCE_MATRIX_VERSION,
        source_name=HPA_NORMAL_TISSUE_SOURCE,
        hpa_version=hpa_version,
        source_path=str(path),
        source_bytes=path.stat().st_size,
        source_sha256=source_digest.hexdigest(),
        verification_status=verification_status,
    )


def derive_tissue_risk_proteins(
    *,
    hpa_frame: Optional[pd.DataFrame] = None,
    provenance: Optional[HPAReferenceProvenance] = None,
    tissue_group_names: Optional[tuple[str, ...]] = None,
    allowed_levels: tuple[str, ...] = DEFAULT_HPA_LEVELS,
    allowed_reliabilities: tuple[str, ...] = DEFAULT_HPA_RELIABILITIES,
) -> TissueRiskAssessment:
    """Apply vaxrank's IHC rule to oncoref-owned HPA and CTA facts.

    A source gene qualifies when any selected cell-type row is Medium or High
    with Enhanced, Supported, or Approved reliability under the default policy.
    Genes in oncoref's full unfiltered CTA candidate universe are then removed
    by source identity only.
    """
    from oncoref.cta import cta_unfiltered_gene_ids
    from oncoref.hpa import hpa_normal_tissue

    policy = resolve_tissue_risk_policy(
        tissue_group_names,
        allowed_levels=allowed_levels,
        allowed_reliabilities=allowed_reliabilities,
    )
    if hpa_frame is None:
        hpa_frame = hpa_normal_tissue()
        if provenance is None:
            provenance = build_hpa_reference_provenance()
    elif provenance is None:
        raise ValueError("Injected HPA data requires explicit provenance")

    required_columns = {
        "Gene", "Gene name", "Tissue", "Cell type", "Level", "Reliability"
    }
    missing_columns = required_columns - set(hpa_frame.columns)
    if missing_columns:
        raise TissueRiskDerivationError(
            f"HPA normal-tissue data is missing columns: {sorted(missing_columns)}"
        )

    cta_gene_ids = frozenset(
        normalize_ensembl_gene_id(gene_id)
        for gene_id in cta_unfiltered_gene_ids()
    )
    tissue_to_groups: dict[str, list[str]] = {}
    for group in policy.tissue_groups:
        for tissue in group.tissues:
            tissue_to_groups.setdefault(tissue.strip().casefold(), []).append(
                group.name
            )

    counts = {
        "total_rows": len(hpa_frame),
        "selected_tissue_rows": 0,
        "missing_gene_id_rows": 0,
        "missing_level_rows": 0,
        "missing_reliability_rows": 0,
        "nonqualifying_level_rows": 0,
        "nonqualifying_reliability_rows": 0,
        "qualifying_rows": 0,
        "cta_source_rows_excluded": 0,
    }
    by_gene: dict[str, dict] = {}
    excluded_cta_gene_ids = set()
    for _, row in hpa_frame.iterrows():
        tissue = _text(row["Tissue"])
        groups = tissue_to_groups.get(tissue.casefold(), ())
        if not groups:
            continue
        counts["selected_tissue_rows"] += 1
        gene_id = (
            "" if pd.isna(row["Gene"])
            else normalize_ensembl_gene_id(row["Gene"])
        )
        level = _text(row["Level"])
        reliability = _text(row["Reliability"])
        if not gene_id:
            counts["missing_gene_id_rows"] += 1
            continue
        if not level:
            counts["missing_level_rows"] += 1
            continue
        if not reliability:
            counts["missing_reliability_rows"] += 1
            continue
        if level not in policy.allowed_levels:
            counts["nonqualifying_level_rows"] += 1
            continue
        if reliability not in policy.allowed_reliabilities:
            counts["nonqualifying_reliability_rows"] += 1
            continue
        counts["qualifying_rows"] += 1
        if gene_id in cta_gene_ids:
            counts["cta_source_rows_excluded"] += 1
            excluded_cta_gene_ids.add(gene_id)
            continue

        gene_name = _text(row["Gene name"])
        slot = by_gene.setdefault(gene_id, {"names": set(), "evidence": set()})
        if gene_name:
            slot["names"].add(gene_name)
        for group in groups:
            slot["evidence"].add(TissueRiskEvidence(
                tissue_group=group,
                hpa_tissue=tissue,
                cell_type=_text(row["Cell type"]),
                level=level,
                reliability=reliability,
            ))

    proteins = []
    for gene_id, slot in by_gene.items():
        if len(slot["names"]) > 1:
            raise TissueRiskDerivationError(
                f"HPA gene {gene_id} has conflicting symbols: "
                f"{sorted(slot['names'])}"
            )
        proteins.append(TissueRiskProtein(
            gene_id=gene_id,
            gene_name=next(iter(slot["names"]), ""),
            evidence=tuple(slot["evidence"]),
        ))

    return TissueRiskAssessment(
        policy=policy,
        proteins=tuple(proteins),
        coverage=HPAEvidenceCoverage(**counts),
        hpa_provenance=provenance,
        cta_unfiltered_gene_count=len(cta_gene_ids),
        cta_unfiltered_gene_ids_sha256=hashlib.sha256(
            "\n".join(sorted(cta_gene_ids)).encode("utf-8")
        ).hexdigest(),
        excluded_cta_source_gene_ids=tuple(excluded_cta_gene_ids),
    )
