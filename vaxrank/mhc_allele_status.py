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

"""Typed cross-species MHC allele status and expression contract
(openvax/vaxrank#921).

Carries both raw quantitative assay measurements (DNA copy-number LOH,
bulk RNA allele-specific expression + tumor purity, and single-cell RNA
tumor-cluster vs immune-control nTPM + resolvable UMIs) and a resolved
categorical presentation status (``retained`` | ``lost`` | ``silenced`` |
``unresolvable`` | ``unknown``).

All allele matching uses ``mhcgnomes`` so human (HLA), canine (DLA),
murine (H-2), swine (SLA), and other vertebrate MHC genotypes work
uniformly without string-prefix hacks.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, field, replace
from numbers import Real
from pathlib import Path
from typing import Any

import mhcgnomes


VALID_ALLELE_STATUSES = ("retained", "lost", "silenced", "unresolvable", "unknown")
VALID_WEIGHT_MODES = ("binary", "continuous_ratio")


def normalize_mhc_allele_name(allele: str) -> str:
    """Return canonical ``mhcgnomes`` allele name across species.

    Falls back to stripped text only when an upstream table carries a
    non-standard label that ``mhcgnomes`` cannot parse (raising on empty).
    """
    text = str(allele or "").strip()
    if not text:
        raise ValueError("MHC allele name must be non-empty")
    try:
        parsed = mhcgnomes.parse(text, raise_on_error=True)
        return parsed.to_string()
    except Exception:
        # Try explicit cross-species fallbacks only via mhcgnomes (never prefix hacks)
        parsed = mhcgnomes.parse(text, raise_on_error=False)
        if parsed is not None:
            return parsed.to_string()
        return text


def _optional_float(value: Any, field_name: str, *, min_val: float | None = None,
                    max_val: float | None = None) -> float | None:
    if value is None or value == "":
        return None
    if isinstance(value, bool):
        raise ValueError(f"{field_name} must be a number, not bool")
    if isinstance(value, str):
        stripped = value.strip().lower()
        if stripped in ("", "na", "n/a", "nan", "none", "null"):
            return None
        try:
            f = float(value)
        except ValueError as exc:
            raise ValueError(f"{field_name} must be numeric; got {value!r}") from exc
    elif isinstance(value, Real):
        f = float(value)
    else:
        raise ValueError(f"{field_name} must be numeric; got {type(value).__name__}")
    if not math.isfinite(f):
        raise ValueError(f"{field_name} must be finite; got {value!r}")
    if min_val is not None and f < min_val:
        raise ValueError(f"{field_name} must be >= {min_val}; got {f}")
    if max_val is not None and f > max_val:
        raise ValueError(f"{field_name} must be <= {max_val}; got {f}")
    return f


def _optional_int(value: Any, field_name: str, *, min_val: int | None = None) -> int | None:
    if value is None or value == "":
        return None
    if isinstance(value, bool):
        raise ValueError(f"{field_name} must be an integer, not bool")
    if isinstance(value, str):
        stripped = value.strip().lower()
        if stripped in ("", "na", "n/a", "nan", "none", "null"):
            return None
        try:
            i = int(value.strip())
        except ValueError as exc:
            raise ValueError(f"{field_name} must be an integer; got {value!r}") from exc
    elif isinstance(value, int):
        i = value
    else:
        raise ValueError(f"{field_name} must be an integer; got {type(value).__name__}")
    if min_val is not None and i < min_val:
        raise ValueError(f"{field_name} must be >= {min_val}; got {i}")
    return i


def _optional_bool(value: Any, field_name: str) -> bool | None:
    if value is None or value == "":
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        low = value.strip().lower()
        if low in ("", "na", "n/a", "none", "null"):
            return None
        if low in ("true", "1", "yes", "y"):
            return True
        if low in ("false", "0", "no", "n"):
            return False
    raise ValueError(f"{field_name} must be a boolean; got {value!r}")


@dataclass(frozen=True)
class MHCAlleleStatusPolicy:
    """Explicit per-trial policy governing MHC allele LOH / silencing resolution
    and ranking-credit adjustment.

    Configurable under ``mhc_allele_status:`` in a per-trial YAML policy.
    """
    mask_lost_alleles: bool = True
    mask_silenced_alleles: bool = True
    mask_unresolvable_alleles: bool = False
    expression_weight_mode: str = "binary"
    dna_loh_max_copy_number: float = 0.5
    min_scrna_resolvable_umis: int = 20
    silenced_tumor_to_immune_ratio: float = 0.1
    silenced_bulk_allele_ratio: float = 0.1
    min_purity_for_bulk_ase_call: float = 0.5

    def __post_init__(self):
        for flag in ("mask_lost_alleles", "mask_silenced_alleles", "mask_unresolvable_alleles"):
            if not isinstance(getattr(self, flag), bool):
                raise ValueError(f"mhc_allele_status.{flag} must be a boolean")
        if self.expression_weight_mode not in VALID_WEIGHT_MODES:
            raise ValueError(
                f"mhc_allele_status.expression_weight_mode must be one of "
                f"{VALID_WEIGHT_MODES}; got {self.expression_weight_mode!r}"
            )
        for name, lo, hi in (
            ("dna_loh_max_copy_number", 0.0, None),
            ("silenced_tumor_to_immune_ratio", 0.0, 1.0),
            ("silenced_bulk_allele_ratio", 0.0, 1.0),
            ("min_purity_for_bulk_ase_call", 0.0, 1.0),
        ):
            val = getattr(self, name)
            if isinstance(val, bool) or not isinstance(val, Real) or not math.isfinite(val):
                raise ValueError(f"mhc_allele_status.{name} must be a finite number")
            if float(val) < lo or (hi is not None and float(val) > hi):
                raise ValueError(f"mhc_allele_status.{name} out of bounds [{lo}, {hi}]")
        if (isinstance(self.min_scrna_resolvable_umis, bool)
                or not isinstance(self.min_scrna_resolvable_umis, int)
                or self.min_scrna_resolvable_umis < 1):
            raise ValueError("mhc_allele_status.min_scrna_resolvable_umis must be an integer >= 1")

    @classmethod
    def resolve(cls, value=None) -> MHCAlleleStatusPolicy:
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if not isinstance(value, dict):
            raise TypeError(
                f"mhc_allele_status policy must be a mapping or MHCAlleleStatusPolicy; "
                f"got {type(value).__name__}"
            )
        defaults = cls().to_dict()
        unknown = set(value) - set(defaults)
        if unknown:
            raise ValueError(
                "Unknown mhc_allele_status policy fields: %s" % ", ".join(sorted(unknown))
            )
        merged = {k: (defaults[k] if value.get(k) is None else value[k]) for k in defaults}
        return cls(**merged)

    def to_dict(self) -> dict[str, Any]:
        return {
            "mask_lost_alleles": self.mask_lost_alleles,
            "mask_silenced_alleles": self.mask_silenced_alleles,
            "mask_unresolvable_alleles": self.mask_unresolvable_alleles,
            "expression_weight_mode": self.expression_weight_mode,
            "dna_loh_max_copy_number": float(self.dna_loh_max_copy_number),
            "min_scrna_resolvable_umis": self.min_scrna_resolvable_umis,
            "silenced_tumor_to_immune_ratio": float(self.silenced_tumor_to_immune_ratio),
            "silenced_bulk_allele_ratio": float(self.silenced_bulk_allele_ratio),
            "min_purity_for_bulk_ase_call": float(self.min_purity_for_bulk_ase_call),
        }


@dataclass(frozen=True)
class MHCAlleleStatus:
    """Per-allele presentation status and multi-assay quantitative evidence."""
    allele: str
    status: str = "unknown"
    # DNA MHC LOH (e.g. LOHHLA / DASH / Polysolver)
    dna_copy_number: float | None = None
    dna_loh: bool | None = None
    dna_loh_pvalue: float | None = None
    # Bulk RNA-seq allele-specific expression (e.g. HLApers / T1K / arcasHLA)
    bulk_rna_tpm: float | None = None
    bulk_rna_allele_ratio: float | None = None
    tumor_purity: float | None = None
    bulk_purity_confounded: bool | None = None
    # Single-cell RNA-seq tumor-cluster vs immune-control allele expression
    scrna_tumor_cluster_ntpm: float | None = None
    scrna_immune_control_ntpm: float | None = None
    scrna_resolvable_umis: int | None = None
    scrna_chemistry: str = ""
    # Provenance
    assay_source: str = ""
    notes: str = ""

    def __post_init__(self):
        canonical = normalize_mhc_allele_name(self.allele)
        object.__setattr__(self, "allele", canonical)
        st = (str(self.status or "unknown").strip().lower()) or "unknown"
        if st not in VALID_ALLELE_STATUSES:
            raise ValueError(
                f"MHCAlleleStatus.status must be one of {VALID_ALLELE_STATUSES}; got {st!r}"
            )
        object.__setattr__(self, "status", st)
        object.__setattr__(
            self, "dna_copy_number",
            _optional_float(self.dna_copy_number, "dna_copy_number", min_val=0.0),
        )
        object.__setattr__(self, "dna_loh", _optional_bool(self.dna_loh, "dna_loh"))
        object.__setattr__(
            self, "dna_loh_pvalue",
            _optional_float(self.dna_loh_pvalue, "dna_loh_pvalue", min_val=0.0, max_val=1.0),
        )
        object.__setattr__(
            self, "bulk_rna_tpm",
            _optional_float(self.bulk_rna_tpm, "bulk_rna_tpm", min_val=0.0),
        )
        object.__setattr__(
            self, "bulk_rna_allele_ratio",
            _optional_float(self.bulk_rna_allele_ratio, "bulk_rna_allele_ratio", min_val=0.0, max_val=1.0),
        )
        object.__setattr__(
            self, "tumor_purity",
            _optional_float(self.tumor_purity, "tumor_purity", min_val=0.0, max_val=1.0),
        )
        object.__setattr__(
            self, "bulk_purity_confounded",
            _optional_bool(self.bulk_purity_confounded, "bulk_purity_confounded"),
        )
        object.__setattr__(
            self, "scrna_tumor_cluster_ntpm",
            _optional_float(self.scrna_tumor_cluster_ntpm, "scrna_tumor_cluster_ntpm", min_val=0.0),
        )
        object.__setattr__(
            self, "scrna_immune_control_ntpm",
            _optional_float(self.scrna_immune_control_ntpm, "scrna_immune_control_ntpm", min_val=0.0),
        )
        object.__setattr__(
            self, "scrna_resolvable_umis",
            _optional_int(self.scrna_resolvable_umis, "scrna_resolvable_umis", min_val=0),
        )
        object.__setattr__(self, "scrna_chemistry", str(self.scrna_chemistry or "").strip())
        object.__setattr__(self, "assay_source", str(self.assay_source or "").strip())
        object.__setattr__(self, "notes", str(self.notes or "").strip())

    def to_dict(self) -> dict[str, Any]:
        return {
            "allele": self.allele,
            "status": self.status,
            "dna_copy_number": self.dna_copy_number,
            "dna_loh": self.dna_loh,
            "dna_loh_pvalue": self.dna_loh_pvalue,
            "bulk_rna_tpm": self.bulk_rna_tpm,
            "bulk_rna_allele_ratio": self.bulk_rna_allele_ratio,
            "tumor_purity": self.tumor_purity,
            "bulk_purity_confounded": self.bulk_purity_confounded,
            "scrna_tumor_cluster_ntpm": self.scrna_tumor_cluster_ntpm,
            "scrna_immune_control_ntpm": self.scrna_immune_control_ntpm,
            "scrna_resolvable_umis": self.scrna_resolvable_umis,
            "scrna_chemistry": self.scrna_chemistry,
            "assay_source": self.assay_source,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class APMStatus:
    """Antigen presentation machinery (APM) gene-level status context."""
    b2m_status: str = ""
    tap1_status: str = ""
    tap2_status: str = ""
    nlrc5_status: str = ""
    jak1_status: str = ""
    jak2_status: str = ""
    notes: str = ""

    def any_recorded(self) -> bool:
        return any(
            bool(getattr(self, k))
            for k in ("b2m_status", "tap1_status", "tap2_status",
                      "nlrc5_status", "jak1_status", "jak2_status", "notes")
        )

    def disrupted_genes(self) -> list[str]:
        disrupted_tokens = ("lost", "loh", "lof", "mutated", "deleted", "silenced", "deficient", "truncating")
        out = []
        for label, attr in (
            ("B2M", "b2m_status"),
            ("TAP1", "tap1_status"),
            ("TAP2", "tap2_status"),
            ("NLRC5", "nlrc5_status"),
            ("JAK1", "jak1_status"),
            ("JAK2", "jak2_status"),
        ):
            val = (getattr(self, attr) or "").strip().lower()
            if val and any(tok in val for tok in disrupted_tokens):
                out.append(label)
        return out

    def to_dict(self) -> dict[str, Any]:
        return {
            "b2m_status": self.b2m_status,
            "tap1_status": self.tap1_status,
            "tap2_status": self.tap2_status,
            "nlrc5_status": self.nlrc5_status,
            "jak1_status": self.jak1_status,
            "jak2_status": self.jak2_status,
            "notes": self.notes,
            "disrupted_genes": self.disrupted_genes(),
        }


@dataclass(frozen=True)
class MHCAlleleStatusBundle:
    """Collection of per-allele MHC status records plus optional APM context."""
    alleles: tuple[MHCAlleleStatus, ...] = field(default_factory=tuple)
    apm: APMStatus | None = None
    source_path: str = ""

    def __post_init__(self):
        object.__setattr__(self, "alleles", tuple(self.alleles))

    def by_allele(self) -> dict[str, MHCAlleleStatus]:
        return {rec.allele: rec for rec in self.alleles}

    def is_empty(self) -> bool:
        return not self.alleles and (self.apm is None or not self.apm.any_recorded())

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": "vaxrank.mhc_allele_status.v1",
            "source_path": self.source_path,
            "alleles": [rec.to_dict() for rec in self.alleles],
            "apm": self.apm.to_dict() if self.apm and self.apm.any_recorded() else None,
        }


@dataclass(frozen=True)
class ResolvedMHCAlleleDecision:
    """Policy-evaluated presentation status and ranking weight for one allele."""
    record: MHCAlleleStatus
    effective_status: str
    ranking_weight: float
    masked: bool
    reason_codes: tuple[str, ...]
    caveats: tuple[str, ...]
    badge: str

    def to_dict(self) -> dict[str, Any]:
        return {
            **self.record.to_dict(),
            "effective_status": self.effective_status,
            "ranking_weight": round(float(self.ranking_weight), 6),
            "masked": self.masked,
            "reason_codes": list(self.reason_codes),
            "caveats": list(self.caveats),
            "badge": self.badge,
        }


_STATUS_BADGE = {
    "retained": "Retained",
    "lost": "MHC-LOH",
    "silenced": "MHC-Silenced",
    "unresolvable": "MHC-Unresolvable",
    "unknown": "Unknown",
}


def resolve_effective_status(
    record: MHCAlleleStatus,
    policy: MHCAlleleStatusPolicy | dict[str, Any] | None = None,
) -> ResolvedMHCAlleleDecision:
    """Derive effective presentation status, caveats, and ranking weight from
    multi-assay evidence under ``policy``.

    Never conflates ``unresolvable`` (e.g. 10x 3' scRNA with low UMIs or
    low-purity bulk RNA) with confirmed ``lost`` or ``silenced``.
    """
    pol = MHCAlleleStatusPolicy.resolve(policy)
    reasons: list[str] = []
    caveats: list[str] = []

    purity_confounded = bool(record.bulk_purity_confounded) or (
        record.tumor_purity is not None
        and record.tumor_purity < pol.min_purity_for_bulk_ase_call
    )
    if purity_confounded and (
        record.bulk_rna_tpm is not None or record.bulk_rna_allele_ratio is not None
    ):
        caveats.append("bulk_rna_purity_confounded")

    scrna_low_umis = (
        record.scrna_resolvable_umis is not None
        and record.scrna_resolvable_umis < pol.min_scrna_resolvable_umis
    )
    if scrna_low_umis:
        caveats.append("scrna_insufficient_resolvable_umis")
    if record.scrna_chemistry and "3p" in record.scrna_chemistry.lower() and scrna_low_umis:
        caveats.append("scrna_3p_exon2_3_coverage_gap")

    # Compute continuous expression ratio when available and unconfounded
    scrna_ratio: float | None = None
    if (not scrna_low_umis
            and record.scrna_tumor_cluster_ntpm is not None
            and record.scrna_immune_control_ntpm is not None
            and record.scrna_immune_control_ntpm > 0):
        scrna_ratio = record.scrna_tumor_cluster_ntpm / record.scrna_immune_control_ntpm

    effective = record.status
    if effective != "unknown":
        reasons.append(f"explicit_status:{effective}")
    else:
        # 1. DNA copy-number LOH (highest specificity for genomic loss)
        if record.dna_loh is True or (
            record.dna_copy_number is not None
            and record.dna_copy_number <= pol.dna_loh_max_copy_number
        ):
            effective = "lost"
            reasons.append("dna_loh_confirmed")
        # 2. Single-cell tumor-cluster expression (purity-independent)
        elif scrna_low_umis:
            effective = "unresolvable"
            reasons.append("scrna_insufficient_resolvable_umis")
        elif scrna_ratio is not None:
            if scrna_ratio <= pol.silenced_tumor_to_immune_ratio:
                effective = "silenced"
                reasons.append("scrna_tumor_cluster_silenced")
            else:
                effective = "retained"
                reasons.append("scrna_tumor_cluster_expressed")
        # 3. Bulk RNA allele-specific expression (subject to purity guard)
        elif record.bulk_rna_allele_ratio is not None:
            if record.bulk_rna_allele_ratio <= pol.silenced_bulk_allele_ratio:
                if purity_confounded:
                    effective = "unresolvable"
                    reasons.append("bulk_ase_low_but_purity_confounded")
                else:
                    effective = "silenced"
                    reasons.append("bulk_rna_ase_silenced")
            else:
                effective = "retained"
                reasons.append("bulk_rna_ase_expressed")
                if purity_confounded:
                    caveats.append("infiltrate_may_mask_tumor_loss")
        elif record.dna_loh is False:
            effective = "retained"
            reasons.append("dna_no_loh")
        else:
            effective = "unknown"
            reasons.append("no_decisive_assay_evidence")

    # Determine ranking weight
    masked = False
    if effective == "lost" and pol.mask_lost_alleles:
        weight = 0.0
        masked = True
    elif effective == "silenced" and pol.mask_silenced_alleles:
        weight = 0.0
        masked = True
    elif effective == "unresolvable" and pol.mask_unresolvable_alleles:
        weight = 0.0
        masked = True
    elif pol.expression_weight_mode == "continuous_ratio":
        if scrna_ratio is not None:
            weight = max(0.0, min(1.0, scrna_ratio))
        elif record.bulk_rna_allele_ratio is not None and not purity_confounded:
            # 0.5 is balanced heterozygous expression -> weight 1.0
            weight = max(0.0, min(1.0, 2.0 * record.bulk_rna_allele_ratio))
        else:
            weight = 1.0
    else:
        weight = 1.0

    return ResolvedMHCAlleleDecision(
        record=record,
        effective_status=effective,
        ranking_weight=weight,
        masked=masked,
        reason_codes=tuple(dict.fromkeys(reasons)),
        caveats=tuple(dict.fromkeys(caveats)),
        badge=_STATUS_BADGE[effective],
    )


def _record_from_mapping(row: dict[str, Any]) -> MHCAlleleStatus:
    allowed = {
        "allele", "status", "dna_copy_number", "dna_loh", "dna_loh_pvalue",
        "bulk_rna_tpm", "bulk_rna_allele_ratio", "tumor_purity",
        "bulk_purity_confounded", "scrna_tumor_cluster_ntpm",
        "scrna_immune_control_ntpm", "scrna_resolvable_umis",
        "scrna_chemistry", "assay_source", "notes",
    }
    kwargs = {k: v for k, v in row.items() if k in allowed}
    if "allele" not in kwargs or not str(kwargs["allele"] or "").strip():
        raise ValueError(f"MHC allele status record missing required 'allele' field: {row!r}")
    return MHCAlleleStatus(**kwargs)


def load_mhc_allele_status(path: str | Path) -> MHCAlleleStatusBundle:
    """Load an ``MHCAlleleStatusBundle`` from JSON or TSV/CSV."""
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(f"MHC allele status file not found: {p}")

    suffix = p.suffix.lower()
    if suffix == ".json":
        with open(p, encoding="utf-8") as handle:
            data = json.load(handle)
        if isinstance(data, list):
            alleles = [_record_from_mapping(item) for item in data]
            return MHCAlleleStatusBundle(alleles=tuple(alleles), source_path=str(p))
        if not isinstance(data, dict):
            raise ValueError(f"Expected JSON object or list in {p}")
        raw_alleles = data.get("alleles")
        if raw_alleles is None:
            raw_alleles = []
        if not isinstance(raw_alleles, list):
            raise ValueError(f"'alleles' in {p} must be a list")
        alleles = [_record_from_mapping(item) for item in raw_alleles]
        apm = None
        if isinstance(data.get("apm"), dict):
            apm_fields = {
                k: str(v or "").strip()
                for k, v in data["apm"].items()
                if k in ("b2m_status", "tap1_status", "tap2_status",
                         "nlrc5_status", "jak1_status", "jak2_status", "notes")
            }
            apm = APMStatus(**apm_fields)
        return MHCAlleleStatusBundle(alleles=tuple(alleles), apm=apm, source_path=str(p))

    delimiter = "\t" if suffix in (".tsv", ".txt") else ","
    alleles = []
    with open(p, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            if not row or not any(str(v or "").strip() for v in row.values()):
                continue
            alleles.append(_record_from_mapping(row))
    return MHCAlleleStatusBundle(alleles=tuple(alleles), source_path=str(p))


def bundle_from_loh_alleles(alleles_input: str | list[str] | tuple[str, ...]) -> MHCAlleleStatusBundle:
    """Create an ``MHCAlleleStatusBundle`` marking listed alleles as ``lost``
    (for the lightweight ``--mhc-loh-alleles`` CLI flag)."""
    if isinstance(alleles_input, str):
        tokens = [t.strip() for t in alleles_input.replace(" ", ",").split(",") if t.strip()]
    else:
        tokens = [str(t).strip() for t in alleles_input if str(t).strip()]
    records = tuple(
        MHCAlleleStatus(
            allele=tok,
            status="lost",
            dna_loh=True,
            assay_source="cli_mhc_loh_alleles",
        )
        for tok in tokens
    )
    return MHCAlleleStatusBundle(alleles=records, source_path="cli:--mhc-loh-alleles")


def merge_mhc_allele_status_bundles(
    base: MHCAlleleStatusBundle | None,
    override: MHCAlleleStatusBundle | None,
) -> MHCAlleleStatusBundle | None:
    """Merge two bundles (e.g. ``--mhc-allele-status`` file + ``--mhc-loh-alleles``).
    Entries in ``override`` take precedence per canonical ``mhcgnomes`` allele."""
    if base is None and override is None:
        return None
    if base is None:
        return override
    if override is None:
        return base
    merged: dict[str, MHCAlleleStatus] = {r.allele: r for r in base.alleles}
    for r in override.alleles:
        merged[r.allele] = r
    apm = override.apm if (override.apm and override.apm.any_recorded()) else base.apm
    src = "; ".join(s for s in (base.source_path, override.source_path) if s)
    return MHCAlleleStatusBundle(alleles=tuple(merged.values()), apm=apm, source_path=src)


def apply_mhc_allele_status_to_ranked(
    ranked_variants_with_vaccine_peptides,
    bundle: MHCAlleleStatusBundle,
    policy: MHCAlleleStatusPolicy | dict[str, Any] | None = None,
) -> tuple[list, dict[str, Any]]:
    """Apply policy-resolved MHC allele status to ranked vaccine peptides.

    Zeroes (or continuously scales under ``continuous_ratio``) positive ranking
    credit on ``CandidateEpitope.per_allele_scores`` for affected alleles while
    preserving all raw ``CandidateEpitope.predictions`` intact. Recomputes
    ``VaccinePeptide`` scores and re-sorts both per-source windows and sources.

    Returns ``(updated_ranked, audit_dict)``.
    """
    pol = MHCAlleleStatusPolicy.resolve(policy)
    decisions_by_allele = {
        rec.allele: resolve_effective_status(rec, pol)
        for rec in bundle.alleles
    }

    affected_epitopes = 0
    total_credit_removed = 0.0
    allele_credit_removed: dict[str, float] = {}
    affected_sources: list[str] = []

    updated_ranked = []
    for source, vaccine_peptides in ranked_variants_with_vaccine_peptides:
        new_vps = []
        source_touched = False
        for vp in vaccine_peptides:
            vp_changed = False
            new_epitopes = []
            for ep in vp.epitopes:
                if not ep.per_allele_scores:
                    new_epitopes.append(ep)
                    continue
                ep_touched = False
                new_scores = dict(ep.per_allele_scores)
                for raw_allele, score in ep.per_allele_scores.items():
                    canonical = normalize_mhc_allele_name(raw_allele)
                    decision = decisions_by_allele.get(canonical)
                    if decision is None or decision.ranking_weight == 1.0:
                        continue
                    adjusted = float(score) * float(decision.ranking_weight)
                    delta = float(score) - adjusted
                    if abs(delta) > 1e-12:
                        new_scores[raw_allele] = adjusted
                        total_credit_removed += delta
                        allele_credit_removed[canonical] = (
                            allele_credit_removed.get(canonical, 0.0) + delta
                        )
                        ep_touched = True
                if ep_touched:
                    new_epitopes.append(replace(ep, per_allele_scores=new_scores))
                    affected_epitopes += 1
                    vp_changed = True
                    source_touched = True
                else:
                    new_epitopes.append(ep)

            if vp_changed:
                new_vps.append(replace(vp, epitopes=new_epitopes))
            else:
                new_vps.append(vp)

        if source_touched:
            new_vps.sort(key=lambda p: p.combined_score, reverse=True)
            desc = getattr(source, "short_description", None) or str(source)
            affected_sources.append(desc)
        updated_ranked.append((source, new_vps))

    # Re-sort sources by top candidate combined_score (descending), preserving
    # relative order on exact ties.
    updated_ranked.sort(
        key=lambda item: (item[1][0].combined_score if item[1] else float("-inf")),
        reverse=True,
    )

    masked_alleles = [
        allele for allele, dec in decisions_by_allele.items() if dec.masked
    ]
    unresolvable_alleles = [
        allele for allele, dec in decisions_by_allele.items()
        if dec.effective_status == "unresolvable"
    ]

    audit = {
        "schema": "vaxrank.mhc_allele_status_audit.v1",
        "source_path": bundle.source_path,
        "policy": pol.to_dict(),
        "alleles": [
            {
                **dec.to_dict(),
                "credit_removed": round(allele_credit_removed.get(allele, 0.0), 6),
            }
            for allele, dec in decisions_by_allele.items()
        ],
        "masked_alleles": masked_alleles,
        "unresolvable_alleles": unresolvable_alleles,
        "apm": bundle.apm.to_dict() if bundle.apm and bundle.apm.any_recorded() else None,
        "affected_epitope_count": affected_epitopes,
        "affected_sources": affected_sources,
        "total_credit_removed": round(total_credit_removed, 6),
    }
    return updated_ranked, audit
