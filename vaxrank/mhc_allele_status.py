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
import re
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
    """Per-allele presentation status and multi-assay quantitative evidence.

    Explicitly separates permanent DNA structural/mutational evidence
    (``dna_copy_number``, ``dna_loh``, ``dna_somatic_variants``) from
    transcriptional RNA expression evidence (bulk ASE and scRNA tumor-cluster
    vs immune-control nTPM).
    """
    allele: str
    status: str = "unknown"
    # DNA MHC LOH & somatic mutations (e.g. LOHHLA / DASH / LILAC / POLYSOLVER)
    dna_copy_number: float | None = None
    dna_loh: bool | None = None
    dna_loh_pvalue: float | None = None
    dna_somatic_variants: tuple[dict[str, Any], ...] = field(default_factory=tuple)
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
        raw_vars = self.dna_somatic_variants or ()
        if isinstance(raw_vars, str):
            raw_vars = json.loads(raw_vars) if raw_vars.strip() else ()
        object.__setattr__(
            self, "dna_somatic_variants",
            tuple(dict(v) for v in raw_vars if isinstance(v, dict)),
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
            "dna_somatic_variants": list(self.dna_somatic_variants),
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


_DISRUPTED_WORDS = frozenset({
    "lost", "loss", "loh", "lof", "mutated", "deleted", "deletion",
    "silenced", "deficient", "truncating", "frameshift", "nonsense", "homdel",
})
_NEGATION_WORDS = frozenset({
    "no", "not", "non", "none", "absent", "without", "intact",
    "wildtype", "wt", "normal", "retained", "unmutated", "undeleted",
    "unmodified", "benign", "silent", "synonymous",
})


def _status_text_is_disrupted(text: str) -> bool:
    """Return True only for affirmative disruption labels, never negated strings
    like ``no_loh``, ``not_deleted``, ``unmutated``, or ``non-truncating``."""
    raw = (text or "").strip().lower()
    if not raw:
        return False
    tokens = [t for t in re.split(r"[^a-z0-9]+", raw) if t]
    if not tokens:
        return False
    if any(tok in ("intact", "wildtype", "wt", "normal", "retained", "unmutated", "undeleted") for tok in tokens):
        return False
    for idx, tok in enumerate(tokens):
        if tok in _DISRUPTED_WORDS:
            if idx > 0 and tokens[idx - 1] in ("no", "not", "non", "without"):
                continue
            return True
    return False


@dataclass(frozen=True)
class APMGeneStatus:
    """Typed per-gene status for an antigen presentation / IFN-gamma pathway gene,
    separating DNA somatic mutations/copy loss from RNA expression."""
    gene: str
    dna_status: str = "unknown"
    rna_status: str = "unknown"
    dna_variants: tuple[dict[str, Any], ...] = field(default_factory=tuple)
    rna_ntpm: float | None = None
    immune_control_ntpm: float | None = None
    assay_source: str = ""
    notes: str = ""

    def __post_init__(self):
        object.__setattr__(self, "gene", str(self.gene or "").strip().upper())
        object.__setattr__(self, "dna_status", str(self.dna_status or "unknown").strip().lower())
        object.__setattr__(self, "rna_status", str(self.rna_status or "unknown").strip().lower())
        object.__setattr__(
            self, "dna_variants",
            tuple(dict(v) for v in (self.dna_variants or ()) if isinstance(v, dict)),
        )
        object.__setattr__(self, "rna_ntpm", _optional_float(self.rna_ntpm, "rna_ntpm", min_val=0.0))
        object.__setattr__(
            self, "immune_control_ntpm",
            _optional_float(self.immune_control_ntpm, "immune_control_ntpm", min_val=0.0),
        )
        object.__setattr__(self, "assay_source", str(self.assay_source or "").strip())
        object.__setattr__(self, "notes", str(self.notes or "").strip())

    @property
    def is_disrupted(self) -> bool:
        return _status_text_is_disrupted(self.dna_status) or _status_text_is_disrupted(self.rna_status)

    def to_dict(self) -> dict[str, Any]:
        return {
            "gene": self.gene,
            "dna_status": self.dna_status,
            "rna_status": self.rna_status,
            "dna_variants": list(self.dna_variants),
            "rna_ntpm": self.rna_ntpm,
            "immune_control_ntpm": self.immune_control_ntpm,
            "is_disrupted": self.is_disrupted,
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
    genes: tuple[APMGeneStatus, ...] = field(default_factory=tuple)

    def __post_init__(self):
        raw_genes = self.genes or ()
        parsed = []
        for g in raw_genes:
            if isinstance(g, APMGeneStatus):
                parsed.append(g)
            elif isinstance(g, dict):
                parsed.append(APMGeneStatus(**{
                    k: v for k, v in g.items()
                    if k in ("gene", "dna_status", "rna_status", "dna_variants",
                             "rna_ntpm", "immune_control_ntpm", "assay_source", "notes")
                }))
        object.__setattr__(self, "genes", tuple(parsed))

    def any_recorded(self) -> bool:
        return bool(self.genes) or any(
            bool(getattr(self, k))
            for k in ("b2m_status", "tap1_status", "tap2_status",
                      "nlrc5_status", "jak1_status", "jak2_status", "notes")
        )

    def disrupted_genes(self) -> list[str]:
        out: list[str] = []
        for label, attr in (
            ("B2M", "b2m_status"),
            ("TAP1", "tap1_status"),
            ("TAP2", "tap2_status"),
            ("NLRC5", "nlrc5_status"),
            ("JAK1", "jak1_status"),
            ("JAK2", "jak2_status"),
        ):
            if _status_text_is_disrupted(getattr(self, attr) or ""):
                out.append(label)
        for gs in self.genes:
            if gs.is_disrupted and gs.gene not in out:
                out.append(gs.gene)
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
            "genes": [g.to_dict() for g in self.genes],
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
    dna_status: str = "unknown"
    rna_status: str = "unknown"
    loss_mechanism: str = "unknown"

    def to_dict(self) -> dict[str, Any]:
        return {
            **self.record.to_dict(),
            "effective_status": self.effective_status,
            "dna_status": self.dna_status,
            "rna_status": self.rna_status,
            "loss_mechanism": self.loss_mechanism,
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

_TRUNCATING_VARIANT_TOKENS = frozenset({
    "truncating", "lof", "frameshift", "nonsense", "prematurestop",
    "stopgained", "stop", "spliceacceptor", "splicedonor", "spliceloss",
    "startloss", "startlost", "exondeletion", "exonicdeletion", "homdel",
    "disruptive",
})
_MISSENSE_VARIANT_TOKENS = frozenset({
    "missense", "substitution", "inframe", "nonsynonymous",
})


def _classify_dna_somatic_variants(
    variants: tuple[dict[str, Any], ...],
) -> tuple[bool, bool]:
    """Return ``(has_truncating_lof, has_missense_or_inframe)`` from recorded
    DNA somatic variant dicts on an MHC allele."""
    has_lof = False
    has_missense = False
    for var in variants:
        if not isinstance(var, dict):
            continue
        if var.get("is_lof") is True or var.get("truncating") is True:
            has_lof = True
            continue
        raw = " ".join(
            str(var.get(k) or "")
            for k in ("impact", "effect_class", "consequence", "effect", "classification")
        ).lower()
        normalized = re.sub(r"[^a-z0-9]+", " ", raw).strip()
        compact_tokens = set(normalized.split()) | {
            t.replace("_", "") for t in re.split(r"[^a-z0-9_]+", raw) if t
        }
        if any(tok in ("silent", "synonymous", "benign", "intronic", "noncoding") for tok in compact_tokens):
            continue
        if "inframe" in normalized.replace(" ", ""):
            has_missense = True
            continue
        if any(tok in _TRUNCATING_VARIANT_TOKENS for tok in compact_tokens):
            has_lof = True
        elif any(tok in _MISSENSE_VARIANT_TOKENS for tok in compact_tokens):
            has_missense = True
    return has_lof, has_missense


def resolve_effective_status(
    record: MHCAlleleStatus,
    policy: MHCAlleleStatusPolicy | dict[str, Any] | None = None,
) -> ResolvedMHCAlleleDecision:
    """Derive effective presentation status, separate DNA vs RNA status,
    caveats, and ranking weight from multi-assay evidence under ``policy``.

    Never conflates ``unresolvable`` (e.g. 10x 3' scRNA with low UMIs or
    low-purity bulk RNA) with confirmed ``lost`` or ``silenced``, and never
    lets low-UMI scRNA short-circuit clean high-purity bulk RNA ASE.
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

    # 1. DNA axis: copy-number LOH and/or somatic truncating LoF variants
    has_dna_lof, has_dna_missense = _classify_dna_somatic_variants(
        record.dna_somatic_variants
    )
    dna_loh_confirmed = record.dna_loh is True or (
        record.dna_copy_number is not None
        and record.dna_copy_number <= pol.dna_loh_max_copy_number
    )
    if dna_loh_confirmed:
        dna_status = "lost"
        if has_dna_lof:
            reasons.append("dna_truncating_lof_variant")
    elif has_dna_lof:
        dna_status = "truncating_lof"
    elif has_dna_missense:
        dna_status = "missense_or_inframe"
        caveats.append("dna_somatic_missense_variant")
    elif record.dna_loh is False or (
        record.dna_copy_number is not None
        and record.dna_copy_number > pol.dna_loh_max_copy_number
    ):
        dna_status = "intact"
    else:
        dna_status = "unknown"

    # 2. RNA axis: purity-independent scRNA first (when resolvable), else bulk RNA ASE
    scrna_ratio: float | None = None
    if (not scrna_low_umis
            and record.scrna_tumor_cluster_ntpm is not None
            and record.scrna_immune_control_ntpm is not None
            and record.scrna_immune_control_ntpm > 0):
        scrna_ratio = record.scrna_tumor_cluster_ntpm / record.scrna_immune_control_ntpm

    rna_reason = ""
    if scrna_ratio is not None:
        if scrna_ratio <= pol.silenced_tumor_to_immune_ratio:
            rna_status = "silenced"
            rna_reason = "scrna_tumor_cluster_silenced"
        else:
            rna_status = "expressed"
            rna_reason = "scrna_tumor_cluster_expressed"
    elif record.bulk_rna_allele_ratio is not None:
        if record.bulk_rna_allele_ratio <= pol.silenced_bulk_allele_ratio:
            if purity_confounded:
                rna_status = "unresolvable"
                rna_reason = "bulk_ase_low_but_purity_confounded"
            else:
                rna_status = "silenced"
                rna_reason = "bulk_rna_ase_silenced"
        else:
            rna_status = "expressed"
            rna_reason = "bulk_rna_ase_expressed"
            if purity_confounded:
                caveats.append("infiltrate_may_mask_tumor_loss")
    elif scrna_low_umis:
        rna_status = "unresolvable"
        rna_reason = "scrna_insufficient_resolvable_umis"
    else:
        rna_status = "unknown"

    # 3. Effective presentation status & loss mechanism
    effective = record.status
    if effective != "unknown":
        reasons.append(f"explicit_status:{effective}")
        if dna_status in ("lost", "truncating_lof") and rna_status == "silenced":
            loss_mechanism = "dna_and_rna"
        elif dna_status == "lost":
            loss_mechanism = "dna_loh"
        elif dna_status == "truncating_lof":
            loss_mechanism = "dna_truncating_lof"
        elif rna_status == "silenced":
            loss_mechanism = "rna_silenced"
        elif effective in ("lost", "silenced"):
            loss_mechanism = "explicit_override"
        elif effective == "retained":
            loss_mechanism = "none"
        else:
            loss_mechanism = effective
    else:
        if dna_status == "lost":
            effective = "lost"
            reasons.append("dna_loh_confirmed")
            if rna_reason:
                reasons.append(rna_reason)
            loss_mechanism = "dna_and_rna" if rna_status == "silenced" else "dna_loh"
        elif dna_status == "truncating_lof":
            effective = "lost"
            reasons.append("dna_truncating_lof_variant")
            if rna_reason:
                reasons.append(rna_reason)
            loss_mechanism = (
                "dna_and_rna" if rna_status == "silenced" else "dna_truncating_lof"
            )
        elif rna_status == "silenced":
            effective = "silenced"
            reasons.append(rna_reason)
            loss_mechanism = "rna_silenced"
        elif rna_status == "expressed":
            effective = "retained"
            reasons.append(rna_reason)
            loss_mechanism = "none"
        elif rna_status == "unresolvable":
            effective = "unresolvable"
            reasons.append(rna_reason)
            loss_mechanism = "unresolvable"
        elif dna_status == "intact":
            effective = "retained"
            reasons.append("dna_no_loh")
            loss_mechanism = "none"
        else:
            effective = "unknown"
            reasons.append("no_decisive_assay_evidence")
            loss_mechanism = "unknown"

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
        dna_status=dna_status,
        rna_status=rna_status,
        loss_mechanism=loss_mechanism,
    )


def _record_from_mapping(row: dict[str, Any]) -> MHCAlleleStatus:
    allowed = {
        "allele", "status", "dna_copy_number", "dna_loh", "dna_loh_pvalue",
        "dna_somatic_variants",
        "bulk_rna_tpm", "bulk_rna_allele_ratio", "tumor_purity",
        "bulk_purity_confounded", "scrna_tumor_cluster_ntpm",
        "scrna_immune_control_ntpm", "scrna_resolvable_umis",
        "scrna_chemistry", "assay_source", "notes",
    }
    kwargs = {k: v for k, v in row.items() if k in allowed}
    if "allele" not in kwargs or not str(kwargs["allele"] or "").strip():
        raise ValueError(f"MHC allele status record missing required 'allele' field: {row!r}")
    return MHCAlleleStatus(**kwargs)


def bundle_from_dict(
    data: dict[str, Any] | list[Any],
    *,
    source_path: str = "",
) -> MHCAlleleStatusBundle:
    """Reconstruct an ``MHCAlleleStatusBundle`` from a parsed JSON dict or list."""
    if isinstance(data, list):
        alleles = [_record_from_mapping(item) for item in data]
        return MHCAlleleStatusBundle(alleles=tuple(alleles), source_path=source_path)
    if not isinstance(data, dict):
        raise ValueError(f"Expected dict or list for MHC allele status bundle, got {type(data).__name__}")
    raw_alleles = data.get("alleles")
    if raw_alleles is None:
        raw_alleles = []
    if not isinstance(raw_alleles, list):
        raise ValueError("'alleles' in MHC allele status bundle must be a list")
    alleles = [_record_from_mapping(item) for item in raw_alleles]
    apm = None
    if isinstance(data.get("apm"), dict):
        apm_fields: dict[str, Any] = {
            k: str(v or "").strip()
            for k, v in data["apm"].items()
            if k in ("b2m_status", "tap1_status", "tap2_status",
                     "nlrc5_status", "jak1_status", "jak2_status", "notes")
        }
        if isinstance(data["apm"].get("genes"), list):
            apm_fields["genes"] = data["apm"]["genes"]
        apm = APMStatus(**apm_fields)
    resolved_src = source_path or str(data.get("source_path") or "")
    return MHCAlleleStatusBundle(alleles=tuple(alleles), apm=apm, source_path=resolved_src)


def load_mhc_allele_status(path: str | Path) -> MHCAlleleStatusBundle:
    """Load an ``MHCAlleleStatusBundle`` from JSON or TSV/CSV."""
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(f"MHC allele status file not found: {p}")

    suffix = p.suffix.lower()
    if suffix == ".json":
        with open(p, encoding="utf-8") as handle:
            data = json.load(handle)
        return bundle_from_dict(data, source_path=str(p))

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


def _merge_apm_status(
    base: APMStatus | None,
    override: APMStatus | None,
) -> APMStatus | None:
    if base is None or not base.any_recorded():
        return override if (override and override.any_recorded()) else None
    if override is None or not override.any_recorded():
        return base
    merged_genes: dict[str, APMGeneStatus] = {g.gene: g for g in base.genes}
    for og in override.genes:
        bg = merged_genes.get(og.gene)
        if bg is None:
            merged_genes[og.gene] = og
        else:
            combined_vars = list(bg.dna_variants)
            for v in og.dna_variants:
                if v not in combined_vars:
                    combined_vars.append(v)
            merged_genes[og.gene] = APMGeneStatus(
                gene=og.gene,
                dna_status=og.dna_status if og.dna_status != "unknown" else bg.dna_status,
                rna_status=og.rna_status if og.rna_status != "unknown" else bg.rna_status,
                dna_variants=tuple(combined_vars),
                rna_ntpm=og.rna_ntpm if og.rna_ntpm is not None else bg.rna_ntpm,
                immune_control_ntpm=(
                    og.immune_control_ntpm
                    if og.immune_control_ntpm is not None
                    else bg.immune_control_ntpm
                ),
                assay_source="; ".join(s for s in (bg.assay_source, og.assay_source) if s),
                notes="; ".join(s for s in (bg.notes, og.notes) if s),
            )
    return APMStatus(
        b2m_status=override.b2m_status or base.b2m_status,
        tap1_status=override.tap1_status or base.tap1_status,
        tap2_status=override.tap2_status or base.tap2_status,
        nlrc5_status=override.nlrc5_status or base.nlrc5_status,
        jak1_status=override.jak1_status or base.jak1_status,
        jak2_status=override.jak2_status or base.jak2_status,
        notes="; ".join(s for s in (base.notes, override.notes) if s),
        genes=tuple(merged_genes.values()),
    )


def _merge_allele_records(
    base: MHCAlleleStatus,
    override: MHCAlleleStatus,
) -> MHCAlleleStatus:
    """Combine complementary multi-assay fields (DNA LOH/variants, bulk RNA ASE,
    scRNA tumor-cluster nTPM) for the same canonical MHC allele, letting
    explicit non-empty values in ``override`` take precedence."""
    combined_vars = list(base.dna_somatic_variants)
    for v in override.dna_somatic_variants:
        if v not in combined_vars:
            combined_vars.append(v)

    def _pick(ov, bv):
        return ov if ov is not None else bv

    status = override.status if override.status != "unknown" else base.status
    assay_src = "; ".join(
        dict.fromkeys(s for s in (base.assay_source, override.assay_source) if s)
    )
    notes = "; ".join(
        dict.fromkeys(s for s in (base.notes, override.notes) if s)
    )
    return MHCAlleleStatus(
        allele=override.allele,
        status=status,
        dna_copy_number=_pick(override.dna_copy_number, base.dna_copy_number),
        dna_loh=_pick(override.dna_loh, base.dna_loh),
        dna_loh_pvalue=_pick(override.dna_loh_pvalue, base.dna_loh_pvalue),
        dna_somatic_variants=tuple(combined_vars),
        bulk_rna_tpm=_pick(override.bulk_rna_tpm, base.bulk_rna_tpm),
        bulk_rna_allele_ratio=_pick(override.bulk_rna_allele_ratio, base.bulk_rna_allele_ratio),
        tumor_purity=_pick(override.tumor_purity, base.tumor_purity),
        bulk_purity_confounded=_pick(override.bulk_purity_confounded, base.bulk_purity_confounded),
        scrna_tumor_cluster_ntpm=_pick(
            override.scrna_tumor_cluster_ntpm, base.scrna_tumor_cluster_ntpm
        ),
        scrna_immune_control_ntpm=_pick(
            override.scrna_immune_control_ntpm, base.scrna_immune_control_ntpm
        ),
        scrna_resolvable_umis=_pick(override.scrna_resolvable_umis, base.scrna_resolvable_umis),
        scrna_chemistry=override.scrna_chemistry or base.scrna_chemistry,
        assay_source=assay_src,
        notes=notes,
    )


def merge_mhc_allele_status_bundles(
    base: MHCAlleleStatusBundle | None,
    override: MHCAlleleStatusBundle | None,
) -> MHCAlleleStatusBundle | None:
    """Merge two bundles (e.g. auto-detected DNA MHC/APM variants +
    ``--mhc-allele-status`` file + ``--mhc-loh-alleles``). Complementary
    multi-assay axes (DNA, bulk RNA, scRNA) are combined per canonical
    ``mhcgnomes`` allele while explicit values in ``override`` win."""
    if base is None and override is None:
        return None
    if base is None:
        return override
    if override is None:
        return base
    merged: dict[str, MHCAlleleStatus] = {r.allele: r for r in base.alleles}
    for r in override.alleles:
        prev = merged.get(r.allele)
        merged[r.allele] = _merge_allele_records(prev, r) if prev is not None else r
    apm = _merge_apm_status(base.apm, override.apm)
    src = "; ".join(s for s in (base.source_path, override.source_path) if s)
    return MHCAlleleStatusBundle(alleles=tuple(merged.values()), apm=apm, source_path=src)


def apply_mhc_allele_status_to_ranked(
    ranked_variants_with_vaccine_peptides,
    bundle: MHCAlleleStatusBundle,
    policy: MHCAlleleStatusPolicy | dict[str, Any] | None = None,
) -> tuple[list, dict[str, Any]]:
    """Apply policy-resolved MHC allele status to ranked vaccine peptides.

    Removes (or continuously scales under ``continuous_ratio``) positive ranking
    credit on ``CandidateEpitope.per_allele_scores`` (and ``ms_scoring`` for
    masked alleles) on target tumor epitopes while preserving all raw
    ``CandidateEpitope.predictions`` intact and leaving autoimmune
    ``self_epitopes`` unmasked (normal tissues retain intact MHC). Recomputes
    ``VaccinePeptide`` scores and re-sorts both per-source windows and sources.

    Returns ``(updated_ranked, audit_dict)``.
    """
    from .ranking import rank_constructs

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
                # Somatic tumor MHC LOH/silencing does NOT affect normal healthy
                # tissues; never mask self-reference or non-targetable epitopes.
                if not ep.is_targetable or ep.occurs_in_self_reference:
                    new_epitopes.append(ep)
                    continue
                if not ep.per_allele_scores and not ep.ms_scoring:
                    new_epitopes.append(ep)
                    continue
                ep_touched = False
                new_scores = dict(ep.per_allele_scores)
                new_ms_scoring = dict(ep.ms_scoring) if ep.ms_scoring else {}
                ep_masked_alleles: list[str] = []
                for raw_allele, score in list(ep.per_allele_scores.items()):
                    canonical = normalize_mhc_allele_name(raw_allele)
                    decision = decisions_by_allele.get(canonical)
                    if decision is None or decision.ranking_weight == 1.0:
                        continue
                    if decision.masked:
                        # Remove key completely so ms_increment cannot leak MS boost
                        # through a zero-valued per_allele_scores entry.
                        new_scores.pop(raw_allele, None)
                        ep_masked_alleles.append(canonical)
                        delta = float(score)
                        if abs(delta) > 1e-12:
                            total_credit_removed += delta
                            allele_credit_removed[canonical] = (
                                allele_credit_removed.get(canonical, 0.0) + delta
                            )
                        ep_touched = True
                    else:
                        adjusted = float(score) * float(decision.ranking_weight)
                        delta = float(score) - adjusted
                        if abs(delta) > 1e-12:
                            new_scores[raw_allele] = adjusted
                            total_credit_removed += delta
                            allele_credit_removed[canonical] = (
                                allele_credit_removed.get(canonical, 0.0) + delta
                            )
                            ep_touched = True

                # Also strip masked alleles from ms_scoring directly
                if new_ms_scoring:
                    for ms_allele in list(new_ms_scoring.keys()):
                        ms_canonical = normalize_mhc_allele_name(ms_allele)
                        ms_dec = decisions_by_allele.get(ms_canonical)
                        if ms_dec is not None and ms_dec.masked:
                            new_ms_scoring.pop(ms_allele, None)
                            if ms_canonical not in ep_masked_alleles:
                                ep_masked_alleles.append(ms_canonical)
                            ep_touched = True

                if ep_touched:
                    old_ms_inc = float(ep.ms_increment)
                    new_ep = replace(
                        ep,
                        per_allele_scores=new_scores,
                        ms_scoring=new_ms_scoring,
                    )
                    ms_delta = old_ms_inc - float(new_ep.ms_increment)
                    if ms_delta > 1e-12:
                        total_credit_removed += ms_delta
                        if len(ep_masked_alleles) == 1:
                            only_a = ep_masked_alleles[0]
                            allele_credit_removed[only_a] = (
                                allele_credit_removed.get(only_a, 0.0) + ms_delta
                            )
                    new_epitopes.append(new_ep)
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
            if len(new_vps) > 1:
                def _vp_window_sort_key(p):
                    try:
                        return (-p.combined_score, p.lexicographic_sort_key())
                    except AttributeError:
                        return (
                            -p.combined_score,
                            -(p.target_epitope_score or 0.0),
                            p.self_epitope_score or 0.0,
                        )
                new_vps.sort(key=_vp_window_sort_key)
            desc = getattr(source, "short_description", None) or str(source)
            affected_sources.append(desc)
        updated_ranked.append((source, new_vps))

    # Re-sort sources using canonical vaxrank source ranking (combined_score,
    # comparable RNA support, target_epitope_score) while preserving tie order.
    updated_ranked = rank_constructs(updated_ranked)

    masked_alleles = [
        allele for allele, dec in decisions_by_allele.items() if dec.masked
    ]
    lost_alleles = [
        allele for allele, dec in decisions_by_allele.items()
        if dec.effective_status == "lost"
    ]
    silenced_alleles = [
        allele for allele, dec in decisions_by_allele.items()
        if dec.effective_status == "silenced"
    ]
    dna_lost_alleles = [
        allele for allele, dec in decisions_by_allele.items()
        if dec.dna_status in ("lost", "truncating_lof")
    ]
    rna_silenced_alleles = [
        allele for allele, dec in decisions_by_allele.items()
        if dec.rna_status == "silenced"
    ]
    unresolvable_alleles = [
        allele for allele, dec in decisions_by_allele.items()
        if dec.effective_status == "unresolvable"
    ]
    disrupted_apm_genes = (
        bundle.apm.disrupted_genes()
        if bundle.apm and bundle.apm.any_recorded()
        else []
    )

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
        "loss_of_mhc": {
            "has_loss": bool(lost_alleles or silenced_alleles),
            "lost_alleles": lost_alleles,
            "silenced_alleles": silenced_alleles,
            "dna_lost_alleles": dna_lost_alleles,
            "rna_silenced_alleles": rna_silenced_alleles,
            "unresolvable_alleles": unresolvable_alleles,
        },
        "loss_of_apm": {
            "has_loss": bool(disrupted_apm_genes),
            "disrupted_genes": disrupted_apm_genes,
        },
        "apm": bundle.apm.to_dict() if bundle.apm and bundle.apm.any_recorded() else None,
        "affected_epitope_count": affected_epitopes,
        "affected_sources": affected_sources,
        "total_credit_removed": round(total_credit_removed, 6),
    }
    return updated_ranked, audit
