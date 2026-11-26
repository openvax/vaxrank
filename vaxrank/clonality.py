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

"""Cross-species, multi-modality AntigenClonality estimation and integration (#4).

Models tumor-cell clonality across **all** vaccine antigen classes without
forcing an SNV-centric definition:

1. **GenomicClonality** (for genomically templated nuclear events: SNVs, indels,
   SVs/fusions, integrated oncoviruses):
   - Pure-Python parsers (zero hard upstream package dependencies) for:
     - PyClone-VI (``*.tsv``)
     - PURPLE somatic VCF (``PURPLE_CN``, ``PURPLE_MACN``, ``PURPLE_VCN``, ``SUBCL``)
       + ``*.purple.purity.tsv``
     - FACETS / facets-suite annotated MAF/TSV (``ccf_expected_copies``, ``tcn.em``,
       ``lcn.em``, ``clonality``)
     - PureCN variants CSV (``CCF``/``CELLFRACTION``, ``ML.C``, ``ML.M.SEGMENT``)
     - Sequenza segments + confints (``*_segments.txt`` + ``*_confints_CP.txt``)
     - LENS / pVACseq clonality columns + canonical Vaxrank TSV/JSON
   - Standalone cross-species **direct Beta-Binomial CCF** estimator from somatic
     VCF / read counts + tumor purity + optional CN segments (with species-aware
     sex-chromosome ploidy).
   - **Mitochondrial boundary**: refuses to apply nuclear diploid CCF math to
     mitochondrial contigs (``MT``/``chrM``/``chrMT``), preserving separate
     ``mt_heteroplasmy_tumor`` / ``mt_heteroplasmy_normal`` with explicit NUMT
     and copy-number caveats.
   - Multi-caller reconciliation on the same reads uses configurable precedence
     (never multiplying correlated callers), flags ``caller_discordance``, and
     widens the credible interval to the cross-caller envelope.

2. **ExpressionalClonality** (for **all** antigen classes when barcoded scRNA
   — especially full-length **long-read scRNA** — is available: SNVs, indels,
   fusions, splice isoforms, ncORFs, CTAs, ERVs, oncoviruses, MT transcripts):
   - Pure-Python parsers for short-read and long-read per-barcode allele/isoform
     outputs (VarTrix MatrixMarket, cellsnp-lite MatrixMarket, IsoQuant
     ``read_info.tsv``, FLAMES ``transcript_count.csv``, and generic per-barcode
     allele/isoform/gene tables), plus direct integration with
     :mod:`vaxrank.scrna_tumor` tumor clusters.
   - Explicitly separates **barcoded single-cell long-read RNA** from unbarcoded
     **bulk long-read RNA** (which provides isoform phasing and bulk expression
     confirmation, not cellular clonality).
   - Models abundance / dropout confounding explicitly via:
     - Allelic dropout power ``q_i = 1 - (1 - p_alt) ** d_i`` on locus-covered
       malignant cells (for mutations, junctions, and isoforms),
     - Zero-truncated Poisson UMI capture power ``1 - exp(-mu_pos)`` (for
       gene-level CTAs, ERVs, and viral transcripts),
     - Hierarchical fallback to **malignant-subcluster pseudobulk clonality**
       when per-cell detection power is too low.

3. **SpatialClonality** & **LongitudinalClonality**:
   - Spatial regions within a timepoint (e.g. multi-region primary resection or
     synchronous metastases) aggregate symmetrically.
   - Directed longitudinal timepoints (``T0 -> T1 -> T2 -> T3``) anchor current
     clonality on ``T_latest`` and classify trajectories as ``persistent``,
     ``emergent_latest``, ``lost_in_latest``, ``expanding``, ``contracting``, or
     ``latest_unassessed``.

4. **Non-double-counting cross-axis combination & asymmetric mixed-class policy**:
   - Allele-/junction-/viral-specific scRNA observes the joint event ``P(G=1, E=1)``
     and is precision-weighted (bounded by DNA) or flags an expression-silencing
     bottleneck — never multiplied by ``C_DNA``.
   - Gene-level scRNA observes marginal ``P(E_gene=1)`` and combines via ``min``
     (or configurable ``product``).
   - Asymmetric policy penalizes *verified subclonality*, never missing assay
     coverage across mixed SNV / SV / CTA / viral runs.
"""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, field
import gzip
import json
import math
from numbers import Real
from pathlib import Path
from typing import Any

DEFAULT_CALLER_PRECEDENCE: tuple[str, ...] = (
    "pyclone_vi",
    "dpclust",
    "mobster",
    "purple",
    "facets",
    "purecn",
    "sequenza",
    "lens",
    "pvacseq",
    "canonical",
    "direct_cn",
    "direct_diploid",
)

MITOCHONDRIAL_CONTIGS = frozenset({"MT", "CHRM", "CHRMT", "M", "MITO", "CHRMITO"})
VALID_CLONALITY_MODES = frozenset({"annotate", "weight", "gate", "tier"})
VALID_GENE_LEVEL_COMBINATIONS = frozenset({"min", "product"})
VALID_MULTI_REGION_AGGREGATIONS = frozenset({"mean", "min"})
VALID_ASSAY_PLATFORMS = frozenset({
    "long_read_ont",
    "long_read_pacbio",
    "short_read_5prime",
    "short_read_3prime",
    "smartseq2",
    "scrna_generic",
    "bulk_long_read",
    "bulk_short_read",
    "unknown",
})
BULK_ASSAY_PLATFORMS = frozenset({"bulk_long_read", "bulk_short_read"})
GENOMICALLY_TEMPLATED_KINDS = frozenset({"mutation", "snv", "indel", "mnv", "sv", "fusion"})


def _open_text(path: str | Path):
    p = Path(path)
    if p.suffix.lower() == ".gz":
        return gzip.open(p, "rt", encoding="utf-8")
    return open(p, encoding="utf-8")


def _safe_float(val: Any) -> float | None:
    if val is None or isinstance(val, bool):
        return None
    if isinstance(val, Real):
        f = float(val)
        return f if math.isfinite(f) else None
    text = str(val).strip()
    if not text or text.lower() in {"na", "nan", "none", "null", ".", ""}:
        return None
    try:
        f = float(text)
    except ValueError:
        return None
    return f if math.isfinite(f) else None


def _safe_int(val: Any) -> int | None:
    f = _safe_float(val)
    if f is None:
        return None
    return int(round(f))


def _clamp01(value: float) -> float:
    return max(0.0, min(1.0, float(value)))


def normalize_contig(contig: str | None) -> str:
    """Normalize chromosome/contig name without mutating mitochondrial detection."""
    text = str(contig or "").strip()
    if not text:
        return ""
    if text.lower().startswith("chr"):
        return text[3:].upper()
    return text.upper()


def is_mitochondrial_contig(contig: str | None) -> bool:
    """Return True if ``contig`` denotes mitochondrial DNA (MT/chrM/chrMT/M)."""
    raw = str(contig or "").strip().upper()
    return raw in MITOCHONDRIAL_CONTIGS or normalize_contig(raw) in MITOCHONDRIAL_CONTIGS


def normal_ploidy_for_contig(contig: str | None, *, sex: str | None = None) -> float:
    """Species-agnostic diploid default with hemizygous male X/Y support."""
    norm = normalize_contig(contig)
    sex_norm = str(sex or "").strip().lower()
    if sex_norm in {"male", "m", "xy"} and norm in {"X", "Y"}:
        return 1.0
    if sex_norm in {"female", "f", "xx"} and norm == "Y":
        return 0.0
    return 2.0


def wilson_interval(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    """95% Wilson score confidence interval for binomial proportion k / n."""
    if n <= 0:
        return (0.0, 1.0)
    k = max(0, min(int(k), int(n)))
    p = k / n
    z2 = z * z
    denom = 1.0 + z2 / n
    center = (p + z2 / (2.0 * n)) / denom
    half = (z * math.sqrt((p * (1.0 - p) + z2 / (4.0 * n)) / n)) / denom
    return (_clamp01(center - half), _clamp01(center + half))


def _solve_zero_truncated_poisson_mean(observed_positive_mean: float) -> float:
    """Solve ``observed_positive_mean = mu / (1 - exp(-mu))`` for ``mu > 0``."""
    y = float(observed_positive_mean)
    if not math.isfinite(y) or y <= 1.0:
        return 1e-4
    # Newton-Raphson on f(mu) = mu / (1 - exp(-mu)) - y = 0
    mu = max(y - 1.0, 1e-3)
    for _ in range(40):
        exp_neg = math.exp(-mu)
        denom = 1.0 - exp_neg
        if denom <= 1e-12:
            break
        val = mu / denom - y
        # Derivative: (1 - exp(-mu) - mu * exp(-mu)) / (1 - exp(-mu))^2
        deriv = (denom - mu * exp_neg) / (denom * denom)
        if abs(deriv) <= 1e-12:
            break
        step = val / deriv
        mu = max(1e-6, mu - step)
        if abs(step) < 1e-9:
            break
    return mu


@dataclass(frozen=True)
class ClonalityPolicy:
    """Per-trial YAML policy governing clonality estimation, combination, and ranking."""

    mode: str = "annotate"
    clonal_threshold: float = 0.85
    subclonal_threshold: float = 0.50
    min_clonal_probability: float = 0.50
    min_clonality: float = 0.30
    weight_exponent: float = 1.0
    subclonal_floor: float = 0.10
    unassessed_weight: float = 1.0
    persistence_bonus: float = 0.0
    caller_discordance_threshold: float = 0.25
    caller_precedence: tuple[str, ...] = DEFAULT_CALLER_PRECEDENCE
    min_scrna_detection_power: float = 0.50
    min_scrna_covered_cells: int = 5
    min_subcluster_cells: int = 5
    min_subcluster_ntpm: float = 1.0
    min_subcluster_alt_reads: int = 2
    min_latest_depth_for_loss: int = 20
    gene_level_combination: str = "min"
    multi_region_aggregation: str = "mean"
    require_measured: bool = False
    by_antigen_kind: dict[str, dict[str, Any]] = field(default_factory=dict)

    def __post_init__(self):
        mode = str(self.mode or "annotate").strip().lower()
        if mode not in VALID_CLONALITY_MODES:
            raise ValueError(
                f"clonality.mode must be one of {sorted(VALID_CLONALITY_MODES)}, got {self.mode!r}"
            )
        object.__setattr__(self, "mode", mode)
        for name in (
            "clonal_threshold", "subclonal_threshold", "min_clonal_probability",
            "min_clonality", "subclonal_floor", "unassessed_weight",
            "caller_discordance_threshold", "min_scrna_detection_power",
        ):
            val = _safe_float(getattr(self, name))
            if val is None or not 0.0 <= val <= 1.0:
                raise ValueError(f"clonality.{name} must be in [0.0, 1.0], got {getattr(self, name)!r}")
            object.__setattr__(self, name, val)
        if self.subclonal_threshold > self.clonal_threshold:
            raise ValueError("clonality.subclonal_threshold must be <= clonal_threshold")
        exp = _safe_float(self.weight_exponent)
        if exp is None or exp < 0.0:
            raise ValueError(f"clonality.weight_exponent must be >= 0, got {self.weight_exponent!r}")
        object.__setattr__(self, "weight_exponent", exp)
        pb = _safe_float(self.persistence_bonus)
        if pb is None or pb < 0.0:
            raise ValueError(f"clonality.persistence_bonus must be >= 0, got {self.persistence_bonus!r}")
        object.__setattr__(self, "persistence_bonus", pb)
        glc = str(self.gene_level_combination or "min").strip().lower()
        if glc not in VALID_GENE_LEVEL_COMBINATIONS:
            raise ValueError(
                f"clonality.gene_level_combination must be one of {sorted(VALID_GENE_LEVEL_COMBINATIONS)}"
            )
        object.__setattr__(self, "gene_level_combination", glc)
        mra = str(self.multi_region_aggregation or "mean").strip().lower()
        if mra not in VALID_MULTI_REGION_AGGREGATIONS:
            raise ValueError(
                f"clonality.multi_region_aggregation must be one of {sorted(VALID_MULTI_REGION_AGGREGATIONS)}"
            )
        object.__setattr__(self, "multi_region_aggregation", mra)
        object.__setattr__(self, "caller_precedence", tuple(str(x).lower() for x in self.caller_precedence))
        object.__setattr__(self, "by_antigen_kind", {
            str(k).lower(): dict(v) for k, v in (self.by_antigen_kind or {}).items()
        })

    def for_antigen_kind(self, antigen_kind: str) -> dict[str, Any]:
        """Return effective per-kind overrides merged with top-level defaults."""
        kind_key = str(antigen_kind or "mutation").strip().lower()
        override = self.by_antigen_kind.get(kind_key, {})
        return {
            "unassessed_weight": _safe_float(override.get("unassessed_weight"))
            if "unassessed_weight" in override else self.unassessed_weight,
            "require_measured": bool(override.get("require_measured", self.require_measured)),
            "min_clonality": _safe_float(override.get("min_clonality"))
            if "min_clonality" in override else self.min_clonality,
            "integrated_driver_prior": _safe_float(override.get("integrated_driver_prior")),
        }

    @classmethod
    def from_dict(cls, raw: dict[str, Any] | None) -> "ClonalityPolicy":
        if not raw:
            return cls()
        data = dict(raw)
        if "caller_precedence" in data and data["caller_precedence"] is not None:
            data["caller_precedence"] = tuple(data["caller_precedence"])
        allowed = {f.name for f in cls.__dataclass_fields__.values()}
        unknown = set(data) - allowed
        if unknown:
            raise ValueError(f"Unknown clonality policy field(s): {sorted(unknown)}")
        return cls(**data)

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["caller_precedence"] = list(self.caller_precedence)
        return d


@dataclass(frozen=True)
class CopyNumberSegment:
    """Allele-specific genomic copy number segment (1-based inclusive)."""

    contig: str
    start: int
    end: int
    cn_total: float
    cn_major: float | None = None
    cn_minor: float | None = None
    cellular_fraction: float | None = None
    source_tool: str = ""

    def contains(self, contig: str, position: int) -> bool:
        return (
            normalize_contig(self.contig) == normalize_contig(contig)
            and self.start <= position <= self.end
        )


@dataclass(frozen=True)
class CallerClonalityEstimate:
    """Single-caller DNA clonality estimate for one variant in one sample."""

    source_tool: str
    ccf: float | None = None
    ci_lower: float | None = None
    ci_upper: float | None = None
    clonal_probability: float | None = None
    multiplicity: float | None = None
    local_cn_total: float | None = None
    local_cn_minor: float | None = None
    cluster_id: str | None = None
    tumor_purity: float | None = None
    depth: int | None = None
    alt_reads: int | None = None
    vaf: float | None = None
    compartment: str = "nuclear"
    mt_heteroplasmy_tumor: float | None = None
    mt_heteroplasmy_normal: float | None = None
    assumed_diploid: bool = False
    sample_id: str = ""
    notes: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["notes"] = list(self.notes)
        return d


@dataclass(frozen=True)
class GenomicClonality:
    """Reconciled DNA-axis clonality for a genomically templated antigen."""

    dna_ccf: float | None = None
    ci_lower: float | None = None
    ci_upper: float | None = None
    clonal_probability: float | None = None
    multiplicity: float | None = None
    local_cn_total: float | None = None
    local_cn_minor: float | None = None
    cluster_id: str | None = None
    tumor_purity: float | None = None
    compartment: str = "nuclear"
    mt_heteroplasmy_tumor: float | None = None
    mt_heteroplasmy_normal: float | None = None
    primary_source: str = ""
    caller_estimates: tuple[CallerClonalityEstimate, ...] = ()
    caller_spread: float | None = None
    caller_discordance: bool = False
    depth: int | None = None
    caveats: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "dna_ccf": self.dna_ccf,
            "ci_lower": self.ci_lower,
            "ci_upper": self.ci_upper,
            "clonal_probability": self.clonal_probability,
            "multiplicity": self.multiplicity,
            "local_cn_total": self.local_cn_total,
            "local_cn_minor": self.local_cn_minor,
            "cluster_id": self.cluster_id,
            "tumor_purity": self.tumor_purity,
            "compartment": self.compartment,
            "mt_heteroplasmy_tumor": self.mt_heteroplasmy_tumor,
            "mt_heteroplasmy_normal": self.mt_heteroplasmy_normal,
            "primary_source": self.primary_source,
            "caller_estimates": [c.to_dict() for c in self.caller_estimates],
            "caller_spread": self.caller_spread,
            "caller_discordance": self.caller_discordance,
            "depth": self.depth,
            "caveats": list(self.caveats),
        }


@dataclass(frozen=True)
class ExpressionalClonality:
    """scRNA expressional clonality across short-read and long-read platforms."""

    assay_platform: str = "scrna_generic"
    measurement_level: str = "gene"  # allele_or_junction | isoform | gene | viral_transcript | mt_transcript | bulk_only
    n_tumor_cells: int = 0
    n_covered_tumor_cells: int = 0
    n_expressing_tumor_cells: int = 0
    locus_coverage_fraction: float | None = None
    detected_cell_fraction: float | None = None
    detected_ci_lower: float | None = None
    detected_ci_upper: float | None = None
    mean_umi_or_depth_positive: float | None = None
    estimated_allelic_fraction: float | None = None
    detection_power: float | None = None
    dropout_adjusted_clonality: float | None = None
    dropout_ci_lower: float | None = None
    dropout_ci_upper: float | None = None
    expressing_subcluster_fraction: float | None = None
    n_powered_subclusters: int = 0
    n_expressing_subclusters: int = 0
    selected_clonality: float | None = None
    selected_ci_lower: float | None = None
    selected_ci_upper: float | None = None
    selected_scale: str = "unassessed"  # per_cell | subcluster_pseudobulk | unassessed
    abundance_confounded: bool = False
    is_bimodal_heterogeneous: bool = False
    bulk_expression_value: float | None = None
    caveats: tuple[str, ...] = ()

    @property
    def is_specific_to_alteration(self) -> bool:
        return self.measurement_level in {
            "allele_or_junction", "isoform", "viral_transcript", "mt_transcript",
        }

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["caveats"] = list(self.caveats)
        return d


@dataclass(frozen=True)
class SampleClonalityObservation:
    """One sample's clonality evidence tagged with timepoint and physical region."""

    sample_id: str
    timepoint_index: int = 0
    timepoint_label: str = "T0"
    region_label: str = "primary"
    genomic: GenomicClonality | None = None
    expressional: ExpressionalClonality | None = None
    depth: int | None = None
    present: bool | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "timepoint_index": self.timepoint_index,
            "timepoint_label": self.timepoint_label,
            "region_label": self.region_label,
            "genomic": self.genomic.to_dict() if self.genomic else None,
            "expressional": self.expressional.to_dict() if self.expressional else None,
            "depth": self.depth,
            "present": self.present,
        }


@dataclass(frozen=True)
class SpatialClonality:
    """Spatial summary across physical regions / synchronous sites at one timepoint."""

    timepoint_index: int
    timepoint_label: str
    regions_sampled: int
    regions_powered: int
    regions_present: int
    region_fraction: float | None
    min_region_clonality: float | None
    mean_region_clonality: float | None
    is_spatially_ubiquitous: bool | None
    region_labels: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["region_labels"] = list(self.region_labels)
        return d


@dataclass(frozen=True)
class LongitudinalClonality:
    """Directed temporal summary across ordered timepoints (T0 -> T1 -> ... -> T_latest)."""

    timepoints_ordered: tuple[str, ...]
    timepoints_powered: int
    timepoints_present: int
    latest_timepoint: str
    present_in_latest: bool | None
    latest_clonality: float | None
    trajectory: str  # persistent | emergent_latest | lost_in_latest | expanding | contracting | latest_unassessed | single_timepoint
    per_timepoint_spatial: tuple[SpatialClonality, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "timepoints_ordered": list(self.timepoints_ordered),
            "timepoints_powered": self.timepoints_powered,
            "timepoints_present": self.timepoints_present,
            "latest_timepoint": self.latest_timepoint,
            "present_in_latest": self.present_in_latest,
            "latest_clonality": self.latest_clonality,
            "trajectory": self.trajectory,
            "per_timepoint_spatial": [s.to_dict() for s in self.per_timepoint_spatial],
        }


@dataclass(frozen=True)
class AntigenClonality:
    """Unified multi-axis tumor-cell clonality record for any antigen class."""

    antigen_id: str
    antigen_kind: str = "mutation"
    genomic: GenomicClonality | None = None
    expressional: ExpressionalClonality | None = None
    spatial: SpatialClonality | None = None
    longitudinal: LongitudinalClonality | None = None
    effective_clonality: float | None = None
    ci_lower: float | None = None
    ci_upper: float | None = None
    clonality_status: str = "unassessed"  # clonal | subclonal | uncertain | lost | unassessed
    combination_rule: str = "unassessed"
    divergences: tuple[str, ...] = ()
    caveats: tuple[str, ...] = ()

    def score_multiplier(self, policy: ClonalityPolicy | None = None) -> float:
        """Compute policy-configured score multiplier (1.0 in annotate mode)."""
        pol = policy or ClonalityPolicy()
        if pol.mode == "annotate":
            return 1.0
        kind_cfg = pol.for_antigen_kind(self.antigen_kind)
        if self.clonality_status == "lost":
            return pol.subclonal_floor
        if self.effective_clonality is None or self.clonality_status in {"unassessed", "uncertain"}:
            return float(kind_cfg["unassessed_weight"])
        base = max(pol.subclonal_floor, _clamp01(self.effective_clonality) ** pol.weight_exponent)
        if (
            pol.persistence_bonus > 0.0
            and self.longitudinal is not None
            and self.longitudinal.trajectory == "persistent"
        ):
            base *= (1.0 + pol.persistence_bonus)
        return base

    def passes_gate(self, policy: ClonalityPolicy | None = None) -> bool:
        """Evaluate hard gate in ``gate`` mode without dropping unassessed antigens by default."""
        pol = policy or ClonalityPolicy()
        if pol.mode != "gate":
            return True
        kind_cfg = pol.for_antigen_kind(self.antigen_kind)
        if self.clonality_status == "lost":
            return False
        if self.effective_clonality is None or self.clonality_status in {"unassessed", "uncertain"}:
            return not bool(kind_cfg["require_measured"])
        min_c = float(kind_cfg["min_clonality"])
        if self.effective_clonality >= min_c:
            return True
        # Guard against dropping shallow-depth / wide-CI candidates whose upper bound
        # is still consistent with clonal expression/occurrence.
        if self.ci_upper is not None and self.ci_upper >= pol.clonal_threshold:
            return True
        return False

    def tier_rank(self) -> int:
        """Lexicographic tier rank for ``tier`` mode (lower is higher priority)."""
        if self.clonality_status == "clonal":
            if self.longitudinal is not None and self.longitudinal.trajectory == "persistent":
                return 0
            return 1
        if self.clonality_status in {"unassessed", "uncertain"}:
            return 2
        if self.clonality_status == "subclonal":
            return 3
        return 4  # lost

    def to_dict(self) -> dict[str, Any]:
        return {
            "antigen_id": self.antigen_id,
            "antigen_kind": self.antigen_kind,
            "genomic": self.genomic.to_dict() if self.genomic else None,
            "expressional": self.expressional.to_dict() if self.expressional else None,
            "spatial": self.spatial.to_dict() if self.spatial else None,
            "longitudinal": self.longitudinal.to_dict() if self.longitudinal else None,
            "effective_clonality": self.effective_clonality,
            "ci_lower": self.ci_lower,
            "ci_upper": self.ci_upper,
            "clonality_status": self.clonality_status,
            "combination_rule": self.combination_rule,
            "divergences": list(self.divergences),
            "caveats": list(self.caveats),
        }


# ---------------------------------------------------------------------------
# Direct Cross-Species Computational Estimators (Zero External Dependencies)
# ---------------------------------------------------------------------------

def estimate_direct_dna_ccf(
    *,
    alt_reads: int,
    total_depth: int,
    tumor_purity: float,
    contig: str = "",
    position: int | None = None,
    cn_total: float | None = None,
    cn_major: float | None = None,
    cn_minor: float | None = None,
    cn_segments: list[CopyNumberSegment] | tuple[CopyNumberSegment, ...] | None = None,
    sex: str | None = None,
    normal_alt_reads: int | None = None,
    normal_depth: int | None = None,
    clonal_threshold: float = 0.85,
    sample_id: str = "",
    grid_points: int = 201,
    seq_error_rate: float = 1e-3,
) -> CallerClonalityEstimate:
    """Direct cross-species Beta-Binomial/Binomial grid CCF estimator.

    Refuses to apply nuclear diploid CCF math to mitochondrial contigs
    (``MT``/``chrM``/``chrMT``), returning separate tumor/normal heteroplasmy
    with ``ccf=None`` and explicit caveats (openvax/vaxrank#4).
    """
    alt = max(0, int(alt_reads))
    depth = max(0, int(total_depth))
    vaf = (alt / depth) if depth > 0 else None

    if is_mitochondrial_contig(contig):
        normal_het = (
            (normal_alt_reads / normal_depth)
            if (normal_alt_reads is not None and normal_depth and normal_depth > 0)
            else None
        )
        return CallerClonalityEstimate(
            source_tool="direct_mt_heteroplasmy",
            ccf=None,
            ci_lower=None,
            ci_upper=None,
            clonal_probability=None,
            multiplicity=None,
            local_cn_total=None,
            local_cn_minor=None,
            tumor_purity=_safe_float(tumor_purity),
            depth=depth,
            alt_reads=alt,
            vaf=vaf,
            compartment="mitochondrial",
            mt_heteroplasmy_tumor=vaf,
            mt_heteroplasmy_normal=normal_het,
            assumed_diploid=False,
            sample_id=sample_id,
            notes=(
                "mitochondrial_non_diploid_ccf_not_applicable",
                "mt_heteroplasmy_confounded_by_copy_number_and_numt_risk",
            ),
        )

    purity = _safe_float(tumor_purity)
    if purity is None or not 0.0 < purity <= 1.0:
        raise ValueError(f"tumor_purity must be in (0.0, 1.0], got {tumor_purity!r}")

    notes: list[str] = []
    assumed_diploid = False
    if cn_total is None and cn_segments and contig and position is not None:
        for seg in cn_segments:
            if seg.contains(contig, position):
                cn_total = seg.cn_total
                cn_major = seg.cn_major if cn_major is None else cn_major
                cn_minor = seg.cn_minor if cn_minor is None else cn_minor
                break

    cn_normal = normal_ploidy_for_contig(contig, sex=sex)
    if cn_total is None:
        cn_total = cn_normal if cn_normal > 0 else 2.0
        cn_major = max(1.0, math.ceil(cn_total / 2.0))
        cn_minor = max(0.0, cn_total - cn_major)
        assumed_diploid = True
        notes.append("assumed_diploid_cn")
    else:
        cn_total = max(0.0, float(cn_total))
        if cn_major is None and cn_minor is not None:
            cn_major = max(0.0, cn_total - float(cn_minor))
        elif cn_major is None:
            cn_major = max(1.0, float(cn_total))
        if cn_minor is None:
            cn_minor = max(0.0, cn_total - float(cn_major))
        if cn_minor == 0.0 and cn_total >= 1.0:
            notes.append("local_loh_or_cnloh")

    if depth == 0 or cn_total <= 0.0:
        notes.append("zero_depth_or_homozygous_deletion")
        return CallerClonalityEstimate(
            source_tool="direct_diploid" if assumed_diploid else "direct_cn",
            ccf=None,
            ci_lower=None,
            ci_upper=None,
            clonal_probability=None,
            local_cn_total=cn_total,
            local_cn_minor=cn_minor,
            tumor_purity=purity,
            depth=depth,
            alt_reads=alt,
            vaf=vaf,
            assumed_diploid=assumed_diploid,
            sample_id=sample_id,
            notes=tuple(notes),
        )

    denom_cn = purity * cn_total + (1.0 - purity) * cn_normal
    max_mult = max(1, int(round(cn_major or cn_total)))
    # Choose valid multiplicity m in {1..max_mult} via empirical Bayes / MAP:
    # prefer m=1 unless observed VAF clearly exceeds 1-copy clonal expectation
    # (prevents overfitting subclonal variants onto amplified multiplicities).
    one_copy_clonal_vaf = (purity * 1.0) / denom_cn
    if max_mult == 1 or (vaf is not None and vaf <= one_copy_clonal_vaf * 1.25):
        candidate_multiplicities = [1]
    else:
        raw_mult = max(1, min(max_mult, int(round((vaf or 0.0) * denom_cn / purity))))
        candidate_multiplicities = sorted({1, raw_mult, min(max_mult, raw_mult + 1)})

    grid = [i / (grid_points - 1) for i in range(grid_points)]
    best_m = 1
    best_log_evidence = -math.inf
    best_post: list[float] = []

    for m in candidate_multiplicities:
        log_liks: list[float] = []
        for ccf in grid:
            expected_p = (purity * m * ccf) / denom_cn
            p_obs = max(seq_error_rate, min(1.0 - seq_error_rate, expected_p))
            ll = alt * math.log(p_obs) + (depth - alt) * math.log(1.0 - p_obs)
            # Mild parsimony prior favoring m=1 when likelihoods tie
            if m > 1:
                ll -= 0.5 * math.log(m)
            log_liks.append(ll)
        max_ll = max(log_liks)
        weights = [math.exp(ll - max_ll) for ll in log_liks]
        total_w = sum(weights)
        log_ev = max_ll + math.log(total_w / len(grid))
        if log_ev > best_log_evidence:
            best_log_evidence = log_ev
            best_m = m
            best_post = [w / total_w for w in weights]

    post_mean = sum(c * w for c, w in zip(grid, best_post))
    clonal_prob = sum(w for c, w in zip(grid, best_post) if c >= clonal_threshold)

    # 95% equal-tailed credible interval
    cdf = 0.0
    ci_lo = 0.0
    ci_hi = 1.0
    found_lo = False
    for c, w in zip(grid, best_post):
        cdf += w
        if not found_lo and cdf >= 0.025:
            ci_lo = c
            found_lo = True
        if cdf >= 0.975:
            ci_hi = c
            break

    if depth < 15 or (ci_hi - ci_lo) > 0.45:
        notes.append("low_depth_uncertain")

    return CallerClonalityEstimate(
        source_tool="direct_diploid" if assumed_diploid else "direct_cn",
        ccf=round(_clamp01(post_mean), 4),
        ci_lower=round(_clamp01(ci_lo), 4),
        ci_upper=round(_clamp01(ci_hi), 4),
        clonal_probability=round(_clamp01(clonal_prob), 4),
        multiplicity=float(best_m),
        local_cn_total=cn_total,
        local_cn_minor=cn_minor,
        tumor_purity=purity,
        depth=depth,
        alt_reads=alt,
        vaf=round(vaf, 4) if vaf is not None else None,
        compartment="nuclear",
        assumed_diploid=assumed_diploid,
        sample_id=sample_id,
        notes=tuple(notes),
    )


def estimate_scrna_gene_expressional_clonality(
    tumor_cell_counts: dict[str, int | float],
    *,
    cell_to_subcluster: dict[str, str] | None = None,
    cell_total_umis: dict[str, int | float] | None = None,
    assay_platform: str = "scrna_generic",
    measurement_level: str = "gene",
    policy: ClonalityPolicy | None = None,
) -> ExpressionalClonality:
    """Estimate scRNA expressional clonality for a gene/CTA/ERV/viral transcript.

    Uses a zero-truncated Poisson capture model to separate genuine bimodal
    subclonal expression from low-abundance technical dropout, and computes
    malignant-subcluster pseudobulk clonality as a high-power fallback.
    """
    pol = policy or ClonalityPolicy()
    if assay_platform in BULK_ASSAY_PLATFORMS:
        return ExpressionalClonality(
            assay_platform=assay_platform,
            measurement_level="bulk_only",
            selected_clonality=None,
            selected_scale="unassessed",
            caveats=("bulk_rna_only_no_single_cell_barcodes",),
        )

    n_tumor = len(tumor_cell_counts)
    if n_tumor == 0:
        return ExpressionalClonality(
            assay_platform=assay_platform,
            measurement_level=measurement_level,
            caveats=("no_tumor_cells",),
        )

    pos_counts = [float(v) for v in tumor_cell_counts.values() if float(v) > 0]
    n_pos = len(pos_counts)
    f_obs = n_pos / n_tumor
    w_lo, w_hi = wilson_interval(n_pos, n_tumor)

    caveats: list[str] = []
    if n_pos == 0:
        mean_pos = 0.0
        power = None
        adj = 0.0
        adj_lo, adj_hi = 0.0, w_hi
        abundance_confounded = True
        is_bimodal = False
        caveats.append("zero_detected_tumor_cells_cannot_exclude_dropout")
    else:
        mean_pos = sum(pos_counts) / n_pos
        mu_pos = _solve_zero_truncated_poisson_mean(mean_pos)
        power = _clamp01(1.0 - math.exp(-mu_pos))
        if power >= 0.05:
            adj = _clamp01(f_obs / power)
            adj_lo = _clamp01(w_lo / power)
            adj_hi = _clamp01(w_hi / power)
        else:
            adj = None
            adj_lo, adj_hi = None, None
        abundance_confounded = power < pol.min_scrna_detection_power
        # High conditional expression (mean UMI >= 3, dropout < 5%) in a minority
        # of tumor cells proves genuine bimodal cellular heterogeneity.
        is_bimodal = (
            mean_pos >= 3.0
            and power >= 0.90
            and f_obs < pol.subclonal_threshold
            and w_hi < pol.clonal_threshold
        )
        if abundance_confounded:
            caveats.append("abundance_confounded_low_per_cell_detection_power")

    # Malignant-subcluster pseudobulk clonality (cell-weighted across powered subclusters)
    sub_frac: float | None = None
    n_powered_sub = 0
    n_expr_sub = 0
    if cell_to_subcluster:
        by_sub: dict[str, list[str]] = {}
        for bc in tumor_cell_counts:
            sub = cell_to_subcluster.get(bc)
            if sub:
                by_sub.setdefault(sub, []).append(bc)
        powered_cells = 0
        expr_cells = 0
        for sub, bcs in by_sub.items():
            if len(bcs) < pol.min_subcluster_cells:
                continue
            n_powered_sub += 1
            powered_cells += len(bcs)
            sub_umis = sum(float(tumor_cell_counts[bc]) for bc in bcs)
            if cell_total_umis:
                tot = sum(float(cell_total_umis.get(bc, 0.0)) for bc in bcs)
                ntpm = (1e6 * sub_umis / tot) if tot > 0 else 0.0
                expressed = ntpm >= pol.min_subcluster_ntpm and sub_umis >= 1.0
            else:
                expressed = sub_umis >= pol.min_subcluster_alt_reads
            if expressed:
                n_expr_sub += 1
                expr_cells += len(bcs)
        if powered_cells > 0:
            sub_frac = expr_cells / powered_cells

    # Hierarchical two-scale selection
    if (
        n_tumor >= pol.min_scrna_covered_cells
        and power is not None
        and power >= pol.min_scrna_detection_power
        and adj is not None
    ):
        selected = round(adj, 4)
        sel_lo = round(adj_lo, 4) if adj_lo is not None else None
        sel_hi = round(adj_hi, 4) if adj_hi is not None else None
        scale = "per_cell"
    elif sub_frac is not None and n_powered_sub >= 1:
        selected = round(_clamp01(sub_frac), 4)
        s_lo, s_hi = wilson_interval(n_expr_sub, n_powered_sub)
        sel_lo, sel_hi = round(s_lo, 4), round(s_hi, 4)
        scale = "subcluster_pseudobulk"
        caveats.append("used_malignant_subcluster_pseudobulk_due_to_per_cell_dropout")
    else:
        selected = None
        sel_lo, sel_hi = None, None
        scale = "unassessed"

    return ExpressionalClonality(
        assay_platform=assay_platform,
        measurement_level=measurement_level,
        n_tumor_cells=n_tumor,
        n_covered_tumor_cells=n_tumor,
        n_expressing_tumor_cells=n_pos,
        locus_coverage_fraction=1.0,
        detected_cell_fraction=round(f_obs, 4),
        detected_ci_lower=round(w_lo, 4),
        detected_ci_upper=round(w_hi, 4),
        mean_umi_or_depth_positive=round(mean_pos, 4) if n_pos > 0 else 0.0,
        detection_power=round(power, 4) if power is not None else None,
        dropout_adjusted_clonality=round(adj, 4) if adj is not None else None,
        dropout_ci_lower=round(adj_lo, 4) if adj_lo is not None else None,
        dropout_ci_upper=round(adj_hi, 4) if adj_hi is not None else None,
        expressing_subcluster_fraction=round(sub_frac, 4) if sub_frac is not None else None,
        n_powered_subclusters=n_powered_sub,
        n_expressing_subclusters=n_expr_sub,
        selected_clonality=selected,
        selected_ci_lower=sel_lo,
        selected_ci_upper=sel_hi,
        selected_scale=scale,
        abundance_confounded=abundance_confounded,
        is_bimodal_heterogeneous=is_bimodal,
        caveats=tuple(caveats),
    )


def estimate_scrna_allele_expressional_clonality(
    tumor_cell_alt_ref: dict[str, tuple[int, int]],
    *,
    n_total_tumor_cells: int | None = None,
    cell_to_subcluster: dict[str, str] | None = None,
    assay_platform: str = "long_read_ont",
    measurement_level: str = "allele_or_junction",
    min_cell_depth: int = 1,
    expected_allelic_fraction: float | None = None,
    policy: ClonalityPolicy | None = None,
) -> ExpressionalClonality:
    """Estimate scRNA expressional clonality at a mutation, splice junction, or isoform locus.

    Works identically for short-read scRNA (VarTrix / cellsnp-lite) and
    full-length long-read scRNA (PacBio MAS-Seq/Kinnex, ONT scNanoGPS/FLAMES/IsoQuant).
    Uses locus-covered malignant cells as denominator and adjusts for per-cell
    allelic dropout ``q_i = 1 - (1 - p_alt) ** d_i``.
    """
    pol = policy or ClonalityPolicy()
    if assay_platform in BULK_ASSAY_PLATFORMS:
        tot_alt = sum(max(0, int(a)) for a, _ in tumor_cell_alt_ref.values())
        tot_ref = sum(max(0, int(r)) for _, r in tumor_cell_alt_ref.values())
        bulk_af = (tot_alt / (tot_alt + tot_ref)) if (tot_alt + tot_ref) > 0 else None
        return ExpressionalClonality(
            assay_platform=assay_platform,
            measurement_level="bulk_only",
            bulk_expression_value=round(bulk_af, 4) if bulk_af is not None else None,
            selected_clonality=None,
            selected_scale="unassessed",
            caveats=("bulk_long_read_or_short_read_no_cell_barcodes",),
        )

    covered = {
        bc: (max(0, int(a)), max(0, int(r)))
        for bc, (a, r) in tumor_cell_alt_ref.items()
        if (max(0, int(a)) + max(0, int(r))) >= min_cell_depth
    }
    n_total = max(int(n_total_tumor_cells or len(tumor_cell_alt_ref)), len(covered))
    n_cov = len(covered)
    locus_cov_frac = (n_cov / n_total) if n_total > 0 else 0.0

    caveats: list[str] = []
    if n_cov == 0:
        caveats.append("zero_locus_covered_tumor_cells")
        return ExpressionalClonality(
            assay_platform=assay_platform,
            measurement_level=measurement_level,
            n_tumor_cells=n_total,
            n_covered_tumor_cells=0,
            n_expressing_tumor_cells=0,
            locus_coverage_fraction=0.0,
            abundance_confounded=True,
            selected_clonality=None,
            selected_scale="unassessed",
            caveats=tuple(caveats),
        )

    pos_cells = {bc: counts for bc, counts in covered.items() if counts[0] >= 1}
    n_pos = len(pos_cells)
    f_cov = n_pos / n_cov
    w_lo, w_hi = wilson_interval(n_pos, n_cov)

    # Estimate allelic expression fraction p_alt among mutant cells:
    # use multi-read positive cells (d >= 2) to avoid conditioning bias at d=1,
    # defaulting to balanced 0.5 (or 0.8 for homozygous/hemizygous/isoform override).
    if expected_allelic_fraction is not None:
        p_alt = max(0.05, min(0.99, float(expected_allelic_fraction)))
    else:
        deep_pos = [(a, a + r) for a, r in pos_cells.values() if (a + r) >= 2]
        if deep_pos:
            sum_a = sum(a for a, _ in deep_pos)
            sum_d = sum(d for _, d in deep_pos)
            # Zero-truncated binomial rough bias adjustment toward 0.5 when shallow
            p_alt = max(0.10, min(0.95, sum_a / sum_d))
        else:
            p_alt = 0.50

    depths = [a + r for a, r in covered.values()]
    per_cell_powers = [1.0 - ((1.0 - p_alt) ** d) for d in depths]
    mean_power = sum(per_cell_powers) / n_cov
    mean_depth_pos = (
        sum(a + r for a, r in pos_cells.values()) / n_pos if n_pos > 0 else 0.0
    )

    adj = _clamp01(f_cov / mean_power) if mean_power >= 0.05 else None
    adj_lo = _clamp01(w_lo / mean_power) if mean_power >= 0.05 else None
    adj_hi = _clamp01(w_hi / mean_power) if mean_power >= 0.05 else None

    coverage_limited = n_cov < pol.min_scrna_covered_cells or mean_power < pol.min_scrna_detection_power
    if n_cov < pol.min_scrna_covered_cells:
        caveats.append("low_covered_tumor_cell_count")
    if mean_power < pol.min_scrna_detection_power:
        caveats.append("low_allelic_detection_power")
    if locus_cov_frac < 0.10:
        caveats.append("sparse_locus_coverage_across_tumor_cells")

    is_bimodal = (
        mean_power >= 0.85
        and n_cov >= pol.min_scrna_covered_cells
        and f_cov < pol.subclonal_threshold
        and w_hi < pol.clonal_threshold
    )

    # Malignant-subcluster pseudobulk pooling of (alt, ref) reads
    sub_frac: float | None = None
    n_powered_sub = 0
    n_expr_sub = 0
    if cell_to_subcluster:
        by_sub: dict[str, list[tuple[int, int]]] = {}
        for bc, counts in covered.items():
            sub = cell_to_subcluster.get(bc)
            if sub:
                by_sub.setdefault(sub, []).append(counts)
        powered_cells = 0
        expr_cells = 0
        for sub, cell_list in by_sub.items():
            tot_alt = sum(a for a, _ in cell_list)
            tot_depth = sum(a + r for a, r in cell_list)
            # Require enough pooled depth so P(0 alt reads | p_alt) < 0.05
            sub_power = 1.0 - ((1.0 - p_alt) ** tot_depth) if tot_depth > 0 else 0.0
            if sub_power < 0.90:
                continue
            n_powered_sub += 1
            powered_cells += len(cell_list)
            if tot_alt >= pol.min_subcluster_alt_reads:
                n_expr_sub += 1
                expr_cells += len(cell_list)
        if powered_cells > 0:
            sub_frac = expr_cells / powered_cells

    if not coverage_limited and adj is not None:
        selected = round(adj, 4)
        sel_lo = round(adj_lo, 4) if adj_lo is not None else None
        sel_hi = round(adj_hi, 4) if adj_hi is not None else None
        scale = "per_cell"
    elif sub_frac is not None and n_powered_sub >= 1:
        selected = round(_clamp01(sub_frac), 4)
        s_lo, s_hi = wilson_interval(n_expr_sub, n_powered_sub)
        sel_lo, sel_hi = round(s_lo, 4), round(s_hi, 4)
        scale = "subcluster_pseudobulk"
        caveats.append("used_malignant_subcluster_pseudobulk_due_to_per_cell_dropout")
    else:
        selected = None
        sel_lo, sel_hi = None, None
        scale = "unassessed"

    return ExpressionalClonality(
        assay_platform=assay_platform,
        measurement_level=measurement_level,
        n_tumor_cells=n_total,
        n_covered_tumor_cells=n_cov,
        n_expressing_tumor_cells=n_pos,
        locus_coverage_fraction=round(locus_cov_frac, 4),
        detected_cell_fraction=round(f_cov, 4),
        detected_ci_lower=round(w_lo, 4),
        detected_ci_upper=round(w_hi, 4),
        mean_umi_or_depth_positive=round(mean_depth_pos, 4),
        estimated_allelic_fraction=round(p_alt, 4),
        detection_power=round(mean_power, 4),
        dropout_adjusted_clonality=round(adj, 4) if adj is not None else None,
        dropout_ci_lower=round(adj_lo, 4) if adj_lo is not None else None,
        dropout_ci_upper=round(adj_hi, 4) if adj_hi is not None else None,
        expressing_subcluster_fraction=round(sub_frac, 4) if sub_frac is not None else None,
        n_powered_subclusters=n_powered_sub,
        n_expressing_subclusters=n_expr_sub,
        selected_clonality=selected,
        selected_ci_lower=sel_lo,
        selected_ci_upper=sel_hi,
        selected_scale=scale,
        abundance_confounded=coverage_limited,
        is_bimodal_heterogeneous=is_bimodal,
        caveats=tuple(caveats),
    )


# ---------------------------------------------------------------------------
# Three-Stage Multi-Caller, Cross-Axis, Spatial & Longitudinal Combiners
# ---------------------------------------------------------------------------

def combine_genomic_callers(
    estimates: list[CallerClonalityEstimate] | tuple[CallerClonalityEstimate, ...],
    policy: ClonalityPolicy | None = None,
) -> GenomicClonality:
    """Reconcile multiple DNA clonality callers on the same tumor reads.

    Selects the primary estimate by configurable precedence (never multiplying
    correlated callers), flags cross-caller discordance, and widens the credible
    interval to the cross-caller envelope when callers disagree.
    """
    pol = policy or ClonalityPolicy()
    ests = tuple(estimates)
    if not ests:
        return GenomicClonality(caveats=("no_dna_clonality_estimates",))

    # Mitochondrial estimates never enter nuclear diploid CCF aggregation.
    if any(e.compartment == "mitochondrial" for e in ests):
        mt = next(e for e in ests if e.compartment == "mitochondrial")
        caveats = list(mt.notes) or ["mitochondrial_non_diploid_ccf_not_applicable"]
        return GenomicClonality(
            dna_ccf=None,
            compartment="mitochondrial",
            mt_heteroplasmy_tumor=mt.mt_heteroplasmy_tumor,
            mt_heteroplasmy_normal=mt.mt_heteroplasmy_normal,
            tumor_purity=mt.tumor_purity,
            primary_source=mt.source_tool,
            caller_estimates=ests,
            depth=mt.depth,
            caveats=tuple(sorted(set(caveats))),
        )

    precedence_map = {name: idx for idx, name in enumerate(pol.caller_precedence)}

    def _sort_key(est: CallerClonalityEstimate):
        return (
            est.ccf is None,
            precedence_map.get(est.source_tool.lower(), len(precedence_map)),
            est.source_tool.lower(),
        )

    ordered = sorted(ests, key=_sort_key)
    primary = ordered[0]
    valid_ccfs = [e.ccf for e in ordered if e.ccf is not None]
    caveats = list(primary.notes)

    spread: float | None = None
    discordant = False
    ci_lo = primary.ci_lower
    ci_hi = primary.ci_upper

    if len(valid_ccfs) >= 2:
        spread = round(max(valid_ccfs) - min(valid_ccfs), 4)
        if spread > pol.caller_discordance_threshold:
            discordant = True
            caveats.append("caller_discordance")
            lowers = [
                e.ci_lower if e.ci_lower is not None else e.ccf
                for e in ordered if e.ccf is not None
            ]
            uppers = [
                e.ci_upper if e.ci_upper is not None else e.ccf
                for e in ordered if e.ccf is not None
            ]
            ci_lo = round(_clamp01(min(lowers)), 4)
            ci_hi = round(_clamp01(max(uppers)), 4)

    max_depth = max((e.depth for e in ordered if e.depth is not None), default=None)
    return GenomicClonality(
        dna_ccf=round(_clamp01(primary.ccf), 4) if primary.ccf is not None else None,
        ci_lower=ci_lo,
        ci_upper=ci_hi,
        clonal_probability=primary.clonal_probability,
        multiplicity=primary.multiplicity,
        local_cn_total=primary.local_cn_total,
        local_cn_minor=primary.local_cn_minor,
        cluster_id=primary.cluster_id,
        tumor_purity=primary.tumor_purity,
        compartment=primary.compartment,
        primary_source=primary.source_tool,
        caller_estimates=tuple(ordered),
        caller_spread=spread,
        caller_discordance=discordant,
        depth=max_depth,
        caveats=tuple(sorted(set(caveats))),
    )


def _classify_status(
    value: float | None,
    ci_lower: float | None,
    ci_upper: float | None,
    *,
    clonal_probability: float | None = None,
    discordant: bool = False,
    abundance_confounded: bool = False,
    policy: ClonalityPolicy,
) -> str:
    if value is None:
        return "unassessed"
    if discordant and ci_lower is not None and ci_upper is not None:
        if ci_lower < policy.subclonal_threshold and ci_upper >= policy.clonal_threshold:
            return "uncertain"
    if (
        value >= policy.clonal_threshold
        and (clonal_probability is None or clonal_probability >= policy.min_clonal_probability)
        and (ci_lower is None or ci_lower >= policy.subclonal_threshold)
    ):
        return "clonal"
    if (
        value < policy.subclonal_threshold
        and (ci_upper is None or ci_upper < policy.clonal_threshold)
        and not abundance_confounded
    ):
        return "subclonal"
    return "uncertain"


def _estimate_variance_from_ci(value: float, lo: float | None, hi: float | None) -> float:
    if lo is not None and hi is not None and hi > lo:
        sigma = (hi - lo) / (2.0 * 1.96)
        return max(sigma * sigma, 1e-4)
    return 0.04  # fallback sd = 0.20 when CI is absent


def combine_sample_clonality(
    antigen_id: str,
    antigen_kind: str = "mutation",
    *,
    genomic: GenomicClonality | None = None,
    expressional: ExpressionalClonality | None = None,
    policy: ClonalityPolicy | None = None,
) -> AntigenClonality:
    """Combine GenomicClonality and ExpressionalClonality for one tumor sample.

    Follows exact probability semantics without double-counting:
    - Allele-/junction-/viral-specific scRNA estimates joint ``P(G=1, E=1) <= P(G=1)``
      -> precision-weighted fusion bounded above by DNA when consistent, or
      expression-silencing bottleneck when scRNA is powered and significantly < DNA.
    - Gene-level scRNA estimates marginal ``P(E_gene=1)`` -> ``min`` (default) or
      ``product`` (configurable) with bottleneck divergence flagged.
    - Non-templated targets (CTA/ERV/splice/isoform) rely on expressional clonality.
    - Viral targets may optionally apply a trial-configured ``integrated_driver_prior``
      explicitly labeled ``prior_assumed`` when unmeasured, overridden by scRNA.
    """
    pol = policy or ClonalityPolicy()
    kind_lower = str(antigen_kind or "mutation").strip().lower()
    kind_cfg = pol.for_antigen_kind(kind_lower)

    caveats: list[str] = []
    divergences: list[str] = []
    if genomic:
        caveats.extend(genomic.caveats)
        if genomic.caller_discordance:
            divergences.append("caller_discordance")
    if expressional:
        caveats.extend(expressional.caveats)

    # Mitochondrial compartment: never use nuclear diploid DNA CCF.
    if genomic is not None and genomic.compartment == "mitochondrial":
        if expressional is not None and expressional.selected_clonality is not None:
            eff = expressional.selected_clonality
            lo, hi = expressional.selected_ci_lower, expressional.selected_ci_upper
            status = _classify_status(
                eff, lo, hi,
                abundance_confounded=expressional.abundance_confounded and expressional.selected_scale != "subcluster_pseudobulk",
                policy=pol,
            )
            return AntigenClonality(
                antigen_id=antigen_id,
                antigen_kind=antigen_kind,
                genomic=genomic,
                expressional=expressional,
                effective_clonality=eff,
                ci_lower=lo,
                ci_upper=hi,
                clonality_status=status,
                combination_rule="mt_scrna_cell_level_only",
                divergences=tuple(sorted(set(divergences))),
                caveats=tuple(sorted(set(caveats))),
            )
        return AntigenClonality(
            antigen_id=antigen_id,
            antigen_kind=antigen_kind,
            genomic=genomic,
            expressional=expressional,
            effective_clonality=None,
            clonality_status="unassessed",
            combination_rule="mt_heteroplasmy_non_diploid",
            divergences=tuple(sorted(set(divergences))),
            caveats=tuple(sorted(set(caveats))),
        )

    dna_val = genomic.dna_ccf if genomic is not None else None
    dna_lo = genomic.ci_lower if genomic is not None else None
    dna_hi = genomic.ci_upper if genomic is not None else None
    expr_val = expressional.selected_clonality if expressional is not None else None
    expr_lo = expressional.selected_ci_lower if expressional is not None else None
    expr_hi = expressional.selected_ci_upper if expressional is not None else None

    eff: float | None = None
    ci_lo: float | None = None
    ci_hi: float | None = None
    rule = "unassessed"
    clonal_prob: float | None = genomic.clonal_probability if genomic else None
    abundance_conf = bool(
        expressional
        and expressional.abundance_confounded
        and expressional.selected_scale != "subcluster_pseudobulk"
    )

    if dna_val is not None and expr_val is not None:
        if expressional is not None and expressional.is_specific_to_alteration:
            # Joint event P(G=1, E=1) <= P(G=1). Check for silencing bottleneck first.
            dna_bound_lo = dna_lo if dna_lo is not None else max(0.0, dna_val - 0.15)
            expr_bound_hi = expr_hi if expr_hi is not None else min(1.0, expr_val + 0.15)
            if (dna_val - expr_val) > 0.25 and expr_bound_hi < dna_bound_lo:
                eff = expr_val
                ci_lo, ci_hi = expr_lo, expr_hi
                clonal_prob = None
                rule = "expression_silencing_bottleneck"
                divergences.append("transcriptional_silencing_or_nmd_suspected")
            else:
                # Precision-weighted fusion without double-counting, bounded by DNA upper CI.
                v_dna = _estimate_variance_from_ci(dna_val, dna_lo, dna_hi)
                v_expr = _estimate_variance_from_ci(expr_val, expr_lo, expr_hi)
                w_dna = 1.0 / v_dna
                w_expr = 1.0 / v_expr
                fused = (w_dna * dna_val + w_expr * expr_val) / (w_dna + w_expr)
                upper_cap = dna_hi if dna_hi is not None else min(1.0, dna_val + 0.10)
                eff = round(_clamp01(min(fused, upper_cap)), 4)
                fused_sd = math.sqrt(1.0 / (w_dna + w_expr))
                ci_lo = round(_clamp01(eff - 1.96 * fused_sd), 4)
                ci_hi = round(_clamp01(min(upper_cap, eff + 1.96 * fused_sd)), 4)
                rule = "precision_weighted_dna_and_allele_scrna"
        else:
            # Gene-level scRNA P(E_gene=1) combined with DNA CCF P(G=1)
            if pol.gene_level_combination == "product":
                eff = round(_clamp01(dna_val * expr_val), 4)
                ci_lo = round(_clamp01((dna_lo or dna_val) * (expr_lo or expr_val)), 4)
                ci_hi = round(_clamp01((dna_hi or dna_val) * (expr_hi or expr_val)), 4)
                rule = "dna_times_gene_level_scrna"
            else:
                eff = round(_clamp01(min(dna_val, expr_val)), 4)
                ci_lo = round(_clamp01(min(dna_lo if dna_lo is not None else dna_val,
                                           expr_lo if expr_lo is not None else expr_val)), 4)
                ci_hi = round(_clamp01(min(dna_hi if dna_hi is not None else dna_val,
                                           expr_hi if expr_hi is not None else expr_val)), 4)
                rule = "min_dna_and_gene_level_scrna"
            if (dna_val - expr_val) > 0.25:
                divergences.append("gene_level_expression_bottleneck")
                clonal_prob = None
    elif dna_val is not None:
        eff = dna_val
        ci_lo, ci_hi = dna_lo, dna_hi
        rule = "dna_only"
        abundance_conf = False
    elif expr_val is not None:
        eff = expr_val
        ci_lo, ci_hi = expr_lo, expr_hi
        rule = (
            "scrna_subcluster_pseudobulk"
            if expressional and expressional.selected_scale == "subcluster_pseudobulk"
            else "scrna_expressional_only"
        )
    elif kind_lower == "viral" and kind_cfg.get("integrated_driver_prior") is not None:
        prior = _clamp01(float(kind_cfg["integrated_driver_prior"]))
        eff = round(prior, 4)
        ci_lo, ci_hi = None, None
        rule = "viral_integrated_driver_prior_assumed"
        caveats.append("prior_assumed_not_measured")
    else:
        if kind_lower in {"sv", "fusion"}:
            caveats.append("sv_or_fusion_no_sv_ccf_or_scrna")
        elif kind_lower in {"cta", "erv", "splice", "isoform"}:
            caveats.append("bulk_rna_only_cannot_resolve_cellular_clonality")

    status = _classify_status(
        eff, ci_lo, ci_hi,
        clonal_probability=clonal_prob,
        discordant=bool(genomic and genomic.caller_discordance),
        abundance_confounded=abundance_conf,
        policy=pol,
    )
    return AntigenClonality(
        antigen_id=antigen_id,
        antigen_kind=antigen_kind,
        genomic=genomic,
        expressional=expressional,
        effective_clonality=eff,
        ci_lower=ci_lo,
        ci_upper=ci_hi,
        clonality_status=status,
        combination_rule=rule,
        divergences=tuple(sorted(set(divergences))),
        caveats=tuple(sorted(set(caveats))),
    )


def combine_spatial_and_longitudinal(
    antigen_id: str,
    antigen_kind: str,
    observations: list[SampleClonalityObservation] | tuple[SampleClonalityObservation, ...],
    policy: ClonalityPolicy | None = None,
) -> AntigenClonality:
    """Aggregate spatial regions within each timepoint, then directed timepoints.

    - Spatial regions at the same timepoint aggregate symmetrically.
    - Directed longitudinal timepoints (``T0 -> T1 -> ... -> T_latest``) anchor
      current clonality on ``T_latest``:
      - ``persistent``: present across multiple timepoints including ``T_latest``
      - ``emergent_latest``: absent earlier, present/clonal in ``T_latest``
        (retains full ``T_latest`` clonality!)
      - ``lost_in_latest``: present earlier, absent at adequate depth in ``T_latest``
        (status becomes ``lost`` with ``effective_clonality=0.0``!)
      - ``expanding`` / ``contracting``: monotonic trend across >=2 timepoints
      - ``latest_unassessed``: ``T_latest`` lacked adequate coverage/power
    """
    pol = policy or ClonalityPolicy()
    if not observations:
        return AntigenClonality(antigen_id=antigen_id, antigen_kind=antigen_kind)

    # Group by (timepoint_index, timepoint_label)
    by_tp: dict[tuple[int, str], list[SampleClonalityObservation]] = {}
    for obs in observations:
        by_tp.setdefault((int(obs.timepoint_index), str(obs.timepoint_label)), []).append(obs)

    ordered_tps = sorted(by_tp.keys(), key=lambda x: (x[0], x[1]))
    per_tp_spatial: list[SpatialClonality] = []
    per_tp_clonality: list[tuple[tuple[int, str], float | None, bool | None, AntigenClonality]] = []

    for tp_key in ordered_tps:
        tp_idx, tp_label = tp_key
        tp_obs = by_tp[tp_key]
        # Group by physical region within this timepoint (multiple sequencing runs
        # on the SAME physical region/timepoint are reconciled, not double-counted).
        by_region: dict[str, list[SampleClonalityObservation]] = {}
        for o in tp_obs:
            by_region.setdefault(o.region_label or o.sample_id, []).append(o)

        region_results: list[tuple[str, AntigenClonality, bool | None, bool]] = []
        for reg_label, reg_obs in sorted(by_region.items()):
            # Merge DNA caller estimates and pick highest-power scRNA across runs of same region.
            all_dna_callers: list[CallerClonalityEstimate] = []
            mt_genomic: GenomicClonality | None = None
            best_expr: ExpressionalClonality | None = None
            max_depth: int | None = None
            explicit_present: list[bool] = []
            for o in reg_obs:
                if o.depth is not None:
                    max_depth = o.depth if max_depth is None else max(max_depth, o.depth)
                if o.present is not None:
                    explicit_present.append(bool(o.present))
                if o.genomic is not None:
                    if o.genomic.depth is not None:
                        max_depth = o.genomic.depth if max_depth is None else max(max_depth, o.genomic.depth)
                    if o.genomic.compartment == "mitochondrial":
                        mt_genomic = o.genomic
                    elif o.genomic.caller_estimates:
                        all_dna_callers.extend(o.genomic.caller_estimates)
                    elif o.genomic.dna_ccf is not None:
                        all_dna_callers.append(CallerClonalityEstimate(
                            source_tool=o.genomic.primary_source or "canonical",
                            ccf=o.genomic.dna_ccf,
                            ci_lower=o.genomic.ci_lower,
                            ci_upper=o.genomic.ci_upper,
                            clonal_probability=o.genomic.clonal_probability,
                            multiplicity=o.genomic.multiplicity,
                            local_cn_total=o.genomic.local_cn_total,
                            local_cn_minor=o.genomic.local_cn_minor,
                            cluster_id=o.genomic.cluster_id,
                            tumor_purity=o.genomic.tumor_purity,
                            depth=o.genomic.depth,
                            sample_id=o.sample_id,
                        ))
                if o.expressional is not None:
                    if best_expr is None or (
                        (o.expressional.detection_power or 0.0, o.expressional.n_covered_tumor_cells)
                        > (best_expr.detection_power or 0.0, best_expr.n_covered_tumor_cells)
                    ):
                        best_expr = o.expressional

            merged_genomic = (
                mt_genomic
                if mt_genomic is not None
                else (combine_genomic_callers(all_dna_callers, policy=pol) if all_dna_callers else None)
            )
            reg_clon = combine_sample_clonality(
                antigen_id, antigen_kind,
                genomic=merged_genomic,
                expressional=best_expr,
                policy=pol,
            )
            if explicit_present:
                reg_present: bool | None = any(explicit_present)
                reg_powered = True
            elif reg_clon.effective_clonality is not None:
                reg_present = reg_clon.effective_clonality > 0.02
                reg_powered = True
            elif max_depth is not None and max_depth >= pol.min_latest_depth_for_loss:
                reg_present = False
                reg_powered = True
            else:
                reg_present = None
                reg_powered = False
            region_results.append((reg_label, reg_clon, reg_present, reg_powered))

        n_sampled = len(region_results)
        powered_regs = [r for r in region_results if r[3]]
        n_powered = len(powered_regs)
        n_present = sum(1 for r in powered_regs if r[2] is True)
        reg_frac = (n_present / n_powered) if n_powered > 0 else None

        clon_vals = [
            (r[1].effective_clonality if r[1].effective_clonality is not None else 0.0)
            for r in powered_regs
            if (r[1].effective_clonality is not None or r[2] is False)
        ]
        min_reg = round(min(clon_vals), 4) if clon_vals else None
        mean_reg = round(sum(clon_vals) / len(clon_vals), 4) if clon_vals else None
        is_ubiq = (n_powered > 0 and n_present == n_powered) if n_powered > 0 else None

        spatial = SpatialClonality(
            timepoint_index=tp_idx,
            timepoint_label=tp_label,
            regions_sampled=n_sampled,
            regions_powered=n_powered,
            regions_present=n_present,
            region_fraction=round(reg_frac, 4) if reg_frac is not None else None,
            min_region_clonality=min_reg,
            mean_region_clonality=mean_reg,
            is_spatially_ubiquitous=is_ubiq,
            region_labels=tuple(r[0] for r in region_results),
        )
        per_tp_spatial.append(spatial)

        # Pick representative region record (highest clonality or first) and apply spatial aggregation.
        rep = max(region_results, key=lambda r: (r[1].effective_clonality is not None, r[1].effective_clonality or -1.0))[1]
        if n_sampled > 1 and clon_vals:
            tp_eff = min_reg if pol.multi_region_aggregation == "min" else mean_reg
        else:
            tp_eff = rep.effective_clonality if rep.effective_clonality is not None else (0.0 if (n_powered > 0 and n_present == 0) else None)
        tp_present: bool | None = (n_present > 0) if n_powered > 0 else None
        per_tp_clonality.append((tp_key, tp_eff, tp_present, rep))

    # Single timepoint -> attach spatial summary directly.
    latest_tp_key, latest_eff, latest_present, latest_rep = per_tp_clonality[-1]
    latest_spatial = per_tp_spatial[-1]

    if len(ordered_tps) == 1:
        status = latest_rep.clonality_status
        caveats = list(latest_rep.caveats)
        rule = latest_rep.combination_rule
        if latest_spatial.regions_sampled > 1 and latest_eff is not None:
            status = _classify_status(latest_eff, latest_rep.ci_lower, latest_rep.ci_upper, policy=pol)
            if latest_spatial.is_spatially_ubiquitous is False:
                caveats.append("spatially_heterogeneous_across_regions")
            rule = f"multi_region_{pol.multi_region_aggregation}_{rule}"
        return AntigenClonality(
            antigen_id=antigen_id,
            antigen_kind=antigen_kind,
            genomic=latest_rep.genomic,
            expressional=latest_rep.expressional,
            spatial=latest_spatial if latest_spatial.regions_sampled > 1 else None,
            longitudinal=None,
            effective_clonality=latest_eff,
            ci_lower=latest_rep.ci_lower,
            ci_upper=latest_rep.ci_upper,
            clonality_status=status,
            combination_rule=rule,
            divergences=latest_rep.divergences,
            caveats=tuple(sorted(set(caveats))),
        )

    # Directed longitudinal trajectory classification across T0 -> ... -> T_latest.
    tp_labels = tuple(k[1] for k in ordered_tps)
    tp_powered = sum(1 for _, _, pres, _ in per_tp_clonality if pres is not None)
    tp_present_count = sum(1 for _, _, pres, _ in per_tp_clonality if pres is True)
    earlier_present = any(pres is True for _, _, pres, _ in per_tp_clonality[:-1])
    earlier_absent = any(pres is False for _, _, pres, _ in per_tp_clonality[:-1])
    measured_series = [eff for _, eff, _, _ in per_tp_clonality if eff is not None]

    if latest_present is None:
        trajectory = "latest_unassessed"
    elif latest_present is False and earlier_present:
        trajectory = "lost_in_latest"
    elif latest_present is True and not earlier_present and earlier_absent:
        trajectory = "emergent_latest"
    elif (
        latest_present is True
        and len(measured_series) >= 2
        and (measured_series[-1] - measured_series[0]) >= 0.25
        and all(b >= a - 0.05 for a, b in zip(measured_series, measured_series[1:]))
    ):
        trajectory = "expanding"
    elif (
        len(measured_series) >= 2
        and (measured_series[0] - measured_series[-1]) >= 0.25
        and all(b <= a + 0.05 for a, b in zip(measured_series, measured_series[1:]))
    ):
        trajectory = "contracting"
    elif latest_present is True and earlier_present:
        trajectory = "persistent"
    else:
        trajectory = "single_timepoint"

    longitudinal = LongitudinalClonality(
        timepoints_ordered=tp_labels,
        timepoints_powered=tp_powered,
        timepoints_present=tp_present_count,
        latest_timepoint=latest_tp_key[1],
        present_in_latest=latest_present,
        latest_clonality=latest_eff,
        trajectory=trajectory,
        per_timepoint_spatial=tuple(per_tp_spatial),
    )

    caveats = list(latest_rep.caveats)
    divergences = list(latest_rep.divergences)
    if trajectory == "lost_in_latest":
        eff = 0.0
        status = "lost"
        rule = "longitudinal_lost_in_latest"
        divergences.append("lost_in_latest_timepoint")
    elif trajectory == "latest_unassessed":
        eff = None
        status = "unassessed"
        rule = "longitudinal_latest_unassessed"
        caveats.append("latest_timepoint_lacks_adequate_depth_or_coverage")
    else:
        eff = latest_eff
        status = _classify_status(eff, latest_rep.ci_lower, latest_rep.ci_upper, policy=pol)
        rule = f"longitudinal_{trajectory}_{latest_rep.combination_rule}"

    return AntigenClonality(
        antigen_id=antigen_id,
        antigen_kind=antigen_kind,
        genomic=latest_rep.genomic,
        expressional=latest_rep.expressional,
        spatial=latest_spatial if latest_spatial.regions_sampled > 1 else None,
        longitudinal=longitudinal,
        effective_clonality=eff,
        ci_lower=0.0 if trajectory == "lost_in_latest" else latest_rep.ci_lower,
        ci_upper=0.0 if trajectory == "lost_in_latest" else latest_rep.ci_upper,
        clonality_status=status,
        combination_rule=rule,
        divergences=tuple(sorted(set(divergences))),
        caveats=tuple(sorted(set(caveats))),
    )


# ---------------------------------------------------------------------------
# Pure-Python Upstream Parsers (Zero External Dependencies)
# ---------------------------------------------------------------------------

def _variant_key(contig: str, pos: int, ref: str, alt: str) -> str:
    return f"{normalize_contig(contig)}:{int(pos)}:{str(ref).upper()}:{str(alt).upper()}"


def load_pyclone_vi(path: str | Path) -> dict[str, list[CallerClonalityEstimate]]:
    """Parse PyClone-VI output TSV (``mutation_id``, ``sample_id``, ``cluster_id``,
    ``cellular_prevalence``, ``cellular_prevalence_std``, ``cluster_assignment_prob``).
    """
    out: dict[str, list[CallerClonalityEstimate]] = {}
    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            mut_id = str(row.get("mutation_id") or "").strip()
            if not mut_id:
                continue
            ccf = _safe_float(row.get("cellular_prevalence"))
            std = _safe_float(row.get("cellular_prevalence_std"))
            ci_lo = _clamp01(ccf - 1.96 * std) if (ccf is not None and std is not None) else None
            ci_hi = _clamp01(ccf + 1.96 * std) if (ccf is not None and std is not None) else None
            cluster_id = str(row.get("cluster_id") or "").strip() or None
            sample_id = str(row.get("sample_id") or "").strip()
            out.setdefault(mut_id, []).append(CallerClonalityEstimate(
                source_tool="pyclone_vi",
                ccf=round(_clamp01(ccf), 4) if ccf is not None else None,
                ci_lower=round(ci_lo, 4) if ci_lo is not None else None,
                ci_upper=round(ci_hi, 4) if ci_hi is not None else None,
                cluster_id=cluster_id,
                sample_id=sample_id,
            ))
    return out


def load_purple_purity_tsv(path: str | Path) -> dict[str, Any]:
    """Parse Hartwig PURPLE ``*.purple.purity.tsv`` global purity/ploidy/WGD metadata."""
    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        row = next(reader, None)
        if not row:
            return {}
        return {
            "purity": _safe_float(row.get("purity")),
            "ploidy": _safe_float(row.get("ploidy")),
            "gender": str(row.get("gender") or row.get("sex") or "").strip().lower() or None,
            "whole_genome_duplication": str(row.get("wholeGenomeDuplication") or "").strip().lower() == "true",
            "qc_status": str(row.get("qcStatus") or "").strip(),
        }


def load_purple_vcf(
    vcf_path: str | Path,
    *,
    purity_tsv_path: str | Path | None = None,
    sample_id: str = "",
) -> dict[str, list[CallerClonalityEstimate]]:
    """Parse Hartwig PURPLE somatic VCF INFO fields (``PURPLE_CN``, ``PURPLE_MACN``,
    ``PURPLE_VCN``, ``SUBCL``, ``PURPLE_AF``).
    """
    purity_meta = load_purple_purity_tsv(purity_tsv_path) if purity_tsv_path else {}
    global_purity = purity_meta.get("purity")
    out: dict[str, list[CallerClonalityEstimate]] = {}
    with _open_text(vcf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos_s, _vid, ref, alt, _qual, _flt, info_str = parts[:8]
            pos = _safe_int(pos_s)
            if pos is None:
                continue
            info: dict[str, str] = {}
            for item in info_str.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info[k] = v
                elif item:
                    info[item] = "true"
            key = _variant_key(chrom, pos, ref, alt)
            if is_mitochondrial_contig(chrom):
                af = _safe_float(info.get("PURPLE_AF") or info.get("AF"))
                out.setdefault(key, []).append(CallerClonalityEstimate(
                    source_tool="purple",
                    ccf=None,
                    compartment="mitochondrial",
                    mt_heteroplasmy_tumor=af,
                    tumor_purity=global_purity,
                    sample_id=sample_id,
                    notes=("mitochondrial_non_diploid_ccf_not_applicable",),
                ))
                continue

            vcn = _safe_float(info.get("PURPLE_VCN"))
            cn = _safe_float(info.get("PURPLE_CN"))
            macn = _safe_float(info.get("PURPLE_MACN"))
            subcl = _safe_float(info.get("SUBCL"))
            p_af = _safe_float(info.get("PURPLE_AF"))
            cn_minor = round(max(0.0, cn - macn), 4) if (cn is not None and macn is not None) else None
            # In PURPLE, PURPLE_VCN is purity-adjusted variant copy number:
            # VCN >= 1.0 corresponds to clonal (multiplicity = round(VCN)), while
            # 0 <= VCN < 1.0 directly reflects CCF on a single copy!
            if vcn is not None:
                ccf = _clamp01(vcn) if vcn < 1.0 else 1.0
                mult = max(1.0, float(round(vcn))) if vcn >= 1.0 else 1.0
            else:
                ccf = None
                mult = None
            clonal_prob = _clamp01(1.0 - subcl) if subcl is not None else None
            out.setdefault(key, []).append(CallerClonalityEstimate(
                source_tool="purple",
                ccf=round(ccf, 4) if ccf is not None else None,
                clonal_probability=round(clonal_prob, 4) if clonal_prob is not None else None,
                multiplicity=mult,
                local_cn_total=cn,
                local_cn_minor=cn_minor,
                tumor_purity=global_purity,
                vaf=p_af,
                sample_id=sample_id,
            ))
    return out


def load_facets_maf(path: str | Path) -> dict[str, list[CallerClonalityEstimate]]:
    """Parse FACETS / facets-suite annotated MAF/TSV (``ccf_expected_copies``,
    ``ccf_expected_copies_lower/upper``, ``tcn.em``, ``lcn.em``, ``clonality``, ``purity``).
    """
    out: dict[str, list[CallerClonalityEstimate]] = {}
    with _open_text(path) as handle:
        # Skip comment lines starting with '#'
        lines = [ln for ln in handle if ln.strip() and not ln.startswith("#")]
    if not lines:
        return out
    reader = csv.DictReader(lines, delimiter="\t")
    for row in reader:
        chrom = row.get("Chromosome") or row.get("chrom") or row.get("chr") or ""
        pos = _safe_int(row.get("Start_Position") or row.get("start") or row.get("pos"))
        ref = row.get("Reference_Allele") or row.get("ref") or ""
        alt = row.get("Tumor_Seq_Allele2") or row.get("alt") or ""
        if not chrom or pos is None:
            continue
        key = _variant_key(chrom, pos, ref, alt)
        sample_id = str(row.get("Tumor_Sample_Barcode") or row.get("sample_id") or "").strip()
        t_alt = _safe_int(row.get("t_alt_count"))
        t_ref = _safe_int(row.get("t_ref_count"))
        depth = (t_alt + t_ref) if (t_alt is not None and t_ref is not None) else None
        if is_mitochondrial_contig(chrom):
            het = (t_alt / depth) if (t_alt is not None and depth and depth > 0) else None
            out.setdefault(key, []).append(CallerClonalityEstimate(
                source_tool="facets",
                ccf=None,
                compartment="mitochondrial",
                mt_heteroplasmy_tumor=het,
                depth=depth,
                alt_reads=t_alt,
                sample_id=sample_id,
                notes=("mitochondrial_non_diploid_ccf_not_applicable",),
            ))
            continue
        ccf = _safe_float(row.get("ccf_expected_copies") or row.get("ccf_Mcopies"))
        lo = _safe_float(row.get("ccf_expected_copies_lower"))
        hi = _safe_float(row.get("ccf_expected_copies_upper"))
        tcn = _safe_float(row.get("tcn.em") or row.get("tcn"))
        lcn = _safe_float(row.get("lcn.em") or row.get("lcn"))
        purity = _safe_float(row.get("purity"))
        clonality_lbl = str(row.get("clonality") or "").strip().upper()
        notes = (f"facets_label:{clonality_lbl}",) if clonality_lbl else ()
        out.setdefault(key, []).append(CallerClonalityEstimate(
            source_tool="facets",
            ccf=round(_clamp01(ccf), 4) if ccf is not None else None,
            ci_lower=round(_clamp01(lo), 4) if lo is not None else None,
            ci_upper=round(_clamp01(hi), 4) if hi is not None else None,
            local_cn_total=tcn,
            local_cn_minor=lcn,
            tumor_purity=purity,
            depth=depth,
            alt_reads=t_alt,
            vaf=(t_alt / depth) if (t_alt is not None and depth and depth > 0) else None,
            sample_id=sample_id,
            notes=notes,
        ))
    return out


def load_purecn_variants(path: str | Path) -> dict[str, list[CallerClonalityEstimate]]:
    """Parse PureCN ``*_variants.csv`` output."""
    out: dict[str, list[CallerClonalityEstimate]] = {}
    with _open_text(path) as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            chrom = row.get("chr") or row.get("Chromosome") or ""
            pos = _safe_int(row.get("start") or row.get("Start_Position"))
            ref = row.get("REF") or row.get("ref") or ""
            alt = row.get("ALT") or row.get("alt") or ""
            if not chrom or pos is None:
                continue
            key = _variant_key(chrom, pos, ref, alt)
            sample_id = str(row.get("Sampleid") or row.get("sample_id") or "").strip()
            depth = _safe_int(row.get("DEPTH"))
            vaf = _safe_float(row.get("AR"))
            if is_mitochondrial_contig(chrom):
                out.setdefault(key, []).append(CallerClonalityEstimate(
                    source_tool="purecn",
                    ccf=None,
                    compartment="mitochondrial",
                    mt_heteroplasmy_tumor=vaf,
                    depth=depth,
                    sample_id=sample_id,
                    notes=("mitochondrial_non_diploid_ccf_not_applicable",),
                ))
                continue
            ccf = _safe_float(row.get("CELLFRACTION") or row.get("CCF") or row.get("ML.CCF"))
            lo = _safe_float(row.get("CELLFRACTION.95.LOWER"))
            hi = _safe_float(row.get("CELLFRACTION.95.UPPER"))
            cn_tot = _safe_float(row.get("ML.C"))
            cn_min = _safe_float(row.get("ML.M.SEGMENT"))
            mult = _safe_float(row.get("ML.M"))
            out.setdefault(key, []).append(CallerClonalityEstimate(
                source_tool="purecn",
                ccf=round(_clamp01(ccf), 4) if ccf is not None else None,
                ci_lower=round(_clamp01(lo), 4) if lo is not None else None,
                ci_upper=round(_clamp01(hi), 4) if hi is not None else None,
                multiplicity=mult,
                local_cn_total=cn_tot,
                local_cn_minor=cn_min,
                depth=depth,
                vaf=vaf,
                sample_id=sample_id,
            ))
    return out


def load_sequenza_segments(
    segments_path: str | Path,
    *,
    confints_path: str | Path | None = None,
) -> tuple[list[CopyNumberSegment], float | None, float | None]:
    """Parse Sequenza ``*_segments.txt`` (and optional ``*_confints_CP.txt``)."""
    segments: list[CopyNumberSegment] = []
    with _open_text(segments_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chrom = str(row.get("chromosome") or row.get("chrom") or row.get("chr") or "").strip().strip('"')
            start = _safe_int(row.get("start.pos") or row.get("start"))
            end = _safe_int(row.get("end.pos") or row.get("end"))
            cnt = _safe_float(row.get("CNt") or row.get("tcn.em") or row.get("cn_total"))
            if not chrom or start is None or end is None or cnt is None:
                continue
            a_cn = _safe_float(row.get("A") or row.get("cn_major"))
            b_cn = _safe_float(row.get("B") or row.get("lcn.em") or row.get("cn_minor"))
            segments.append(CopyNumberSegment(
                contig=chrom,
                start=start,
                end=end,
                cn_total=cnt,
                cn_major=a_cn,
                cn_minor=b_cn,
                source_tool="sequenza",
            ))

    purity: float | None = None
    ploidy: float | None = None
    if confints_path:
        with _open_text(confints_path) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
            # Sequenza confints_CP.txt row 0 is lower CI, row 1 is point estimate, row 2 is upper CI.
            target = rows[1] if len(rows) >= 2 else (rows[0] if rows else None)
            if target:
                purity = _safe_float(target.get("cellularity"))
                ploidy = _safe_float(target.get("ploidy.estimate"))
    return segments, purity, ploidy


def _read_matrix_market_int(path: str | Path) -> dict[tuple[int, int], int]:
    """Pure-Python 1-based MatrixMarket integer coordinate matrix reader."""
    entries: dict[tuple[int, int], int] = {}
    with _open_text(path) as handle:
        dims_seen = False
        for line in handle:
            s = line.strip()
            if not s or s.startswith("%"):
                continue
            if not dims_seen:
                dims_seen = True
                continue
            parts = s.split()
            if len(parts) < 3:
                continue
            r, c = int(parts[0]) - 1, int(parts[1]) - 1
            v = int(round(float(parts[2])))
            if v != 0:
                entries[(r, c)] = v
    return entries


def _read_vcf_variant_keys(vcf_path: str | Path) -> list[str]:
    keys: list[str] = []
    with _open_text(vcf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 5:
                continue
            pos = _safe_int(parts[1])
            if pos is not None:
                keys.append(_variant_key(parts[0], pos, parts[3], parts[4]))
    return keys


def load_vartrix_matrices(
    alt_mtx: str | Path,
    ref_mtx: str | Path,
    barcodes_path: str | Path,
    variants_vcf: str | Path,
    *,
    tumor_barcodes: set[str] | None = None,
) -> dict[str, dict[str, tuple[int, int]]]:
    """Parse VarTrix ``alt.mtx`` / ``ref.mtx`` + ``barcodes.tsv`` + ``variants.vcf``.

    Returns ``variant_key -> {barcode: (alt_count, ref_count)}`` filtered to
    ``tumor_barcodes`` when provided.
    """
    with _open_text(barcodes_path) as handle:
        barcodes = [ln.strip().split("\t")[0] for ln in handle if ln.strip()]
    variants = _read_vcf_variant_keys(variants_vcf)
    alts = _read_matrix_market_int(alt_mtx)
    refs = _read_matrix_market_int(ref_mtx)

    coords = set(alts.keys()) | set(refs.keys())
    result: dict[str, dict[str, tuple[int, int]]] = {v: {} for v in variants}
    for r, c in coords:
        if r < 0 or r >= len(variants) or c < 0 or c >= len(barcodes):
            continue
        bc = barcodes[c]
        if tumor_barcodes is not None and bc not in tumor_barcodes:
            continue
        result[variants[r]][bc] = (alts.get((r, c), 0), refs.get((r, c), 0))
    return result


def load_cellsnp_lite_dir(
    directory: str | Path,
    *,
    tumor_barcodes: set[str] | None = None,
) -> dict[str, dict[str, tuple[int, int]]]:
    """Parse cellsnp-lite directory (``cellSNP.tag.AD.mtx``, ``cellSNP.tag.DP.mtx``,
    ``cellSNP.samples.tsv``, ``cellSNP.base.vcf`` / ``.vcf.gz``).
    """
    d = Path(directory)
    ad_path = d / "cellSNP.tag.AD.mtx"
    dp_path = d / "cellSNP.tag.DP.mtx"
    samples_path = d / "cellSNP.samples.tsv"
    vcf_path = d / "cellSNP.base.vcf.gz"
    if not vcf_path.exists():
        vcf_path = d / "cellSNP.base.vcf"

    with _open_text(samples_path) as handle:
        barcodes = [ln.strip().split("\t")[0] for ln in handle if ln.strip()]
    variants = _read_vcf_variant_keys(vcf_path)
    ads = _read_matrix_market_int(ad_path)
    dps = _read_matrix_market_int(dp_path)

    coords = set(ads.keys()) | set(dps.keys())
    result: dict[str, dict[str, tuple[int, int]]] = {v: {} for v in variants}
    for r, c in coords:
        if r < 0 or r >= len(variants) or c < 0 or c >= len(barcodes):
            continue
        bc = barcodes[c]
        if tumor_barcodes is not None and bc not in tumor_barcodes:
            continue
        alt_c = ads.get((r, c), 0)
        dp_c = dps.get((r, c), alt_c)
        ref_c = max(0, dp_c - alt_c)
        result[variants[r]][bc] = (alt_c, ref_c)
    return result


def load_isoquant_read_info(
    path: str | Path,
    *,
    tumor_barcodes: set[str] | None = None,
) -> tuple[dict[str, dict[str, tuple[int, int]]], bool]:
    """Parse IsoQuant ``read_info.tsv[.gz]`` into per-isoform ``{barcode: (iso_umis, other_gene_umis)}``.

    Returns ``(isoform_counts, is_single_cell)``. When ``barcode`` is ``'.'``
    across all rows (IsoQuant bulk mode), ``is_single_cell`` is ``False`` so
    bulk long-read RNA is never confused with single-cell expressional clonality.
    """
    # Track unique (gene, isoform, barcode, umi_or_read_id)
    gene_iso_barcodes: dict[str, dict[str, set[str]]] = {}  # gene -> {barcode -> set of molecule ids}
    iso_barcodes: dict[str, dict[str, set[str]]] = {}       # isoform -> {barcode -> set of molecule ids}
    iso_to_gene: dict[str, str] = {}
    saw_real_barcode = False

    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            bc = str(row.get("barcode") or "").strip()
            gene_id = str(row.get("gene_id") or "").strip()
            iso_id = str(row.get("isoform_id") or "").strip()
            if not gene_id or not iso_id or iso_id == ".":
                continue
            if bc and bc != ".":
                saw_real_barcode = True
                if tumor_barcodes is not None and bc not in tumor_barcodes:
                    continue
            else:
                bc = "__BULK__"
            umi = str(row.get("umi") or "").strip()
            mol_id = umi if (umi and umi != ".") else str(row.get("read_id") or "")
            iso_to_gene[iso_id] = gene_id
            gene_iso_barcodes.setdefault(gene_id, {}).setdefault(bc, set()).add(mol_id)
            iso_barcodes.setdefault(iso_id, {}).setdefault(bc, set()).add(mol_id)

    result: dict[str, dict[str, tuple[int, int]]] = {}
    for iso_id, bc_map in iso_barcodes.items():
        gene_id = iso_to_gene[iso_id]
        gene_bcs = gene_iso_barcodes.get(gene_id, {})
        all_bcs = set(gene_bcs.keys()) | set(bc_map.keys())
        iso_dict: dict[str, tuple[int, int]] = {}
        for bc in all_bcs:
            iso_c = len(bc_map.get(bc, ()))
            tot_c = len(gene_bcs.get(bc, ()))
            other_c = max(0, tot_c - iso_c)
            iso_dict[bc] = (iso_c, other_c)
        result[iso_id] = iso_dict
    return result, saw_real_barcode


def load_flames_transcript_counts(
    path: str | Path,
    *,
    tumor_barcodes: set[str] | None = None,
) -> dict[str, dict[str, int]]:
    """Parse FLAMES ``transcript_count.csv[.gz]`` (rows = isoforms, cols = barcodes)."""
    result: dict[str, dict[str, int]] = {}
    with _open_text(path) as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if not header or len(header) < 2:
            return result
        # FLAMES header is either [transcript_id, gene_id, BC1, BC2, ...] or [transcript_id, BC1, BC2, ...]
        start_col = 2 if len(header) >= 3 and header[1].lower() in {"gene_id", "gene", "gene_name"} else 1
        barcodes = [h.strip() for h in header[start_col:]]
        for row in reader:
            if len(row) <= start_col:
                continue
            iso_id = row[0].strip()
            if not iso_id:
                continue
            counts: dict[str, int] = {}
            for bc, raw_val in zip(barcodes, row[start_col:]):
                if tumor_barcodes is not None and bc not in tumor_barcodes:
                    continue
                c = _safe_int(raw_val)
                counts[bc] = max(0, c or 0)
            result[iso_id] = counts
    return result


def load_canonical_clonality_table(
    path: str | Path,
    policy: ClonalityPolicy | None = None,
) -> dict[str, AntigenClonality]:
    """Load Vaxrank's canonical TSV/CSV/JSON multi-sample clonality table.

    Supports both single-sample and multi-region / longitudinal rows with
    columns: ``antigen_id``, ``antigen_kind``, ``sample_id``, ``timepoint_index``,
    ``timepoint_label``, ``region_label``, ``source_tool``, ``dna_ccf``,
    ``ci_lower``, ``ci_upper``, ``clonal_probability``, ``multiplicity``,
    ``local_cn_total``, ``local_cn_minor``, ``cluster_id``, ``tumor_purity``,
    ``depth``, ``alt_reads``, ``contig``, ``compartment``,
    ``mt_heteroplasmy_tumor``, ``mt_heteroplasmy_normal``, ``present``,
    ``scrna_clonality``, ``scrna_power``, ``scrna_platform``, ``scrna_level``.
    """
    p = Path(path)
    if p.suffix.lower() == ".json":
        payload = json.loads(p.read_text(encoding="utf-8"))
        rows = payload if isinstance(payload, list) else payload.get("records", [])
    else:
        with _open_text(p) as handle:
            sample_line = handle.readline()
            handle.seek(0)
            delim = "\t" if "\t" in sample_line else ","
            rows = list(csv.DictReader(handle, delimiter=delim))

    pol = policy or ClonalityPolicy()
    grouped: dict[str, tuple[str, list[SampleClonalityObservation]]] = {}
    for row in rows:
        antigen_id = str(row.get("antigen_id") or row.get("mutation_id") or "").strip()
        if not antigen_id:
            continue
        antigen_kind = str(row.get("antigen_kind") or "mutation").strip()
        sample_id = str(row.get("sample_id") or "sample_0").strip()
        tp_idx = _safe_int(row.get("timepoint_index"))
        if tp_idx is None:
            tp_idx = 0
        tp_label = str(row.get("timepoint_label") or f"T{tp_idx}").strip()
        region_label = str(row.get("region_label") or "primary").strip()
        contig = str(row.get("contig") or "").strip()
        compartment = str(row.get("compartment") or ("mitochondrial" if is_mitochondrial_contig(contig) else "nuclear")).strip().lower()
        depth = _safe_int(row.get("depth"))
        alt_reads = _safe_int(row.get("alt_reads"))
        ccf = _safe_float(row.get("dna_ccf") if "dna_ccf" in row else row.get("ccf"))
        genomic: GenomicClonality | None = None
        if compartment == "mitochondrial":
            genomic = GenomicClonality(
                dna_ccf=None,
                compartment="mitochondrial",
                mt_heteroplasmy_tumor=_safe_float(row.get("mt_heteroplasmy_tumor")),
                mt_heteroplasmy_normal=_safe_float(row.get("mt_heteroplasmy_normal")),
                depth=depth,
                primary_source=str(row.get("source_tool") or "canonical"),
                caveats=("mitochondrial_non_diploid_ccf_not_applicable",),
            )
        elif ccf is not None:
            est = CallerClonalityEstimate(
                source_tool=str(row.get("source_tool") or "canonical").strip().lower(),
                ccf=_clamp01(ccf),
                ci_lower=_safe_float(row.get("ci_lower")),
                ci_upper=_safe_float(row.get("ci_upper")),
                clonal_probability=_safe_float(row.get("clonal_probability")),
                multiplicity=_safe_float(row.get("multiplicity")),
                local_cn_total=_safe_float(row.get("local_cn_total")),
                local_cn_minor=_safe_float(row.get("local_cn_minor")),
                cluster_id=str(row.get("cluster_id") or "").strip() or None,
                tumor_purity=_safe_float(row.get("tumor_purity")),
                depth=depth,
                alt_reads=alt_reads,
                sample_id=sample_id,
            )
            genomic = combine_genomic_callers([est], policy=pol)

        expr_val = _safe_float(row.get("scrna_clonality"))
        expressional: ExpressionalClonality | None = None
        if expr_val is not None:
            power = _safe_float(row.get("scrna_power"))
            expressional = ExpressionalClonality(
                assay_platform=str(row.get("scrna_platform") or "scrna_generic").strip().lower(),
                measurement_level=str(row.get("scrna_level") or "allele_or_junction").strip().lower(),
                detection_power=power,
                dropout_adjusted_clonality=_clamp01(expr_val),
                selected_clonality=_clamp01(expr_val),
                selected_scale="per_cell",
                abundance_confounded=bool(power is not None and power < pol.min_scrna_detection_power),
            )

        raw_pres = row.get("present")
        present: bool | None = None
        if raw_pres is not None and str(raw_pres).strip() != "":
            present = str(raw_pres).strip().lower() in {"true", "1", "yes"}

        obs = SampleClonalityObservation(
            sample_id=sample_id,
            timepoint_index=tp_idx,
            timepoint_label=tp_label,
            region_label=region_label,
            genomic=genomic,
            expressional=expressional,
            depth=depth,
            present=present,
        )
        if antigen_id not in grouped:
            grouped[antigen_id] = (antigen_kind, [])
        grouped[antigen_id][1].append(obs)

    return {
        aid: combine_spatial_and_longitudinal(aid, akind, obs_list, policy=pol)
        for aid, (akind, obs_list) in grouped.items()
    }


def detect_and_load_clonality_file(
    path: str | Path,
    *,
    purity_tsv_path: str | Path | None = None,
    policy: ClonalityPolicy | None = None,
) -> dict[str, AntigenClonality]:
    """Auto-detect upstream clonality format (PyClone-VI, PURPLE VCF, FACETS MAF,
    PureCN CSV, or canonical Vaxrank JSON/TSV) and return ``antigen_id -> AntigenClonality``.
    """
    p = Path(path)
    pol = policy or ClonalityPolicy()
    name_lower = p.name.lower()
    if name_lower.endswith(".json"):
        return load_canonical_clonality_table(p, policy=pol)
    if name_lower.endswith(".vcf") or name_lower.endswith(".vcf.gz"):
        by_key = load_purple_vcf(p, purity_tsv_path=purity_tsv_path)
        return {
            k: combine_sample_clonality(k, "mutation", genomic=combine_genomic_callers(v, pol), policy=pol)
            for k, v in by_key.items()
        }

    with _open_text(p) as handle:
        header_line = ""
        for line in handle:
            if line.strip() and not line.startswith("#"):
                header_line = line.strip()
                break
    cols = {c.strip() for c in (header_line.split("\t") if "\t" in header_line else header_line.split(","))}
    if {"mutation_id", "cellular_prevalence"} <= cols:
        by_key = load_pyclone_vi(p)
        return {
            k: combine_sample_clonality(k, "mutation", genomic=combine_genomic_callers(v, pol), policy=pol)
            for k, v in by_key.items()
        }
    if "ccf_expected_copies" in cols or "ccf_Mcopies" in cols:
        by_key = load_facets_maf(p)
        return {
            k: combine_sample_clonality(k, "mutation", genomic=combine_genomic_callers(v, pol), policy=pol)
            for k, v in by_key.items()
        }
    if "CELLFRACTION" in cols or "ML.C" in cols or "POSTERIOR.SOMATIC" in cols:
        by_key = load_purecn_variants(p)
        return {
            k: combine_sample_clonality(k, "mutation", genomic=combine_genomic_callers(v, pol), policy=pol)
            for k, v in by_key.items()
        }
    return load_canonical_clonality_table(p, policy=pol)


@dataclass
class ClonalityBundle:
    """Container of AntigenClonality records indexed by antigen/variant/gene ID."""

    policy: ClonalityPolicy = field(default_factory=ClonalityPolicy)
    records: dict[str, AntigenClonality] = field(default_factory=dict)

    def lookup(
        self,
        *,
        antigen_id: str | None = None,
        variant: Any = None,
        gene_name: str | None = None,
    ) -> AntigenClonality | None:
        candidates: list[str] = []
        if antigen_id:
            candidates.append(str(antigen_id))
        if variant is not None:
            contig = getattr(variant, "contig", None)
            start = getattr(variant, "start", None)
            ref = getattr(variant, "ref", None)
            alt = getattr(variant, "alt", None)
            if contig and start is not None and ref is not None and alt is not None:
                candidates.append(_variant_key(contig, start, ref, alt))
                candidates.append(f"{contig}:{start}:{ref}:{alt}")
                candidates.append(f"chr{normalize_contig(contig)}:{start}:{ref}:{alt}")
            short = getattr(variant, "short_description", None)
            if short:
                candidates.append(str(short))
        if gene_name:
            candidates.append(str(gene_name))
            candidates.append(str(gene_name).upper())
        for key in candidates:
            if key in self.records:
                return self.records[key]
        return None

    def summary(self) -> dict[str, Any]:
        counts: dict[str, int] = {}
        trajectories: dict[str, int] = {}
        for rec in self.records.values():
            counts[rec.clonality_status] = counts.get(rec.clonality_status, 0) + 1
            if rec.longitudinal is not None:
                t = rec.longitudinal.trajectory
                trajectories[t] = trajectories.get(t, 0) + 1
        return {
            "policy_mode": self.policy.mode,
            "n_antigens": len(self.records),
            "status_counts": dict(sorted(counts.items())),
            "longitudinal_trajectories": dict(sorted(trajectories.items())),
        }

    def to_dict(self) -> dict[str, Any]:
        return {
            "policy": self.policy.to_dict(),
            "summary": self.summary(),
            "records": {k: v.to_dict() for k, v in sorted(self.records.items())},
        }


def merge_clonality_records(
    record_maps: list[dict[str, AntigenClonality]],
    policy: ClonalityPolicy | None = None,
) -> dict[str, AntigenClonality]:
    """Merge multiple ``antigen_id -> AntigenClonality`` maps, reconciling
    multiple DNA callers by precedence & discordance envelope rather than
    silently overwriting earlier inputs.
    """
    pol = policy or ClonalityPolicy()
    grouped: dict[str, list[AntigenClonality]] = {}
    for rmap in record_maps:
        for aid, rec in (rmap or {}).items():
            grouped.setdefault(aid, []).append(rec)

    merged: dict[str, AntigenClonality] = {}
    for aid, recs in grouped.items():
        if len(recs) == 1:
            merged[aid] = recs[0]
            continue
        # If any record carries longitudinal trajectories, preserve it if no other
        # record adds callers; otherwise reconcile DNA callers + pick highest-power scRNA.
        all_callers: list[CallerClonalityEstimate] = []
        mt_genomic: GenomicClonality | None = None
        best_expr: ExpressionalClonality | None = None
        akind = recs[0].antigen_kind
        for r in recs:
            if r.genomic is not None:
                if r.genomic.compartment == "mitochondrial":
                    mt_genomic = r.genomic
                elif r.genomic.caller_estimates:
                    all_callers.extend(r.genomic.caller_estimates)
                elif r.genomic.dna_ccf is not None:
                    all_callers.append(CallerClonalityEstimate(
                        source_tool=r.genomic.primary_source or "canonical",
                        ccf=r.genomic.dna_ccf,
                        ci_lower=r.genomic.ci_lower,
                        ci_upper=r.genomic.ci_upper,
                        clonal_probability=r.genomic.clonal_probability,
                        multiplicity=r.genomic.multiplicity,
                        local_cn_total=r.genomic.local_cn_total,
                        local_cn_minor=r.genomic.local_cn_minor,
                        cluster_id=r.genomic.cluster_id,
                        tumor_purity=r.genomic.tumor_purity,
                        depth=r.genomic.depth,
                    ))
            if r.expressional is not None:
                if best_expr is None or (
                    (r.expressional.detection_power or 0.0, r.expressional.n_covered_tumor_cells)
                    > (best_expr.detection_power or 0.0, best_expr.n_covered_tumor_cells)
                ):
                    best_expr = r.expressional
        genomic = (
            mt_genomic
            if mt_genomic is not None
            else (combine_genomic_callers(all_callers, policy=pol) if all_callers else None)
        )
        merged[aid] = combine_sample_clonality(
            aid, akind, genomic=genomic, expressional=best_expr, policy=pol,
        )
    return merged


def apply_clonality_to_ranked(
    ranked_variants_with_vaccine_peptides: list[tuple[Any, list[Any]]],
    bundle: ClonalityBundle,
) -> tuple[list[tuple[Any, list[Any]]], dict[str, Any]]:
    """Apply ``ClonalityBundle`` to ranked ``(source, vaccine_peptides)`` entries.

    - In ``annotate`` mode: preserves exact ordering and scores; returns full
      audit summary.
    - In ``gate`` mode: drops only candidates failing ``clonality.passes_gate(policy)``
      (proven subclonal or lost in latest timepoint; never drops unassessed
      SVs/CTAs/viruses unless ``require_measured=True``).
    - In ``weight`` mode: re-sorts sources by ``top_score * score_multiplier(policy)``.
    - In ``tier`` mode: lexicographically ranks by ``(tier_rank, -top_score)``.
    """
    pol = bundle.policy
    audit: dict[str, Any] = {
        "mode": pol.mode,
        "dropped": [],
        "annotated": {},
    }
    if pol.mode == "annotate" or not ranked_variants_with_vaccine_peptides:
        for source, _peptides in ranked_variants_with_vaccine_peptides:
            sid = str(getattr(source, "short_description", None) or source)
            rec = bundle.lookup(antigen_id=sid, variant=source)
            if rec is not None:
                audit["annotated"][sid] = rec.to_dict()
        return list(ranked_variants_with_vaccine_peptides), audit

    kept: list[tuple[float, int, int, Any, list[Any]]] = []
    for idx, (source, peptides) in enumerate(ranked_variants_with_vaccine_peptides):
        sid = str(getattr(source, "short_description", None) or source)
        gene = getattr(source, "gene_name", None)
        rec = bundle.lookup(antigen_id=sid, variant=source, gene_name=gene)
        if rec is None:
            akind = str(getattr(source, "antigen_kind", None) or "mutation")
            rec = AntigenClonality(antigen_id=sid, antigen_kind=akind)
        audit["annotated"][sid] = rec.to_dict()

        if pol.mode == "gate" and not rec.passes_gate(pol):
            audit["dropped"].append({"antigen_id": sid, "clonality_status": rec.clonality_status})
            continue

        top_score = float(getattr(peptides[0], "combined_score", 0.0) or 0.0) if peptides else 0.0
        weighted = top_score * rec.score_multiplier(pol)
        tier = rec.tier_rank() if pol.mode == "tier" else 0
        kept.append((weighted, tier, idx, source, peptides))

    if pol.mode == "weight":
        kept.sort(key=lambda item: (-item[0], item[2]))
    elif pol.mode == "tier":
        kept.sort(key=lambda item: (item[1], -item[0], item[2]))

    return [(src, peps) for _, _, _, src, peps in kept], audit
