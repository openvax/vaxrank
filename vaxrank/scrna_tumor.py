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

"""Cross-species single-cell RNA-seq tumor-cluster identification, per-cluster
gene nTPM quantification, and allele-specific MHC expression (#936).

Why single-cell matters for personalized vaccine design:
- Bulk tumor RNA-seq mixes malignant cells with CD45+ immune infiltrate, which
  expresses high levels of intact MHC Class I and APM genes and can mask
  tumor-cell MHC silencing or distort target gene expression.
- Identifying the malignant cluster(s) and an internal immune-control cluster
  (PTPRC/CD45+ T/B/myeloid cells) enables:
  1. Pure **tumor-cluster gene nTPM** (``1e6 * gene_umi / total_umi`` within
     tumor cells) for every gene, including APM/IFN-gamma pathway genes.
  2. Purity-independent **tumor vs immune-control allele-specific MHC nTPM**
     from per-barcode allele UMI matrices (e.g. ``scHLAcount``, ``T1K-sc``,
     or allele-tagged 10x 5' / Smart-seq2 BAMs), with strict UMI resolvability
     guards so 10x 3' exon-2/3 coverage gaps are labeled ``unresolvable``
     rather than falsely ``silenced``.

Works cross-species (human, canine, murine, swine) using case-insensitive
ortholog lineage marker panels and ``mhcgnomes`` allele normalization.
"""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, field
import gzip
import json
import logging
import math
from numbers import Real
from pathlib import Path
from typing import Any

from .apm_variants import _APM_ALIAS_TO_CANONICAL
from .mhc_allele_status import (
    APMGeneStatus,
    APMStatus,
    MHCAlleleStatus,
    MHCAlleleStatusBundle,
    normalize_mhc_allele_name,
)

logger = logging.getLogger(__name__)

# Cross-species conserved immune & stromal lineage marker sets (uppercase; matched
# case-insensitively so mouse Ptprc/Cd3d, human PTPRC/CD3D, dog PTPRC/CD3D match).
IMMUNE_LINEAGE_MARKERS = frozenset({
    "PTPRC",  # CD45 pan-leukocyte
    "CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "CD4",  # T cells
    "MS4A1", "CD79A", "CD79B", "CD19",              # B cells
    "CD68", "CD14", "LYZ", "C1QA", "C1QB", "TYROBP", "AIF1",  # Myeloid / macrophage
    "NKG7", "GNLY", "KLRD1", "NCAM1",               # NK / cytotoxic
})
STROMAL_LINEAGE_MARKERS = frozenset({
    "PECAM1", "VWF", "CDH5", "KDR",                 # Endothelial
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "ACTA2", "FAP", "PDGFRB",  # Fibroblast/CAF
})


def _open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def _safe_float(val: Any) -> float | None:
    if val is None or isinstance(val, bool):
        return None
    if isinstance(val, Real):
        f = float(val)
        return f if math.isfinite(f) else None
    text = str(val).strip()
    if not text:
        return None
    try:
        f = float(text)
    except ValueError:
        return None
    return f if math.isfinite(f) else None


@dataclass(frozen=True)
class SCRNATumorClusterPolicy:
    """Explicit per-trial YAML policy for scRNA tumor-cluster extraction,
    per-cluster gene nTPM, and per-cluster MHC allele ASE."""
    min_cells_per_cluster: int = 5
    min_somatic_variant_cells_for_tumor: int = 3
    min_somatic_variant_fraction_for_tumor: float = 0.05
    max_immune_marker_ntpm_for_tumor: float = 1500.0
    min_immune_marker_ntpm_for_control: float = 1000.0
    max_stromal_marker_ntpm_for_tumor: float = 3000.0
    silenced_apm_tumor_to_immune_ratio: float = 0.1
    min_scrna_resolvable_umis: int = 20
    chemistry: str = "unknown"
    tumor_markers: tuple[str, ...] = field(default_factory=tuple)
    explicit_tumor_clusters: tuple[str, ...] = field(default_factory=tuple)
    explicit_immune_clusters: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self):
        for name, lo in (
            ("min_cells_per_cluster", 1),
            ("min_somatic_variant_cells_for_tumor", 1),
            ("min_scrna_resolvable_umis", 0),
        ):
            val = getattr(self, name)
            if isinstance(val, bool) or not isinstance(val, int) or val < lo:
                raise ValueError(f"scrna_tumor.{name} must be an integer >= {lo}")
        for name in (
            "min_somatic_variant_fraction_for_tumor",
            "silenced_apm_tumor_to_immune_ratio",
        ):
            v = _safe_float(getattr(self, name))
            if v is None or not 0.0 <= v <= 1.0:
                raise ValueError(f"scrna_tumor.{name} must be in [0.0, 1.0]")
            object.__setattr__(self, name, v)
        for name in (
            "max_immune_marker_ntpm_for_tumor",
            "min_immune_marker_ntpm_for_control",
            "max_stromal_marker_ntpm_for_tumor",
        ):
            v = _safe_float(getattr(self, name))
            if v is None or v < 0.0:
                raise ValueError(f"scrna_tumor.{name} must be >= 0.0")
            object.__setattr__(self, name, v)
        object.__setattr__(self, "chemistry", str(self.chemistry or "unknown").strip())
        object.__setattr__(
            self, "tumor_markers",
            tuple(str(m).strip().upper() for m in (self.tumor_markers or ()) if str(m).strip()),
        )
        object.__setattr__(
            self, "explicit_tumor_clusters",
            tuple(str(c).strip() for c in (self.explicit_tumor_clusters or ()) if str(c).strip()),
        )
        object.__setattr__(
            self, "explicit_immune_clusters",
            tuple(str(c).strip() for c in (self.explicit_immune_clusters or ()) if str(c).strip()),
        )

    @classmethod
    def resolve(
        cls,
        value: SCRNATumorClusterPolicy | dict[str, Any] | None = None,
    ) -> SCRNATumorClusterPolicy:
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if not isinstance(value, dict):
            raise ValueError(
                f"scrna_tumor must be a dict or SCRNATumorClusterPolicy, got {type(value).__name__}"
            )
        defaults = cls().to_dict()
        unknown = set(value) - set(defaults)
        if unknown:
            raise ValueError(
                f"Unknown scrna_tumor keys: {', '.join(sorted(unknown))}. "
                f"Valid keys: {', '.join(sorted(defaults))}"
            )
        merged = {k: (defaults[k] if value.get(k) is None else value[k]) for k in defaults}
        return cls(**merged)

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["tumor_markers"] = list(self.tumor_markers)
        d["explicit_tumor_clusters"] = list(self.explicit_tumor_clusters)
        d["explicit_immune_clusters"] = list(self.explicit_immune_clusters)
        return d


@dataclass(frozen=True)
class SCRNACountMatrix:
    """Compact in-memory representation of a single-cell gene x barcode count table."""
    genes: tuple[str, ...]
    barcodes: tuple[str, ...]
    # Sparse per-barcode counts: barcode -> {gene_upper: umi_count}
    counts_by_barcode: dict[str, dict[str, float]]
    # Preserve original gene capitalization keyed by uppercase symbol
    display_gene_name: dict[str, str] = field(default_factory=dict)


def _find_first_existing(directory: Path, candidates: tuple[str, ...]) -> Path | None:
    for name in candidates:
        p = directory / name
        if p.is_file():
            return p
    return None


def _load_10x_mtx_dir(directory: Path) -> SCRNACountMatrix:
    mtx_p = _find_first_existing(directory, ("matrix.mtx.gz", "matrix.mtx"))
    bar_p = _find_first_existing(directory, ("barcodes.tsv.gz", "barcodes.tsv"))
    feat_p = _find_first_existing(
        directory, ("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv")
    )
    if not (mtx_p and bar_p and feat_p):
        raise FileNotFoundError(
            f"10x MTX directory {directory} must contain matrix.mtx[.gz], "
            "barcodes.tsv[.gz], and features.tsv[.gz] (or genes.tsv[.gz])."
        )
    with _open_text(bar_p) as fh:
        barcodes = [ln.strip().split("\t")[0] for ln in fh if ln.strip()]
    genes: list[str] = []
    display_map: dict[str, str] = {}
    with _open_text(feat_p) as fh:
        for ln in fh:
            s = ln.strip()
            if not s:
                continue
            cols = s.split("\t")
            gname = cols[1].strip() if len(cols) >= 2 and cols[1].strip() else cols[0].strip()
            genes.append(gname)
            display_map.setdefault(gname.upper(), gname)

    counts_by_barcode: dict[str, dict[str, float]] = {bc: {} for bc in barcodes}
    with _open_text(mtx_p) as fh:
        header_seen = False
        for ln in fh:
            s = ln.strip()
            if not s or s.startswith("%"):
                continue
            if not header_seen:
                header_seen = True
                continue
            parts = s.split()
            if len(parts) < 3:
                continue
            row_idx = int(parts[0]) - 1
            col_idx = int(parts[1]) - 1
            val = float(parts[2])
            if val <= 0:
                continue
            if 0 <= row_idx < len(genes) and 0 <= col_idx < len(barcodes):
                bc = barcodes[col_idx]
                gu = genes[row_idx].upper()
                counts_by_barcode[bc][gu] = counts_by_barcode[bc].get(gu, 0.0) + val

    return SCRNACountMatrix(
        genes=tuple(genes),
        barcodes=tuple(barcodes),
        counts_by_barcode=counts_by_barcode,
        display_gene_name=display_map,
    )


def load_scrna_counts(path: str | Path) -> SCRNACountMatrix:
    """Load scRNA UMI counts from a 10x MTX directory, TSV/CSV matrix, or
    optional ``.h5ad`` file."""
    p = Path(path).expanduser().resolve()
    if p.is_dir():
        return _load_10x_mtx_dir(p)
    if not p.is_file():
        raise FileNotFoundError(f"scRNA counts input not found: {p}")

    if p.suffix.lower() == ".h5ad":
        try:
            import anndata
        except ImportError as exc:
            raise ImportError(
                "Reading .h5ad requires optional dependency 'anndata' (`pip install anndata`). "
                "Alternatively, export 10x MTX or a gene-by-barcode TSV/CSV."
            ) from exc
        adata = anndata.read_h5ad(str(p))
        genes = [str(g) for g in adata.var_names]
        barcodes = [str(b) for b in adata.obs_names]
        display_map = {g.upper(): g for g in genes}
        counts_by_barcode: dict[str, dict[str, float]] = {bc: {} for bc in barcodes}
        X = adata.X
        if hasattr(X, "tocoo"):
            coo = X.tocoo()
            for r, c, v in zip(coo.row, coo.col, coo.data):
                if v > 0:
                    bc = barcodes[int(r)]
                    gu = genes[int(c)].upper()
                    counts_by_barcode[bc][gu] = counts_by_barcode[bc].get(gu, 0.0) + float(v)
        else:
            for r, bc in enumerate(barcodes):
                row_vals = X[r]
                for c, v in enumerate(row_vals):
                    fv = float(v)
                    if fv > 0:
                        gu = genes[c].upper()
                        counts_by_barcode[bc][gu] = counts_by_barcode[bc].get(gu, 0.0) + fv
        return SCRNACountMatrix(
            genes=tuple(genes),
            barcodes=tuple(barcodes),
            counts_by_barcode=counts_by_barcode,
            display_gene_name=display_map,
        )

    # Wide TSV/CSV: either genes-as-rows (first col named gene/gene_name/symbol)
    # or barcodes-as-rows (first col named barcode/cell/cell_id).
    delimiter = "\t" if p.suffix.lower() in (".tsv", ".txt") or p.name.endswith(".tsv.gz") else ","
    with _open_text(p) as fh:
        reader = csv.reader(fh, delimiter=delimiter)
        header = next(reader, None)
        if not header or len(header) < 2:
            raise ValueError(f"Empty or malformed scRNA counts table: {p}")
        first_col = header[0].strip().lower()
        rows = [row for row in reader if row and any(c.strip() for c in row)]

    if first_col in ("barcode", "cell", "cell_id", "cell_barcode"):
        genes = [h.strip() for h in header[1:]]
        display_map = {g.upper(): g for g in genes}
        barcodes: list[str] = []
        counts_by_barcode = {}
        for row in rows:
            bc = row[0].strip()
            barcodes.append(bc)
            cdict: dict[str, float] = {}
            for g, raw_v in zip(genes, row[1:]):
                fv = _safe_float(raw_v)
                if fv is not None and fv > 0:
                    gu = g.upper()
                    cdict[gu] = cdict.get(gu, 0.0) + fv
            counts_by_barcode[bc] = cdict
        return SCRNACountMatrix(
            genes=tuple(genes),
            barcodes=tuple(barcodes),
            counts_by_barcode=counts_by_barcode,
            display_gene_name=display_map,
        )

    # Default: genes as rows, barcodes as columns
    barcodes = [h.strip() for h in header[1:]]
    genes = []
    display_map = {}
    counts_by_barcode = {bc: {} for bc in barcodes}
    for row in rows:
        gname = row[0].strip()
        if not gname:
            continue
        genes.append(gname)
        gu = gname.upper()
        display_map.setdefault(gu, gname)
        for bc, raw_v in zip(barcodes, row[1:]):
            fv = _safe_float(raw_v)
            if fv is not None and fv > 0:
                counts_by_barcode[bc][gu] = counts_by_barcode[bc].get(gu, 0.0) + fv
    return SCRNACountMatrix(
        genes=tuple(genes),
        barcodes=tuple(barcodes),
        counts_by_barcode=counts_by_barcode,
        display_gene_name=display_map,
    )


def load_cluster_assignments(path: str | Path) -> dict[str, str]:
    """Load ``barcode -> cluster_id`` mapping from TSV/CSV."""
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(f"Cluster assignments file not found: {p}")
    delimiter = "\t" if p.suffix.lower() in (".tsv", ".txt") or p.name.endswith(".tsv.gz") else ","
    mapping: dict[str, str] = {}
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        for row in reader:
            row_l = {str(k).strip().lower(): str(v or "").strip() for k, v in row.items() if k is not None}
            bc = row_l.get("barcode") or row_l.get("cell") or row_l.get("cell_id") or ""
            cl = (
                row_l.get("cluster")
                or row_l.get("cluster_id")
                or row_l.get("leiden")
                or row_l.get("seurat_clusters")
                or row_l.get("cell_type")
                or ""
            )
            if bc and cl:
                mapping[bc] = cl
    return mapping


def load_somatic_barcode_calls(path: str | Path) -> dict[str, tuple[int, int]]:
    """Load per-barcode somatic variant read/UMI support (e.g. from VarTrix or
    cellsnp-lite) as ``{barcode: (alt_umis, ref_umis)}``."""
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(f"Somatic barcode calls file not found: {p}")
    delimiter = "\t" if p.suffix.lower() in (".tsv", ".txt") or p.name.endswith(".tsv.gz") else ","
    calls: dict[str, list[int]] = {}
    with _open_text(p) as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        for row in reader:
            row_l = {str(k).strip().lower(): str(v or "").strip() for k, v in row.items() if k is not None}
            bc = row_l.get("barcode") or row_l.get("cell") or row_l.get("cell_id") or ""
            if not bc:
                continue
            alt = int(_safe_float(row_l.get("alt_count") or row_l.get("alt_umis") or row_l.get("alt")) or 0)
            ref = int(_safe_float(row_l.get("ref_count") or row_l.get("ref_umis") or row_l.get("ref")) or 0)
            cur = calls.setdefault(bc, [0, 0])
            cur[0] += max(0, alt)
            cur[1] += max(0, ref)
    return {bc: (v[0], v[1]) for bc, v in calls.items()}


def load_scrna_mhc_allele_umis(path: str | Path) -> dict[str, dict[str, float]]:
    """Load per-barcode MHC allele UMI counts from a ``scHLAcount`` directory
    (``count_matrix.mtx`` + ``alleles.tsv`` + ``barcodes.tsv``) OR a TSV/CSV
    table (long ``barcode,allele,umi_count`` or wide ``allele x barcode``).

    Returns ``{barcode: {canonical_mhc_allele: umi_count}}``.
    """
    p = Path(path).expanduser().resolve()
    if p.is_dir():
        mtx_p = _find_first_existing(p, ("count_matrix.mtx.gz", "count_matrix.mtx", "matrix.mtx.gz", "matrix.mtx"))
        bar_p = _find_first_existing(p, ("barcodes.tsv.gz", "barcodes.tsv"))
        all_p = _find_first_existing(p, ("alleles.tsv.gz", "alleles.tsv", "labels.tsv.gz", "labels.tsv"))
        if not (mtx_p and bar_p and all_p):
            raise FileNotFoundError(
                f"scHLAcount directory {p} must contain count_matrix.mtx[.gz], "
                "barcodes.tsv[.gz], and alleles.tsv[.gz]."
            )
        with _open_text(bar_p) as fh:
            barcodes = [ln.strip().split("\t")[0] for ln in fh if ln.strip()]
        with _open_text(all_p) as fh:
            alleles = [normalize_mhc_allele_name(ln.strip().split("\t")[0]) for ln in fh if ln.strip()]
        out: dict[str, dict[str, float]] = {}
        with _open_text(mtx_p) as fh:
            header_seen = False
            for ln in fh:
                s = ln.strip()
                if not s or s.startswith("%"):
                    continue
                if not header_seen:
                    header_seen = True
                    continue
                parts = s.split()
                if len(parts) < 3:
                    continue
                r, c, v = int(parts[0]) - 1, int(parts[1]) - 1, float(parts[2])
                if v > 0 and 0 <= r < len(alleles) and 0 <= c < len(barcodes):
                    bc = barcodes[c]
                    a = alleles[r]
                    bdict = out.setdefault(bc, {})
                    bdict[a] = bdict.get(a, 0.0) + v
        return out

    if not p.is_file():
        raise FileNotFoundError(f"scRNA MHC allele UMI file not found: {p}")

    delimiter = "\t" if p.suffix.lower() in (".tsv", ".txt") or p.name.endswith(".tsv.gz") else ","
    with _open_text(p) as fh:
        reader = csv.reader(fh, delimiter=delimiter)
        header = next(reader, None)
        if not header or len(header) < 2:
            return {}
        hl = [h.strip().lower() for h in header]
        rows = [r for r in reader if r and any(c.strip() for c in r)]

    out = {}
    if "barcode" in hl and "allele" in hl:
        bc_i = hl.index("barcode")
        al_i = hl.index("allele")
        cnt_i = next(
            (hl.index(c) for c in ("umi_count", "umis", "count", "umi") if c in hl),
            2 if len(hl) > 2 else -1,
        )
        for r in rows:
            if len(r) <= max(bc_i, al_i, cnt_i):
                continue
            bc = r[bc_i].strip()
            al = normalize_mhc_allele_name(r[al_i].strip())
            v = _safe_float(r[cnt_i])
            if bc and al and v is not None and v > 0:
                bdict = out.setdefault(bc, {})
                bdict[al] = bdict.get(al, 0.0) + v
        return out

    # Wide format: rows are alleles, columns are barcodes
    barcodes = [h.strip() for h in header[1:]]
    for r in rows:
        al = normalize_mhc_allele_name(r[0].strip())
        if not al:
            continue
        for bc, raw_v in zip(barcodes, r[1:]):
            v = _safe_float(raw_v)
            if v is not None and v > 0:
                bdict = out.setdefault(bc, {})
                bdict[al] = bdict.get(al, 0.0) + v
    return out


@dataclass(frozen=True)
class ClusterSummary:
    cluster_id: str
    classification: str  # tumor | immune_control | stromal | other
    n_cells: int
    total_umis: float
    immune_marker_ntpm: float
    stromal_marker_ntpm: float
    tumor_marker_ntpm: float
    somatic_variant_cells: int
    somatic_variant_fraction: float
    reason: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class SCRNATumorAnalysisResult:
    """Output of single-cell tumor-cluster extraction, per-cluster gene nTPM,
    and per-cluster MHC allele ASE."""
    clusters: tuple[ClusterSummary, ...]
    tumor_cluster_ids: tuple[str, ...]
    immune_control_cluster_ids: tuple[str, ...]
    tumor_gene_ntpm: dict[str, float]
    immune_control_gene_ntpm: dict[str, float]
    per_cluster_gene_ntpm: dict[str, dict[str, float]]
    mhc_allele_status_bundle: MHCAlleleStatusBundle
    policy: SCRNATumorClusterPolicy

    def to_mhc_allele_status_bundle(self) -> MHCAlleleStatusBundle:
        return self.mhc_allele_status_bundle

    def write_tumor_cluster_ntpm_tsv(self, path: str | Path) -> str:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        all_genes = sorted(set(self.tumor_gene_ntpm) | set(self.immune_control_gene_ntpm))
        with open(p, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow([
                "gene",
                "tumor_cluster_ntpm",
                "immune_control_ntpm",
                "tumor_to_immune_ratio",
            ])
            for g in all_genes:
                t_val = self.tumor_gene_ntpm.get(g, 0.0)
                i_val = self.immune_control_gene_ntpm.get(g, 0.0)
                ratio = f"{t_val / i_val:.6f}" if i_val > 0 else ""
                writer.writerow([g, f"{t_val:.4f}", f"{i_val:.4f}", ratio])
        return str(p)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": "vaxrank.scrna_tumor_analysis.v1",
            "policy": self.policy.to_dict(),
            "tumor_cluster_ids": list(self.tumor_cluster_ids),
            "immune_control_cluster_ids": list(self.immune_control_cluster_ids),
            "clusters": [c.to_dict() for c in self.clusters],
            "gene_count": len(self.tumor_gene_ntpm),
            "mhc_allele_status": self.mhc_allele_status_bundle.to_dict(),
        }

    def write_summary_json(self, path: str | Path) -> str:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        with open(p, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh, indent=2, allow_nan=False)
        return str(p)


def analyze_scrna_tumor_clusters(
    counts: SCRNACountMatrix,
    clusters: dict[str, str],
    *,
    somatic_barcodes: dict[str, tuple[int, int]] | None = None,
    mhc_allele_umis: dict[str, dict[str, float]] | None = None,
    policy: SCRNATumorClusterPolicy | dict[str, Any] | None = None,
) -> SCRNATumorAnalysisResult:
    """Identify tumor vs immune-control clusters from scRNA counts, compute
    per-cluster and aggregated tumor-cluster gene nTPM, derive APM gene RNA
    expression status, and aggregate per-cluster MHC allele ASE UMIs."""
    pol = SCRNATumorClusterPolicy.resolve(policy)
    somatic_barcodes = somatic_barcodes or {}
    mhc_allele_umis = mhc_allele_umis or {}

    # Group barcodes by cluster
    barcodes_by_cluster: dict[str, list[str]] = {}
    for bc in counts.barcodes:
        cl = clusters.get(bc)
        if cl:
            barcodes_by_cluster.setdefault(cl, []).append(bc)

    per_cluster_gene_sums: dict[str, dict[str, float]] = {}
    per_cluster_total_umis: dict[str, float] = {}
    per_cluster_gene_ntpm: dict[str, dict[str, float]] = {}
    summaries: list[ClusterSummary] = []

    for cl_id, bcs in sorted(barcodes_by_cluster.items()):
        gene_sums: dict[str, float] = {}
        total_umis = 0.0
        som_cells = 0
        for bc in bcs:
            cdict = counts.counts_by_barcode.get(bc, {})
            for gu, v in cdict.items():
                gene_sums[gu] = gene_sums.get(gu, 0.0) + v
                total_umis += v
            alt_c, _ = somatic_barcodes.get(bc, (0, 0))
            if alt_c > 0:
                som_cells += 1

        per_cluster_gene_sums[cl_id] = gene_sums
        per_cluster_total_umis[cl_id] = total_umis
        ntpm_map: dict[str, float] = {}
        if total_umis > 0:
            scale = 1e6 / total_umis
            for gu, v in gene_sums.items():
                ntpm_map[gu] = round(v * scale, 4)
        per_cluster_gene_ntpm[cl_id] = ntpm_map

        immune_ntpm = sum(ntpm_map.get(m, 0.0) for m in IMMUNE_LINEAGE_MARKERS)
        stromal_ntpm = sum(ntpm_map.get(m, 0.0) for m in STROMAL_LINEAGE_MARKERS)
        tumor_m_ntpm = sum(ntpm_map.get(m, 0.0) for m in pol.tumor_markers)
        som_frac = som_cells / len(bcs) if bcs else 0.0

        # Deterministic, transparent classification hierarchy
        if cl_id in pol.explicit_tumor_clusters:
            cls_label, reason = "tumor", "explicit_tumor_cluster"
        elif cl_id in pol.explicit_immune_clusters:
            cls_label, reason = "immune_control", "explicit_immune_cluster"
        elif len(bcs) < pol.min_cells_per_cluster:
            cls_label, reason = "other", f"too_few_cells ({len(bcs)} < {pol.min_cells_per_cluster})"
        elif (
            somatic_barcodes
            and som_cells >= pol.min_somatic_variant_cells_for_tumor
            and som_frac >= pol.min_somatic_variant_fraction_for_tumor
            and immune_ntpm <= pol.max_immune_marker_ntpm_for_tumor
        ):
            cls_label, reason = "tumor", f"somatic_variant_enriched ({som_cells}/{len(bcs)} cells)"
        elif immune_ntpm >= pol.min_immune_marker_ntpm_for_control and immune_ntpm >= stromal_ntpm:
            cls_label, reason = "immune_control", f"immune_markers_high (nTPM={immune_ntpm:.1f})"
        elif stromal_ntpm > pol.max_stromal_marker_ntpm_for_tumor and stromal_ntpm > tumor_m_ntpm:
            cls_label, reason = "stromal", f"stromal_markers_high (nTPM={stromal_ntpm:.1f})"
        elif (
            immune_ntpm <= pol.max_immune_marker_ntpm_for_tumor
            and stromal_ntpm <= pol.max_stromal_marker_ntpm_for_tumor
            and (not pol.tumor_markers or tumor_m_ntpm > 0)
        ):
            cls_label, reason = "tumor", "non_immune_non_stromal_malignant_heuristic"
        else:
            cls_label, reason = "other", "unclassified_lineage"

        summaries.append(
            ClusterSummary(
                cluster_id=cl_id,
                classification=cls_label,
                n_cells=len(bcs),
                total_umis=round(total_umis, 2),
                immune_marker_ntpm=round(immune_ntpm, 2),
                stromal_marker_ntpm=round(stromal_ntpm, 2),
                tumor_marker_ntpm=round(tumor_m_ntpm, 2),
                somatic_variant_cells=som_cells,
                somatic_variant_fraction=round(som_frac, 4),
                reason=reason,
            )
        )

    tumor_ids = tuple(s.cluster_id for s in summaries if s.classification == "tumor")
    immune_ids = tuple(s.cluster_id for s in summaries if s.classification == "immune_control")

    def _aggregate_ntpm(cl_ids: tuple[str, ...]) -> tuple[dict[str, float], float]:
        tot = sum(per_cluster_total_umis.get(cid, 0.0) for cid in cl_ids)
        gene_tot: dict[str, float] = {}
        for cid in cl_ids:
            for gu, v in per_cluster_gene_sums.get(cid, {}).items():
                gene_tot[gu] = gene_tot.get(gu, 0.0) + v
        if tot <= 0:
            return {}, 0.0
        scale = 1e6 / tot
        return {
            counts.display_gene_name.get(gu, gu): round(v * scale, 4)
            for gu, v in gene_tot.items()
        }, tot

    tumor_gene_ntpm, tumor_total_umis = _aggregate_ntpm(tumor_ids)
    immune_gene_ntpm, immune_total_umis = _aggregate_ntpm(immune_ids)

    # Uppercase lookup helpers for APM gene RNA status
    tumor_ntpm_u = {g.upper(): v for g, v in tumor_gene_ntpm.items()}
    immune_ntpm_u = {g.upper(): v for g, v in immune_gene_ntpm.items()}

    # Build APMGeneStatus (strictly RNA expression axis: dna_status="unknown")
    apm_gene_records: list[APMGeneStatus] = []
    seen_apm: set[str] = set()
    for gu in sorted(set(tumor_ntpm_u) | set(immune_ntpm_u)):
        if gu not in _APM_ALIAS_TO_CANONICAL:
            continue
        canon_gene, _ = _APM_ALIAS_TO_CANONICAL[gu]
        if canon_gene in seen_apm:
            continue
        seen_apm.add(canon_gene)
        t_v = tumor_ntpm_u.get(gu, 0.0)
        i_v = immune_ntpm_u.get(gu, 0.0)
        if tumor_total_umis > 0 and immune_total_umis > 0 and i_v > 0:
            rna_status = (
                "silenced"
                if (t_v / i_v) <= pol.silenced_apm_tumor_to_immune_ratio
                else "expressed"
            )
        elif tumor_total_umis > 0:
            rna_status = "expressed" if t_v > 0 else "unknown"
        else:
            rna_status = "unknown"
        apm_gene_records.append(
            APMGeneStatus(
                gene=canon_gene,
                dna_status="unknown",
                rna_status=rna_status,
                rna_ntpm=t_v if tumor_total_umis > 0 else None,
                immune_control_ntpm=i_v if immune_total_umis > 0 else None,
                assay_source="scrna_tumor_cluster_ntpm",
            )
        )

    # Aggregate per-barcode MHC allele UMIs across tumor vs immune-control clusters
    tumor_bc_set = {bc for cid in tumor_ids for bc in barcodes_by_cluster.get(cid, ())}
    immune_bc_set = {bc for cid in immune_ids for bc in barcodes_by_cluster.get(cid, ())}

    all_alleles: set[str] = set()
    tumor_allele_umis: dict[str, float] = {}
    immune_allele_umis: dict[str, float] = {}
    for bc, adict in mhc_allele_umis.items():
        for raw_a, u in adict.items():
            a = normalize_mhc_allele_name(raw_a)
            all_alleles.add(a)
            if bc in tumor_bc_set:
                tumor_allele_umis[a] = tumor_allele_umis.get(a, 0.0) + u
            if bc in immune_bc_set:
                immune_allele_umis[a] = immune_allele_umis.get(a, 0.0) + u

    mhc_records: list[MHCAlleleStatus] = []
    for allele in sorted(all_alleles):
        t_umi = tumor_allele_umis.get(allele, 0.0)
        i_umi = immune_allele_umis.get(allele, 0.0)
        resolvable = int(round(t_umi + i_umi))
        t_ntpm = round((t_umi / tumor_total_umis) * 1e6, 4) if tumor_total_umis > 0 else None
        i_ntpm = round((i_umi / immune_total_umis) * 1e6, 4) if immune_total_umis > 0 else None
        mhc_records.append(
            MHCAlleleStatus(
                allele=allele,
                status="unknown",
                scrna_tumor_cluster_ntpm=t_ntpm,
                scrna_immune_control_ntpm=i_ntpm,
                scrna_resolvable_umis=resolvable,
                scrna_chemistry=pol.chemistry,
                assay_source="scrna_tumor_cluster_mhc_ase",
                notes=(
                    f"tumor_clusters={','.join(tumor_ids) or 'none'}; "
                    f"immune_clusters={','.join(immune_ids) or 'none'}; "
                    f"tumor_umis={int(round(t_umi))}; immune_umis={int(round(i_umi))}"
                ),
            )
        )

    apm_status = APMStatus(genes=tuple(apm_gene_records)) if apm_gene_records else None
    bundle = MHCAlleleStatusBundle(
        alleles=tuple(mhc_records),
        apm=apm_status,
        source_path="scrna:tumor_cluster_analysis",
    )
    return SCRNATumorAnalysisResult(
        clusters=tuple(summaries),
        tumor_cluster_ids=tumor_ids,
        immune_control_cluster_ids=immune_ids,
        tumor_gene_ntpm=tumor_gene_ntpm,
        immune_control_gene_ntpm=immune_gene_ntpm,
        per_cluster_gene_ntpm=per_cluster_gene_ntpm,
        mhc_allele_status_bundle=bundle,
        policy=pol,
    )
