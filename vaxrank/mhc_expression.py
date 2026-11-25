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

"""Cross-species bulk RNA-seq allele-specific MHC expression and upstream
MHC LOH/mutation tool ingestion (#936).

Supports BOTH:
1. **Upstream companion tool ingestors** (normalized cross-species via
   ``mhcgnomes`` across human HLA, canine DLA, murine H-2, swine SLA, etc.,
   strictly keeping DNA LOH/mutations separate from RNA expression):
   - Bulk RNA ASE: ``T1K``, ``arcasHLA``, ``HLApers``
   - DNA copy-number LOH & somatic MHC mutations: ``LOHHLA``, ``DASH``,
     ``LILAC``, ``POLYSOLVER``
2. **Direct computational bulk RNA-seq allele-specific MHC estimators** (no
   external HLA-only binaries required; works for human, dog, mouse, swine):
   - Personalized patient MHC allele FASTA discriminating k-mer read quantifier
     over bulk RNA-seq FASTQ/FASTA reads (with tumor-purity caveat tracking)
   - Optional ``isovar``-assisted local RNA assembly around polymorphic MHC
     exon variants distinguishing heterozygous alleles at an MHC locus
"""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
import gzip
import json
import logging
import math
from numbers import Real
from pathlib import Path
import re
from typing import Any, Iterable

import mhcgnomes

from .mhc_allele_status import (
    MHCAlleleStatus,
    MHCAlleleStatusBundle,
    normalize_mhc_allele_name,
)

logger = logging.getLogger(__name__)

_DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_DNA_COMPLEMENT)[::-1]


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
    if not text or text.lower() in ("na", "nan", "none", "null", ".", "-"):
        return None
    try:
        f = float(text)
    except ValueError:
        return None
    return f if math.isfinite(f) else None


def _allele_locus_key(allele: str) -> str:
    """Group canonical MHC alleles by ``(species_prefix, gene_name)`` via
    ``mhcgnomes`` so within-locus ASE ratios compare homologous alleles
    (e.g. HLA-A vs HLA-A, DLA-88 vs DLA-88, H2-K vs H2-K, SLA-1 vs SLA-1)."""
    canon = normalize_mhc_allele_name(allele)
    try:
        parsed = mhcgnomes.parse(canon, raise_on_error=False)
    except Exception:
        parsed = None
    if parsed is not None:
        sp = getattr(parsed, "species_prefix", None) or ""
        gene = getattr(parsed, "gene_name", None) or ""
        if sp and gene:
            return f"{sp}-{gene}"
    # Fallback prefix before '*'
    return canon.split("*", 1)[0]


@dataclass(frozen=True)
class MHCExpressionPolicy:
    """Explicit per-trial YAML policy for bulk MHC ASE & companion tool ingestion."""
    kmer_size: int = 21
    min_discriminating_kmers: int = 5
    min_locus_supporting_reads: int = 10
    default_tumor_purity: float | None = None
    lohhla_pvalue_threshold: float = 0.01
    dash_deletion_threshold: float = 0.5
    lilac_min_copy_number: float = 0.5

    def __post_init__(self):
        for name, lo in (
            ("kmer_size", 9),
            ("min_discriminating_kmers", 1),
            ("min_locus_supporting_reads", 1),
        ):
            val = getattr(self, name)
            if isinstance(val, bool) or not isinstance(val, int) or val < lo:
                raise ValueError(f"mhc_expression.{name} must be an integer >= {lo}")
        if self.default_tumor_purity is not None:
            p = _safe_float(self.default_tumor_purity)
            if p is None or not 0.0 <= p <= 1.0:
                raise ValueError("mhc_expression.default_tumor_purity must be in [0.0, 1.0]")
            object.__setattr__(self, "default_tumor_purity", p)
        for name in ("lohhla_pvalue_threshold", "dash_deletion_threshold"):
            v = _safe_float(getattr(self, name))
            if v is None or not 0.0 <= v <= 1.0:
                raise ValueError(f"mhc_expression.{name} must be in [0.0, 1.0]")
            object.__setattr__(self, name, v)
        cn = _safe_float(self.lilac_min_copy_number)
        if cn is None or cn < 0.0:
            raise ValueError("mhc_expression.lilac_min_copy_number must be >= 0.0")
        object.__setattr__(self, "lilac_min_copy_number", cn)

    @classmethod
    def resolve(cls, value: MHCExpressionPolicy | dict[str, Any] | None = None) -> MHCExpressionPolicy:
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if not isinstance(value, dict):
            raise ValueError(f"mhc_expression must be a dict or MHCExpressionPolicy, got {type(value).__name__}")
        defaults = cls().to_dict()
        unknown = set(value) - set(defaults)
        if unknown:
            raise ValueError(
                f"Unknown mhc_expression keys: {', '.join(sorted(unknown))}. "
                f"Valid keys: {', '.join(sorted(defaults))}"
            )
        merged = {k: (defaults[k] if value.get(k) is None else value[k]) for k in defaults}
        return cls(**merged)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _bundle_from_allele_abundances(
    abundances: dict[str, float],
    *,
    assay_source: str,
    source_path: str,
    tumor_purity: float | None = None,
    notes_by_allele: dict[str, str] | None = None,
) -> MHCAlleleStatusBundle:
    """Convert canonical ``{allele: tpm_or_abundance}`` into an
    ``MHCAlleleStatusBundle`` with within-locus ``bulk_rna_allele_ratio``."""
    notes_by_allele = notes_by_allele or {}
    by_locus: dict[str, list[tuple[str, float]]] = {}
    for raw_allele, abund in abundances.items():
        canon = normalize_mhc_allele_name(raw_allele)
        locus = _allele_locus_key(canon)
        by_locus.setdefault(locus, []).append((canon, max(0.0, float(abund))))

    records: list[MHCAlleleStatus] = []
    for locus, entries in sorted(by_locus.items()):
        locus_total = sum(v for _, v in entries)
        for canon, abund in entries:
            ratio = (abund / locus_total) if locus_total > 0 else None
            records.append(
                MHCAlleleStatus(
                    allele=canon,
                    status="unknown",
                    bulk_rna_tpm=round(abund, 4),
                    bulk_rna_allele_ratio=(round(ratio, 6) if ratio is not None else None),
                    tumor_purity=tumor_purity,
                    assay_source=assay_source,
                    notes=notes_by_allele.get(canon, f"locus={locus}; bulk_rna_ase"),
                )
            )
    return MHCAlleleStatusBundle(alleles=tuple(records), source_path=source_path)


# ---------------------------------------------------------------------------
# 1. Upstream companion tool ingestors (cross-species via mhcgnomes)
# ---------------------------------------------------------------------------


def load_t1k_mhc_expression(
    path: str | Path,
    *,
    tumor_purity: float | None = None,
    policy: MHCExpressionPolicy | dict[str, Any] | None = None,
) -> MHCAlleleStatusBundle:
    """Parse T1K ``*_genotype.tsv`` (8-column headerless or headered) or
    ``*_expression.tsv`` into an ``MHCAlleleStatusBundle`` (RNA axis only).

    Supports human HLA, canine DLA, murine H-2, swine SLA, etc.
    """
    pol = MHCExpressionPolicy.resolve(policy)
    purity = tumor_purity if tumor_purity is not None else pol.default_tumor_purity
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(f"T1K file not found: {p}")

    abundances: dict[str, float] = {}
    with _open_text(p) as handle:
        lines = [ln.strip() for ln in handle if ln.strip() and not ln.startswith("#")]

    if not lines:
        return MHCAlleleStatusBundle(source_path=str(p))

    first_cols = lines[0].split("\t")
    header_lower = [c.strip().lower() for c in first_cols]
    has_header = any(tok in header_lower for tok in ("allele", "allele1", "gene", "abundance", "tpm"))

    if has_header:
        reader = csv.DictReader(lines, delimiter="\t")
        for row in reader:
            row_l = {str(k).strip().lower(): str(v or "").strip() for k, v in row.items() if k is not None}
            if "allele1" in row_l:
                for a_col, ab_col in (("allele1", "abund1"), ("allele2", "abund2")):
                    a_name = row_l.get(a_col, "")
                    if not a_name or a_name == ".":
                        continue
                    ab = _safe_float(row_l.get(ab_col))
                    if ab is not None:
                        canon = normalize_mhc_allele_name(a_name)
                        abundances[canon] = abundances.get(canon, 0.0) + ab
            else:
                a_name = row_l.get("allele") or row_l.get("target_id") or ""
                ab = _safe_float(row_l.get("tpm") or row_l.get("abundance") or row_l.get("count"))
                if a_name and a_name != "." and ab is not None:
                    canon = normalize_mhc_allele_name(a_name)
                    abundances[canon] = abundances.get(canon, 0.0) + ab
    else:
        # Standard headerless T1K *_genotype.tsv format:
        # gene  num_alleles  allele1  abund1  qual1  allele2  abund2  qual2
        for line in lines:
            cols = line.split("\t")
            if len(cols) >= 5:
                for a_idx, ab_idx in ((2, 3), (5, 6)):
                    if len(cols) > ab_idx:
                        a_name = cols[a_idx].strip()
                        ab = _safe_float(cols[ab_idx])
                        if a_name and a_name != "." and ab is not None:
                            canon = normalize_mhc_allele_name(a_name)
                            abundances[canon] = abundances.get(canon, 0.0) + ab
            elif len(cols) >= 2:
                a_name = cols[0].strip()
                ab = _safe_float(cols[1])
                if a_name and a_name != "." and ab is not None:
                    canon = normalize_mhc_allele_name(a_name)
                    abundances[canon] = abundances.get(canon, 0.0) + ab

    return _bundle_from_allele_abundances(
        abundances,
        assay_source="T1K",
        source_path=str(p),
        tumor_purity=purity,
    )


def load_arcashla_mhc_expression(
    path: str | Path,
    *,
    tumor_purity: float | None = None,
    policy: MHCExpressionPolicy | dict[str, Any] | None = None,
) -> MHCAlleleStatusBundle:
    """Parse arcasHLA ``*.quant.tsv`` (or CSV/JSON allele quantification) into
    an ``MHCAlleleStatusBundle`` (RNA axis only)."""
    pol = MHCExpressionPolicy.resolve(policy)
    purity = tumor_purity if tumor_purity is not None else pol.default_tumor_purity
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(f"arcasHLA file not found: {p}")

    abundances: dict[str, float] = {}
    if p.suffix.lower() == ".json":
        with _open_text(p) as handle:
            data = json.load(handle)
        if isinstance(data, dict):
            for k, v in data.items():
                ab = _safe_float(v)
                if ab is not None:
                    abundances[normalize_mhc_allele_name(k)] = ab
    else:
        delimiter = "\t" if p.suffix.lower() in (".tsv", ".txt") else ","
        with _open_text(p) as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            for row in reader:
                row_l = {str(k).strip().lower(): str(v or "").strip() for k, v in row.items() if k is not None}
                a_name = row_l.get("allele") or row_l.get("transcript") or row_l.get("target_id") or ""
                ab = _safe_float(
                    row_l.get("tpm")
                    or row_l.get("abundance")
                    or row_l.get("est_counts")
                    or row_l.get("count")
                )
                if a_name and ab is not None:
                    canon = normalize_mhc_allele_name(a_name)
                    abundances[canon] = abundances.get(canon, 0.0) + ab

    return _bundle_from_allele_abundances(
        abundances,
        assay_source="arcasHLA",
        source_path=str(p),
        tumor_purity=purity,
    )


def load_hlapers_mhc_expression(
    path: str | Path,
    *,
    tumor_purity: float | None = None,
    policy: MHCExpressionPolicy | dict[str, Any] | None = None,
) -> MHCAlleleStatusBundle:
    """Parse HLApers allele quantification TSV/CSV into an
    ``MHCAlleleStatusBundle`` (RNA axis only)."""
    pol = MHCExpressionPolicy.resolve(policy)
    purity = tumor_purity if tumor_purity is not None else pol.default_tumor_purity
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(f"HLApers file not found: {p}")

    delimiter = "\t" if p.suffix.lower() in (".tsv", ".txt") else ","
    abundances: dict[str, float] = {}
    with _open_text(p) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            row_l = {str(k).strip().lower(): str(v or "").strip() for k, v in row.items() if k is not None}
            a_name = (
                row_l.get("allele")
                or row_l.get("target_id")
                or row_l.get("name")
                or ""
            )
            ab = _safe_float(
                row_l.get("tpm")
                or row_l.get("est_counts")
                or row_l.get("numreads")
                or row_l.get("count")
            )
            if a_name and ab is not None:
                canon = normalize_mhc_allele_name(a_name)
                abundances[canon] = abundances.get(canon, 0.0) + ab

    return _bundle_from_allele_abundances(
        abundances,
        assay_source="HLApers",
        source_path=str(p),
        tumor_purity=purity,
    )


def _normalize_lohhla_allele_token(token: str) -> str:
    """Convert LOHHLA ``hla_a_02_01_01_01`` or ``HLA-A*02:01`` into canonical
    ``mhcgnomes`` form."""
    text = str(token or "").strip()
    if not text:
        return ""
    m = re.match(r"^(hla|dla|sla|h2)_([a-z0-9]+)_(\d+)_(\d+)", text, re.I)
    if m:
        sp, gene, f1, f2 = m.group(1).upper(), m.group(2).upper(), m.group(3), m.group(4)
        text = f"{sp}-{gene}*{f1}:{f2}"
    return normalize_mhc_allele_name(text)


def load_lohhla_dna_loh(
    path: str | Path,
    *,
    policy: MHCExpressionPolicy | dict[str, Any] | None = None,
) -> MHCAlleleStatusBundle:
    """Parse LOHHLA ``*.HLAlossPrediction*`` TSV into an ``MHCAlleleStatusBundle``
    (strictly DNA copy-number LOH axis)."""
    pol = MHCExpressionPolicy.resolve(policy)
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(f"LOHHLA file not found: {p}")

    delimiter = "\t" if p.suffix.lower() in (".tsv", ".txt", ".xls") else ","
    records_by_allele: dict[str, MHCAlleleStatus] = {}
    with _open_text(p) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            row_l = {str(k).strip().lower(): str(v or "").strip() for k, v in row.items() if k is not None}
            a1 = row_l.get("hla_a_type1") or row_l.get("allele1") or row_l.get("hla_type1") or ""
            a2 = row_l.get("hla_a_type2") or row_l.get("allele2") or row_l.get("hla_type2") or ""
            cn1 = _safe_float(
                row_l.get("hla_type1copynum_withbafbin")
                or row_l.get("hla_type1copynum_withbaf")
                or row_l.get("cn1")
            )
            cn2 = _safe_float(
                row_l.get("hla_type2copynum_withbafbin")
                or row_l.get("hla_type2copynum_withbaf")
                or row_l.get("cn2")
            )
            pval = _safe_float(
                row_l.get("pval_unique")
                or row_l.get("pval_total")
                or row_l.get("pvalue")
            )
            for raw_a, cn in ((a1, cn1), (a2, cn2)):
                if not raw_a:
                    continue
                canon = _normalize_lohhla_allele_token(raw_a)
                clamped_cn = max(0.0, cn) if cn is not None else None
                is_loh = None
                if clamped_cn is not None:
                    sig = True if pval is None else (pval <= pol.lohhla_pvalue_threshold)
                    is_loh = bool(clamped_cn < 0.5 and sig)
                records_by_allele[canon] = MHCAlleleStatus(
                    allele=canon,
                    status="unknown",
                    dna_copy_number=clamped_cn,
                    dna_loh=is_loh,
                    dna_loh_pvalue=pval,
                    assay_source="LOHHLA",
                )
    return MHCAlleleStatusBundle(alleles=tuple(records_by_allele.values()), source_path=str(p))


def load_dash_dna_loh(
    path: str | Path,
    *,
    policy: MHCExpressionPolicy | dict[str, Any] | None = None,
) -> MHCAlleleStatusBundle:
    """Parse DASH allele-specific MHC deletion output TSV/CSV into an
    ``MHCAlleleStatusBundle`` (DNA axis)."""
    pol = MHCExpressionPolicy.resolve(policy)
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(f"DASH file not found: {p}")

    delimiter = "\t" if p.suffix.lower() in (".tsv", ".txt") else ","
    records: list[MHCAlleleStatus] = []
    with _open_text(p) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            row_l = {str(k).strip().lower(): str(v or "").strip() for k, v in row.items() if k is not None}
            raw_a = row_l.get("allele") or row_l.get("hla_allele") or ""
            if not raw_a:
                continue
            canon = _normalize_lohhla_allele_token(raw_a)
            pred_del = row_l.get("predicted_deletion") or row_l.get("dash_deletion") or ""
            del_prob = _safe_float(row_l.get("deletion_probability") or row_l.get("dash_probability"))
            cn = _safe_float(row_l.get("copy_number") or row_l.get("tumor_copy_number"))
            dna_loh = None
            if pred_del.lower() in ("true", "1", "yes", "deleted", "loh"):
                dna_loh = True
            elif pred_del.lower() in ("false", "0", "no", "intact", "retained"):
                dna_loh = False
            elif del_prob is not None:
                dna_loh = del_prob >= pol.dash_deletion_threshold
            records.append(
                MHCAlleleStatus(
                    allele=canon,
                    status="unknown",
                    dna_copy_number=(max(0.0, cn) if cn is not None else None),
                    dna_loh=dna_loh,
                    assay_source="DASH",
                    notes=(f"deletion_prob={del_prob:.3f}" if del_prob is not None else ""),
                )
            )
    return MHCAlleleStatusBundle(alleles=tuple(records), source_path=str(p))


def load_lilac_mhc_status(
    path: str | Path,
    *,
    policy: MHCExpressionPolicy | dict[str, Any] | None = None,
) -> MHCAlleleStatusBundle:
    """Parse Hartwig LILAC ``*.lilac.tsv`` into an ``MHCAlleleStatusBundle``
    capturing BOTH DNA copy-number LOH and somatic MHC LoF/missense mutations."""
    pol = MHCExpressionPolicy.resolve(policy)
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(f"LILAC file not found: {p}")

    delimiter = "\t" if p.suffix.lower() in (".tsv", ".txt") else ","
    records: list[MHCAlleleStatus] = []
    with _open_text(p) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            row_l = {str(k).strip().lower(): str(v or "").strip() for k, v in row.items() if k is not None}
            raw_a = row_l.get("allele") or ""
            if not raw_a:
                continue
            canon = normalize_mhc_allele_name(raw_a)
            cn = _safe_float(row_l.get("tumorcopynumber") or row_l.get("tumor_copy_number"))
            dna_loh = (cn <= pol.lilac_min_copy_number) if cn is not None else None
            somatic_vars: list[dict[str, Any]] = []
            lof_n = _safe_float(row_l.get("somaticnonsenseorframeshift")) or 0.0
            splice_n = _safe_float(row_l.get("somaticsplice")) or 0.0
            mis_n = _safe_float(row_l.get("somaticmissense")) or 0.0
            indel_n = _safe_float(row_l.get("somaticinframeindel")) or 0.0
            if lof_n > 0:
                somatic_vars.append({"consequence": "truncating_lof", "is_lof": True, "source": "LILAC", "count": lof_n})
            if splice_n > 0:
                somatic_vars.append({"consequence": "splice_donor", "is_lof": True, "source": "LILAC", "count": splice_n})
            if mis_n > 0:
                somatic_vars.append({"consequence": "missense", "is_lof": False, "source": "LILAC", "count": mis_n})
            if indel_n > 0:
                somatic_vars.append({"consequence": "inframe_indel", "is_lof": False, "source": "LILAC", "count": indel_n})
            records.append(
                MHCAlleleStatus(
                    allele=canon,
                    status="unknown",
                    dna_copy_number=(max(0.0, cn) if cn is not None else None),
                    dna_loh=dna_loh,
                    dna_somatic_variants=tuple(somatic_vars),
                    assay_source="LILAC",
                )
            )
    return MHCAlleleStatusBundle(alleles=tuple(records), source_path=str(p))


def load_polysolver_mutations(path: str | Path) -> MHCAlleleStatusBundle:
    """Parse POLYSOLVER somatic MHC mutation output TSV/CSV into an
    ``MHCAlleleStatusBundle`` (DNA somatic variant axis)."""
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(f"POLYSOLVER file not found: {p}")

    delimiter = "\t" if p.suffix.lower() in (".tsv", ".txt") else ","
    vars_by_allele: dict[str, list[dict[str, Any]]] = {}
    with _open_text(p) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            row_l = {str(k).strip().lower(): str(v or "").strip() for k, v in row.items() if k is not None}
            raw_a = row_l.get("allele") or row_l.get("hla_allele") or row_l.get("chr") or ""
            if not raw_a:
                continue
            canon = _normalize_lohhla_allele_token(raw_a)
            consequence = (
                row_l.get("consequence")
                or row_l.get("effect")
                or row_l.get("variant_classification")
                or row_l.get("type")
                or "missense"
            )
            hgvs_p = row_l.get("hgvs_p") or row_l.get("protein_change") or row_l.get("mutation") or ""
            vars_by_allele.setdefault(canon, []).append({
                "consequence": consequence,
                "hgvs_p": hgvs_p,
                "source": "POLYSOLVER",
            })
    records = tuple(
        MHCAlleleStatus(
            allele=canon,
            status="unknown",
            dna_somatic_variants=tuple(vlist),
            assay_source="POLYSOLVER",
        )
        for canon, vlist in sorted(vars_by_allele.items())
    )
    return MHCAlleleStatusBundle(alleles=records, source_path=str(p))


def load_mhc_status_from_tool(
    path: str | Path,
    tool: str,
    *,
    tumor_purity: float | None = None,
    policy: MHCExpressionPolicy | dict[str, Any] | None = None,
) -> MHCAlleleStatusBundle:
    """Unified dispatcher for companion DNA LOH / somatic MHC & bulk RNA ASE tools."""
    key = str(tool or "").strip().lower().replace("-", "_")
    if key == "t1k":
        return load_t1k_mhc_expression(path, tumor_purity=tumor_purity, policy=policy)
    if key in ("arcashla", "arcas_hla"):
        return load_arcashla_mhc_expression(path, tumor_purity=tumor_purity, policy=policy)
    if key == "hlapers":
        return load_hlapers_mhc_expression(path, tumor_purity=tumor_purity, policy=policy)
    if key == "lohhla":
        return load_lohhla_dna_loh(path, policy=policy)
    if key == "dash":
        return load_dash_dna_loh(path, policy=policy)
    if key == "lilac":
        return load_lilac_mhc_status(path, policy=policy)
    if key == "polysolver":
        return load_polysolver_mutations(path)
    raise ValueError(
        f"Unknown MHC status companion tool {tool!r}; expected one of: "
        "t1k, arcashla, hlapers, lohhla, dash, lilac, polysolver."
    )


# ---------------------------------------------------------------------------
# 2. Direct computational cross-species bulk RNA-seq MHC ASE estimator
#    (discriminating k-mer quantification on personalized patient MHC FASTA)
# ---------------------------------------------------------------------------


def _read_fasta_sequences(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    current = None
    with _open_text(path) as handle:
        for line in handle:
            s = line.strip()
            if not s:
                continue
            if s.startswith(">"):
                header = s[1:].strip().split()[0]
                current = normalize_mhc_allele_name(header)
                seqs[current] = []
            elif current is not None:
                seqs[current].append(s.upper().replace("U", "T"))
    return {k: "".join(v) for k, v in seqs.items() if v}


def _iter_read_sequences(path: Path) -> Iterable[str]:
    """Stream nucleotide read sequences from FASTQ/FASTQ.gz or FASTA/FASTA.gz."""
    with _open_text(path) as handle:
        first = ""
        for ln in handle:
            if ln.strip():
                first = ln.strip()
                break
        if not first:
            return
        if first.startswith("@"):
            # FASTQ format: 4 lines per read
            seq_line = next(handle, "").strip()
            if seq_line:
                yield seq_line.upper().replace("U", "T")
            next(handle, None)
            next(handle, None)
            while True:
                hdr = next(handle, None)
                if hdr is None:
                    break
                seq_line = next(handle, "").strip()
                next(handle, None)
                next(handle, None)
                if seq_line:
                    yield seq_line.upper().replace("U", "T")
        elif first.startswith(">"):
            buf: list[str] = []
            for ln in handle:
                s = ln.strip()
                if not s:
                    continue
                if s.startswith(">"):
                    if buf:
                        yield "".join(buf)
                        buf = []
                else:
                    buf.append(s.upper().replace("U", "T"))
            if buf:
                yield "".join(buf)
        else:
            # Plain one-sequence-per-line fallback
            yield first.upper().replace("U", "T")
            for ln in handle:
                s = ln.strip()
                if s:
                    yield s.upper().replace("U", "T")


def estimate_bulk_mhc_ase_from_reads(
    allele_fasta_path: str | Path,
    reads_paths: str | Path | Iterable[str | Path],
    *,
    tumor_purity: float | None = None,
    policy: MHCExpressionPolicy | dict[str, Any] | None = None,
) -> MHCAlleleStatusBundle:
    """Directly estimate cross-species bulk RNA-seq allele-specific MHC
    expression by counting unambiguous allele-discriminating k-mers across
    patient MHC allele reference sequences (IPD-IMGT/HLA or IPD-MHC FASTA).

    Works identically for human HLA, canine DLA, murine H-2, swine SLA, etc.
    Always preserves tumor-purity confounding guards (never claims tumor
    silencing on low-purity bulk RNA without caveats).
    """
    pol = MHCExpressionPolicy.resolve(policy)
    purity = tumor_purity if tumor_purity is not None else pol.default_tumor_purity
    fasta_p = Path(allele_fasta_path).expanduser().resolve()
    if not fasta_p.is_file():
        raise FileNotFoundError(f"Patient MHC allele FASTA not found: {fasta_p}")

    allele_seqs = _read_fasta_sequences(fasta_p)
    if not allele_seqs:
        raise ValueError(f"No allele sequences found in {fasta_p}")

    k = pol.kmer_size
    # Build canonical k-mer -> set of alleles containing that k-mer (fwd or rc)
    kmer_to_alleles: dict[str, set[str]] = {}
    for allele, seq in allele_seqs.items():
        if len(seq) < k:
            continue
        for i in range(len(seq) - k + 1):
            km = seq[i:i + k]
            if "N" in km:
                continue
            rc = _revcomp(km)
            canon_km = km if km <= rc else rc
            kmer_to_alleles.setdefault(canon_km, set()).add(allele)

    # Retain only allele-discriminating k-mers (unique to a single patient allele)
    unique_kmer_owner: dict[str, str] = {}
    disc_kmer_counts: dict[str, int] = {a: 0 for a in allele_seqs}
    for km, owners in kmer_to_alleles.items():
        if len(owners) == 1:
            owner = next(iter(owners))
            unique_kmer_owner[km] = owner
            disc_kmer_counts[owner] += 1

    if isinstance(reads_paths, (str, Path)):
        read_files = [Path(reads_paths).expanduser().resolve()]
    else:
        read_files = [Path(rp).expanduser().resolve() for rp in reads_paths]
    for rf in read_files:
        if not rf.is_file():
            raise FileNotFoundError(f"Bulk RNA reads file not found: {rf}")

    allele_read_counts: dict[str, int] = {a: 0 for a in allele_seqs}
    total_reads = 0
    for rf in read_files:
        for seq in _iter_read_sequences(rf):
            total_reads += 1
            if len(seq) < k:
                continue
            matched_owners: set[str] = set()
            for i in range(len(seq) - k + 1):
                km = seq[i:i + k]
                rc = _revcomp(km)
                canon_km = km if km <= rc else rc
                owner = unique_kmer_owner.get(canon_km)
                if owner is not None:
                    matched_owners.add(owner)
                    if len(matched_owners) > 1:
                        break
            if len(matched_owners) == 1:
                only = next(iter(matched_owners))
                allele_read_counts[only] += 1

    # Group by MHC gene locus and compute normalized within-locus ASE ratios
    by_locus: dict[str, list[str]] = {}
    for allele in allele_seqs:
        by_locus.setdefault(_allele_locus_key(allele), []).append(allele)

    records: list[MHCAlleleStatus] = []
    for locus, alleles_in_locus in sorted(by_locus.items()):
        locus_reads = sum(allele_read_counts[a] for a in alleles_in_locus)
        norm_expr = {
            a: (
                allele_read_counts[a] / max(1, disc_kmer_counts[a])
                if disc_kmer_counts[a] >= pol.min_discriminating_kmers
                else 0.0
            )
            for a in alleles_in_locus
        }
        locus_norm_total = sum(norm_expr.values())

        for a in alleles_in_locus:
            n_disc = disc_kmer_counts[a]
            n_reads = allele_read_counts[a]
            # RPM-like proxy scaled by discriminating k-mer length in kb
            tpm_proxy = (
                round((n_reads / max(1, total_reads)) * 1e6 / max(0.05, n_disc / 1000.0), 3)
                if total_reads > 0 and n_disc >= pol.min_discriminating_kmers
                else 0.0
            )
            if n_disc < pol.min_discriminating_kmers:
                records.append(
                    MHCAlleleStatus(
                        allele=a,
                        status="unresolvable",
                        bulk_rna_tpm=None,
                        bulk_rna_allele_ratio=None,
                        tumor_purity=purity,
                        assay_source="vaxrank_kmer_bulk_ase",
                        notes=(
                            f"locus={locus}; insufficient_discriminating_kmers "
                            f"({n_disc} < {pol.min_discriminating_kmers})"
                        ),
                    )
                )
            elif locus_reads < pol.min_locus_supporting_reads:
                records.append(
                    MHCAlleleStatus(
                        allele=a,
                        status="unresolvable",
                        bulk_rna_tpm=tpm_proxy,
                        bulk_rna_allele_ratio=None,
                        tumor_purity=purity,
                        assay_source="vaxrank_kmer_bulk_ase",
                        notes=(
                            f"locus={locus}; insufficient_locus_reads "
                            f"({locus_reads} < {pol.min_locus_supporting_reads}); "
                            f"reads={n_reads}, disc_kmers={n_disc}"
                        ),
                    )
                )
            else:
                ratio = norm_expr[a] / locus_norm_total if locus_norm_total > 0 else 0.0
                records.append(
                    MHCAlleleStatus(
                        allele=a,
                        status="unknown",
                        bulk_rna_tpm=tpm_proxy,
                        bulk_rna_allele_ratio=round(ratio, 6),
                        tumor_purity=purity,
                        assay_source="vaxrank_kmer_bulk_ase",
                        notes=f"locus={locus}; reads={n_reads}; disc_kmers={n_disc}",
                    )
                )

    return MHCAlleleStatusBundle(
        alleles=tuple(records),
        source_path=f"kmer_ase:{fasta_p.name}",
    )


# ---------------------------------------------------------------------------
# 3. Isovar-assisted local RNA assembly around polymorphic MHC exon variants
# ---------------------------------------------------------------------------


def estimate_bulk_mhc_ase_from_isovar(
    allele_variant_map: dict[str, dict[str, str]],
    isovar_results: Iterable[Any],
    *,
    tumor_purity: float | None = None,
    policy: MHCExpressionPolicy | dict[str, Any] | None = None,
) -> MHCAlleleStatusBundle:
    """Aggregate ``isovar`` local RNA assembly read counts across polymorphic
    MHC exon variants distinguishing a patient's heterozygous MHC alleles.

    Parameters
    ----------
    allele_variant_map : dict[str, dict[str, str]]
        Maps ``variant.short_description`` -> ``{"ref_allele": ..., "alt_allele": ...}``
        identifying which patient MHC allele carries the reference vs alternate
        nucleotide at each polymorphic MHC exon site (cross-species).
    isovar_results : iterable of IsovarResult
        Results from ``isovar.run_isovar(...)`` on the tumor RNA BAM.
    """
    pol = MHCExpressionPolicy.resolve(policy)
    purity = tumor_purity if tumor_purity is not None else pol.default_tumor_purity

    allele_reads: dict[str, int] = {}
    marker_counts: dict[str, int] = {}
    for res in isovar_results or ():
        v = getattr(res, "variant", None)
        desc = str(getattr(v, "short_description", None) or v or "")
        mapping = allele_variant_map.get(desc)
        if not mapping:
            continue
        ref_a = normalize_mhc_allele_name(mapping["ref_allele"])
        alt_a = normalize_mhc_allele_name(mapping["alt_allele"])
        n_ref = getattr(res, "num_ref_supporting_reads", None)
        if n_ref is None:
            n_ref = getattr(res, "n_ref_reads", 0)
        n_alt = getattr(res, "num_alt_supporting_reads", None)
        if n_alt is None:
            n_alt = getattr(res, "n_alt_reads", 0)
        allele_reads[ref_a] = allele_reads.get(ref_a, 0) + int(n_ref or 0)
        allele_reads[alt_a] = allele_reads.get(alt_a, 0) + int(n_alt or 0)
        marker_counts[ref_a] = marker_counts.get(ref_a, 0) + 1
        marker_counts[alt_a] = marker_counts.get(alt_a, 0) + 1

    by_locus: dict[str, list[str]] = {}
    for a in allele_reads:
        by_locus.setdefault(_allele_locus_key(a), []).append(a)

    records: list[MHCAlleleStatus] = []
    for locus, alleles_in_locus in sorted(by_locus.items()):
        locus_total = sum(allele_reads[a] for a in alleles_in_locus)
        for a in alleles_in_locus:
            cnt = allele_reads[a]
            if locus_total < pol.min_locus_supporting_reads:
                records.append(
                    MHCAlleleStatus(
                        allele=a,
                        status="unresolvable",
                        tumor_purity=purity,
                        assay_source="isovar_mhc_local_assembly",
                        notes=(
                            f"locus={locus}; insufficient_locus_reads "
                            f"({locus_total} < {pol.min_locus_supporting_reads}); "
                            f"reads={cnt}; markers={marker_counts.get(a, 0)}"
                        ),
                    )
                )
            else:
                records.append(
                    MHCAlleleStatus(
                        allele=a,
                        status="unknown",
                        bulk_rna_allele_ratio=round(cnt / locus_total, 6),
                        tumor_purity=purity,
                        assay_source="isovar_mhc_local_assembly",
                        notes=f"locus={locus}; reads={cnt}; markers={marker_counts.get(a, 0)}",
                    )
                )
    return MHCAlleleStatusBundle(
        alleles=tuple(records),
        source_path="isovar:mhc_local_assembly",
    )


def run_isovar_mhc_ase(
    allele_variant_map: dict[str, dict[str, str]],
    variants: Any,
    alignment_file: Any,
    *,
    tumor_purity: float | None = None,
    policy: MHCExpressionPolicy | dict[str, Any] | None = None,
    **isovar_kwargs,
) -> MHCAlleleStatusBundle:
    """Invoke ``isovar.run_isovar`` on polymorphic MHC exon variants and a
    tumor RNA BAM, then quantify allele-specific MHC expression."""
    from isovar import run_isovar

    results = run_isovar(
        variants=variants,
        alignment_file=alignment_file,
        **isovar_kwargs,
    )
    return estimate_bulk_mhc_ase_from_isovar(
        allele_variant_map,
        results,
        tumor_purity=tumor_purity,
        policy=policy,
    )
