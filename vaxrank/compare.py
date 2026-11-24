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

"""Structured run-to-run, policy-to-policy, and RNA-vs-DNA comparison reports
(``vaxrank compare``; openvax/vaxrank#920)."""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from numbers import Real
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class ComparePolicy:
    """Explicit per-trial policy governing ``vaxrank compare`` thresholds.

    Configurable under ``compare:`` in a per-trial YAML policy or via CLI flags.
    """
    top_k_values: tuple[int, ...] = (5, 10, 20)
    primary_top_k: int = 10
    score_relative_change_threshold: float = 0.25
    expression_fold_change_threshold: float = 1.5
    rank_shift_highlight_threshold: int = 3
    include_pdf: bool = False

    def __post_init__(self):
        if isinstance(self.top_k_values, (str, bytes)):
            raise ValueError("compare.top_k_values must be a sequence of integers >= 1")
        raw_ks = tuple(self.top_k_values)
        if not raw_ks or any(
            isinstance(k, bool) or not isinstance(k, int) or k < 1 for k in raw_ks
        ):
            raise ValueError("compare.top_k_values must be a non-empty sequence of integers >= 1")
        object.__setattr__(self, "top_k_values", tuple(sorted(set(raw_ks))))
        if (isinstance(self.primary_top_k, bool)
                or not isinstance(self.primary_top_k, int)
                or self.primary_top_k < 1):
            raise ValueError("compare.primary_top_k must be an integer >= 1")
        if self.primary_top_k not in self.top_k_values:
            object.__setattr__(
                self,
                "top_k_values",
                tuple(sorted(set(self.top_k_values) | {self.primary_top_k})),
            )
        for name in ("score_relative_change_threshold",):
            val = getattr(self, name)
            if (isinstance(val, bool)
                    or not isinstance(val, Real)
                    or not math.isfinite(val)
                    or float(val) < 0.0):
                raise ValueError(f"compare.{name} must be a finite number >= 0")
        if (isinstance(self.expression_fold_change_threshold, bool)
                or not isinstance(self.expression_fold_change_threshold, Real)
                or not math.isfinite(self.expression_fold_change_threshold)
                or float(self.expression_fold_change_threshold) <= 1.0):
            raise ValueError("compare.expression_fold_change_threshold must be > 1.0")
        if (isinstance(self.rank_shift_highlight_threshold, bool)
                or not isinstance(self.rank_shift_highlight_threshold, int)
                or self.rank_shift_highlight_threshold < 1):
            raise ValueError("compare.rank_shift_highlight_threshold must be an integer >= 1")
        if not isinstance(self.include_pdf, bool):
            raise ValueError("compare.include_pdf must be a boolean")

    @classmethod
    def resolve(cls, value=None) -> ComparePolicy:
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if not isinstance(value, dict):
            raise TypeError(
                f"compare policy must be a mapping or ComparePolicy; got {type(value).__name__}"
            )
        defaults = cls().to_dict()
        unknown = set(value) - set(defaults)
        if unknown:
            raise ValueError("Unknown compare policy fields: %s" % ", ".join(sorted(unknown)))
        merged = {k: (defaults[k] if value.get(k) is None else value[k]) for k in defaults}
        if isinstance(merged["top_k_values"], list):
            merged["top_k_values"] = tuple(merged["top_k_values"])
        return cls(**merged)

    def to_dict(self) -> dict[str, Any]:
        return {
            "top_k_values": list(self.top_k_values),
            "primary_top_k": self.primary_top_k,
            "score_relative_change_threshold": float(self.score_relative_change_threshold),
            "expression_fold_change_threshold": float(self.expression_fold_change_threshold),
            "rank_shift_highlight_threshold": self.rank_shift_highlight_threshold,
            "include_pdf": self.include_pdf,
        }


def _parse_float(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        f = float(value)
        return f if math.isfinite(f) else None
    text = str(value).strip()
    if not text or text.lower() in ("nan", "none", "null", "na", "n/a", ""):
        return None
    try:
        f = float(text)
    except ValueError:
        return None
    return f if math.isfinite(f) else None


def _parse_int(value: Any) -> int | None:
    f = _parse_float(value)
    if f is None:
        return None
    return int(round(f))


def _target_key_from_row(row: dict[str, Any]) -> str:
    chrom = str(row.get("chr") or "").strip()
    pos = _parse_int(row.get("pos"))
    ref = str(row.get("ref") or "").strip().upper()
    alt = str(row.get("alt") or "").strip().upper()
    if chrom and pos is not None and (ref or alt):
        norm_chr = chrom[3:] if chrom.lower().startswith("chr") else chrom
        return f"{norm_chr}:{pos}:{ref}>{alt}"
    source_id = str(row.get("source_identifier") or row.get("variant") or "").strip()
    if source_id:
        return source_id
    gene = str(row.get("gene_name") or "").strip()
    aa = str(row.get("amino_acids") or row.get("peptide") or "").strip()
    return f"{gene}:{aa}" if gene else aa


def _load_predictions_tsv(tsv_path: Path) -> tuple[list[dict[str, Any]], dict[str, Any] | None]:
    """Load target rankings and optional embedded policy metadata from a TSV
    prediction table (``candidate_predictions.tsv`` / ``neoepitopes.tsv``)."""
    policy_meta: dict[str, Any] | None = None
    data_lines: list[str] = []
    with open(tsv_path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                # Inspect Topiary/vaxrank JSON comment headers for policy metadata
                stripped = line.lstrip("#").strip()
                if stripped.startswith("{") and ("selection_policy" in stripped or "run_configuration" in stripped):
                    try:
                        parsed = json.loads(stripped)
                        if isinstance(parsed, dict):
                            policy_meta = parsed.get("run_configuration") or parsed
                    except json.JSONDecodeError:
                        pass
                continue
            if line.strip():
                data_lines.append(line)

    if not data_lines:
        return [], policy_meta

    reader = csv.DictReader(data_lines, delimiter="\t")
    by_key: dict[str, dict[str, Any]] = {}
    for row in reader:
        key = (
            str(row.get("source_identifier") or row.get("variant") or row.get("source_sequence_name") or "").strip()
            or _target_key_from_row(row)
        )
        if not key:
            continue
        score = _parse_float(
            row.get("combined_score")
            or row.get("score")
            or row.get("target_epitope_score")
            or row.get("value")
        )
        expr = _parse_float(
            row.get("expression_score") or row.get("Variant allele RNA read count")
        )
        gene = str(row.get("gene_name") or row.get("gene") or "").strip()
        aa = str(row.get("amino_acids") or row.get("peptide") or "").strip()
        existing = by_key.get(key)
        if existing is None or (
            score is not None and (existing["combined_score"] is None or score > existing["combined_score"])
        ):
            by_key[key] = {
                "target_key": key,
                "source_rank": 0,
                "source_identifier": key,
                "gene_name": gene,
                "source_kind": str(row.get("source_kind") or "prediction_tsv").strip(),
                "variant_types": str(row.get("variant_types") or "").strip(),
                "rna_product": str(row.get("rna_product") or "").strip(),
                "amino_acids": aa,
                "combined_score": score,
                "expression_score": expr,
                "target_epitope_score": _parse_float(row.get("target_epitope_score") or row.get("score")),
                "biological_evidence": "",
                "cta_targeting": "",
            }

    records = sorted(
        by_key.values(),
        key=lambda r: (
            -(r["combined_score"] if r["combined_score"] is not None else float("-inf")),
            r["target_key"],
        ),
    )
    for rank_idx, rec in enumerate(records, start=1):
        rec["source_rank"] = rank_idx
        rec["rank"] = rank_idx
    return records, policy_meta


def _load_ranked_csv(csv_path: Path) -> list[dict[str, Any]]:
    records_by_key: dict[str, dict[str, Any]] = {}
    ordered_keys: list[str] = []
    with open(csv_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for idx, row in enumerate(reader, start=1):
            window_rank = _parse_int(row.get("window_rank"))
            if window_rank is not None and window_rank > 1:
                # Primary per-source record is window_rank == 1
                continue
            key = _target_key_from_row(row)
            if not key or key in records_by_key:
                continue
            source_rank = _parse_int(row.get("source_rank"))
            if source_rank is None:
                source_rank = len(ordered_keys) + 1
            rec = {
                "target_key": key,
                "source_rank": source_rank,
                "source_identifier": str(row.get("source_identifier") or key).strip(),
                "gene_name": str(row.get("gene_name") or "").strip(),
                "source_kind": str(row.get("source_kind") or "").strip(),
                "variant_types": str(row.get("variant_types") or "").strip(),
                "rna_product": str(row.get("rna_product") or "").strip(),
                "amino_acids": str(row.get("amino_acids") or "").strip(),
                "combined_score": _parse_float(row.get("combined_score")),
                "expression_score": _parse_float(row.get("expression_score")),
                "target_epitope_score": _parse_float(row.get("target_epitope_score")),
                "biological_evidence": str(row.get("biological_evidence") or "").strip(),
                "cta_targeting": str(row.get("cta_targeting") or "").strip(),
            }
            records_by_key[key] = rec
            ordered_keys.append(key)
    records = [records_by_key[k] for k in ordered_keys]
    records.sort(key=lambda r: (r["source_rank"], r["target_key"]))
    for rank_idx, rec in enumerate(records, start=1):
        rec["rank"] = rank_idx
    return records


def _load_manifest(manifest_path: Path) -> list[dict[str, Any]]:
    with open(manifest_path, encoding="utf-8") as handle:
        data = json.load(handle)
    if isinstance(data, list):
        return [item for item in data if isinstance(item, dict)]
    if isinstance(data, dict) and isinstance(data.get("constructs"), list):
        return [item for item in data["constructs"] if isinstance(item, dict)]
    return []


def _find_ranked_csv(root: Path) -> Path | None:
    for name in (
        "ranked_candidates.csv",
        "vaccine-peptides.csv",
        "vaccine_peptides.csv",
        "ranked-candidates.csv",
    ):
        cand = root / name
        if cand.is_file():
            return cand
    # Fallback: scan top-level CSV files for a vaxrank ranked schema
    for cand in sorted(root.glob("*.csv")):
        try:
            with open(cand, newline="", encoding="utf-8") as handle:
                header = next(csv.reader(handle), [])
        except OSError:
            continue
        if "amino_acids" in header and ("source_rank" in header or "combined_score" in header):
            return cand
    return None


def _find_selection_policy(root: Path) -> dict[str, Any] | None:
    for rel in ("selection_policy.json", "policy/selection_policy.json"):
        cand = root / rel
        if cand.is_file():
            try:
                with open(cand, encoding="utf-8") as handle:
                    data = json.load(handle)
                if isinstance(data, dict):
                    return data
            except (OSError, json.JSONDecodeError):
                return None
    return None


def load_run_bundle(path: str | Path, label: str | None = None) -> dict[str, Any]:
    """Load ranked candidates, emitted constructs, and policy metadata from a
    vaxrank output directory or a single ranked CSV / manifest JSON file.
    """
    p = Path(path).expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(f"Comparison input path does not exist: {p}")

    resolved_label = label or p.name
    ranked: list[dict[str, Any]] = []
    peptide_manifest: list[dict[str, Any]] = []
    mrna_manifest: list[dict[str, Any]] = []
    policy_meta: dict[str, Any] | None = None
    files_loaded: dict[str, str] = {}

    if p.is_file():
        suffix = p.suffix.lower()
        if suffix == ".csv":
            ranked = _load_ranked_csv(p)
            files_loaded["ranked_csv"] = str(p)
        elif suffix == ".tsv":
            ranked, policy_meta = _load_predictions_tsv(p)
            files_loaded["predictions_tsv"] = str(p)
        elif suffix == ".json":
            entries = _load_manifest(p)
            for item in entries:
                mod = item.get("modality")
                if mod == "mrna":
                    mrna_manifest.append(item)
                else:
                    peptide_manifest.append(item)
            files_loaded["manifest_json"] = str(p)
        else:
            raise ValueError(
                f"Unsupported single-file comparison input {p}; pass a vaxrank output "
                f"directory, ranked_candidates.csv, candidate_predictions.tsv, or manifest.json."
            )
    else:
        csv_path = _find_ranked_csv(p)
        if csv_path is not None:
            ranked = _load_ranked_csv(csv_path)
            files_loaded["ranked_csv"] = str(csv_path)

        tsv_path = p / "candidate_predictions.tsv"
        if tsv_path.is_file():
            tsv_ranked, tsv_policy = _load_predictions_tsv(tsv_path)
            files_loaded["predictions_tsv"] = str(tsv_path)
            if not ranked and tsv_ranked:
                ranked = tsv_ranked
            if policy_meta is None and tsv_policy is not None:
                policy_meta = tsv_policy

        pep_path = p / "peptide" / "manifest.json"
        if pep_path.is_file():
            peptide_manifest = _load_manifest(pep_path)
            files_loaded["peptide_manifest"] = str(pep_path)

        mrna_path = p / "mrna" / "manifest.json"
        if mrna_path.is_file():
            mrna_manifest = _load_manifest(mrna_path)
            files_loaded["mrna_manifest"] = str(mrna_path)

        top_manifest = p / "manifest.json"
        if top_manifest.is_file() and not (peptide_manifest or mrna_manifest):
            entries = _load_manifest(top_manifest)
            for item in entries:
                if item.get("modality") == "mrna":
                    mrna_manifest.append(item)
                else:
                    peptide_manifest.append(item)
            files_loaded["manifest_json"] = str(top_manifest)

        found_policy = _find_selection_policy(p)
        if found_policy is not None:
            policy_meta = found_policy
            files_loaded["selection_policy"] = str(
                (p / "selection_policy.json")
                if (p / "selection_policy.json").is_file()
                else (p / "policy" / "selection_policy.json")
            )

    if not ranked and not peptide_manifest and not mrna_manifest:
        raise ValueError(
            f"No ranked candidates CSV or construct manifest.json found under {p}"
        )

    return {
        "label": resolved_label,
        "path": str(p),
        "ranked": ranked,
        "peptide_manifest": peptide_manifest,
        "mrna_manifest": mrna_manifest,
        "policy": policy_meta,
        "files_loaded": files_loaded,
    }


def _average_ranks(values: list[float]) -> list[float]:
    n = len(values)
    indexed = sorted(range(n), key=lambda i: values[i])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and values[indexed[j + 1]] == values[indexed[i]]:
            j += 1
        avg_rank = 0.5 * (i + j) + 1.0
        for k in range(i, j + 1):
            ranks[indexed[k]] = avg_rank
        i = j + 1
    return ranks


def _spearman_rho(xs: list[float], ys: list[float]) -> float | None:
    n = len(xs)
    if n < 2:
        return None
    rx = _average_ranks(xs)
    ry = _average_ranks(ys)
    mx = sum(rx) / n
    my = sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    if dx == 0.0 or dy == 0.0:
        return None
    return max(-1.0, min(1.0, num / (dx * dy)))


def _kendall_tau_b(xs: list[float], ys: list[float]) -> float | None:
    n = len(xs)
    if n < 2:
        return None
    concordant = 0
    discordant = 0
    ties_x = 0
    ties_y = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            dx = xs[i] - xs[j]
            dy = ys[i] - ys[j]
            if dx == 0 and dy == 0:
                continue
            if dx == 0:
                ties_x += 1
            elif dy == 0:
                ties_y += 1
            elif (dx > 0 and dy > 0) or (dx < 0 and dy < 0):
                concordant += 1
            else:
                discordant += 1
    denom = math.sqrt(
        (concordant + discordant + ties_x) * (concordant + discordant + ties_y)
    )
    if denom == 0.0:
        return None
    return max(-1.0, min(1.0, (concordant - discordant) / denom))


def _is_proxy_expression(rec: dict[str, Any] | None) -> bool:
    if not rec:
        return False
    ev = (rec.get("biological_evidence") or "").lower()
    return "proxy" in ev or "dna-only" in ev or "dna_only" in ev


def _has_safety_flag(rec: dict[str, Any] | None) -> bool:
    if not rec:
        return False
    text = f"{rec.get('biological_evidence', '')} {rec.get('cta_targeting', '')}".lower()
    return "unsafe" in text or "self_overlap" in text or "excluded" in text


def _attribute_target_change(
    rec_a: dict[str, Any] | None,
    rec_b: dict[str, Any] | None,
    *,
    policy: ComparePolicy,
    policies_differ: bool,
) -> list[str]:
    reasons: list[str] = []

    if rec_a is None and rec_b is not None:
        if _is_proxy_expression(rec_b):
            reasons.append("gained_via_dna_proxy_expression")
        elif rec_b.get("rna_product") or (rec_b.get("expression_score") or 0) > 0:
            reasons.append("gained_with_rna_or_filter_inclusion")
        else:
            reasons.append("only_in_run_b")
        return reasons

    if rec_b is None and rec_a is not None:
        if _has_safety_flag(rec_a):
            reasons.append("dropped_safety_or_cta_exclusion")
        elif _is_proxy_expression(rec_a):
            reasons.append("dropped_without_rna_support")
        else:
            reasons.append("only_in_run_a")
        return reasons

    assert rec_a is not None and rec_b is not None
    seq_a = rec_a.get("amino_acids") or ""
    seq_b = rec_b.get("amino_acids") or ""
    if seq_a and seq_b and seq_a != seq_b:
        reasons.append("peptide_window_shift")

    proxy_a = _is_proxy_expression(rec_a)
    proxy_b = _is_proxy_expression(rec_b)
    if proxy_a != proxy_b:
        reasons.append("expression_source_change_rna_vs_proxy")

    expr_a = rec_a.get("expression_score")
    expr_b = rec_b.get("expression_score")
    if expr_a is not None and expr_b is not None:
        if (expr_a == 0) != (expr_b == 0):
            reasons.append("expression_support_shift")
        elif expr_a > 0 and expr_b > 0:
            fold = max(expr_a / expr_b, expr_b / expr_a)
            if fold >= policy.expression_fold_change_threshold:
                reasons.append("expression_tpm_or_read_shift")
    elif (expr_a is None) != (expr_b is None):
        reasons.append("expression_availability_change")

    epi_a = rec_a.get("target_epitope_score")
    epi_b = rec_b.get("target_epitope_score")
    if epi_a is not None and epi_b is not None:
        denom = max(abs(epi_a), abs(epi_b), 1e-9)
        if abs(epi_a - epi_b) / denom >= policy.score_relative_change_threshold:
            reasons.append("epitope_score_shift")

    score_a = rec_a.get("combined_score")
    score_b = rec_b.get("combined_score")
    if score_a is not None and score_b is not None:
        denom = max(abs(score_a), abs(score_b), 1e-9)
        if abs(score_a - score_b) / denom >= policy.score_relative_change_threshold:
            if "epitope_score_shift" not in reasons:
                reasons.append("policy_or_combined_score_shift")

    if policies_differ and not reasons and rec_a["rank"] != rec_b["rank"]:
        reasons.append("policy_change")

    if not reasons:
        if rec_a["rank"] != rec_b["rank"]:
            reasons.append("rank_displacement")
        else:
            reasons.append("unchanged")
    return reasons


def _extract_peptide_slots(manifest: list[dict[str, Any]]) -> list[dict[str, Any]]:
    slots = []
    for entry in manifest:
        seq = str(entry.get("sequence") or "").strip()
        names = list(entry.get("antigen_names") or [])
        label = names[0] if names else str(entry.get("name") or "")
        slots.append({
            "construct_name": str(entry.get("name") or ""),
            "antigen_label": label,
            "sequence": seq,
            "pool": entry.get("pool"),
        })
    return slots


def _extract_mrna_slots(manifest: list[dict[str, Any]]) -> list[dict[str, Any]]:
    slots = []
    for entry in manifest:
        cname = str(entry.get("name") or "")
        antigens = entry.get("antigens") or []
        if antigens:
            for ag in antigens:
                slots.append({
                    "construct_name": cname,
                    "antigen_label": str(ag.get("name") or ""),
                    "sequence": str(ag.get("aa") or ""),
                    "genes": list(ag.get("gene_names") or []),
                })
        else:
            for name in entry.get("antigen_names") or []:
                slots.append({
                    "construct_name": cname,
                    "antigen_label": str(name),
                    "sequence": "",
                    "genes": [],
                })
    return slots


def _sequences_overlap(seq_a: str, seq_b: str, min_overlap: int = 8) -> bool:
    if not seq_a or not seq_b:
        return False
    if seq_a in seq_b or seq_b in seq_a:
        return True
    k = min(min_overlap, len(seq_a), len(seq_b))
    if k < 4:
        return False
    kmers_a = {seq_a[i:i + k] for i in range(len(seq_a) - k + 1)}
    return any(seq_b[i:i + k] in kmers_a for i in range(len(seq_b) - k + 1))


def _diff_modality_slots(
    slots_a: list[dict[str, Any]],
    slots_b: list[dict[str, Any]],
    *,
    construct_count_a: int,
    construct_count_b: int,
) -> dict[str, Any]:
    if not slots_a and not slots_b:
        return {"available": False}

    # 1-to-1 exact sequence matching preserves exact slot accounting when
    # duplicate sequences appear across constructs.
    exact_matches = []
    matched_a_exact: set[int] = set()
    matched_b_exact: set[int] = set()
    for idx_a, sa in enumerate(slots_a):
        seq = sa["sequence"]
        if not seq:
            continue
        for idx_b, sb in enumerate(slots_b):
            if idx_b in matched_b_exact:
                continue
            if sb["sequence"] == seq:
                matched_a_exact.add(idx_a)
                matched_b_exact.add(idx_b)
                exact_matches.append({
                    "sequence": seq,
                    "antigen_a": sa["antigen_label"],
                    "antigen_b": sb["antigen_label"],
                    "construct_a": sa["construct_name"],
                    "construct_b": sb["construct_name"],
                    **({"pool_a": sa.get("pool"), "pool_b": sb.get("pool")}
                       if (sa.get("pool") is not None or sb.get("pool") is not None)
                       else {}),
                })
                break

    rem_a = [sa for idx_a, sa in enumerate(slots_a) if idx_a not in matched_a_exact]
    rem_b = [sb for idx_b, sb in enumerate(slots_b) if idx_b not in matched_b_exact]

    shifted_windows = []
    matched_b_indices: set[int] = set()
    dropped = []
    for sa in rem_a:
        match_idx = None
        for idx_b, sb in enumerate(rem_b):
            if idx_b in matched_b_indices:
                continue
            same_label = bool(sa["antigen_label"] and sa["antigen_label"] == sb["antigen_label"])
            overlap = _sequences_overlap(sa["sequence"], sb["sequence"])
            if same_label or overlap:
                match_idx = idx_b
                break
        if match_idx is not None:
            matched_b_indices.add(match_idx)
            sb = rem_b[match_idx]
            shifted_windows.append({
                "antigen_a": sa["antigen_label"],
                "antigen_b": sb["antigen_label"],
                "sequence_a": sa["sequence"],
                "sequence_b": sb["sequence"],
                "construct_a": sa["construct_name"],
                "construct_b": sb["construct_name"],
            })
        else:
            dropped.append(sa)

    added = [sb for idx_b, sb in enumerate(rem_b) if idx_b not in matched_b_indices]

    return {
        "available": True,
        "construct_count_a": construct_count_a,
        "construct_count_b": construct_count_b,
        "slot_count_a": len(slots_a),
        "slot_count_b": len(slots_b),
        "exact_match_count": len(exact_matches),
        "shifted_window_count": len(shifted_windows),
        "dropped_count": len(dropped),
        "added_count": len(added),
        "exact_matches": exact_matches,
        "shifted_windows": shifted_windows,
        "dropped_from_a": dropped,
        "added_in_b": added,
    }


def _policy_summary(policy_meta: dict[str, Any] | None) -> dict[str, Any]:
    if not policy_meta:
        return {}
    out: dict[str, Any] = {}
    for key in ("name", "policy_name", "digest", "sha256", "policy_sha256"):
        if key in policy_meta and policy_meta[key] is not None:
            out[key] = policy_meta[key]
    vp = policy_meta.get("vaccine_peptides") or {}
    if isinstance(vp, dict) and vp.get("combined_score_expr"):
        out["combined_score_expr"] = vp["combined_score_expr"]
    return out


def compare_runs(
    run_a: str | Path | dict[str, Any],
    run_b: str | Path | dict[str, Any],
    *,
    policy: ComparePolicy | dict[str, Any] | None = None,
    label_a: str | None = None,
    label_b: str | None = None,
) -> dict[str, Any]:
    """Compare two vaxrank runs and return a structured comparison summary."""
    resolved_policy = ComparePolicy.resolve(policy)
    bundle_a = (
        run_a if isinstance(run_a, dict) and "ranked" in run_a
        else load_run_bundle(run_a, label=label_a)
    )
    bundle_b = (
        run_b if isinstance(run_b, dict) and "ranked" in run_b
        else load_run_bundle(run_b, label=label_b)
    )
    if label_a:
        bundle_a = dict(bundle_a, label=label_a)
    if label_b:
        bundle_b = dict(bundle_b, label=label_b)

    ranked_a = bundle_a["ranked"]
    ranked_b = bundle_b["ranked"]
    by_key_a = {r["target_key"]: r for r in ranked_a}
    by_key_b = {r["target_key"]: r for r in ranked_b}

    keys_a = [r["target_key"] for r in ranked_a]
    keys_b = [r["target_key"] for r in ranked_b]
    set_a = set(keys_a)
    set_b = set(keys_b)
    shared_keys = [k for k in keys_a if k in set_b]
    only_a_keys = [k for k in keys_a if k not in set_b]
    only_b_keys = [k for k in keys_b if k not in set_a]

    shared_ranks_a = [float(by_key_a[k]["rank"]) for k in shared_keys]
    shared_ranks_b = [float(by_key_b[k]["rank"]) for k in shared_keys]
    spearman = _spearman_rho(shared_ranks_a, shared_ranks_b)
    kendall = _kendall_tau_b(shared_ranks_a, shared_ranks_b)

    pol_a = _policy_summary(bundle_a.get("policy"))
    pol_b = _policy_summary(bundle_b.get("policy"))
    policies_differ = bool(pol_a and pol_b and pol_a != pol_b)

    # Top-K overlap statistics
    top_k_stats = {}
    for k in resolved_policy.top_k_values:
        tk_a = set(keys_a[:k])
        tk_b = set(keys_b[:k])
        overlap = tk_a & tk_b
        union = tk_a | tk_b
        exact_seq = sum(
            1 for key in overlap
            if (by_key_a[key]["amino_acids"] and
                by_key_a[key]["amino_acids"] == by_key_b[key]["amino_acids"])
        )
        shifted_seq = sum(
            1 for key in overlap
            if (by_key_a[key]["amino_acids"] and by_key_b[key]["amino_acids"] and
                by_key_a[key]["amino_acids"] != by_key_b[key]["amino_acids"])
        )
        top_k_stats[str(k)] = {
            "k": k,
            "count_a": len(tk_a),
            "count_b": len(tk_b),
            "overlap_count": len(overlap),
            "jaccard": round(len(overlap) / len(union), 6) if union else None,
            "gained_in_b_count": len(tk_b - tk_a),
            "lost_from_a_count": len(tk_a - tk_b),
            "exact_sequence_matches": exact_seq,
            "shifted_window_matches": shifted_seq,
        }

    # Build unified per-target comparison table
    pk = resolved_policy.primary_top_k
    pk_a = set(keys_a[:pk])
    pk_b = set(keys_b[:pk])
    all_keys = list(dict.fromkeys([*keys_a, *keys_b]))

    target_rows = []
    gained_top_k = []
    lost_top_k = []
    notable_movers = []
    shifted_window_targets = []
    attribution_counts: dict[str, int] = {}

    for key in all_keys:
        ra = by_key_a.get(key)
        rb = by_key_b.get(key)
        rank_a = ra["rank"] if ra else None
        rank_b = rb["rank"] if rb else None
        rank_delta = (rank_a - rank_b) if (rank_a is not None and rank_b is not None) else None
        in_pk_a = key in pk_a
        in_pk_b = key in pk_b
        if in_pk_a and in_pk_b:
            top_k_status = f"retained_top_{pk}"
        elif in_pk_b and not in_pk_a:
            top_k_status = f"gained_in_top_{pk}"
        elif in_pk_a and not in_pk_b:
            top_k_status = f"lost_from_top_{pk}"
        else:
            top_k_status = "outside_top_k"

        status = (
            "shared" if (ra and rb)
            else ("only_in_a" if ra else "only_in_b")
        )
        seq_a = ra["amino_acids"] if ra else ""
        seq_b = rb["amino_acids"] if rb else ""
        if ra and rb:
            seq_match = "exact" if (seq_a == seq_b) else "shifted"
        else:
            seq_match = "n/a"

        reasons = _attribute_target_change(
            ra, rb, policy=resolved_policy, policies_differ=policies_differ
        )
        if in_pk_a or in_pk_b:
            for reason in reasons:
                attribution_counts[reason] = attribution_counts.get(reason, 0) + 1

        row = {
            "target_key": key,
            "source_identifier": (rb or ra)["source_identifier"],
            "gene_name": (rb or ra)["gene_name"],
            "source_kind": (rb or ra)["source_kind"],
            "status": status,
            "top_k_status": top_k_status,
            "rank_a": rank_a,
            "rank_b": rank_b,
            "rank_delta": rank_delta,
            "score_a": ra["combined_score"] if ra else None,
            "score_b": rb["combined_score"] if rb else None,
            "expression_a": ra["expression_score"] if ra else None,
            "expression_b": rb["expression_score"] if rb else None,
            "epitope_score_a": ra["target_epitope_score"] if ra else None,
            "epitope_score_b": rb["target_epitope_score"] if rb else None,
            "amino_acids_a": seq_a,
            "amino_acids_b": seq_b,
            "sequence_match": seq_match,
            "attribution": reasons,
        }
        target_rows.append(row)
        if seq_match == "shifted":
            shifted_window_targets.append(row)
        if in_pk_b and not in_pk_a:
            gained_top_k.append(row)
        elif in_pk_a and not in_pk_b:
            lost_top_k.append(row)
        elif (rank_delta is not None
              and abs(rank_delta) >= resolved_policy.rank_shift_highlight_threshold):
            notable_movers.append(row)

    gained_top_k.sort(key=lambda r: (r["rank_b"] or 10**9, r["target_key"]))
    lost_top_k.sort(key=lambda r: (r["rank_a"] or 10**9, r["target_key"]))
    notable_movers.sort(key=lambda r: (-abs(r["rank_delta"] or 0), r["target_key"]))
    shifted_window_targets.sort(key=lambda r: (r["rank_b"] or 10**9, r["target_key"]))

    # Construct diffs
    peptide_diff = _diff_modality_slots(
        _extract_peptide_slots(bundle_a["peptide_manifest"]),
        _extract_peptide_slots(bundle_b["peptide_manifest"]),
        construct_count_a=len(bundle_a["peptide_manifest"]),
        construct_count_b=len(bundle_b["peptide_manifest"]),
    )
    mrna_diff = _diff_modality_slots(
        _extract_mrna_slots(bundle_a["mrna_manifest"]),
        _extract_mrna_slots(bundle_b["mrna_manifest"]),
        construct_count_a=len(bundle_a["mrna_manifest"]),
        construct_count_b=len(bundle_b["mrna_manifest"]),
    )

    return {
        "schema": "vaxrank.comparison.v1",
        "policy": resolved_policy.to_dict(),
        "run_a": {
            "label": bundle_a["label"],
            "path": bundle_a["path"],
            "ranked_count": len(ranked_a),
            "files_loaded": bundle_a["files_loaded"],
            "policy": pol_a,
        },
        "run_b": {
            "label": bundle_b["label"],
            "path": bundle_b["path"],
            "ranked_count": len(ranked_b),
            "files_loaded": bundle_b["files_loaded"],
            "policy": pol_b,
        },
        "policies_differ": policies_differ,
        "rank_preservation": {
            "total_unique_targets": len(all_keys),
            "shared_target_count": len(shared_keys),
            "only_in_a_count": len(only_a_keys),
            "only_in_b_count": len(only_b_keys),
            "spearman_rho": round(spearman, 6) if spearman is not None else None,
            "kendall_tau_b": round(kendall, 6) if kendall is not None else None,
            "top_k": top_k_stats,
        },
        "primary_top_k": pk,
        "attribution_summary": attribution_counts,
        "gained_in_top_k": gained_top_k,
        "lost_from_top_k": lost_top_k,
        "notable_rank_movers": notable_movers,
        "shifted_window_targets": shifted_window_targets,
        "construct_diffs": {
            "peptide": peptide_diff,
            "mrna": mrna_diff,
        },
        "targets": target_rows,
    }


def _fmt_num(value: float | None, digits: int = 3) -> str:
    if value is None:
        return "n/a"
    return f"{value:.{digits}f}"


def render_comparison_ascii(comparison: dict[str, Any]) -> str:
    """Render a concise ASCII diff report for two vaxrank runs."""
    ra = comparison["run_a"]
    rb = comparison["run_b"]
    rp = comparison["rank_preservation"]
    pk = comparison["primary_top_k"]

    lines = [
        "=== Vaxrank Run Comparison ===",
        f"Run A: {ra['label']} ({ra['ranked_count']} ranked targets) — {ra['path']}",
        f"Run B: {rb['label']} ({rb['ranked_count']} ranked targets) — {rb['path']}",
        "",
        "--- Rank Preservation ---",
        f"Shared targets: {rp['shared_target_count']} / {rp['total_unique_targets']} "
        f"(only in A: {rp['only_in_a_count']}, only in B: {rp['only_in_b_count']})",
        f"Spearman rho: {_fmt_num(rp['spearman_rho'])} | "
        f"Kendall tau-b: {_fmt_num(rp['kendall_tau_b'])}",
    ]
    for k_str, stats in rp["top_k"].items():
        lines.append(
            f"  Top-{k_str}: overlap {stats['overlap_count']} "
            f"(Jaccard {_fmt_num(stats['jaccard'])}; +{stats['gained_in_b_count']} gained, "
            f"-{stats['lost_from_a_count']} lost; "
            f"exact seq {stats['exact_sequence_matches']}, "
            f"shifted window {stats['shifted_window_matches']})"
        )

    lines.append("")
    lines.append(f"--- Gained in Top-{pk} ({rb['label']}) [{len(comparison['gained_in_top_k'])}] ---")
    if comparison["gained_in_top_k"]:
        for row in comparison["gained_in_top_k"]:
            lines.append(
                f"  #{row['rank_b']} (was {row['rank_a'] or 'absent'}) "
                f"{row['gene_name'] or row['source_identifier']} ({row['target_key']}): "
                f"score {_fmt_num(row['score_a'])} -> {_fmt_num(row['score_b'])} "
                f"[{', '.join(row['attribution'])}]"
            )
    else:
        lines.append("  (none)")

    lines.append("")
    lines.append(f"--- Lost from Top-{pk} ({ra['label']}) [{len(comparison['lost_from_top_k'])}] ---")
    if comparison["lost_from_top_k"]:
        for row in comparison["lost_from_top_k"]:
            lines.append(
                f"  was #{row['rank_a']} (now {row['rank_b'] or 'absent'}) "
                f"{row['gene_name'] or row['source_identifier']} ({row['target_key']}): "
                f"score {_fmt_num(row['score_a'])} -> {_fmt_num(row['score_b'])} "
                f"[{', '.join(row['attribution'])}]"
            )
    else:
        lines.append("  (none)")

    if comparison.get("shifted_window_targets"):
        lines.append("")
        lines.append(
            f"--- Shifted Peptide Windows [{len(comparison['shifted_window_targets'])}] ---"
        )
        for row in comparison["shifted_window_targets"]:
            lines.append(
                f"  #{row['rank_a']} -> #{row['rank_b']} "
                f"{row['gene_name'] or row['source_identifier']}: "
                f"{row['amino_acids_a']} -> {row['amino_acids_b']} "
                f"[{', '.join(row['attribution'])}]"
            )

    for mod in ("peptide", "mrna"):
        diff = comparison["construct_diffs"][mod]
        if not diff.get("available"):
            continue
        lines.append("")
        lines.append(f"--- Emitted {mod.upper()} Construct Diff ---")
        lines.append(
            f"Constructs: {diff['construct_count_a']} -> {diff['construct_count_b']} | "
            f"Antigen slots: {diff['slot_count_a']} -> {diff['slot_count_b']}"
        )
        lines.append(
            f"Exact sequence matches: {diff['exact_match_count']} | "
            f"Shifted windows: {diff['shifted_window_count']} | "
            f"Added: {diff['added_count']} | Dropped: {diff['dropped_count']}"
        )
        for item in diff.get("shifted_windows", ()):
            lines.append(
                f"  [shifted] {item['antigen_a']}: {item['sequence_a']} -> {item['sequence_b']}"
            )
        for item in diff.get("added_in_b", ()):
            lines.append(
                f"  [added]   {item['antigen_label']} ({item['construct_name']}): {item['sequence']}"
            )
        for item in diff.get("dropped_from_a", ()):
            lines.append(
                f"  [dropped] {item['antigen_label']} ({item['construct_name']}): {item['sequence']}"
            )
    return "\n".join(lines) + "\n"


def render_comparison_html(comparison: dict[str, Any]) -> str:
    """Render a self-contained HTML comparison report via Jinja2 template."""
    from jinja2 import Environment, PackageLoader

    env = Environment(
        loader=PackageLoader("vaxrank", "templates"),
        autoescape=False,
    )
    template = env.get_template("compare_report.html")
    return template.render(
        comparison=comparison,
        run_a=comparison["run_a"],
        run_b=comparison["run_b"],
        rp=comparison["rank_preservation"],
        primary_top_k=comparison["primary_top_k"],
        policies_differ=comparison.get("policies_differ", False),
        fmt_num=_fmt_num,
    )


def write_comparison_outputs(
    comparison: dict[str, Any],
    output_dir: str | Path,
    *,
    policy: ComparePolicy | dict[str, Any] | None = None,
) -> dict[str, str]:
    """Write structured JSON, CSV, ASCII, HTML, and optional PDF comparison outputs."""
    resolved_policy = ComparePolicy.resolve(
        policy if policy is not None else comparison.get("policy")
    )
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    written: dict[str, str] = {}

    json_path = out / "comparison_summary.json"
    with open(json_path, "w", encoding="utf-8") as handle:
        json.dump(comparison, handle, indent=2, allow_nan=False)
    written["json"] = str(json_path)

    csv_path = out / "comparable_targets.csv"
    fieldnames = [
        "target_key", "source_identifier", "gene_name", "source_kind",
        "status", "top_k_status", "rank_a", "rank_b", "rank_delta",
        "score_a", "score_b", "expression_a", "expression_b",
        "epitope_score_a", "epitope_score_b",
        "amino_acids_a", "amino_acids_b", "sequence_match", "attribution",
    ]
    with open(csv_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in comparison["targets"]:
            flat = dict(row)
            flat["attribution"] = ";".join(row["attribution"])
            writer.writerow(flat)
    written["csv"] = str(csv_path)

    txt_path = out / "comparison_report.txt"
    txt_path.write_text(render_comparison_ascii(comparison), encoding="utf-8")
    written["ascii"] = str(txt_path)

    html_content = render_comparison_html(comparison)
    html_path = out / "comparison_report.html"
    html_path.write_text(html_content, encoding="utf-8")
    written["html"] = str(html_path)

    if resolved_policy.include_pdf:
        pdf_path = out / "comparison_report.pdf"
        import weasyprint
        weasyprint.HTML(string=html_content).write_pdf(str(pdf_path))
        written["pdf"] = str(pdf_path)

    return written
