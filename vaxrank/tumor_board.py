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

"""AI-ready ``tumor_board_dossier.json``, deterministic safety-gate checks,
and audited ``--selection-overrides`` replay for human Immunogenomics Tumor
Boards and automated AI reviewers (openvax/vaxrank#922)."""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from numbers import Real
from pathlib import Path
from typing import Any

from .mhc_allele_status import normalize_mhc_allele_name
from .version import __version__


DOSSIER_SCHEMA = "vaxrank.tumor_board_dossier.v1"
OVERRIDES_SCHEMA = "vaxrank.tumor_board_overrides.v1"
VALID_DECISIONS = ("Accept", "Review", "Reject")
VALID_REVIEWER_KINDS = ("human", "ai")


class TumorBoardValidationError(ValueError):
    """Raised when ``tumor_board_overrides.json`` violates integrity,
    citation, or deterministic safety-gate rules."""


@dataclass(frozen=True)
class TumorBoardPolicy:
    """Explicit per-trial policy for deterministic hard safety gates, soft
    review flags, and override validation.

    Configurable under ``tumor_board:`` in a per-trial YAML policy.
    """
    forbid_override_exact_unsafe_self: bool = True
    forbid_override_all_mhc_lost: bool = True
    max_normal_alt_vaf: float = 0.02
    min_retained_mhc_score: float = 0.0
    require_dossier_sha256_match: bool = True
    require_citation_keys: bool = True
    require_rationale: bool = True
    max_gravy_review_threshold: float = 2.0
    min_rna_alt_reads_review_threshold: int = 3
    flag_dna_only_proxy_for_review: bool = True

    def __post_init__(self):
        for flag in (
            "forbid_override_exact_unsafe_self",
            "forbid_override_all_mhc_lost",
            "require_dossier_sha256_match",
            "require_citation_keys",
            "require_rationale",
            "flag_dna_only_proxy_for_review",
        ):
            if not isinstance(getattr(self, flag), bool):
                raise ValueError(f"tumor_board.{flag} must be a boolean")
        for name, lo, hi in (
            ("max_normal_alt_vaf", 0.0, 1.0),
            ("min_retained_mhc_score", 0.0, None),
            ("max_gravy_review_threshold", -5.0, 5.0),
        ):
            val = getattr(self, name)
            if isinstance(val, bool) or not isinstance(val, Real) or not math.isfinite(val):
                raise ValueError(f"tumor_board.{name} must be a finite number")
            if float(val) < lo or (hi is not None and float(val) > hi):
                raise ValueError(f"tumor_board.{name} out of bounds [{lo}, {hi}]")
        if (isinstance(self.min_rna_alt_reads_review_threshold, bool)
                or not isinstance(self.min_rna_alt_reads_review_threshold, int)
                or self.min_rna_alt_reads_review_threshold < 0):
            raise ValueError("tumor_board.min_rna_alt_reads_review_threshold must be an integer >= 0")

    @classmethod
    def resolve(cls, value=None) -> TumorBoardPolicy:
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if not isinstance(value, dict):
            raise TypeError(
                f"tumor_board policy must be a mapping or TumorBoardPolicy; got {type(value).__name__}"
            )
        defaults = cls().to_dict()
        unknown = set(value) - set(defaults)
        if unknown:
            raise ValueError("Unknown tumor_board policy fields: %s" % ", ".join(sorted(unknown)))
        merged = {k: (defaults[k] if value.get(k) is None else value[k]) for k in defaults}
        return cls(**merged)

    def to_dict(self) -> dict[str, Any]:
        return {
            "forbid_override_exact_unsafe_self": self.forbid_override_exact_unsafe_self,
            "forbid_override_all_mhc_lost": self.forbid_override_all_mhc_lost,
            "max_normal_alt_vaf": float(self.max_normal_alt_vaf),
            "min_retained_mhc_score": float(self.min_retained_mhc_score),
            "require_dossier_sha256_match": self.require_dossier_sha256_match,
            "require_citation_keys": self.require_citation_keys,
            "require_rationale": self.require_rationale,
            "max_gravy_review_threshold": float(self.max_gravy_review_threshold),
            "min_rna_alt_reads_review_threshold": self.min_rna_alt_reads_review_threshold,
            "flag_dna_only_proxy_for_review": self.flag_dna_only_proxy_for_review,
        }


@dataclass(frozen=True)
class TumorBoardDecision:
    """One audited human or AI decision for a candidate in the dossier."""
    candidate_id: str
    decision: str
    rationale: str = ""
    confidence: float | None = None
    cited_evidence_keys: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self):
        cid = str(self.candidate_id or "").strip()
        if not cid:
            raise TumorBoardValidationError("Override decision requires a non-empty candidate_id")
        object.__setattr__(self, "candidate_id", cid)
        if self.decision not in VALID_DECISIONS:
            raise TumorBoardValidationError(
                f"Decision for {cid!r} must be one of {VALID_DECISIONS}; got {self.decision!r}"
            )
        if self.confidence is not None:
            if (isinstance(self.confidence, bool)
                    or not isinstance(self.confidence, Real)
                    or not math.isfinite(self.confidence)
                    or not 0.0 <= float(self.confidence) <= 1.0):
                raise TumorBoardValidationError(
                    f"Confidence for {cid!r} must be in [0.0, 1.0]; got {self.confidence!r}"
                )
            object.__setattr__(self, "confidence", float(self.confidence))
        keys = tuple(str(k).strip() for k in (self.cited_evidence_keys or ()) if str(k).strip())
        object.__setattr__(self, "cited_evidence_keys", keys)
        object.__setattr__(self, "rationale", str(self.rationale or "").strip())

    def to_dict(self) -> dict[str, Any]:
        return {
            "candidate_id": self.candidate_id,
            "decision": self.decision,
            "confidence": self.confidence,
            "cited_evidence_keys": list(self.cited_evidence_keys),
            "rationale": self.rationale,
        }


@dataclass(frozen=True)
class TumorBoardOverrides:
    """Shared override bundle produced by either the interactive HTML review UI
    or an automated AI reviewer."""
    dossier_sha256: str
    reviewer_id: str
    reviewer_kind: str = "human"
    decisions: tuple[TumorBoardDecision, ...] = field(default_factory=tuple)
    policy_sha256: str = ""
    notes: str = ""
    schema: str = OVERRIDES_SCHEMA

    def __post_init__(self):
        if self.schema != OVERRIDES_SCHEMA:
            raise TumorBoardValidationError(
                f"Unsupported overrides schema {self.schema!r}; expected {OVERRIDES_SCHEMA!r}"
            )
        rid = str(self.reviewer_id or "").strip()
        if not rid:
            raise TumorBoardValidationError("Overrides bundle requires a non-empty reviewer_id")
        object.__setattr__(self, "reviewer_id", rid)
        rk = str(self.reviewer_kind or "human").strip().lower()
        if rk not in VALID_REVIEWER_KINDS:
            raise TumorBoardValidationError(
                f"reviewer_kind must be one of {VALID_REVIEWER_KINDS}; got {rk!r}"
            )
        object.__setattr__(self, "reviewer_kind", rk)
        object.__setattr__(self, "dossier_sha256", str(self.dossier_sha256 or "").strip())
        object.__setattr__(self, "decisions", tuple(self.decisions))

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "dossier_sha256": self.dossier_sha256,
            "reviewer_id": self.reviewer_id,
            "reviewer_kind": self.reviewer_kind,
            "policy_sha256": self.policy_sha256,
            "notes": self.notes,
            "decisions": [d.to_dict() for d in self.decisions],
        }


def _canonical_sha256(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _finite_or_none(val: Any) -> float | None:
    if val is None or isinstance(val, bool):
        return None
    if isinstance(val, Real):
        f = float(val)
        return round(f, 6) if math.isfinite(f) else None
    return None


def _flatten_evidence_keys(layers: dict[str, Any]) -> list[str]:
    keys: list[str] = []
    for layer_name, body in layers.items():
        if isinstance(body, dict):
            keys.append(layer_name)
            for k in body:
                keys.append(f"{layer_name}.{k}")
    return sorted(keys)


def _candidate_identifier(source: Any, vp: Any, rank: int) -> str:
    desc = getattr(source, "short_description", None)
    if desc:
        return str(desc)
    antigen = getattr(vp, "antigen", None)
    if antigen is not None and getattr(antigen, "display_identifier", None):
        return str(antigen.display_identifier)
    gene = getattr(antigen, "gene_name", None) if antigen is not None else None
    return f"C{rank:03d}:{gene or 'candidate'}"


def build_tumor_board_dossier(
    ranked_variants_with_vaccine_peptides,
    *,
    patient_id: str = "",
    mhc_alleles: list[str] | tuple[str, ...] | None = None,
    mhc_allele_status: dict[str, Any] | None = None,
    dna_vaf_by_variant: dict | None = None,
    products: dict[str, list] | None = None,
    policy: TumorBoardPolicy | dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Assemble the machine-readable, 6-layer ``tumor_board_dossier.json``
    bundle with deterministic hard safety gates and soft review flags.
    """
    pol = TumorBoardPolicy.resolve(policy)
    dna_vaf_by_variant = dna_vaf_by_variant or {}
    products = products or {}

    status_by_allele: dict[str, dict[str, Any]] = {}
    if isinstance(mhc_allele_status, dict):
        for row in mhc_allele_status.get("alleles") or ():
            if isinstance(row, dict) and row.get("allele"):
                status_by_allele[normalize_mhc_allele_name(row["allele"])] = row
    apm_disrupted = (
        list((mhc_allele_status or {}).get("apm", {}).get("disrupted_genes", ()))
        if isinstance((mhc_allele_status or {}).get("apm"), dict)
        else []
    )

    # Map emitted product names & peptide pools by amino-acid sequence
    products_by_aa: dict[str, list[str]] = {}
    pool_by_aa: dict[str, int] = {}
    for modality, constructs in products.items():
        for construct in constructs or ():
            if modality == "peptide":
                seq = getattr(construct, "sequence", "")
                if seq:
                    products_by_aa.setdefault(seq, []).append(construct.name)
                    pool = (getattr(construct, "components", None) or {}).get("pool")
                    if pool is not None:
                        pool_by_aa[seq] = pool
            elif modality == "mrna":
                for ag in getattr(construct, "antigens", ()) or ():
                    seq = ag.get("aa") if isinstance(ag, dict) else ""
                    if seq:
                        products_by_aa.setdefault(seq, []).append(construct.name)

    normalized_alleles = [
        normalize_mhc_allele_name(a) for a in (mhc_alleles or ()) if str(a or "").strip()
    ]

    candidates = []
    seen_ids: dict[str, int] = {}
    hard_failed_count = 0
    review_flagged_count = 0

    for rank_idx, (source, vaccine_peptides) in enumerate(
        ranked_variants_with_vaccine_peptides, start=1
    ):
        if not vaccine_peptides:
            continue
        vp = vaccine_peptides[0]
        base_id = _candidate_identifier(source, vp, rank_idx)
        seen_ids[base_id] = seen_ids.get(base_id, 0) + 1
        candidate_id = (
            base_id if seen_ids[base_id] == 1 else f"{base_id}#{seen_ids[base_id]}"
        )

        fragment = getattr(vp, "mutant_protein_fragment", None)
        antigen = getattr(vp, "antigen", None)
        gene_name = (
            getattr(antigen, "gene_name", None)
            or getattr(fragment, "gene_name", None)
            or (getattr(source, "gene_names", None) or [""])[0]
            or ""
        )
        source_kind = getattr(antigen, "kind", None) or ("variant" if fragment is not None else "antigen")
        aa = vp.amino_acids

        # Layer 1: Genomic / read-artifact telemetry
        vaf_info = None
        try:
            vaf_info = dna_vaf_by_variant.get(source)
        except TypeError:
            vaf_info = None
        if vaf_info is None:
            vaf_info = dna_vaf_by_variant.get(base_id) or {}
        if not isinstance(vaf_info, dict):
            vaf_info = {"tumor_vaf": _finite_or_none(vaf_info)}
        normal_alt_vaf = _finite_or_none(vaf_info.get("normal_alt_vaf") or vaf_info.get("normal_vaf"))
        genomic_telemetry = {
            "chr": getattr(source, "contig", None),
            "pos": getattr(source, "original_start", None),
            "ref": getattr(source, "original_ref", None),
            "alt": getattr(source, "original_alt", None),
            "tumor_vaf": _finite_or_none(vaf_info.get("tumor_vaf") or vaf_info.get("vaf")),
            "tumor_depth": vaf_info.get("tumor_depth"),
            "normal_alt_vaf": normal_alt_vaf,
            "normal_depth": vaf_info.get("normal_depth"),
            "strand_bias_pvalue": _finite_or_none(vaf_info.get("strand_bias_pvalue")),
            "read_end_bias": _finite_or_none(vaf_info.get("read_end_bias")),
            "soft_clip_fraction": _finite_or_none(vaf_info.get("soft_clip_fraction")),
            "homopolymer_context": vaf_info.get("homopolymer_context"),
        }

        # Layer 2: Transcriptomic / isoform integrity
        n_rna_alt = getattr(fragment, "n_alt_reads", None) if fragment is not None else None
        n_rna_ref = getattr(fragment, "n_ref_reads", None) if fragment is not None else None
        n_rna_sup = (
            getattr(fragment, "n_alt_reads_supporting_protein_sequence", None)
            if fragment is not None else None
        )
        proxy_tpm = getattr(fragment, "proxy_expression_tpm", None) if fragment is not None else None
        proxy_cancer = getattr(fragment, "proxy_expression_cancer_code", "") if fragment is not None else ""
        is_dna_proxy = bool(proxy_tpm is not None or (fragment is not None and (n_rna_alt or 0) == 0 and proxy_cancer))
        transcriptomic_integrity = {
            "n_rna_alt_reads": n_rna_alt,
            "n_rna_ref_reads": n_rna_ref,
            "n_rna_alt_supporting_sequence": n_rna_sup,
            "transcript_ids": list(getattr(antigen, "transcript_ids", ()) or ()),
            "is_dna_only_proxy": is_dna_proxy,
            "proxy_expression_tpm": _finite_or_none(proxy_tpm),
            "proxy_expression_cancer_code": proxy_cancer or None,
        }

        # Layer 3: Tumor evolution & MHC presentation competence
        binding_alleles_set: set[str] = set()
        raw_positive_alleles: set[str] = set()
        retained_positive_alleles: set[str] = set()
        for ep in vp.target_epitopes:
            for allele, sc in (ep.per_allele_scores or {}).items():
                norm_a = normalize_mhc_allele_name(allele)
                binding_alleles_set.add(norm_a)
                if sc and float(sc) > 0:
                    retained_positive_alleles.add(norm_a)
            for pred in ep.predictions_flat():
                if pred.allele:
                    norm_a = normalize_mhc_allele_name(pred.allele)
                    binding_alleles_set.add(norm_a)
                    if (pred.percentile_rank is not None and pred.percentile_rank <= 2.0) or (
                        pred.score is not None and pred.score > 0
                    ):
                        raw_positive_alleles.add(norm_a)

        allele_statuses_here = {
            a: {
                "effective_status": status_by_allele[a]["effective_status"],
                "badge": status_by_allele[a]["badge"],
                "masked": status_by_allele[a]["masked"],
                "caveats": status_by_allele[a]["caveats"],
            }
            for a in sorted(binding_alleles_set)
            if a in status_by_allele
        }
        all_mhc_lost_or_silenced = bool(
            raw_positive_alleles
            and not retained_positive_alleles
            and all(
                status_by_allele.get(a, {}).get("masked", False)
                for a in raw_positive_alleles
            )
        )
        tumor_evolution_and_mhc = {
            "source_species": getattr(antigen, "species", None),
            "mhc_alleles_binding": sorted(binding_alleles_set),
            "mhc_allele_status": allele_statuses_here,
            "all_restricting_mhc_lost_or_silenced": all_mhc_lost_or_silenced,
            "apm_disrupted_genes": apm_disrupted,
        }

        # Layer 4: Immunogenicity & epitope geometry
        top_epitopes = []
        best_rank = None
        best_allele = None
        ms_supported = False
        for ep in vp.target_epitopes[:5]:
            if getattr(ep, "ms_increment", 0.0) or getattr(ep, "ms_scoring", None):
                ms_supported = True
            ep_best_rank = None
            ep_best_allele = None
            for pred in ep.predictions_flat():
                pr = _finite_or_none(pred.percentile_rank)
                if pr is not None and (ep_best_rank is None or pr < ep_best_rank):
                    ep_best_rank = pr
                    ep_best_allele = pred.allele
            if ep_best_rank is not None and (best_rank is None or ep_best_rank < best_rank):
                best_rank = ep_best_rank
                best_allele = ep_best_allele
            top_epitopes.append({
                "sequence": ep.sequence,
                "offset": ep.offset,
                "overlaps_mutation": bool(ep.overlaps_mutation),
                "best_percentile_rank": ep_best_rank,
                "best_allele": ep_best_allele,
                "score": _finite_or_none(sum((ep.per_allele_scores or {}).values())),
            })
        immunogenicity_and_geometry = {
            "target_epitope_count": len(vp.target_epitopes),
            "best_percentile_rank": best_rank,
            "best_allele": best_allele,
            "ms_supported": ms_supported,
            "top_epitopes": top_epitopes,
        }

        # Layer 5: Safety & autoimmunity
        exact_self_count = len(vp.self_epitopes)
        targetable_aa = (
            sum(iv.end - iv.start for iv in antigen.targetable_mask.intervals)
            if antigen is not None and getattr(antigen, "targetable_mask", None) is not None
            else len(aa)
        )
        normal_obs = getattr(antigen, "normal_observations", ()) if antigen is not None else ()
        safety_and_autoimmunity = {
            "exact_self_overlap_count": exact_self_count,
            "targetable_aa_count": targetable_aa,
            "self_epitope_score": _finite_or_none(vp.self_epitope_score),
            "normal_observations_count": len(normal_obs or ()),
            "only_exact_self_epitopes": bool(exact_self_count > 0 and len(vp.target_epitopes) == 0),
        }

        # Layer 6: Construct & manufacturability
        mfg = vp.manufacturability_scores
        max_gravy = _finite_or_none(getattr(mfg, "max_7mer_gravy_score", None))
        construct_and_manufacturability = {
            "length_aa": len(aa),
            "max_7mer_gravy": max_gravy,
            "cterm_7mer_gravy": _finite_or_none(getattr(mfg, "cterm_7mer_gravy_score", None)),
            "cysteine_count": getattr(mfg, "cysteine_count", None),
            "difficult_n_terminal": bool(getattr(mfg, "difficult_n_terminal_residue", False)),
            "aspartate_proline_bonds": getattr(mfg, "aspartate_proline_bond_count", None),
            "emitted_products": products_by_aa.get(aa, []),
            "peptide_pool": pool_by_aa.get(aa),
        }

        layers = {
            "genomic_telemetry": genomic_telemetry,
            "transcriptomic_integrity": transcriptomic_integrity,
            "tumor_evolution_and_mhc": tumor_evolution_and_mhc,
            "immunogenicity_and_geometry": immunogenicity_and_geometry,
            "safety_and_autoimmunity": safety_and_autoimmunity,
            "construct_and_manufacturability": construct_and_manufacturability,
        }

        # Deterministic hard safety/technical gates
        hard_gates: list[dict[str, str]] = []
        if pol.forbid_override_exact_unsafe_self and safety_and_autoimmunity["only_exact_self_epitopes"]:
            hard_gates.append({
                "gate": "exact_unsafe_self_only",
                "detail": "Candidate has exact self-proteome epitopes and zero targetable epitopes.",
                "evidence_key": "safety_and_autoimmunity.only_exact_self_epitopes",
            })
        if normal_alt_vaf is not None and normal_alt_vaf > pol.max_normal_alt_vaf:
            hard_gates.append({
                "gate": "normal_alt_leakage",
                "detail": (
                    f"Normal alt VAF {normal_alt_vaf:.4f} exceeds policy maximum "
                    f"{pol.max_normal_alt_vaf:.4f}."
                ),
                "evidence_key": "genomic_telemetry.normal_alt_vaf",
            })
        if pol.forbid_override_all_mhc_lost and all_mhc_lost_or_silenced:
            hard_gates.append({
                "gate": "all_restricting_mhc_lost_or_silenced",
                "detail": "All binding MHC alleles for this candidate are confirmed lost or silenced.",
                "evidence_key": "tumor_evolution_and_mhc.all_restricting_mhc_lost_or_silenced",
            })
        if (pol.min_retained_mhc_score > 0.0
                and (vp.target_epitope_score or 0.0) < pol.min_retained_mhc_score):
            hard_gates.append({
                "gate": "below_min_retained_mhc_score",
                "detail": (
                    f"Retained MHC epitope score {vp.target_epitope_score:.3f} < "
                    f"policy minimum {pol.min_retained_mhc_score:.3f}."
                ),
                "evidence_key": "immunogenicity_and_geometry.target_epitope_count",
            })

        # Nuanced soft review flags (trade-offs for human ITB or AI review)
        review_flags: list[dict[str, str]] = []
        if pol.flag_dna_only_proxy_for_review and is_dna_proxy:
            review_flags.append({
                "flag": "dna_only_proxy_expression",
                "detail": "Expression relies on cohort proxy (no direct tumor RNA reads).",
                "evidence_key": "transcriptomic_integrity.is_dna_only_proxy",
            })
        if (not is_dna_proxy
                and n_rna_alt is not None
                and n_rna_alt < pol.min_rna_alt_reads_review_threshold):
            review_flags.append({
                "flag": "low_rna_alt_reads",
                "detail": (
                    f"Tumor RNA alt reads ({n_rna_alt}) below review threshold "
                    f"({pol.min_rna_alt_reads_review_threshold})."
                ),
                "evidence_key": "transcriptomic_integrity.n_rna_alt_reads",
            })
        if max_gravy is not None and max_gravy > pol.max_gravy_review_threshold:
            review_flags.append({
                "flag": "high_hydrophobicity",
                "detail": (
                    f"Max 7-mer GRAVY ({max_gravy:.2f}) exceeds review threshold "
                    f"({pol.max_gravy_review_threshold:.2f})."
                ),
                "evidence_key": "construct_and_manufacturability.max_7mer_gravy",
            })
        if apm_disrupted:
            review_flags.append({
                "flag": "apm_pathway_disruption",
                "detail": f"APM disruption reported in {', '.join(apm_disrupted)}.",
                "evidence_key": "tumor_evolution_and_mhc.apm_disrupted_genes",
            })
        if any(
            st.get("effective_status") == "unresolvable"
            for st in allele_statuses_here.values()
        ):
            review_flags.append({
                "flag": "mhc_allele_unresolvable",
                "detail": "One or more restricting MHC alleles have unresolvable expression status.",
                "evidence_key": "tumor_evolution_and_mhc.mhc_allele_status",
            })

        if hard_gates:
            hard_failed_count += 1
            default_rec = "Reject"
        elif review_flags:
            review_flagged_count += 1
            default_rec = "Review"
        else:
            default_rec = "Accept"

        candidates.append({
            "candidate_id": candidate_id,
            "rank": rank_idx,
            "gene_name": gene_name,
            "source_kind": source_kind,
            "source_identifier": base_id,
            "amino_acids": aa,
            "combined_score": _finite_or_none(vp.combined_score),
            "target_epitope_score": _finite_or_none(vp.target_epitope_score),
            "expression_score": _finite_or_none(getattr(vp, "expression_score", None)),
            "default_recommendation": default_rec,
            "hard_gate_failures": hard_gates,
            "review_flags": review_flags,
            "evidence_keys": _flatten_evidence_keys(layers),
            **layers,
        })

    body = {
        "schema": DOSSIER_SCHEMA,
        "vaxrank_version": __version__,
        "patient_id": str(patient_id or ""),
        "mhc_alleles": normalized_alleles,
        "mhc_allele_status": mhc_allele_status,
        "policy": pol.to_dict(),
        "candidate_count": len(candidates),
        "hard_gate_failed_count": hard_failed_count,
        "review_flagged_count": review_flagged_count,
        "candidates": candidates,
    }
    dossier_sha256 = _canonical_sha256(body)
    return {"dossier_sha256": dossier_sha256, **body}


def write_tumor_board_dossier(dossier: dict[str, Any], path: str | Path) -> str:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as handle:
        json.dump(dossier, handle, indent=2, allow_nan=False)
    return str(p)


def load_tumor_board_overrides(path: str | Path) -> TumorBoardOverrides:
    """Load and parse a ``tumor_board_overrides.json`` bundle (human or AI)."""
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(f"Selection overrides file not found: {p}")
    with open(p, encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise TumorBoardValidationError(f"Expected JSON object in {p}")
    raw_decisions = data.get("decisions")
    if not isinstance(raw_decisions, list):
        raise TumorBoardValidationError(f"'decisions' must be a list in {p}")
    decisions = tuple(
        TumorBoardDecision(
            candidate_id=item.get("candidate_id", ""),
            decision=item.get("decision", ""),
            rationale=item.get("rationale", ""),
            confidence=item.get("confidence"),
            cited_evidence_keys=tuple(item.get("cited_evidence_keys") or ()),
        )
        for item in raw_decisions
    )
    return TumorBoardOverrides(
        schema=data.get("schema", OVERRIDES_SCHEMA),
        dossier_sha256=data.get("dossier_sha256", ""),
        reviewer_id=data.get("reviewer_id", ""),
        reviewer_kind=data.get("reviewer_kind", "human"),
        policy_sha256=data.get("policy_sha256", ""),
        notes=data.get("notes", ""),
        decisions=decisions,
    )


def validate_and_apply_overrides(
    ranked_variants_with_vaccine_peptides,
    dossier: dict[str, Any],
    overrides: TumorBoardOverrides,
    *,
    policy: TumorBoardPolicy | dict[str, Any] | None = None,
) -> tuple[list, dict[str, Any]]:
    """Deterministically validate ``overrides`` against ``dossier`` and
    ``policy``, blocking hard safety-gate bypasses, then apply decisions to
    ``ranked_variants_with_vaccine_peptides``.

    Returns ``(filtered_ranked, override_audit)``.
    """
    pol = TumorBoardPolicy.resolve(policy if policy is not None else dossier.get("policy"))

    expected_sha = dossier.get("dossier_sha256", "")
    if pol.require_dossier_sha256_match:
        if not overrides.dossier_sha256 or overrides.dossier_sha256 != expected_sha:
            raise TumorBoardValidationError(
                f"Override dossier_sha256 mismatch: got {overrides.dossier_sha256!r}, "
                f"expected {expected_sha!r}. Re-export overrides from the current dossier."
            )

    by_id = {c["candidate_id"]: c for c in dossier.get("candidates", ())}
    decisions_by_id: dict[str, TumorBoardDecision] = {}

    for dec in overrides.decisions:
        if dec.candidate_id in decisions_by_id:
            raise TumorBoardValidationError(
                f"Duplicate override decision for candidate_id {dec.candidate_id!r}"
            )
        cand = by_id.get(dec.candidate_id)
        if cand is None:
            raise TumorBoardValidationError(
                f"Unknown candidate_id {dec.candidate_id!r} in overrides; "
                f"must match a candidate in the dossier."
            )

        # Hard safety gates cannot be bypassed by Accept
        hard_gates = cand.get("hard_gate_failures") or []
        if hard_gates and dec.decision == "Accept":
            gate_names = ", ".join(g["gate"] for g in hard_gates)
            raise TumorBoardValidationError(
                f"Cannot override candidate {dec.candidate_id!r} to 'Accept': "
                f"failed deterministic hard gate(s): {gate_names}."
            )

        # Validate cited evidence keys
        valid_keys = set(cand.get("evidence_keys") or ())
        invalid_keys = [k for k in dec.cited_evidence_keys if k not in valid_keys]
        if invalid_keys:
            raise TumorBoardValidationError(
                f"Candidate {dec.candidate_id!r} cites unknown evidence key(s): "
                f"{', '.join(invalid_keys)}"
            )
        if pol.require_citation_keys and not dec.cited_evidence_keys:
            raise TumorBoardValidationError(
                f"Candidate {dec.candidate_id!r} decision {dec.decision!r} must cite "
                f"at least one valid dossier evidence key."
            )
        if pol.require_rationale and (
            dec.decision == "Reject" or dec.decision != cand.get("default_recommendation")
        ) and not dec.rationale:
            raise TumorBoardValidationError(
                f"Candidate {dec.candidate_id!r} decision {dec.decision!r} requires a non-empty rationale."
            )
        decisions_by_id[dec.candidate_id] = dec

    # Map ranked items to dossier candidate_ids in the same deterministic order
    seen_ids: dict[str, int] = {}
    filtered_ranked = []
    rejected_ids: list[str] = []
    accepted_ids: list[str] = []
    review_ids: list[str] = []

    for rank_idx, (source, vaccine_peptides) in enumerate(
        ranked_variants_with_vaccine_peptides, start=1
    ):
        if not vaccine_peptides:
            continue
        vp = vaccine_peptides[0]
        base_id = _candidate_identifier(source, vp, rank_idx)
        seen_ids[base_id] = seen_ids.get(base_id, 0) + 1
        cid = base_id if seen_ids[base_id] == 1 else f"{base_id}#{seen_ids[base_id]}"
        dec = decisions_by_id.get(cid)
        if dec is not None and dec.decision == "Reject":
            rejected_ids.append(cid)
            continue
        if dec is not None and dec.decision == "Accept":
            accepted_ids.append(cid)
        elif dec is not None and dec.decision == "Review":
            review_ids.append(cid)
        filtered_ranked.append((source, vaccine_peptides))

    audit = {
        "schema": "vaxrank.tumor_board_override_audit.v1",
        "dossier_sha256": expected_sha,
        "reviewer_id": overrides.reviewer_id,
        "reviewer_kind": overrides.reviewer_kind,
        "notes": overrides.notes,
        "policy": pol.to_dict(),
        "total_decisions": len(decisions_by_id),
        "accepted_ids": accepted_ids,
        "review_ids": review_ids,
        "rejected_ids": rejected_ids,
        "decisions": [d.to_dict() for d in overrides.decisions],
    }
    return filtered_ranked, audit
