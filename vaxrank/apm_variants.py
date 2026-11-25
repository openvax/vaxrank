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

"""Cross-species DNA somatic MHC and antigen-presentation machinery (APM)
variant and loss-of-function detector (#935).

Scans somatic variants across species (human, canine, murine, swine, NHP, etc.)
in:
- Core Class I APM genes (B2M, TAP1, TAP2, TAPBP, NLRC5, CALR, CANX, PDIA3,
  ERAP1, ERAP2, PSMB8, PSMB9, PSMB10)
- IFN-gamma signaling / MHC transcriptional regulation genes (JAK1, JAK2,
  IFNGR1, IFNGR2, STAT1, IRF1)
- Classical MHC Class I genes (HLA-A/B/C, DLA-88/12/64/79, H2-K1/D1/L,
  SLA-1/2/3, plus any classical MHC gene recognized by ``mhcgnomes``)

Classifies each variant's DNA functional impact cleanly into:
- ``truncating_lof`` (frameshift, premature stop, start-loss, exon-loss,
  canonical splice-donor/acceptor)
- ``missense_or_inframe`` (substitution, in-frame indel, stop-loss; flagged for
  review, never mislabeled as confirmed LoF)
- ``silent_or_noncoding`` (synonymous/intronic/UTR/intergenic; never flags APM
  disruption)

Strictly separates DNA somatic mutations from RNA expression (optionally
attaching Isovar RNA alt-read confirmation when a tumor RNA BAM was run) and
feeds structured ``APMGeneStatus`` / ``MHCAlleleStatus`` records into reports
and ``tumor_board_dossier.json``.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import logging
import re
from typing import Any, Iterable

import mhcgnomes

from .mhc_allele_status import (
    APMGeneStatus,
    APMStatus,
    MHCAlleleStatus,
    MHCAlleleStatusBundle,
    normalize_mhc_allele_name,
)
from .varcode_effects import select_varcode_effect_outcome

logger = logging.getLogger(__name__)

# Canonical symbol -> (pathway_role, aliases)
_APM_AND_IFNG_GENE_TABLE: dict[str, tuple[str, tuple[str, ...]]] = {
    # Core Class I antigen presentation machinery
    "B2M": ("apm", ("B2M",)),
    "TAP1": ("apm", ("TAP1", "ABCB2", "RING4")),
    "TAP2": ("apm", ("TAP2", "ABCB3", "RING11")),
    "TAPBP": ("apm", ("TAPBP", "TAPASIN", "TAPA")),
    "NLRC5": ("apm", ("NLRC5", "CITA")),
    "CALR": ("apm", ("CALR", "CRT")),
    "CANX": ("apm", ("CANX", "CNX")),
    "PDIA3": ("apm", ("PDIA3", "ERP57", "GRP58")),
    "ERAP1": ("apm", ("ERAP1", "ARTS1", "PILSAP")),
    "ERAP2": ("apm", ("ERAP2", "LRAP")),
    "PSMB8": ("apm", ("PSMB8", "LMP7", "RING10")),
    "PSMB9": ("apm", ("PSMB9", "LMP2", "RING12")),
    "PSMB10": ("apm", ("PSMB10", "MECL1")),
    # IFN-gamma signaling / MHC Class I transcriptional axis
    "JAK1": ("ifng", ("JAK1",)),
    "JAK2": ("ifng", ("JAK2",)),
    "IFNGR1": ("ifng", ("IFNGR1",)),
    "IFNGR2": ("ifng", ("IFNGR2",)),
    "STAT1": ("ifng", ("STAT1",)),
    "IRF1": ("ifng", ("IRF1",)),
}

_APM_ALIAS_TO_CANONICAL: dict[str, tuple[str, str]] = {}
for _canon, (_role, _aliases) in _APM_AND_IFNG_GENE_TABLE.items():
    for _alias in _aliases:
        _APM_ALIAS_TO_CANONICAL[_alias.upper()] = (_canon, _role)

# Known cross-species classical MHC Class I gene symbols (normalized uppercase).
# Map variant genomic symbol -> canonical mhcgnomes (species_prefix, gene_name).
_CLASSICAL_MHC_GENE_ALIASES: dict[str, tuple[str, str]] = {
    # Human
    "HLA-A": ("HLA", "A"),
    "HLA-B": ("HLA", "B"),
    "HLA-C": ("HLA", "C"),
    # Canine (Canis lupus familiaris)
    "DLA-88": ("DLA", "88"),
    "DLA88": ("DLA", "88"),
    "DLA-12": ("DLA", "12"),
    "DLA12": ("DLA", "12"),
    "DLA-64": ("DLA", "64"),
    "DLA64": ("DLA", "64"),
    "DLA-79": ("DLA", "79"),
    "DLA79": ("DLA", "79"),
    # Murine (Mus musculus) — Ensembl names H2-K1, H2-D1, H2-L
    "H2-K1": ("H2", "K"),
    "H2-K": ("H2", "K"),
    "H-2-K": ("H2", "K"),
    "H-2-K1": ("H2", "K"),
    "H2-D1": ("H2", "D"),
    "H2-D": ("H2", "D"),
    "H-2-D": ("H2", "D"),
    "H-2-D1": ("H2", "D"),
    "H2-L": ("H2", "L"),
    "H-2-L": ("H2", "L"),
    # Swine (Sus scrofa)
    "SLA-1": ("SLA", "1"),
    "SLA-2": ("SLA", "2"),
    "SLA-3": ("SLA", "3"),
    # Non-human primates (rhesus / cynomolgus macaque)
    "MAMU-A": ("Mamu", "A"),
    "MAMU-B": ("Mamu", "B"),
    "MAFA-A": ("Mafa", "A"),
    "MAFA-B": ("Mafa", "B"),
}

_TRUNCATING_EFFECT_NAMES = frozenset({
    "frameshift",
    "frameshifttruncation",
    "prematurestop",
    "stopgained",
    "nonsense",
    "startloss",
    "startlost",
    "exonloss",
    "exonicdeletion",
    "splicedonor",
    "spliceacceptor",
    "exonicsplicesite",
    "truncating",
    "truncatinglof",
    "lof",
    "homdel",
})
_MISSENSE_EFFECT_NAMES = frozenset({
    "substitution",
    "complexsubstitution",
    "missense",
    "missensevariant",
    "insertion",
    "deletion",
    "inframeinsertion",
    "inframedeletion",
    "inframe",
    "stoploss",
    "stoplost",
    "alternatestartcodon",
    "incompletetranscript",
})
_SILENT_EFFECT_NAMES = frozenset({
    "silent",
    "synonymous",
    "synonymousvariant",
    "intronic",
    "intronvariant",
    "intergenic",
    "fiveprimeutr",
    "threeprimeutr",
    "noncodingtranscript",
    "intronicsplicesite",
    "crypticexon",
    "downstream",
    "upstream",
})

_LEGACY_APM_ATTR = {
    "B2M": "b2m_status",
    "TAP1": "tap1_status",
    "TAP2": "tap2_status",
    "NLRC5": "nlrc5_status",
    "JAK1": "jak1_status",
    "JAK2": "jak2_status",
}


@dataclass(frozen=True)
class APMVariantPolicy:
    """Explicit per-trial YAML policy for DNA somatic MHC & APM variant detection."""
    enabled: bool = True
    flag_missense_as_review: bool = True
    auto_mask_mhc_gene_lof_when_homozygous_locus: bool = True
    min_rna_alt_reads_confirmed: int = 1
    extra_apm_genes: tuple[str, ...] = field(default_factory=tuple)
    extra_mhc_genes: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self):
        for flag in (
            "enabled",
            "flag_missense_as_review",
            "auto_mask_mhc_gene_lof_when_homozygous_locus",
        ):
            if not isinstance(getattr(self, flag), bool):
                raise ValueError(f"apm_variants.{flag} must be a boolean")
        if (
            isinstance(self.min_rna_alt_reads_confirmed, bool)
            or not isinstance(self.min_rna_alt_reads_confirmed, int)
            or self.min_rna_alt_reads_confirmed < 0
        ):
            raise ValueError("apm_variants.min_rna_alt_reads_confirmed must be an integer >= 0")
        object.__setattr__(
            self, "extra_apm_genes",
            tuple(str(g).strip().upper() for g in (self.extra_apm_genes or ()) if str(g).strip()),
        )
        object.__setattr__(
            self, "extra_mhc_genes",
            tuple(str(g).strip().upper() for g in (self.extra_mhc_genes or ()) if str(g).strip()),
        )

    @classmethod
    def resolve(cls, value: APMVariantPolicy | dict[str, Any] | None = None) -> APMVariantPolicy:
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if not isinstance(value, dict):
            raise ValueError(f"apm_variants must be a dict or APMVariantPolicy, got {type(value).__name__}")
        defaults = cls().to_dict()
        unknown = set(value) - set(defaults)
        if unknown:
            raise ValueError(
                f"Unknown apm_variants keys: {', '.join(sorted(unknown))}. "
                f"Valid keys: {', '.join(sorted(defaults))}"
            )
        merged = {k: (defaults[k] if value.get(k) is None else value[k]) for k in defaults}
        return cls(**merged)

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["extra_apm_genes"] = list(self.extra_apm_genes)
        d["extra_mhc_genes"] = list(self.extra_mhc_genes)
        return d


def classify_gene_role(
    gene_name: str,
    policy: APMVariantPolicy | None = None,
) -> tuple[str | None, str | None, tuple[str, str] | None]:
    """Return ``(canonical_gene, pathway_role, mhc_locus_pair)`` if ``gene_name``
    belongs to the cross-species APM / IFN-gamma / classical MHC panel.

    Returns ``(None, None, None)`` for unrelated genes.
    """
    raw = str(gene_name or "").strip()
    if not raw:
        return None, None, None
    upper = raw.upper()
    pol = APMVariantPolicy.resolve(policy)

    if upper in _APM_ALIAS_TO_CANONICAL:
        canon, role = _APM_ALIAS_TO_CANONICAL[upper]
        return canon, role, None
    if upper in pol.extra_apm_genes:
        return upper, "apm", None
    if upper in _CLASSICAL_MHC_GENE_ALIASES:
        sp_prefix, gene_part = _CLASSICAL_MHC_GENE_ALIASES[upper]
        return f"{sp_prefix}-{gene_part}", "mhc_class_i", (sp_prefix, gene_part)
    if upper in pol.extra_mhc_genes:
        return upper, "mhc_class_i", None

    # Cross-species mhcgnomes fallback for any classical MHC gene symbol
    # (e.g. BoLA-1, Eqca-1, Patr-A, RT1-A, etc.)
    if "-" in raw or upper.startswith(("HLA", "DLA", "SLA", "H2", "MAMU", "MAFA", "BOLA", "EQCA", "PATR")):
        try:
            parsed = mhcgnomes.parse(raw, raise_on_error=False)
        except Exception:
            parsed = None
        if isinstance(parsed, mhcgnomes.Gene) and getattr(parsed, "is_class1", False):
            sp_prefix = getattr(parsed, "species_prefix", "") or ""
            g_name = getattr(parsed, "name", "") or ""
            if sp_prefix and g_name:
                return f"{sp_prefix}-{g_name}", "mhc_class_i", (sp_prefix, g_name)
    return None, None, None


def classify_variant_effect_impact(effect: Any) -> tuple[str, str, str]:
    """Classify a ``varcode`` effect (or stub/dict) into
    ``(impact_class, effect_type, effect_description)``.

    ``impact_class`` is strictly one of:
    - ``"truncating_lof"``
    - ``"missense_or_inframe"``
    - ``"silent_or_noncoding"``
    """
    if effect is None:
        return "silent_or_noncoding", "None", ""

    if isinstance(effect, dict):
        eff_type = str(effect.get("effect_type") or effect.get("consequence") or effect.get("effect") or "")
        eff_desc = str(effect.get("short_description") or effect.get("hgvs_p") or eff_type)
    else:
        eff_type = type(effect).__name__
        eff_desc = str(getattr(effect, "short_description", "") or eff_type)

    token = re.sub(r"[^a-z0-9]+", "", eff_type.lower())
    desc_token = re.sub(r"[^a-z0-9]+", "", eff_desc.lower())

    if token in _SILENT_EFFECT_NAMES:
        return "silent_or_noncoding", eff_type, eff_desc
    if token in _TRUNCATING_EFFECT_NAMES:
        return "truncating_lof", eff_type, eff_desc
    if token in _MISSENSE_EFFECT_NAMES:
        return "missense_or_inframe", eff_type, eff_desc

    # Fallback inspection on description / type substrings (negation-safe)
    if any(k in token or k in desc_token for k in ("silent", "synonymous", "intron", "utr", "intergenic", "noncoding")):
        return "silent_or_noncoding", eff_type, eff_desc
    if any(k in token or k in desc_token for k in ("frameshift", "prematurestop", "stopgain", "nonsense", "startloss", "exonloss", "splicedonor", "spliceacceptor")):
        return "truncating_lof", eff_type, eff_desc
    if re.search(r"p\.[A-Za-z]+\d+\*$", eff_desc) or re.search(r"p\.[A-Za-z]+\d+fs", eff_desc, re.I):
        return "truncating_lof", eff_type, eff_desc
    if any(k in token or k in desc_token for k in ("missense", "substitution", "inframe", "insertion", "deletion", "stoploss")):
        return "missense_or_inframe", eff_type, eff_desc
    return "silent_or_noncoding", eff_type, eff_desc


def _match_patient_alleles_for_mhc_locus(
    mhc_locus: tuple[str, str] | None,
    patient_mhc_alleles: Iterable[str] | None,
) -> tuple[str, ...]:
    """Return canonical patient MHC alleles belonging to ``(species_prefix, gene_name)``."""
    if not mhc_locus or not patient_mhc_alleles:
        return ()
    target_sp, target_gene = mhc_locus[0].upper(), mhc_locus[1].upper()
    matched: list[str] = []
    for raw in patient_mhc_alleles:
        if not str(raw or "").strip():
            continue
        canon = normalize_mhc_allele_name(raw)
        try:
            parsed = mhcgnomes.parse(canon, raise_on_error=False)
        except Exception:
            parsed = None
        if parsed is None:
            continue
        sp = (getattr(parsed, "species_prefix", "") or "").upper()
        gene = (getattr(parsed, "gene_name", "") or "").upper()
        if sp == target_sp and gene == target_gene and canon not in matched:
            matched.append(canon)
    return tuple(matched)


@dataclass(frozen=True)
class SomaticAPMorMHCVariant:
    """One somatic DNA variant detected in an APM, IFN-gamma, or classical MHC gene."""
    gene: str
    raw_gene: str
    pathway_role: str
    variant_id: str
    impact_class: str
    effect_type: str
    effect_description: str
    dna_vaf: float | None = None
    rna_alt_reads: int | None = None
    rna_ref_reads: int | None = None
    rna_confirmed: bool | None = None
    matched_patient_alleles: tuple[str, ...] = field(default_factory=tuple)
    phasing_status: str = "not_applicable"
    species: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "gene": self.gene,
            "raw_gene": self.raw_gene,
            "pathway_role": self.pathway_role,
            "variant_id": self.variant_id,
            "impact_class": self.impact_class,
            "consequence": self.impact_class,
            "is_lof": self.impact_class == "truncating_lof",
            "effect_type": self.effect_type,
            "effect_description": self.effect_description,
            "dna_vaf": self.dna_vaf,
            "rna_alt_reads": self.rna_alt_reads,
            "rna_ref_reads": self.rna_ref_reads,
            "rna_confirmed": self.rna_confirmed,
            "matched_patient_alleles": list(self.matched_patient_alleles),
            "phasing_status": self.phasing_status,
            "species": self.species,
        }


@dataclass(frozen=True)
class APMVariantDetectionResult:
    """Summary of cross-species somatic DNA variants detected in MHC & APM genes."""
    variants: tuple[SomaticAPMorMHCVariant, ...] = field(default_factory=tuple)
    silent_variants: tuple[SomaticAPMorMHCVariant, ...] = field(default_factory=tuple)
    policy: APMVariantPolicy = field(default_factory=APMVariantPolicy)

    @property
    def has_findings(self) -> bool:
        return bool(self.variants)

    @property
    def has_disruptive_findings(self) -> bool:
        return any(v.impact_class == "truncating_lof" for v in self.variants)

    def to_mhc_allele_status_bundle(
        self,
        *,
        source_path: str = "auto:dna_somatic_mhc_apm_detector",
    ) -> MHCAlleleStatusBundle:
        """Convert detected DNA somatic MHC & APM variants into an
        ``MHCAlleleStatusBundle`` strictly keeping DNA status separate from RNA
        expression status."""
        by_gene: dict[str, list[SomaticAPMorMHCVariant]] = {}
        allele_variants: dict[str, list[dict[str, Any]]] = {}

        for v in self.variants:
            by_gene.setdefault(v.gene, []).append(v)
            if v.pathway_role == "mhc_class_i" and v.matched_patient_alleles:
                if (
                    len(v.matched_patient_alleles) == 1
                    and self.policy.auto_mask_mhc_gene_lof_when_homozygous_locus
                ):
                    only_allele = v.matched_patient_alleles[0]
                    allele_variants.setdefault(only_allele, []).append(v.to_dict())

        apm_genes: list[APMGeneStatus] = []
        legacy_fields: dict[str, str] = {}

        for gene, items in sorted(by_gene.items()):
            has_lof = any(x.impact_class == "truncating_lof" for x in items)
            dna_status = "truncating_lof" if has_lof else "missense_or_inframe"
            notes_parts: list[str] = []
            for x in items:
                tag = f"{x.variant_id} ({x.effect_description or x.effect_type})"
                if x.rna_confirmed is True:
                    tag += f" [RNA confirmed: {x.rna_alt_reads} alt reads]"
                elif x.rna_confirmed is False:
                    tag += " [0 tumor RNA alt reads]"
                if x.phasing_status == "unphased_heterozygous_mhc_locus":
                    tag += (
                        f" [unphased across heterozygous alleles: "
                        f"{', '.join(x.matched_patient_alleles)}]"
                    )
                notes_parts.append(tag)
            apm_genes.append(
                APMGeneStatus(
                    gene=gene,
                    dna_status=dna_status,
                    rna_status="unknown",
                    dna_variants=tuple(x.to_dict() for x in items),
                    assay_source="dna_somatic_variant_detector",
                    notes="; ".join(notes_parts),
                )
            )
            if gene in _LEGACY_APM_ATTR:
                legacy_fields[_LEGACY_APM_ATTR[gene]] = (
                    "truncating_lof" if has_lof else "missense_variant"
                )

        mhc_records = tuple(
            MHCAlleleStatus(
                allele=allele,
                status="unknown",
                dna_somatic_variants=tuple(vdicts),
                assay_source="dna_somatic_variant_detector",
                notes="Auto-detected DNA somatic variant on homozygous MHC locus",
            )
            for allele, vdicts in sorted(allele_variants.items())
        )
        apm_obj = (
            APMStatus(genes=tuple(apm_genes), **legacy_fields)
            if apm_genes
            else None
        )
        return MHCAlleleStatusBundle(
            alleles=mhc_records,
            apm=apm_obj,
            source_path=source_path,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": "vaxrank.apm_variant_detection.v1",
            "policy": self.policy.to_dict(),
            "has_findings": self.has_findings,
            "has_disruptive_findings": self.has_disruptive_findings,
            "variants": [v.to_dict() for v in self.variants],
            "silent_variant_count": len(self.silent_variants),
        }


def _extract_variant_genes_and_effect(variant: Any, isovar_result: Any = None):
    """Extract ``(gene_names, top_effect, species)`` safely from a ``varcode.Variant``
    or stub without crashing on unannotated contigs."""
    species = ""
    genome = getattr(variant, "genome", None)
    if genome is not None:
        species = str(getattr(genome, "species", "") or "")

    effect = getattr(isovar_result, "predicted_effect", None) if isovar_result is not None else None
    if effect is None and hasattr(variant, "effects"):
        try:
            effects = variant.effects()
            if effects is not None:
                effect = effects.top_priority_effect()
                effect = select_varcode_effect_outcome(effect)
        except Exception:
            effect = None

    gene_names: list[str] = []
    if effect is not None:
        eg = getattr(effect, "gene_name", None)
        if eg:
            gene_names.append(str(eg))
        gene_obj = getattr(effect, "gene", None)
        if gene_obj is not None and getattr(gene_obj, "name", None):
            gn = str(gene_obj.name)
            if gn not in gene_names:
                gene_names.append(gn)
    try:
        for gn in (getattr(variant, "gene_names", None) or ()):
            if gn and str(gn) not in gene_names:
                gene_names.append(str(gn))
    except Exception:
        pass
    return gene_names, effect, species


def detect_dna_mhc_and_apm_variants(
    variants: Iterable[Any],
    *,
    isovar_results: Iterable[Any] | None = None,
    patient_mhc_alleles: Iterable[str] | None = None,
    dna_vaf_by_variant: dict[Any, Any] | None = None,
    policy: APMVariantPolicy | dict[str, Any] | None = None,
) -> APMVariantDetectionResult:
    """Scan somatic variants for cross-species DNA mutations in APM, IFN-gamma,
    and classical MHC Class I genes.

    Parameters
    ----------
    variants : iterable of varcode.Variant (or compatible objects)
        All loaded somatic variants (including variants that may later fail
        vaccine-peptide filters — e.g. a frameshift in B2M that yields no
        binders MUST still be caught!).
    isovar_results : optional iterable of IsovarResult
        Optional RNA-backed results to attach ``rna_alt_reads`` / ``rna_confirmed``.
    patient_mhc_alleles : optional iterable of str
        Patient MHC genotype for homozygous vs heterozygous MHC locus phasing.
    dna_vaf_by_variant : optional dict
        DNA VAF map keyed by variant or ``variant.short_description``.
    policy : APMVariantPolicy or dict or None
        Explicit per-trial policy.
    """
    pol = APMVariantPolicy.resolve(policy)
    if not pol.enabled:
        return APMVariantDetectionResult(policy=pol)

    dna_vaf_by_variant = dna_vaf_by_variant or {}
    isovar_by_var: dict[Any, Any] = {}
    isovar_by_desc: dict[str, Any] = {}
    for res in isovar_results or ():
        v = getattr(res, "variant", None)
        if v is not None:
            try:
                isovar_by_var[v] = res
            except TypeError:
                pass
            desc = getattr(v, "short_description", None)
            if desc:
                isovar_by_desc[str(desc)] = res

    findings: list[SomaticAPMorMHCVariant] = []
    silent: list[SomaticAPMorMHCVariant] = []

    for variant in variants or ():
        desc = str(getattr(variant, "short_description", None) or variant)
        iso_res = None
        try:
            iso_res = isovar_by_var.get(variant)
        except TypeError:
            iso_res = None
        if iso_res is None:
            iso_res = isovar_by_desc.get(desc)

        gene_names, effect, species = _extract_variant_genes_and_effect(variant, iso_res)
        matched_gene_info = None
        raw_gene = ""
        for gn in gene_names:
            canon, role, locus = classify_gene_role(gn, pol)
            if canon is not None and role is not None:
                matched_gene_info = (canon, role, locus)
                raw_gene = gn
                break
        if matched_gene_info is None:
            continue

        canon_gene, role, mhc_locus = matched_gene_info
        impact_class, effect_type, effect_desc = classify_variant_effect_impact(effect)

        # DNA VAF lookup
        vaf_val = None
        raw_vaf = None
        try:
            raw_vaf = dna_vaf_by_variant.get(variant)
        except TypeError:
            raw_vaf = None
        if raw_vaf is None:
            raw_vaf = dna_vaf_by_variant.get(desc)
        if isinstance(raw_vaf, dict):
            v = raw_vaf.get("tumor_vaf")
            if v is None:
                v = raw_vaf.get("vaf")
            vaf_val = float(v) if isinstance(v, (int, float)) else None
        elif isinstance(raw_vaf, (int, float)):
            vaf_val = float(raw_vaf)

        # Optional RNA alt-read confirmation from Isovar (kept strictly separate from DNA status)
        rna_alt = None
        rna_ref = None
        rna_confirmed = None
        if iso_res is not None:
            alt_c = getattr(iso_res, "num_alt_supporting_reads", None)
            if alt_c is None:
                alt_c = getattr(iso_res, "n_alt_reads", None)
            ref_c = getattr(iso_res, "num_ref_supporting_reads", None)
            if ref_c is None:
                ref_c = getattr(iso_res, "n_ref_reads", None)
            if alt_c is not None:
                rna_alt = int(alt_c)
                rna_confirmed = rna_alt >= pol.min_rna_alt_reads_confirmed
            if ref_c is not None:
                rna_ref = int(ref_c)

        matched_alleles = (
            _match_patient_alleles_for_mhc_locus(mhc_locus, patient_mhc_alleles)
            if role == "mhc_class_i"
            else ()
        )
        if role == "mhc_class_i":
            if len(matched_alleles) == 1:
                phasing = "homozygous_mhc_locus"
            elif len(matched_alleles) > 1:
                phasing = "unphased_heterozygous_mhc_locus"
            else:
                phasing = "mhc_locus_not_in_genotype"
        else:
            phasing = "not_applicable"

        record = SomaticAPMorMHCVariant(
            gene=canon_gene,
            raw_gene=raw_gene,
            pathway_role=role,
            variant_id=desc,
            impact_class=impact_class,
            effect_type=effect_type,
            effect_description=effect_desc,
            dna_vaf=vaf_val,
            rna_alt_reads=rna_alt,
            rna_ref_reads=rna_ref,
            rna_confirmed=rna_confirmed,
            matched_patient_alleles=matched_alleles,
            phasing_status=phasing,
            species=species,
        )
        if impact_class == "silent_or_noncoding":
            silent.append(record)
        else:
            findings.append(record)

    return APMVariantDetectionResult(
        variants=tuple(findings),
        silent_variants=tuple(silent),
        policy=pol,
    )
