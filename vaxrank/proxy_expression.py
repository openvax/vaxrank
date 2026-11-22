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

"""Whole-gene cohort proxy expression from OncoRef for DNA-only ranking.

Cohort reference expression (e.g. TCGA/Treehouse median clean TPM for ``SKCM``)
is a prior over whether a gene is expressed in a given human cancer type when
patient tumor RNA-seq is unavailable. It is *never* treated as observed patient
RNA-seq reads: provenance stays explicitly labeled as cohort proxy expression
and ``n_alt_reads`` remains 0.
"""

from dataclasses import dataclass
import math
from typing import Optional

HUMAN_TAXON_ID = 9606


def _optional_float(value) -> Optional[float]:
    if value is None:
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _optional_int(value) -> Optional[int]:
    parsed = _optional_float(value)
    return int(round(parsed)) if parsed is not None else None


def _normalize_gene_id(gene_id: Optional[str]) -> str:
    if not gene_id or not isinstance(gene_id, str):
        return ""
    cleaned = gene_id.strip()
    if cleaned.startswith("ENSG"):
        cleaned = cleaned.split(".", 1)[0]
    return cleaned.upper()


def _normalize_symbol(symbol: Optional[str]) -> str:
    if not symbol or not isinstance(symbol, str):
        return ""
    return symbol.strip().upper()


@dataclass(frozen=True)
class ProxyGeneExpression:
    """Whole-gene cohort expression summary for one gene from OncoRef."""

    tpm: float
    cancer_code: str
    q1_tpm: Optional[float] = None
    q3_tpm: Optional[float] = None
    n_samples: Optional[int] = None
    source_cohort: str = ""
    normalization: str = "tpm_clean"
    data_version: str = ""
    matched_by: str = ""


class OncoRefProxyExpressionLookup:
    """Fast whole-gene median TPM lookup for a single OncoRef human cancer type.

    Parameters
    ----------
    cancer_type : str
        Canonical OncoRef cancer code (e.g. ``SKCM``, ``BRCA``, ``SARC_OS``)
        or common-name alias (e.g. ``melanoma``, ``breast``, ``osteosarcoma``).
    """

    def __init__(self, cancer_type: str):
        if not isinstance(cancer_type, str) or not cancer_type.strip():
            raise ValueError("--proxy-expression-cancer-type must be a non-empty cancer code or name")
        try:
            import oncoref
        except ImportError as error:
            raise ValueError(
                "--proxy-expression-cancer-type requires the `oncoref` package to be installed"
            ) from error

        self.requested_cancer_type = cancer_type.strip()
        self.cancer_code = oncoref.resolve_cancer_type(self.requested_cancer_type)
        df = oncoref.cancer_reference_expression(self.cancer_code)
        if df is None or df.empty:
            raise ValueError(
                f"OncoRef has no whole-gene reference expression cohort for cancer type "
                f"{self.cancer_code!r} (requested as {self.requested_cancer_type!r})."
            )

        meta = {}
        if hasattr(oncoref, "cancer_reference_expression_source_metadata"):
            try:
                meta = oncoref.cancer_reference_expression_source_metadata(self.cancer_code) or {}
            except Exception:
                meta = {}
        self.source_cohort = str(
            meta.get("source_cohort")
            or (df["source_cohort"].dropna().iloc[0] if "source_cohort" in df and not df["source_cohort"].dropna().empty else "")
        )
        self.source_project = str(
            meta.get("source_project")
            or (df["source_project"].dropna().iloc[0] if "source_project" in df and not df["source_project"].dropna().empty else "")
        )
        self.n_reference_samples = _optional_int(
            meta.get("n_reference_samples")
            if meta.get("n_reference_samples") is not None
            else (df["n_samples"].dropna().iloc[0] if "n_samples" in df and not df["n_samples"].dropna().empty else None)
        )
        self.data_version = str(
            df["data_version"].dropna().iloc[0]
            if "data_version" in df and not df["data_version"].dropna().empty
            else ""
        )
        self.oncoref_version = str(getattr(oncoref, "__version__", "") or "")
        self.source_label = self._build_source_label()

        self._by_gene_id: dict[str, ProxyGeneExpression] = {}
        self._by_symbol: dict[str, ProxyGeneExpression] = {}

        has_gene_id = "Ensembl_Gene_ID" in df.columns
        has_symbol = "Symbol" in df.columns
        has_q1 = "q1" in df.columns
        has_q3 = "q3" in df.columns
        has_n = "n_samples" in df.columns
        has_norm = "normalization" in df.columns
        has_cohort = "source_cohort" in df.columns
        has_ver = "data_version" in df.columns

        for row in df.itertuples(index=False):
            tpm = _optional_float(getattr(row, "expression", None))
            if tpm is None:
                continue
            tpm = max(0.0, tpm)
            q1 = _optional_float(getattr(row, "q1", None)) if has_q1 else None
            q3 = _optional_float(getattr(row, "q3", None)) if has_q3 else None
            n_samples = _optional_int(getattr(row, "n_samples", None)) if has_n else self.n_reference_samples
            norm = str(getattr(row, "normalization", "") or "tpm_clean") if has_norm else "tpm_clean"
            cohort = str(getattr(row, "source_cohort", "") or self.source_cohort) if has_cohort else self.source_cohort
            ver = str(getattr(row, "data_version", "") or self.data_version) if has_ver else self.data_version

            gene_id = _normalize_gene_id(getattr(row, "Ensembl_Gene_ID", None)) if has_gene_id else ""
            symbol = _normalize_symbol(getattr(row, "Symbol", None)) if has_symbol else ""

            if gene_id:
                existing = self._by_gene_id.get(gene_id)
                if existing is None or tpm > existing.tpm:
                    self._by_gene_id[gene_id] = ProxyGeneExpression(
                        tpm=tpm,
                        cancer_code=self.cancer_code,
                        q1_tpm=q1,
                        q3_tpm=q3,
                        n_samples=n_samples,
                        source_cohort=cohort,
                        normalization=norm,
                        data_version=ver,
                        matched_by="ensembl_gene_id",
                    )
            if symbol:
                existing = self._by_symbol.get(symbol)
                if existing is None or tpm > existing.tpm:
                    self._by_symbol[symbol] = ProxyGeneExpression(
                        tpm=tpm,
                        cancer_code=self.cancer_code,
                        q1_tpm=q1,
                        q3_tpm=q3,
                        n_samples=n_samples,
                        source_cohort=cohort,
                        normalization=norm,
                        data_version=ver,
                        matched_by="symbol",
                    )

    def _build_source_label(self) -> str:
        parts = ["OncoRef"]
        if self.source_project:
            parts.append(self.source_project)
        elif self.source_cohort:
            parts.append(self.source_cohort)
        if self.n_reference_samples:
            parts.append(f"n={self.n_reference_samples}")
        return " · ".join(parts)

    def lookup(
        self,
        *,
        gene_id: Optional[str] = None,
        gene_name: Optional[str] = None,
    ) -> Optional[ProxyGeneExpression]:
        """Look up whole-gene cohort median clean TPM by Ensembl gene ID, then symbol."""
        norm_id = _normalize_gene_id(gene_id)
        if norm_id and norm_id in self._by_gene_id:
            return self._by_gene_id[norm_id]
        norm_sym = _normalize_symbol(gene_name)
        if norm_sym and norm_sym in self._by_symbol:
            return self._by_symbol[norm_sym]
        return None


def validate_human_host_for_proxy_expression(host_taxon_id: Optional[int], cancer_type: str) -> None:
    """Reject non-human hosts when human OncoRef cohort expression is requested."""
    if host_taxon_id is not None and int(host_taxon_id) != HUMAN_TAXON_ID:
        raise ValueError(
            f"--proxy-expression-cancer-type {cancer_type!r} uses human OncoRef cancer cohorts "
            f"(NCBI taxon {HUMAN_TAXON_ID}), which cannot be applied to host taxon {host_taxon_id}."
        )


def _transcript_gene_identity(transcript) -> tuple[Optional[str], Optional[str]]:
    if transcript is None:
        return None, None
    gene_id = getattr(transcript, "gene_id", None) or getattr(getattr(transcript, "gene", None), "id", None)
    gene_name = getattr(transcript, "gene_name", None) or getattr(getattr(transcript, "gene", None), "name", None)
    return (str(gene_id) if gene_id else None), (str(gene_name) if gene_name else None)


def _resolve_transcript_by_id(transcript_id: Optional[str], transcripts_by_id: dict, genome) -> object:
    if not transcript_id:
        return None
    if transcript_id in transcripts_by_id:
        return transcripts_by_id[transcript_id]
    if genome is not None and hasattr(genome, "transcript_by_id"):
        try:
            return genome.transcript_by_id(transcript_id)
        except Exception:
            return None
    return None


def resolve_proxy_expression_gene(fragment) -> tuple[Optional[str], Optional[str], bool]:
    """Resolve ``(gene_id, gene_name, is_upstream_complex)`` for cohort proxy lookup.

    For structural / complex variants (e.g. gene fusions, translocations,
    multi-segment rearrangements), transcription of the chimeric product is
    driven by the 5' gene *upstream of the breakpoint* rather than the 3'
    downstream partner. Resolution priority:

    1. Explicit 5' (upstream-of-breakpoint) transcript from
       ``dna_effect_selection['five_prime_transcript_id']`` or
       ``effect.five_prime_transcript``.
    2. First (5'-most, upstream-of-breakpoint) transcript segment in
       ``effect.mutant_transcript.reference_segments`` for multi-segment
       rearrangements.
    3. If ``dna_effect_selection['three_prime_transcript_id']`` matches the
       primary transcript on a two-transcript fusion without an explicit 5' ID,
       select the partner transcript as the 5' upstream gene.
    4. Primary transcript (``supporting_reference_transcripts[0]`` /
       ``effect.transcript``) and ``fragment.gene_name``.
    """
    transcripts = tuple(getattr(fragment, "supporting_reference_transcripts", ()) or ())
    transcripts_by_id = {
        str(getattr(t, "id", None) or getattr(t, "transcript_id", "")): t
        for t in transcripts
        if getattr(t, "id", None) or getattr(t, "transcript_id", None)
    }
    variant = getattr(fragment, "variant", None)
    genome = getattr(variant, "genome", None)
    selection = getattr(fragment, "dna_effect_selection", None) or {}

    # 1. Explicit 5' (upstream-of-breakpoint) transcript from selection provenance
    five_prime_id = selection.get("five_prime_transcript_id")
    if five_prime_id:
        tx = _resolve_transcript_by_id(str(five_prime_id), transcripts_by_id, genome)
        gene_id, gene_name = _transcript_gene_identity(tx)
        if gene_id or gene_name:
            return gene_id, gene_name, True

    # Inspect predicted effect for complex / structural variants when available.
    effect = None
    try:
        effect = fragment.predicted_effect() if hasattr(fragment, "predicted_effect") else None
    except Exception:
        effect = None

    five_prime_tx = getattr(effect, "five_prime_transcript", None)
    if five_prime_tx is not None:
        gene_id, gene_name = _transcript_gene_identity(five_prime_tx)
        if gene_id or gene_name:
            return gene_id, gene_name, True

    # 2. First (5'-most / upstream-of-breakpoint) segment in multi-segment MutantTranscript
    model = getattr(effect, "mutant_transcript", None)
    segments = getattr(model, "reference_segments", None) or ()
    if segments:
        first_source = getattr(segments[0], "source", None)
        if isinstance(first_source, str):
            first_tx = _resolve_transcript_by_id(first_source, transcripts_by_id, genome)
        else:
            first_tx = first_source
        gene_id, gene_name = _transcript_gene_identity(first_tx)
        if gene_id or gene_name:
            is_complex = len(segments) > 1 or bool(selection.get("partner_transcript_id"))
            return gene_id, gene_name, is_complex

    # 3. If primary transcript is explicitly the 3' partner, pick the other transcript
    three_prime_id = selection.get("three_prime_transcript_id")
    primary_id = selection.get("transcript_id")
    partner_id = selection.get("partner_transcript_id")
    if three_prime_id and primary_id and three_prime_id == primary_id and partner_id:
        tx = _resolve_transcript_by_id(str(partner_id), transcripts_by_id, genome)
        gene_id, gene_name = _transcript_gene_identity(tx)
        if gene_id or gene_name:
            return gene_id, gene_name, True

    # 4. Default: primary transcript / fragment gene
    primary_tx = transcripts[0] if transcripts else getattr(effect, "transcript", None)
    gene_id, gene_name = _transcript_gene_identity(primary_tx)
    gene_name = gene_name or (getattr(fragment, "gene_name", None) or None)
    return gene_id, gene_name, False
