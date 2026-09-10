# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Lossless processing evidence, not peptide ranking or a safety verdict.

Use mhctools' input/model/site contract. Bond b splits sequence[:b] from
sequence[b:]; neither endpoint is an internal peptide bond.
"""

from dataclasses import dataclass, field
import math
from numbers import Real
from typing import Optional

from mhctools.cleavage import CleavageInput, CleavageModel, CleavageResult, CleavageSite
from serializable import DataclassSerializable


@dataclass(frozen=True)
class CleavageProfile(DataclassSerializable):
    """Exact input/model observations with execution and coverage provenance.

    Unreturned bonds remain unassessed. Pepsickle's complete residue-indexed
    vector is retained separately: the final entry is an endpoint sentinel,
    not evidence of resistance. Padded bonds used absent sequence context.
    Native model scores and qualitative motif decisions remain distinct.
    """

    peptide: CleavageInput
    model: CleavageModel
    status: str
    sites: tuple[CleavageSite, ...] = field(default_factory=tuple)
    reason_codes: tuple[str, ...] = field(default_factory=tuple)
    error_message: str = ""
    backend_version: str = ""
    settings: tuple[tuple[str, str], ...] = field(default_factory=tuple)
    model_asset_sha256: str = ""
    raw_residue_scores: tuple[float, ...] = field(default_factory=tuple)
    padded_bonds: tuple[int, ...] = field(default_factory=tuple)

    def __post_init__(self):
        if not isinstance(self.peptide, CleavageInput) or not isinstance(self.model, CleavageModel):
            raise ValueError("Cleavage profile requires typed mhctools input and model")
        if any(not isinstance(value, str) for value in (
                self.error_message, self.backend_version, self.model_asset_sha256)):
            raise ValueError("Cleavage execution provenance must be text")
        if self.status not in {"observations_returned", "no_observations", "unassessed", "failed"}:
            raise ValueError("Unknown cleavage profile status")
        for name in ("sites", "reason_codes", "raw_residue_scores", "padded_bonds"):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        if any(not isinstance(reason, str) or not reason for reason in self.reason_codes):
            raise ValueError("Cleavage reasons must be nonempty text")
        settings = tuple(tuple(pair) for pair in self.settings)
        if any(len(pair) != 2 or not all(isinstance(v, str) for v in pair) for pair in settings):
            raise ValueError("Cleavage settings require immutable string key/value pairs")
        if len({pair[0] for pair in settings}) != len(settings):
            raise ValueError("Duplicate cleavage setting")
        object.__setattr__(self, "settings", tuple(sorted(settings)))
        if any(not isinstance(site, CleavageSite) for site in self.sites):
            raise ValueError("Cleavage sites must be typed mhctools observations")
        try:
            hash((self.peptide, self.model, self.sites))
        except TypeError as error:
            raise ValueError("Cleavage input/model/site metadata must be immutable") from error
        # Upstream owns bond bounds, duplicate sites and native score semantics.
        CleavageResult(self.peptide, self.model, self.sites)
        if bool(self.sites) != (self.status == "observations_returned"):
            raise ValueError("Cleavage status disagrees with returned sites")
        if self.status != "observations_returned" and not self.reason_codes:
            raise ValueError("Missing cleavage observations require an explanation")
        if self.status == "failed" and not self.error_message:
            raise ValueError("Failed cleavage assessment requires its error")
        if self.model_asset_sha256 and (
                len(self.model_asset_sha256) != 64
                or any(c not in "0123456789abcdef" for c in self.model_asset_sha256)):
            raise ValueError("Invalid cleavage model asset SHA-256")
        if self.raw_residue_scores:
            scores = _residue_scores(self.raw_residue_scores, len(self.peptide.sequence))
            object.__setattr__(self, "raw_residue_scores", scores)
            if self.status not in {"observations_returned", "no_observations"}:
                raise ValueError("Unassessed/failed profiles cannot claim residue scores")
            expected = tuple((b, scores[b - 1]) for b in range(1, len(scores)))
            if tuple((site.bond, site.score) for site in self.sites) != expected:
                raise ValueError("Bond scores disagree with complete residue output")
        if any(type(b) is not int or not 0 < b < len(self.peptide.sequence)
               for b in self.padded_bonds):
            raise ValueError("Padded positions must name internal peptide bonds")
        if len(set(self.padded_bonds)) != len(self.padded_bonds):
            raise ValueError("Duplicate padded bond")
        if not set(self.padded_bonds) <= {site.bond for site in self.sites}:
            raise ValueError("Padding annotations require observed sites")

    @property
    def unassessed_bonds(self):
        observed = {site.bond for site in self.sites}
        return tuple(b for b in range(1, len(self.peptide.sequence)) if b not in observed)

    @property
    def terminal_output(self):
        return self.raw_residue_scores[-1] if self.raw_residue_scores else None

    def interval_evidence(self, start, end):
        """Every internal/boundary observation for an exact half-open interval.

        A zero-width target junction is a boundary, not an internal ligand cut.
        """
        if type(start) is not int or type(end) is not int or not (
                0 <= start <= end <= len(self.peptide.sequence)):
            raise ValueError("Invalid cleavage overlay interval")
        by_bond = {site.bond: site for site in self.sites}
        def boundary_status(position):
            if position in {0, len(self.peptide.sequence)}:
                return "sequence_endpoint"
            return "observed" if position in by_bond else "unassessed"
        return CleavageIntervalEvidence(
            start, end,
            tuple(by_bond[b] for b in range(start + 1, end) if b in by_bond),
            tuple(b for b in range(start + 1, end) if b not in by_bond),
            by_bond.get(start), by_bond.get(end),
            boundary_status(start), boundary_status(end))


@dataclass(frozen=True)
class CleavageIntervalEvidence(DataclassSerializable):
    """Exact bonds relative to one ligand or target interval; no summary score."""

    start: int
    end: int
    internal_sites: tuple[CleavageSite, ...]
    unassessed_internal_bonds: tuple[int, ...]
    n_boundary: Optional[CleavageSite]
    c_boundary: Optional[CleavageSite]
    n_boundary_status: str
    c_boundary_status: str

    def __post_init__(self):
        if type(self.start) is not int or type(self.end) is not int or not 0 <= self.start <= self.end:
            raise ValueError("Invalid cleavage interval coordinates")
        object.__setattr__(self, "internal_sites", tuple(self.internal_sites))
        object.__setattr__(self, "unassessed_internal_bonds", tuple(self.unassessed_internal_bonds))
        if any(not isinstance(site, CleavageSite) for site in self.internal_sites):
            raise ValueError("Internal sites must be typed cleavage observations")
        observed = [site.bond for site in self.internal_sites]
        missing = list(self.unassessed_internal_bonds)
        if any(type(b) is not int for b in missing) or (
                len(set(observed + missing)) != len(observed + missing)
                or set(observed + missing) != set(range(self.start + 1, self.end))):
            raise ValueError("Internal bond coverage must be complete and disjoint")
        for position, site, status in (
                (self.start, self.n_boundary, self.n_boundary_status),
                (self.end, self.c_boundary, self.c_boundary_status)):
            if status not in {"observed", "unassessed", "sequence_endpoint"}:
                raise ValueError("Unknown boundary coverage status")
            if (site is not None) != (status == "observed") or (
                    site is not None and (not isinstance(site, CleavageSite) or site.bond != position)):
                raise ValueError("Boundary evidence disagrees with coordinates/status")


def _residue_scores(values, length):
    scores = tuple(values)
    if len(scores) != length:
        raise ValueError("Cleavage output must have exactly one value per input residue")
    if any(isinstance(v, bool) or not isinstance(v, Real)
           or not math.isfinite(v) or not 0 <= v <= 1 for v in scores):
        raise ValueError("Cleavage residue outputs must be finite probabilities in [0, 1]")
    return tuple(float(v) for v in scores)


def profile_from_residue_scores(peptide, model, scores, *, upstream_context=0,
                               downstream_context=0, **provenance):
    """Adapt processing output without counting its terminal sentinel as a cut.

    Context sizes refer to residues before/after the scored residue (P1), not
    the bond. Pepsickle's epitope model uses 8 before P1 and 8 after P1.
    """
    scores = _residue_scores(scores, len(peptide.sequence))
    if any(type(n) is not int or n < 0 for n in (upstream_context, downstream_context)):
        raise ValueError("Context sizes must be nonnegative integers")
    sites = tuple(CleavageSite(b, "scored", "residue_indexed_model_output", scores[b - 1])
                  for b in range(1, len(scores)))
    padded = tuple(b for b in range(1, len(scores))
                   if b - 1 < upstream_context or len(scores) - b < downstream_context)
    reasons = ("terminal_context_padded",) if padded else ()
    if not sites:
        reasons += ("sequence_has_no_internal_bonds",)
    return CleavageProfile(
        peptide, model, "observations_returned" if sites else "no_observations",
        sites=sites, raw_residue_scores=scores, padded_bonds=padded,
        reason_codes=reasons, **provenance)


def profile_from_cleavage_result(result, **provenance):
    """Retain a canonical enzyme result, including limited bond coverage."""
    if not isinstance(result, CleavageResult):
        raise ValueError("Expected a mhctools CleavageResult")
    if result.unsupported_reason is not None:
        status, reasons = "unassessed", (result.unsupported_reason,)
    elif result.sites:
        status, reasons = "observations_returned", ()
    else:
        status, reasons = "no_observations", ("model_returned_no_bond_observations",)
    return CleavageProfile(result.peptide, result.model, status,
                          sites=result.sites, reason_codes=reasons, **provenance)
