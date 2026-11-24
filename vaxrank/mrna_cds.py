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

"""Synthesis-safe mRNA CDS optimization, sequence QC, and LinearDesign/MFE hooks.

Extends DnaChisel back-translation with explicit, per-trial YAML-configurable
manufacturing and translational constraints (``mrna.cds_optimization:``):

- Sliding-window GC bounds (preventing Twist/IDT synthesis complexity failures)
- Homopolymer avoidance (e.g. runs >= 6 nt)
- Cryptic T7 RNA polymerase class-II terminator (``ATCTGTT``) and restriction
  site avoidance (named IVT linearization/Golden Gate enzymes or raw motifs)
- Uridine depletion (minimizing U / UU content to reduce TLR7/8 sensing and
  +1 ribosomal frameshifting)
- Optional ``LinearDesign`` backend hook for joint MFE + CAI optimization
- Always-on CDS sequence QC telemetry (GC%, min/max window GC%, U%, UU count,
  max homopolymer, T7 terminator count, CAI, 5' initiation-window GC%, and
  optional ViennaRNA MFE) recorded in ``manifest.json`` and reports.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields
import logging
import math
import shutil
import subprocess
from typing import Optional

from Bio.Seq import Seq
from dnachisel import (
    AvoidChanges,
    AvoidPattern,
    CodonOptimize,
    DnaOptimizationProblem,
    EnforceGCContent,
    EnforceTranslation,
    Location,
)
from dnachisel.DnaOptimizationProblem import NoSolutionError
from dnachisel.biotools import reverse_translate

logger = logging.getLogger(__name__)

T7_CLASS_II_TERMINATOR = "ATCTGTT"

# Common IVT linearization and Golden Gate assembly recognition sites (5' -> 3').
RESTRICTION_SITES: dict[str, str] = {
    "BsaI": "GGTCTC",
    "BbsI": "GAAGAC",
    "BsmBI": "CGTCTC",
    "SapI": "GCTCTTC",
    "BspQI": "GCTCTTC",
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "XhoI": "CTCGAG",
    "NotI": "GCGGCCGC",
    "SpeI": "ACTAGT",
    "XbaI": "TCTAGA",
    "PmeI": "GTTTAAAC",
    "NheI": "GCTAGC",
    "KpnI": "GGTACC",
    "SbfI": "CCTGCAGG",
    "PacI": "TTAATTAA",
}


def resolve_restriction_motif(token: str) -> str:
    """Resolve an enzyme name (case-insensitive) or raw DNA motif."""
    cleaned = str(token or "").strip()
    if not cleaned:
        raise ValueError("Empty restriction site / pattern in avoid_restriction_sites")
    for name, motif in RESTRICTION_SITES.items():
        if name.lower() == cleaned.lower():
            return motif
    upper = cleaned.upper().replace("U", "T")
    if all(ch in "ACGT" for ch in upper):
        return upper
    raise ValueError(
        f"Unknown restriction enzyme or non-ACGT motif {token!r}. "
        f"Known enzymes: {', '.join(sorted(RESTRICTION_SITES))}"
    )


def _check_prob(name: str, value: Optional[float]) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"MRNACDSOptimizationPolicy.{name} must be finite; got {value!r}")
    if not (0.0 <= float(value) <= 1.0):
        raise ValueError(f"MRNACDSOptimizationPolicy.{name} must be in [0, 1]; got {value!r}")
    return float(value)


def _check_nonneg_float(name: str, value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError(f"MRNACDSOptimizationPolicy.{name} must be finite and >= 0; got {value!r}")
    return float(value)


def _check_pos_int(name: str, value: Optional[int], *, allow_none: bool = False) -> Optional[int]:
    if value is None and allow_none:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"MRNACDSOptimizationPolicy.{name} must be a positive integer; got {value!r}")
    return value


@dataclass(frozen=True)
class MRNACDSOptimizationPolicy:
    """Explicit, serializable per-trial policy for mRNA CDS optimization & QC.

    Lives under ``mrna.cds_optimization:`` in YAML configs. Default values leave
    historical DnaChisel behavior untouched so existing designs replay with
    identical sequences and checksums; enable synthesis-safe rules via YAML,
    :meth:`synthesis_safe`, or ``--mrna-cds-synthesis-safe``.
    """

    optimizer: str = "dnachisel"  # 'dnachisel' | 'lineardesign'
    lineardesign_executable: str = "lineardesign"
    lineardesign_lambda: float = 0.0
    gc_min: Optional[float] = None
    gc_max: Optional[float] = None
    gc_window_nt: int = 50
    max_homopolymer_length: Optional[int] = None
    avoid_t7_terminator: bool = False
    avoid_restriction_sites: tuple[str, ...] = ()
    uridine_depletion: bool = False
    uridine_depletion_weight: float = 1.0
    avoid_uu_dinucleotides: bool = False
    compute_mfe: bool = False
    initiation_window_nt: int = 40

    def __post_init__(self):
        opt = str(self.optimizer or "").strip().lower()
        if opt not in ("dnachisel", "lineardesign"):
            raise ValueError(
                f"MRNACDSOptimizationPolicy.optimizer must be 'dnachisel' or 'lineardesign'; got {self.optimizer!r}"
            )
        object.__setattr__(self, "optimizer", opt)
        _check_nonneg_float("lineardesign_lambda", self.lineardesign_lambda)
        gmin = _check_prob("gc_min", self.gc_min)
        gmax = _check_prob("gc_max", self.gc_max)
        if gmin is not None and gmax is not None and gmin > gmax:
            raise ValueError(
                f"MRNACDSOptimizationPolicy.gc_min ({gmin}) must be <= gc_max ({gmax})"
            )
        _check_pos_int("gc_window_nt", self.gc_window_nt)
        _check_pos_int("max_homopolymer_length", self.max_homopolymer_length, allow_none=True)
        _check_pos_int("initiation_window_nt", self.initiation_window_nt)
        _check_nonneg_float("uridine_depletion_weight", self.uridine_depletion_weight)
        sites = tuple(str(s).strip() for s in (self.avoid_restriction_sites or ()) if str(s).strip())
        for s in sites:
            resolve_restriction_motif(s)
        object.__setattr__(self, "avoid_restriction_sites", sites)
        for flag in ("avoid_t7_terminator", "uridine_depletion", "avoid_uu_dinucleotides", "compute_mfe"):
            if not isinstance(getattr(self, flag), bool):
                raise ValueError(
                    f"MRNACDSOptimizationPolicy.{flag} must be a boolean; got {getattr(self, flag)!r}"
                )

    @classmethod
    def synthesis_safe(cls, **overrides) -> "MRNACDSOptimizationPolicy":
        """Preset enforcing DNA-template synthesis bounds, T7 avoidance, and U-depletion."""
        defaults = dict(
            gc_min=0.35,
            gc_max=0.68,
            gc_window_nt=50,
            max_homopolymer_length=5,
            avoid_t7_terminator=True,
            uridine_depletion=True,
            uridine_depletion_weight=1.0,
        )
        defaults.update(overrides)
        return cls(**defaults)

    @classmethod
    def resolve(cls, value=None) -> "MRNACDSOptimizationPolicy":
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if not isinstance(value, dict):
            raise TypeError(
                f"mrna.cds_optimization must be a mapping or MRNACDSOptimizationPolicy; got {type(value).__name__}"
            )
        valid = {f.name for f in fields(cls)}
        unknown = sorted(set(value) - valid)
        if unknown:
            raise ValueError(
                f"Unknown mrna.cds_optimization keys: {', '.join(unknown)}"
            )
        cleaned = {}
        for k, v in value.items():
            if v is None:
                continue
            if k == "avoid_restriction_sites" and isinstance(v, str):
                cleaned[k] = tuple(tok.strip() for tok in v.split(",") if tok.strip())
            elif k == "avoid_restriction_sites":
                cleaned[k] = tuple(v)
            else:
                cleaned[k] = v
        return cls(**cleaned)

    def has_active_constraints(self) -> bool:
        return (
            self.optimizer != "dnachisel"
            or self.gc_min is not None
            or self.gc_max is not None
            or self.max_homopolymer_length is not None
            or self.avoid_t7_terminator
            or bool(self.avoid_restriction_sites)
            or self.uridine_depletion
            or self.avoid_uu_dinucleotides
        )

    def to_dict(self) -> dict:
        data = asdict(self)
        data["avoid_restriction_sites"] = list(self.avoid_restriction_sites)
        return data


def _reverse_complement(dna: str) -> str:
    return str(Seq(dna).reverse_complement())


def _run_lineardesign_segment(aa: str, policy: MRNACDSOptimizationPolicy) -> str:
    """Invoke external LinearDesign CLI on one protein segment and verify translation."""
    exe = shutil.which(policy.lineardesign_executable)
    if exe is None:
        raise RuntimeError(
            f"LinearDesign executable {policy.lineardesign_executable!r} not found on PATH. "
            "Install LinearDesign or set mrna.cds_optimization.lineardesign_executable to its path."
        )
    cmd = [exe, "--lambda", str(policy.lineardesign_lambda)]
    proc = subprocess.run(
        cmd,
        input=aa + "\n",
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"LinearDesign exited with status {proc.returncode}: {proc.stderr.strip() or proc.stdout.strip()}"
        )
    seq_dna = None
    for line in proc.stdout.splitlines():
        text = line.strip()
        if not text:
            continue
        if ":" in text and text.lower().startswith("mrna sequence"):
            candidate = text.split(":", 1)[1].strip().upper().replace("U", "T")
        else:
            candidate = text.upper().replace("U", "T")
        if len(candidate) == len(aa) * 3 and all(ch in "ACGT" for ch in candidate):
            seq_dna = candidate
            break
    if seq_dna is None:
        raise RuntimeError(
            f"Could not parse a {len(aa) * 3}-nt coding sequence from LinearDesign output"
        )
    translated = str(Seq(seq_dna).translate())
    if translated != aa:
        raise RuntimeError(
            "LinearDesign output did not translate back to the requested amino-acid sequence"
        )
    return seq_dna


def optimize_cds_dna(
    amino_acids: str,
    *,
    species: str,
    method: str = "use_best_codon",
    avoid_patterns: tuple[str, ...] = (),
    frozen_segments: tuple = (),
    policy: Optional[MRNACDSOptimizationPolicy] = None,
    reverse_translate_fn=reverse_translate,
) -> str:
    """Optimize coding DNA under DnaChisel (with synthesis-safe rules) or LinearDesign."""
    if not amino_acids:
        return ""
    resolved = MRNACDSOptimizationPolicy.resolve(policy)
    frozen_segments = sorted(frozen_segments)

    if resolved.optimizer == "lineardesign":
        parts = []
        cursor = 0
        for start_aa, end_aa, blessed_dna in frozen_segments:
            if start_aa < cursor:
                raise ValueError("frozen_segments must not overlap")
            if (end_aa - start_aa) * 3 != len(blessed_dna):
                raise ValueError(
                    "blessed DNA length %d does not match AA window %d:%d"
                    % (len(blessed_dna), start_aa, end_aa)
                )
            prefix_aa = amino_acids[cursor:start_aa]
            if prefix_aa:
                parts.append(_run_lineardesign_segment(prefix_aa, resolved))
            parts.append(blessed_dna)
            cursor = end_aa
        suffix_aa = amino_acids[cursor:]
        if suffix_aa:
            parts.append(_run_lineardesign_segment(suffix_aa, resolved))
        return "".join(parts)

    initial_parts = []
    frozen_nt_ranges = []
    cursor = 0
    nt_offset = 0
    for start_aa, end_aa, blessed_dna in frozen_segments:
        if start_aa < cursor:
            raise ValueError("frozen_segments must not overlap")
        if (end_aa - start_aa) * 3 != len(blessed_dna):
            raise ValueError(
                "blessed DNA length %d does not match AA window %d:%d"
                % (len(blessed_dna), start_aa, end_aa)
            )
        prefix_aa = amino_acids[cursor:start_aa]
        if prefix_aa:
            prefix_dna = reverse_translate_fn(prefix_aa)
            initial_parts.append(prefix_dna)
            nt_offset += len(prefix_dna)
        initial_parts.append(blessed_dna)
        frozen_nt_ranges.append((nt_offset, nt_offset + len(blessed_dna)))
        nt_offset += len(blessed_dna)
        cursor = end_aa
    suffix_aa = amino_acids[cursor:]
    if suffix_aa:
        initial_parts.append(reverse_translate_fn(suffix_aa))
    initial = "".join(initial_parts)

    hard_constraints = [EnforceTranslation()]
    for start_nt, end_nt in frozen_nt_ranges:
        hard_constraints.append(AvoidChanges(location=Location(start_nt, end_nt)))
    for pattern in avoid_patterns:
        hard_constraints.append(AvoidPattern(pattern))

    extra_constraints = []
    if resolved.gc_min is not None or resolved.gc_max is not None:
        win = min(resolved.gc_window_nt, len(initial))
        extra_constraints.append(EnforceGCContent(
            mini=resolved.gc_min if resolved.gc_min is not None else 0.0,
            maxi=resolved.gc_max if resolved.gc_max is not None else 1.0,
            window=win,
        ))
    if resolved.max_homopolymer_length is not None:
        run_len = resolved.max_homopolymer_length + 1
        for base in ("A", "T", "C", "G"):
            extra_constraints.append(AvoidPattern(base * run_len))
    if resolved.avoid_t7_terminator:
        extra_constraints.append(AvoidPattern(T7_CLASS_II_TERMINATOR))
    for site in resolved.avoid_restriction_sites:
        motif = resolve_restriction_motif(site)
        extra_constraints.append(AvoidPattern(motif))

    objectives = [CodonOptimize(species=species, method=method)]
    if resolved.uridine_depletion and resolved.uridine_depletion_weight > 0:
        objectives.append(AvoidPattern("T", boost=resolved.uridine_depletion_weight))
    if resolved.avoid_uu_dinucleotides:
        objectives.append(AvoidPattern("TT", boost=max(resolved.uridine_depletion_weight, 1.0)))

    problem = DnaOptimizationProblem(
        sequence=initial,
        constraints=hard_constraints + extra_constraints,
        objectives=objectives,
        logger=None,
    )
    try:
        problem.resolve_constraints()
        problem.optimize()
        return problem.sequence
    except NoSolutionError:
        if not extra_constraints:
            raise
        logger.warning(
            "Strict synthesis-safe CDS constraints could not all be satisfied simultaneously "
            "(likely due to a low-complexity amino-acid run); relaxing synthesis constraints "
            "to high-priority soft objectives."
        )
        soft_objectives = [
            *(spec.copy_with_changes(boost=10.0) if hasattr(spec, "copy_with_changes") else spec
              for spec in extra_constraints),
            *objectives,
        ]
        relaxed = DnaOptimizationProblem(
            sequence=initial,
            constraints=hard_constraints,
            objectives=soft_objectives,
            logger=None,
        )
        relaxed.resolve_constraints()
        relaxed.optimize()
        return relaxed.sequence


def _compute_cai(coding_dna: str, species: Optional[str]) -> Optional[float]:
    if not coding_dna or not species or len(coding_dna) < 3:
        return None
    try:
        import python_codon_tables as pct
        table = pct.get_codons_table(species)
    except Exception:
        return None

    # Build per-codon relative adaptiveness w_c = f_c / max_synonymous(f).
    weights: dict[str, float] = {}
    for aa, codon_freqs in table.items():
        if aa == "*" or not codon_freqs:
            continue
        max_freq = max(codon_freqs.values())
        if max_freq <= 0:
            continue
        for codon, freq in codon_freqs.items():
            weights[codon.upper().replace("U", "T")] = max(freq / max_freq, 1e-6)

    log_sum = 0.0
    n_codons = 0
    seq = coding_dna.upper().replace("U", "T")
    for pos in range(0, len(seq) - (len(seq) % 3), 3):
        codon = seq[pos:pos + 3]
        w = weights.get(codon)
        if w is not None:
            log_sum += math.log(w)
            n_codons += 1
    if n_codons == 0:
        return None
    return round(math.exp(log_sum / n_codons), 4)


def _window_gc_bounds(seq: str, window: int) -> tuple[Optional[float], Optional[float]]:
    if not seq:
        return None, None
    w = min(max(1, window), len(seq))
    gc = sum(1 for ch in seq[:w] if ch in ("G", "C"))
    min_gc = gc / w
    max_gc = min_gc
    for i in range(w, len(seq)):
        if seq[i] in ("G", "C"):
            gc += 1
        if seq[i - w] in ("G", "C"):
            gc -= 1
        frac = gc / w
        if frac < min_gc:
            min_gc = frac
        if frac > max_gc:
            max_gc = frac
    return round(min_gc, 4), round(max_gc, 4)


def _max_homopolymer(seq: str) -> tuple[int, str]:
    if not seq:
        return 0, ""
    best_len = 1
    best_base = seq[0]
    cur_len = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            cur_len += 1
            if cur_len > best_len:
                best_len = cur_len
                best_base = seq[i]
        else:
            cur_len = 1
    return best_len, best_base


def _fold_mfe(dna_or_rna: str) -> tuple[Optional[float], str]:
    try:
        import RNA  # ViennaRNA optional dependency
    except ImportError:
        return None, "viennarna_not_installed"
    rna = dna_or_rna.upper().replace("T", "U")
    _, mfe = RNA.fold(rna)
    return round(float(mfe), 2), "viennarna"


def compute_mrna_sequence_qc(
    cds_nt: str,
    *,
    species: Optional[str] = None,
    policy: Optional[MRNACDSOptimizationPolicy] = None,
) -> dict:
    """Compute synthesis and translational sequence QC metrics for an mRNA CDS."""
    resolved = MRNACDSOptimizationPolicy.resolve(policy)
    seq = str(cds_nt or "").upper().replace("U", "T")
    n = len(seq)
    if n == 0:
        return {
            "cds_length_nt": 0,
            "gc_fraction": None,
            "gc_window_nt": resolved.gc_window_nt,
            "min_window_gc": None,
            "max_window_gc": None,
            "uridine_fraction": None,
            "uu_dinucleotide_count": 0,
            "max_homopolymer_length": 0,
            "homopolymer_base": "",
            "t7_terminator_count": 0,
            "cai": None,
            "initiation_window_nt": resolved.initiation_window_nt,
            "initiation_window_gc": None,
            "cds_mfe_kcal_mol": None,
            "initiation_window_mfe_kcal_mol": None,
            "mfe_status": "not_requested",
            "synthesis_flags": [],
        }

    gc_frac = round((seq.count("G") + seq.count("C")) / n, 4)
    u_frac = round(seq.count("T") / n, 4)
    uu_count = sum(1 for i in range(n - 1) if seq[i] == "T" and seq[i + 1] == "T")
    min_wgc, max_wgc = _window_gc_bounds(seq, resolved.gc_window_nt)
    max_hp, hp_base = _max_homopolymer(seq)
    t7_count = seq.count(T7_CLASS_II_TERMINATOR) + seq.count(_reverse_complement(T7_CLASS_II_TERMINATOR))
    init_slice = seq[:resolved.initiation_window_nt]
    init_gc = round((init_slice.count("G") + init_slice.count("C")) / len(init_slice), 4)
    cai = _compute_cai(seq, species)

    cds_mfe = None
    init_mfe = None
    mfe_status = "not_requested"
    if resolved.compute_mfe:
        cds_mfe, mfe_status = _fold_mfe(seq)
        if mfe_status == "viennarna":
            init_mfe, _ = _fold_mfe(init_slice)

    flags = []
    if resolved.gc_min is not None and min_wgc is not None and min_wgc < resolved.gc_min:
        flags.append(f"min_{resolved.gc_window_nt}nt_gc={min_wgc:.2f}<{resolved.gc_min:.2f}")
    if resolved.gc_max is not None and max_wgc is not None and max_wgc > resolved.gc_max:
        flags.append(f"max_{resolved.gc_window_nt}nt_gc={max_wgc:.2f}>{resolved.gc_max:.2f}")
    if resolved.max_homopolymer_length is not None and max_hp > resolved.max_homopolymer_length:
        flags.append(f"homopolymer_{hp_base}x{max_hp}>{resolved.max_homopolymer_length}")
    if resolved.avoid_t7_terminator and t7_count > 0:
        flags.append(f"t7_terminator_sites={t7_count}")
    for site in resolved.avoid_restriction_sites:
        motif = resolve_restriction_motif(site)
        rc = _reverse_complement(motif)
        hits = seq.count(motif) + (seq.count(rc) if rc != motif else 0)
        if hits > 0:
            flags.append(f"restriction_site_{site}={hits}")

    return {
        "cds_length_nt": n,
        "gc_fraction": gc_frac,
        "gc_window_nt": min(resolved.gc_window_nt, n),
        "min_window_gc": min_wgc,
        "max_window_gc": max_wgc,
        "uridine_fraction": u_frac,
        "uu_dinucleotide_count": uu_count,
        "max_homopolymer_length": max_hp,
        "homopolymer_base": hp_base,
        "t7_terminator_count": t7_count,
        "cai": cai,
        "initiation_window_nt": len(init_slice),
        "initiation_window_gc": init_gc,
        "cds_mfe_kcal_mol": cds_mfe,
        "initiation_window_mfe_kcal_mol": init_mfe,
        "mfe_status": mfe_status,
        "synthesis_flags": flags,
    }
