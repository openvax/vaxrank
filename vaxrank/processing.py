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

"""Proteasomal cleavage credibility tagging for MHC ligand predictions.

For each predicted MHC ligand, attach three scores. All three carry
the ``pepsickle_`` prefix so the predictor is unambiguous in CSVs,
log lines, and report columns; a future per-position cleavage
predictor (NetChop, PAProC, …) can land alongside without colliding.

  pepsickle_c_term_cleavage_prob    pepsickle's probability the
                                    proteasome cuts at the ligand's
                                    C-terminus (clean release)
  pepsickle_max_internal_cut_prob   peak cleavage probability strictly
                                    inside the ligand (high → ligand
                                    is destroyed before reaching MHC)
  pepsickle_processing_score        composite ``c_term * (1 - max_internal)``;
                                    1.0 = ideal release, 0.0 = no
                                    clean release OR near-certain
                                    destruction

The annotations are purely additive — vaccine ranking is unaffected.
Reports surface the three columns when at least one prediction in
the per-VaccinePeptide list has been annotated.

Pepsickle (Weeder et al., Bioinformatics 2021) runs in an isolated
subprocess (issue #266: torch's libomp clashes with the parent's
pandas/numpy/pyarrow OpenMP runtime on macOS, segfaulting the
process). The subprocess imports only torch + pepsickle, no pandas;
single libomp loads, no clash. Linux installs work either way.

Issue: openvax/vaxrank#249.
"""

import logging

logger = logging.getLogger(__name__)


# ----------------------------------------------------------------------------
# Per-peptide arithmetic on a precomputed per-position cleavage array
# ----------------------------------------------------------------------------

def _per_position_processing(seq_probs, start, length):
    """Slice a per-position cleavage array down to one peptide's
    ``(pepsickle_c_term_cleavage_prob, pepsickle_max_internal_cut_prob)``.

    ``seq_probs[i]`` is "would the proteasome cut C-terminally to
    residue *i*?". For a peptide spanning ``[start, start+length)``:
      - C-terminal cut = ``seq_probs[start+length-1]`` (release at end)
      - Max internal cut = ``max(seq_probs[start : start+length-1])``
        (positions strictly inside the peptide; destruction)

    Returns ``(c_term, max_internal)`` floats in [0, 1], or
    ``(None, None)`` if the peptide is out of range. Single-residue
    peptides have no "inside" — max_internal is 0.0.
    """
    if length < 1 or start < 0 or start + length > len(seq_probs):
        return None, None
    c_term = float(seq_probs[start + length - 1])
    if length >= 2:
        max_internal = float(max(seq_probs[start:start + length - 1]))
    else:
        max_internal = 0.0
    return c_term, max_internal


def _composite_processing_score(c_term, max_internal):
    """``c_term * (1 - max_internal)``. Range [0, 1].

    1.0 = clean cut at end + no internal cut.
    0.0 = no clean release, or near-certain internal destruction.

    ``1 - max_internal`` is a conservative approximation for "no
    internal cut anywhere"; the strictly correct expression is the
    joint product ``Π(1 - p_i)`` over all internal positions.
    Proteasomes cut roughly once per substrate molecule, so the
    single highest-probability cut dominates and this heuristic is
    a reasonable proxy. Same formula at any peptide length.

    Returns ``None`` when either input is None.
    """
    if c_term is None or max_internal is None:
        return None
    return float(c_term) * (1.0 - float(max_internal))


def _closest_occurrence(source, peptide, declared_offset):
    """Offset of ``peptide`` in ``source`` closest to ``declared_offset``,
    or ``None`` if not found.

    Bias toward the declared offset so repeated peptides (homopolymer
    tracts, short ligands appearing in multiple loops) don't snap to
    position 0 when the upstream loader recorded a later occurrence.
    """
    if not peptide or not source:
        return None
    best = None
    best_drift = None
    pos = source.find(peptide)
    while pos >= 0:
        drift = abs(pos - declared_offset)
        if best_drift is None or drift < best_drift:
            best, best_drift = pos, drift
        pos = source.find(peptide, pos + 1)
    return best


# ----------------------------------------------------------------------------
# Pepsickle invocation
# ----------------------------------------------------------------------------

# Subprocess script: reads {sequences, human_only, threshold} JSON
# from stdin, writes {sequence: [cleavage_probs]} JSON to stdout.
# Lives in a *clean* Python interpreter (only torch+pepsickle+json),
# so torch's libomp doesn't clash with the parent's pandas/numpy
# libomp on macOS. Issue #266.
_SUBPROCESS_SCRIPT = r"""
import sys, json
data = json.loads(sys.stdin.read())
try:
    from mhctools import Pepsickle
except Exception as e:
    sys.stderr.write("import_error: %s\n" % e)
    sys.exit(2)
try:
    p = Pepsickle(human_only=data["human_only"], threshold=data["threshold"])
except Exception as e:
    sys.stderr.write("instantiate_error: %s\n" % e)
    sys.exit(3)
out = {}
for seq in data["sequences"]:
    try:
        out[seq] = [float(x) for x in p.cleavage_probs(seq)]
    except Exception as e:
        sys.stderr.write("inference_error %r: %s\n" % (seq[:40], e))
sys.stdout.write(json.dumps(out))
"""


def _cleavage_probs_via_subprocess(sequences, human_only=False, threshold=0.5,
                                   timeout=600):
    """Score ``sequences`` in a fresh Python subprocess and return
    ``{sequence: [per-position cleavage probabilities]}``.

    Subprocess isolation is the fix for issue #266 (macOS libomp
    duplication between torch and pandas). Returns ``{}`` and logs
    a warning on any subprocess-level failure (pepsickle missing,
    timeout, unparseable output) — never raises, so a vaxrank run
    continues without annotations.

    Per-source inference failures are silently dropped from the
    returned dict; the caller treats missing keys as "no signal."

    Parameters
    ----------
    sequences : iterable of str
        Source sequences to score. Deduplicated before sending.
    human_only, threshold : forwarded to ``mhctools.Pepsickle``.
    timeout : seconds before SIGTERM. Default 600 — pepsickle is
        ~1s/source on CPU, so even thousands of sources finish well
        within this cap.
    """
    import json
    import subprocess
    import sys

    seqs = sorted({s for s in sequences if s})
    if not seqs:
        return {}
    payload = json.dumps({
        'sequences': seqs,
        'human_only': bool(human_only),
        'threshold': float(threshold),
    })
    try:
        result = subprocess.run(
            [sys.executable, '-c', _SUBPROCESS_SCRIPT],
            input=payload, capture_output=True, text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        logger.warning(
            "Pepsickle subprocess timed out after %ds (%d sources); "
            "skipping proteasome-cleavage annotation.",
            timeout, len(seqs))
        return {}
    except Exception as e:
        logger.warning(
            "Failed to launch pepsickle subprocess (%s); skipping "
            "proteasome-cleavage annotation.", e)
        logger.debug("Subprocess launch traceback:", exc_info=True)
        return {}
    if result.returncode != 0:
        first_err = (result.stderr or "").strip().split("\n")[0]
        logger.warning(
            "Pepsickle subprocess exited with code %d (%s); skipping "
            "proteasome-cleavage annotation. Common causes: pepsickle "
            "not installed (`pip install pepsickle`), or torch model "
            "load failure.", result.returncode, first_err[:200])
        return {}
    try:
        return json.loads(result.stdout)
    except json.JSONDecodeError as e:
        logger.warning(
            "Pepsickle subprocess produced unparseable output (%s); "
            "skipping. First 200 chars of stdout: %r",
            e, (result.stdout or "")[:200])
        return {}


# ----------------------------------------------------------------------------
# Public entry point
# ----------------------------------------------------------------------------

def annotate_processing(predictions, predictor=None,
                        human_only=False, threshold=0.5):
    """Attach pepsickle credibility scores to each EpitopePrediction
    in place.

    Parameters
    ----------
    predictions : iterable of EpitopePrediction
        Mutated in place: each gets ``pepsickle_c_term_cleavage_prob``
        / ``pepsickle_max_internal_cut_prob`` /
        ``pepsickle_processing_score`` populated when annotation
        succeeds. Predictions with no usable source sequence are
        passed through unchanged.
    predictor : optional, object with a ``cleavage_probs(sequence) ->
        list[float]`` method
        Test seam. When None (default), pepsickle runs in a fresh
        subprocess (see :func:`_cleavage_probs_via_subprocess`) to
        avoid the macOS libomp clash. Tests inject an in-process
        stub.
    human_only, threshold : forwarded to pepsickle on the subprocess
        path. Ignored when ``predictor`` is supplied (caller is
        responsible for configuring their own predictor).

    Returns
    -------
    int
        Number of predictions successfully annotated.
    """
    predictions_list = list(predictions)
    if not predictions_list:
        return 0

    # Group by source so we score each unique sequence exactly once.
    by_source = {}
    for p in predictions_list:
        source = getattr(p, 'source_sequence', None) or ''
        if not source:
            continue
        by_source.setdefault(source, []).append(p)
    if not by_source:
        return 0

    # Resolve cleavage probabilities for every unique source.
    if predictor is None:
        probs_by_source = _cleavage_probs_via_subprocess(
            by_source.keys(),
            human_only=human_only, threshold=threshold)
    else:
        probs_by_source = {}
        for source in by_source:
            try:
                probs_by_source[source] = predictor.cleavage_probs(source)
            except Exception as e:
                logger.warning(
                    "Predictor failed for source (len=%d, prefix=%r): "
                    "%s — skipping its %d predictions.",
                    len(source), source[:20], e, len(by_source[source]))

    # Apply per-peptide annotations using the precomputed arrays.
    n_annotated = 0
    for source, preds in by_source.items():
        seq_probs = probs_by_source.get(source)
        if not seq_probs or len(seq_probs) < len(source):
            continue
        for p in preds:
            peptide = p.peptide_sequence or ''
            if not peptide:
                continue
            offset = _resolve_peptide_offset(source, peptide, p)
            if offset is None:
                continue
            c_term, max_internal = _per_position_processing(
                seq_probs, offset, len(peptide))
            if c_term is None:
                continue
            p.pepsickle_c_term_cleavage_prob = c_term
            p.pepsickle_max_internal_cut_prob = max_internal
            p.pepsickle_processing_score = _composite_processing_score(
                c_term, max_internal)
            n_annotated += 1

    if n_annotated:
        logger.info(
            "Pepsickle credibility tagging: annotated %d / %d "
            "EpitopePrediction(s) across %d unique source sequence(s).",
            n_annotated, len(predictions_list), len(by_source))
    return n_annotated


def _resolve_peptide_offset(source, peptide, prediction):
    """Locate the peptide's offset within its source.

    Trust the prediction's declared ``offset`` first; re-locate via
    closest-substring search when the declared offset doesn't match.
    Warn when re-location moves the offset by more than the
    drift threshold (3aa absolute, 5% of source length, whichever is
    larger) — that's a sign the upstream loader's offset accounting
    is broken, not a 0/1-based indexing typo.

    Returns the offset (int) or None when the peptide isn't in the source.
    """
    declared = getattr(prediction, 'offset', None) or 0
    length = len(peptide)
    if (declared + length <= len(source)
            and source[declared:declared + length] == peptide):
        return declared
    found = _closest_occurrence(source, peptide, declared)
    if found is None:
        return None
    drift = abs(found - declared)
    threshold = max(3, int(len(source) * 0.05))
    if drift > threshold:
        logger.warning(
            "Peptide %r re-located in source (len=%d) by %d positions "
            "(declared=%d, found=%d, threshold=%d). Annotation proceeds "
            "with the located position; check the upstream loader's "
            "offset accounting.",
            peptide, len(source), drift, declared, found, threshold)
    return found
