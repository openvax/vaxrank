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
  pepsickle_processing_score        composite ``sqrt(c_term *
                                    (1 - max_internal))`` — the
                                    geometric mean of the two factors;
                                    1.0 = ideal release, 0.0 = no
                                    clean release OR near-certain
                                    destruction. Geometric mean
                                    rather than the raw product so
                                    a balanced (0.6, 0.6) row scores
                                    ~0.6 instead of 0.36.

The annotations are purely additive — vaccine ranking is unaffected.
Reports surface the three columns when at least one prediction in
the per-VaccinePeptide list has been annotated.

Component extractors come from ``mhctools.processing_predictor``
(``ProcessingPredictor.c_term_prob`` / ``max_internal_prob``), so
the slicing semantics live in one place. The composite-score
formula is local: vaxrank uses the geometric mean of the two
factors rather than the raw product mhctools defaults to (which
is also offered there as ``score_cterm_anti_max_internal``).
Geometric mean is a gentler penalty when both factors are
mid-range, and matches how downstream readers interpret the
score as a "credibility tag" rather than a strict joint
probability.

Pepsickle (Weeder et al., Bioinformatics 2021) inference runs in an
isolated subprocess to avoid the macOS libomp clash between torch and
the parent's pandas/numpy/pyarrow OpenMP runtime (issue #266). As of
mhctools 3.13.3 (openvax/mhctools#201) that isolation is built in:
``Pepsickle(isolate_subprocess=True)`` does the right thing,
batching unique sequences per call. Vaxrank used to ship its own
launcher + subprocess module here; now it just sets the kwarg.

Issue: openvax/vaxrank#249.
"""

import logging

from mhctools.processing_predictor import ProcessingPredictor

logger = logging.getLogger(__name__)


# ----------------------------------------------------------------------------
# Per-peptide arithmetic on a precomputed per-position cleavage array
# ----------------------------------------------------------------------------
#
# The component extractors and the composite ``c_term * (1 - max_internal)``
# formula live in ``mhctools.processing_predictor``; we just import them.
# Locally we only need to defend against out-of-range slicing (peptide
# extending past the source) — ``ProcessingPredictor.c_term_prob`` would
# IndexError, so we range-check first and return ``(None, None)`` to mean
# "no signal" upstream.

def _component_probs(seq_probs, start, length):
    """Return ``(c_term, max_internal)`` for a peptide at
    ``seq_probs[start:start+length]``, or ``(None, None)`` when the
    span doesn't fit. Both components come straight from the public
    mhctools helpers — we just guard the array bounds.
    """
    if length < 1 or start < 0 or start + length > len(seq_probs):
        return None, None
    c_term = float(ProcessingPredictor.c_term_prob(seq_probs, start, length))
    max_internal = float(
        ProcessingPredictor.max_internal_prob(seq_probs, start, length))
    return c_term, max_internal


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
# Pepsickle invocation (delegates to mhctools' built-in subprocess isolation)
# ----------------------------------------------------------------------------

def _load_default_predictor(human_only=False, threshold=0.5):
    """Construct an ``mhctools.Pepsickle`` with subprocess isolation
    on, or return ``None`` if mhctools / pepsickle isn't installed.

    ``isolate_subprocess=True`` (mhctools 3.13.3+) routes
    ``cleavage_probs`` calls through a fresh interpreter so torch's
    libomp doesn't clash with the parent's pandas / numpy / pyarrow
    libomp on macOS (issue #266).
    """
    try:
        from mhctools import Pepsickle
    except ImportError as e:
        logger.warning(
            "mhctools.Pepsickle is not importable (%s); skipping "
            "proteasome-cleavage annotation. Install with `pip install "
            "mhctools[pepsickle]` or pass --no-processing-aware-annotation.",
            e)
        return None
    try:
        return Pepsickle(
            human_only=bool(human_only),
            threshold=float(threshold),
            isolate_subprocess=True)
    except Exception as e:
        logger.warning(
            "Failed to construct mhctools.Pepsickle (%s); skipping "
            "proteasome-cleavage annotation.", e)
        logger.debug("Pepsickle construction traceback:", exc_info=True)
        return None


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
        Test seam. When None (default), constructs
        ``mhctools.Pepsickle(isolate_subprocess=True)`` so torch's
        libomp doesn't clash with pandas / numpy on macOS. Tests
        inject an in-process stub to avoid the subprocess startup
        cost.
    human_only, threshold : forwarded to ``mhctools.Pepsickle`` when
        constructing the default predictor. Ignored when ``predictor``
        is supplied (caller is responsible for configuring their own
        predictor).

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

    # Resolve the predictor, then ask it for per-position cleavage
    # probabilities on each unique source sequence. The default path
    # uses mhctools' built-in subprocess isolation; the test-seam
    # path uses whatever object the caller passed in.
    if predictor is None:
        predictor = _load_default_predictor(
            human_only=human_only, threshold=threshold)
        if predictor is None:
            return 0

    # Batch every unique source into ONE predictor call when the
    # predictor exposes ``cleavage_probs_many`` (mhctools 3.13.3+).
    # The per-source path used to spawn a fresh subprocess per
    # sequence — ~1-2s startup × N sources turned a 30s run into 15
    # minutes on Pt02. Stub predictors in tests only implement
    # ``cleavage_probs``; fall back to the per-source loop for those.
    sources = list(by_source)
    probs_by_source = {}
    cleavage_probs_many = getattr(predictor, 'cleavage_probs_many', None)
    if cleavage_probs_many is not None:
        logger.info(
            "Pepsickle: scoring %d unique source sequence(s) in a "
            "single isolated-subprocess batch...", len(sources))
        try:
            probs_by_source = dict(cleavage_probs_many(sources))
        except Exception as e:
            logger.warning(
                "Batch predictor call failed (%s); skipping "
                "proteasome-cleavage annotation. Run with "
                "--no-processing-aware-annotation to suppress.", e)
            logger.debug("Batch predictor traceback:", exc_info=True)
            probs_by_source = {}
    else:
        for source in sources:
            try:
                probs_by_source[source] = predictor.cleavage_probs(source)
            except Exception as e:
                logger.warning(
                    "Predictor failed for source (len=%d, prefix=%r): "
                    "%s — skipping its %d predictions.",
                    len(source), source[:20], e, len(by_source[source]))

    # Apply per-peptide annotations using the precomputed arrays.
    n_annotated = 0
    # Track peptides that aren't substrings of their pep_context source
    # — likely a LENS upstream issue where the ``peptide`` column was
    # built from a different isoform / annotation than ``pep_context``
    # (see Pt02 analysis: the peptide column matches the canonical
    # mutant translation, but pep_context carries an extra residue
    # change that doesn't exist in the canonical reference protein).
    # Skip such rows rather than fabricate an offset; aggregate the
    # count and emit a single WARN at the end with one example for
    # upstream-bug filing.
    n_skipped_peptide_not_in_context = 0
    skipped_examples = []
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
                n_skipped_peptide_not_in_context += 1
                if len(skipped_examples) < 3:
                    skipped_examples.append((peptide, source))
                continue
            c_term, max_internal = _component_probs(
                seq_probs, offset, len(peptide))
            if c_term is None:
                continue
            p.pepsickle_c_term_cleavage_prob = c_term
            p.pepsickle_max_internal_cut_prob = max_internal
            # Composite score is the **geometric mean** of the two
            # factors — ``sqrt(c_term * (1 - max_internal))``. Range
            # [0, 1]; 1 = ideal release, 0 = no clean release OR
            # near-certain internal destruction. Geometric mean
            # penalizes both factors symmetrically and is more
            # forgiving than the raw product when one factor is
            # mid-range, which better matches the "credibility tag"
            # reading: a peptide with c_term=0.6 and (1 -
            # max_internal)=0.6 should score ~0.6, not 0.36.
            import math
            anti_max = max(0.0, 1.0 - max_internal)
            p.pepsickle_processing_score = math.sqrt(c_term * anti_max)
            n_annotated += 1

    if n_annotated:
        logger.info(
            "Pepsickle credibility tagging: annotated %d / %d "
            "EpitopePrediction(s) across %d unique source sequence(s).",
            n_annotated, len(predictions_list), len(by_source))
    if n_skipped_peptide_not_in_context:
        ex_peptide, ex_source = skipped_examples[0]
        logger.warning(
            "Pepsickle credibility tagging: skipped %d / %d "
            "EpitopePrediction(s) because the peptide is not a "
            "substring of its pep_context source — peptide and "
            "pep_context were built from different isoforms / "
            "annotation snapshots. Example: peptide=%r not found "
            "in pep_context=%r.",
            n_skipped_peptide_not_in_context, len(predictions_list),
            ex_peptide, ex_source)
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
