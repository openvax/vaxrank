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

mhcflurry-presentation already includes an antigen-processing prior
conditioned on flanking residues, so a per-k-mer "presentation score"
implicitly captures whether *this specific 9-mer* would be cleaved
out and presented. What it can't tell us is whether the proteasome
would cut *inside* the ligand and destroy it before MHC, or whether
the C-terminal cut at the ligand's boundary is clean.

Pepsickle (Weeder et al., Bioinformatics 2021, doi:10.1093/bioinformatics/btab628)
predicts per-position cleavage probabilities. We use it to
**credibility-tag** existing MHC ligand predictions: each
``EpitopePrediction`` gains three optional fields exposing the per-
position pepsickle signal at the ligand's location in its source
sequence. The annotations don't change ranking by default — they're
extra context for clinical review and downstream filtering.

Issue: openvax/vaxrank#249.
"""

import logging

logger = logging.getLogger(__name__)


def _per_position_processing(seq_probs, start, length):
    """Compute (c_term, max_internal) cleavage probabilities for a
    peptide at offset ``start`` within ``seq_probs``.

    seq_probs is the per-position cleavage probability array for the
    *source sequence* — pepsickle returns one float per residue,
    representing "would the proteasome cut C-terminally to this
    residue?". For a peptide spanning ``[start, start+length)``:

      - C-terminal cut probability = seq_probs[start + length - 1]
      - Max internal cut probability = max of seq_probs[start ..
        start+length-2] (positions strictly inside the peptide;
        last position is the C-term cut, not an internal one)

    Returns ``(c_term, max_internal)`` as floats in [0, 1].
    """
    if length < 1 or start < 0 or start + length > len(seq_probs):
        # Out-of-range request — return a "no signal" tuple.
        return None, None
    c_term = float(seq_probs[start + length - 1])
    if length >= 2:
        internal_window = seq_probs[start:start + length - 1]
        max_internal = float(max(internal_window))
    else:
        # Single-residue peptide has no "inside" — degenerate but
        # don't crash.
        max_internal = 0.0
    return c_term, max_internal


def _composite_processing_score(c_term, max_internal):
    """Composite credibility: ``c_term * (1 - max_internal)``.

    Range [0, 1]: 1.0 = ideal release (probable clean cut at C-term,
    no internal-cut risk); 0.0 = either no clean C-term cut or a
    near-certain internal cut destroys the ligand.

    Returns ``None`` when either input is None (not annotated).
    """
    if c_term is None or max_internal is None:
        return None
    return float(c_term) * (1.0 - float(max_internal))


def _load_default_predictor():
    """Load pepsickle on demand. Returns None if unavailable so the
    caller can degrade gracefully instead of failing the whole run.
    """
    try:
        from mhctools import Pepsickle
    except Exception as e:
        logger.warning(
            "Could not import mhctools.Pepsickle (%s); skipping "
            "proteasome-cleavage credibility tagging. Install "
            "pepsickle (`pip install pepsickle`) to enable.", e)
        return None
    try:
        return Pepsickle()
    except Exception as e:
        logger.warning(
            "Could not instantiate Pepsickle (%s); skipping "
            "proteasome-cleavage credibility tagging.", e)
        return None


def annotate_processing(predictions, predictor=None):
    """Attach pepsickle credibility scores to each EpitopePrediction.

    Groups predictions by ``source_sequence`` and runs the proteasome
    predictor *once per unique source*; the per-position cleavage
    array is then sliced for each peptide based on its ``offset`` /
    ``peptide_sequence`` length.

    Predictions that lack a usable source sequence (empty string,
    too short to contain the peptide, or peptide not findable at the
    declared offset) are passed through unchanged — the credibility
    fields stay None.

    Parameters
    ----------
    predictions : iterable of EpitopePrediction
        Mutated in place; each gets ``c_term_cleavage_prob`` /
        ``max_internal_cut_prob`` / ``processing_score`` set when the
        annotation succeeds.
    predictor : object with a ``cleavage_probs(sequence) -> list[float]``
        method, or ``None``. When ``None``, ``mhctools.Pepsickle`` is
        loaded lazily; if pepsickle isn't installed we log a warning
        and return without modifying any prediction.

    Returns
    -------
    int
        Number of predictions successfully annotated.
    """
    predictions_list = list(predictions)
    if not predictions_list:
        return 0

    if predictor is None:
        predictor = _load_default_predictor()
        if predictor is None:
            return 0

    # Group by source sequence so we run the predictor once per unique
    # SLP / fragment / construct, not once per peptide.
    by_source = {}
    for p in predictions_list:
        source = getattr(p, 'source_sequence', None) or ''
        if not source:
            continue
        by_source.setdefault(source, []).append(p)

    if not by_source:
        return 0

    n_annotated = 0
    for source, preds in by_source.items():
        try:
            seq_probs = predictor.cleavage_probs(source)
        except Exception as e:
            logger.warning(
                "pepsickle prediction failed for source sequence "
                "(len=%d, prefix=%r); skipping its %d predictions. "
                "Error: %s", len(source), source[:20], len(preds), e)
            continue
        if not seq_probs or len(seq_probs) < len(source):
            logger.debug(
                "pepsickle returned %d probs for source of length %d; "
                "shape mismatch — skipping.",
                len(seq_probs) if seq_probs else 0, len(source))
            continue
        for p in preds:
            peptide = p.peptide_sequence or ''
            if not peptide:
                continue
            # Trust the declared offset first; fall back to a string
            # search in the source. Either path can be wrong if the
            # source was truncated relative to the original protein,
            # so re-locate when the declared offset doesn't match.
            offset = getattr(p, 'offset', None) or 0
            length = len(peptide)
            if (offset + length > len(source) or
                    source[offset:offset + length] != peptide):
                # Re-find by substring search; keeps the annotation
                # robust against off-by-one offset drift.
                idx = source.find(peptide)
                if idx < 0:
                    continue
                offset = idx
            c_term, max_internal = _per_position_processing(
                seq_probs, offset, length)
            if c_term is None:
                continue
            p.c_term_cleavage_prob = c_term
            p.max_internal_cut_prob = max_internal
            p.processing_score = _composite_processing_score(
                c_term, max_internal)
            n_annotated += 1

    if n_annotated:
        logger.info(
            "Pepsickle credibility tagging: annotated %d / %d "
            "EpitopePrediction(s) (%d unique source sequences).",
            n_annotated, len(predictions_list), len(by_source))
    return n_annotated
