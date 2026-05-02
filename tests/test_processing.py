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

"""Pepsickle proteasome-cleavage credibility tagging (issue #249).

Pins:
- Annotation is purely additive — never mutates ranking.
- Continuous (not binary) per-position scores attach to each
  EpitopePrediction.
- Re-locates peptides in the source sequence when declared offset is
  off by one (off-by-one drift in offsets is common in upstream
  loaders).
- Composite ``pepsickle_processing_score`` = c_term × (1 − max_internal).
- Single pepsickle pass per unique source sequence (not per peptide).
"""

from vaxrank.epitope_prediction import EpitopePrediction
from vaxrank.processing import (
    _composite_processing_score,
    _per_position_processing,
    annotate_processing,
)


def _ep(peptide, source, offset, ic50=100.0, allele="HLA-A*02:01"):
    """Build a minimal EpitopePrediction for tests."""
    return EpitopePrediction(
        allele=allele,
        peptide_sequence=peptide,
        wt_peptide_sequence="",
        ic50=ic50,
        wt_ic50=10000.0,
        percentile_rank=ic50 / 100.0,
        prediction_method_name="stub",
        overlaps_mutation=True,
        source_sequence=source,
        offset=offset,
        occurs_in_reference=False,
    )


class StubPepsickle:
    """In-memory stub returning configurable per-position cleavage
    probabilities. Tracks how many times ``cleavage_probs`` is called
    so tests can pin "one pass per unique source"."""
    def __init__(self, probs_by_source):
        self.probs_by_source = probs_by_source
        self.call_count = 0
        self.calls = []

    def cleavage_probs(self, sequence):
        self.call_count += 1
        self.calls.append(sequence)
        return self.probs_by_source.get(sequence, [0.0] * len(sequence))


# ---- Helpers (unit-level) -------------------------------------------------

def test_per_position_processing_basic():
    # Source: 10-mer; peptide at offset 2, length 5 (positions 2..6).
    # Last position of the peptide is index 6.
    # Internal positions are 2..5 (4 positions strictly inside).
    seq_probs = [0.01, 0.05, 0.10, 0.50, 0.90, 0.20, 0.95, 0.30, 0.40, 0.05]
    c_term, max_internal = _per_position_processing(
        seq_probs, start=2, length=5)
    assert c_term == 0.95
    assert max_internal == 0.90  # peak of [0.10, 0.50, 0.90, 0.20]


def test_per_position_processing_out_of_range():
    seq_probs = [0.1, 0.2, 0.3]
    # Peptide extends past source — return None tuple.
    c_term, max_internal = _per_position_processing(
        seq_probs, start=2, length=5)
    assert c_term is None
    assert max_internal is None


def test_per_position_processing_single_residue():
    """Degenerate: a 1-residue peptide has no internal positions.
    Returns max_internal=0.0 (no signal) and c_term from the single
    position."""
    c_term, max_internal = _per_position_processing(
        [0.1, 0.7, 0.3], start=1, length=1)
    assert c_term == 0.7
    assert max_internal == 0.0


def test_composite_processing_score():
    # c_term=1.0, max_internal=0.0 → 1.0 (ideal)
    assert _composite_processing_score(1.0, 0.0) == 1.0
    # c_term=0.9, max_internal=0.5 → 0.45
    assert abs(_composite_processing_score(0.9, 0.5) - 0.45) < 1e-9
    # c_term=0.0 → 0 regardless of internal
    assert _composite_processing_score(0.0, 0.0) == 0.0
    # max_internal=1.0 → 0 (ligand always destroyed)
    assert _composite_processing_score(0.9, 1.0) == 0.0
    # None inputs → None
    assert _composite_processing_score(None, 0.5) is None
    assert _composite_processing_score(0.5, None) is None


# ---- Annotation (integration with EpitopePrediction) ---------------------

def test_annotate_processing_attaches_continuous_scores():
    """Each EpitopePrediction gets c_term / max_internal / processing
    set as floats. None of them should be left as the default (None)
    for a successfully-annotated prediction."""
    source = "AAAAKLMNPVAAAA"  # 14 aa
    # Peptide KLMNPV at offset 4, length 6
    probs = [0.0, 0.0, 0.0, 0.05,   # source positions 0-3
             0.10, 0.20, 0.05, 0.50, 0.30, 0.85,   # peptide positions 4-9
             0.0, 0.0, 0.0, 0.0]
    pred = _ep("KLMNPV", source, offset=4)
    n = annotate_processing(
        [pred], predictor=StubPepsickle({source: probs}))
    assert n == 1
    # C-term = probs at index 9 = 0.85 (clean release)
    assert abs(pred.pepsickle_c_term_cleavage_prob - 0.85) < 1e-9
    # max internal = max of probs[4..8] = max(0.10, 0.20, 0.05, 0.50, 0.30) = 0.50
    assert abs(pred.pepsickle_max_internal_cut_prob - 0.50) < 1e-9
    # processing = 0.85 * (1 - 0.50) = 0.425
    assert abs(pred.pepsickle_processing_score - 0.425) < 1e-9


def test_annotate_processing_one_pass_per_unique_source():
    """Critical performance pin: the predictor runs once per unique
    source_sequence, not once per peptide."""
    source_a = "AAAAAAAAAAAAA"
    source_b = "BBBBBBBBBBBBB"
    preds = [
        _ep("AAAAA", source_a, 2),
        _ep("AAAAAA", source_a, 3),  # same source
        _ep("AAAAAA", source_a, 1),  # same source again
        _ep("BBBBB", source_b, 0),
    ]
    stub = StubPepsickle({
        source_a: [0.1] * len(source_a),
        source_b: [0.2] * len(source_b),
    })
    annotate_processing(preds, predictor=stub)
    assert stub.call_count == 2, (
        "Expected one pepsickle call per unique source sequence "
        "(2 sources × 1 call); got %d." % stub.call_count)
    assert sorted(stub.calls) == sorted([source_a, source_b])


def test_annotate_processing_skips_predictions_without_source():
    """Predictions with empty source_sequence are passed through
    untouched (annotation requires a source to slice probs against)."""
    pred = _ep("KLMNPV", source="", offset=0)
    n = annotate_processing(
        [pred], predictor=StubPepsickle({}))
    assert n == 0
    assert pred.pepsickle_c_term_cleavage_prob is None


def test_annotate_processing_relocates_peptide_when_offset_off():
    """Real-world LENS / pVACseq loaders sometimes record the offset
    at the wrong index (off-by-one). The annotator should re-find the
    peptide in the source by substring search rather than failing
    silently."""
    source = "PADKLMNPVKVS"
    # Declared offset is wrong (says 1, peptide is actually at index 3)
    pred = _ep("KLMNPV", source, offset=1)
    probs = [0.0, 0.0, 0.0,   # 0-2
             0.10, 0.20, 0.05, 0.50, 0.30, 0.95,   # peptide at 3-8
             0.0, 0.0, 0.0]
    n = annotate_processing(
        [pred], predictor=StubPepsickle({source: probs}))
    assert n == 1
    # C-term should pick up probs[8] = 0.95 (re-located, not probs[1+5]=0.05)
    assert abs(pred.pepsickle_c_term_cleavage_prob - 0.95) < 1e-9


def test_annotate_processing_does_not_touch_ranking_score():
    """The annotation must not alter the prediction's existing scoring
    fields — ranking is unchanged whether or not pepsickle ran."""
    source = "AAAAKLMNPVAAAA"
    pred = _ep("KLMNPV", source, offset=4, ic50=42.0)
    pre_ic50 = pred.ic50
    pre_rank = pred.percentile_rank
    pre_logistic = pred.logistic_epitope_score()
    annotate_processing(
        [pred],
        predictor=StubPepsickle({source: [0.1] * len(source)}))
    # Annotation fields are now set
    assert pred.pepsickle_c_term_cleavage_prob is not None
    # Ranking-driving fields untouched
    assert pred.ic50 == pre_ic50
    assert pred.percentile_rank == pre_rank
    assert pred.logistic_epitope_score() == pre_logistic


def test_annotate_processing_predictor_failure_degrades_gracefully():
    """A pepsickle exception on one source shouldn't kill annotation
    for predictions from other sources."""
    class FlakyPredictor:
        def __init__(self):
            self.call_count = 0
        def cleavage_probs(self, sequence):
            self.call_count += 1
            if "FAIL" in sequence:
                raise RuntimeError("synthetic failure")
            return [0.5] * len(sequence)
    pred_ok = _ep("AAAAA", "AAAAAAAAAA", offset=2)
    pred_fail = _ep("FFFFF", "FFAILFFFFFF", offset=0)
    n = annotate_processing(
        [pred_ok, pred_fail], predictor=FlakyPredictor())
    # Only the OK one annotated; the failing source skipped, no crash.
    assert n == 1
    assert pred_ok.pepsickle_c_term_cleavage_prob is not None
    assert pred_fail.pepsickle_c_term_cleavage_prob is None


def test_annotate_processing_empty_input_returns_zero():
    n = annotate_processing([], predictor=StubPepsickle({}))
    assert n == 0


# ---- Report integration --------------------------------------------------

def test_epitope_data_surfaces_processing_columns_when_annotated():
    """report.TemplateDataCreator._epitope_data adds three extra
    columns (Pepsickle C-term cut, Pepsickle max internal cut,
    Pepsickle processing score) when ``include_processing=True``.
    Default ``include_processing=False`` keeps the original 6-column
    shape so unannotated reports don't change. Predictor name is in
    the column header so a future per-position predictor (NetChop, …)
    can land alongside without ambiguity."""
    from collections import OrderedDict

    source = "AAAAKLMNPVAAAA"
    pred = _ep("KLMNPV", source, offset=4)

    from vaxrank.report import TemplateDataCreator
    creator = TemplateDataCreator.__new__(TemplateDataCreator)
    # Default: 6-column legacy shape.
    pre = creator._epitope_data(pred)
    assert isinstance(pre, OrderedDict)
    assert 'Pepsickle C-term cut' not in pre
    assert 'Pepsickle processing score' not in pre

    # After annotation, ``include_processing=True`` surfaces the columns.
    annotate_processing(
        [pred],
        predictor=StubPepsickle({source: [0.1] * 9 + [0.85] + [0.0] * 4}))
    post = creator._epitope_data(pred, include_processing=True)
    assert 'Pepsickle C-term cut' in post
    assert 'Pepsickle max internal cut' in post
    assert 'Pepsickle processing score' in post


# ---- entry_point integration --------------------------------------------

def test_processing_aware_annotation_default_on():
    """On by default — pepsickle runs in an isolated subprocess so the
    macOS libomp clash with pandas/numpy doesn't crash the parent."""
    from vaxrank.cli.arg_parser import parse_vaxrank_args
    args = parse_vaxrank_args(['--input-lens', '/dev/null'])
    assert args.processing_aware_annotation is True


def test_processing_aware_annotation_opt_out():
    from vaxrank.cli.arg_parser import parse_vaxrank_args
    args = parse_vaxrank_args([
        '--input-lens', '/dev/null',
        '--no-processing-aware-annotation',
    ])
    assert args.processing_aware_annotation is False


def test_pepsickle_cli_param_passthrough():
    """--pepsickle-human-only and --pepsickle-threshold flags get
    threaded through to the predictor constructor."""
    from vaxrank.cli.arg_parser import parse_vaxrank_args
    args = parse_vaxrank_args([
        '--input-lens', '/dev/null',
        '--pepsickle-human-only',
        '--pepsickle-threshold', '0.7',
    ])
    assert args.pepsickle_human_only is True
    assert args.pepsickle_threshold == 0.7


# ---- Review-fix coverage: report header consistency, peptide re-location

def test_epitope_data_header_consistent_when_some_predictions_unannotated():
    """Mixed annotated/unannotated predictions in the same VP epitope
    list must produce consistent column keys when ``include_processing``
    is True. Pre-fix: header (from first epitope) had 6 keys but later
    annotated rows had 9 — table malformed. Now: caller passes
    ``include_processing=True`` for any list with ANY annotated
    prediction; every row gets all 9 keys with '—' for unannotated."""
    from vaxrank.report import TemplateDataCreator

    creator = TemplateDataCreator.__new__(TemplateDataCreator)
    source = "AAAAKLMNPVAAAA"
    annotated = _ep("KLMNPV", source, offset=4)
    unannotated = _ep("AAAAA", source, offset=0)
    annotate_processing(
        [annotated],
        predictor=StubPepsickle(
            {source: [0.1] * 9 + [0.85] + [0.0] * 4}))
    # Verify mixed state: only one of the two has the field set
    assert annotated.pepsickle_c_term_cleavage_prob is not None
    assert unannotated.pepsickle_c_term_cleavage_prob is None

    # When the caller turns include_processing on, BOTH rows have
    # the same key set — table renders cleanly.
    row_a = creator._epitope_data(annotated, include_processing=True)
    row_u = creator._epitope_data(unannotated, include_processing=True)
    assert list(row_a.keys()) == list(row_u.keys()), (
        "Annotated and unannotated rows should share identical "
        "column keys when include_processing=True; got "
        "annotated=%s vs unannotated=%s" % (
            list(row_a.keys()), list(row_u.keys())))
    # Unannotated row uses the '—' placeholder
    assert row_u['Pepsickle C-term cut'] == '—'
    assert row_u['Pepsickle max internal cut'] == '—'
    assert row_u['Pepsickle processing score'] == '—'


def test_re_location_picks_closest_to_declared_offset():
    """When the peptide appears more than once in the source (a
    repeated motif), re-location should pick the occurrence closest
    to the declared offset rather than always the first."""
    # Source has 'AAAAA' at offsets 0 and 8. Declared offset says
    # the peptide is the second occurrence; re-location should snap
    # to position 8, not position 0.
    source = "AAAAA" + "BBB" + "AAAAA"  # length 13
    pred = _ep("AAAAA", source, offset=7)  # off-by-one near the second occurrence
    # probs distinguish position 0 (low cleavage) from position 12 (high)
    probs = [0.0, 0.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0, 0.95]  # cleavage at last position only
    annotate_processing(
        [pred], predictor=StubPepsickle({source: probs}))
    # If re-location snapped to position 0 (first occurrence), c_term
    # would be probs[4] = 0.0; if it correctly snapped to position 8,
    # c_term = probs[12] = 0.95.
    assert pred.pepsickle_c_term_cleavage_prob == 0.95, (
        "Re-location should pick the closest occurrence to declared "
        "offset (8 → 12), not the first occurrence (0 → 4); got "
        "c_term=%s" % pred.pepsickle_c_term_cleavage_prob)


def test_re_location_warns_on_large_offset_drift(caplog):
    """If re-location moves the offset by more than 3 positions, the
    upstream loader's offset accounting is probably wrong — warn."""
    import logging
    source = "PADDING_PADDING_KLMNPV_PADDING"
    # Declared offset is way off — peptide is at index 16, not 2.
    pred = _ep("KLMNPV", source, offset=2)
    with caplog.at_level(logging.WARNING):
        annotate_processing(
            [pred],
            predictor=StubPepsickle({source: [0.5] * len(source)}))
    assert any(
        "re-located" in r.message and "positions" in r.message
        for r in caplog.records), \
        "Expected a 're-located by N positions' warning for >3aa drift"


def test_dedup_by_content_when_duplicate_objects(monkeypatch):
    """The CLI annotation dispatcher dedups by both id() AND a
    content key (peptide, allele, source, offset). A future loader
    that copies an EpitopePrediction into a VP would produce two
    distinct objects with the same content; only one should be
    annotated. Uses pytest's monkeypatch fixture so cleanup is
    automatic even if the test fails."""
    from types import SimpleNamespace

    from vaxrank.cli.entry_point import (
        _annotate_predictions_with_processing,
    )
    import vaxrank.processing as proc_mod

    # Capture the real function BEFORE patching so the wrapper can
    # call into it without recursing through the patched name.
    real_annotate = proc_mod.annotate_processing

    # Two distinct EpitopePrediction objects with identical content.
    source = "AAAAKLMNPVAAAA"
    pred_a = _ep("KLMNPV", source, offset=4)
    pred_b = _ep("KLMNPV", source, offset=4)  # different object, same content
    assert id(pred_a) != id(pred_b)

    fragment = SimpleNamespace(amino_acids=source, gene_name='G')
    vp = SimpleNamespace(
        mutant_protein_fragment=fragment,
        mutant_epitope_predictions=[pred_a])
    ranked = [(SimpleNamespace(), [vp])]
    lens_predictions = [pred_b]

    captured = {}
    stub = StubPepsickle({source: [0.1] * 14})

    def _capture(pred_list, predictor=None, human_only=False, threshold=0.5):
        captured['list'] = list(pred_list)
        return real_annotate(pred_list, predictor=stub)

    monkeypatch.setattr(proc_mod, 'annotate_processing', _capture)
    _annotate_predictions_with_processing(ranked, lens_predictions)

    assert 'list' in captured
    deduped = captured['list']
    assert len(deduped) == 1, (
        "Content-key dedup should collapse two distinct objects with "
        "identical (peptide, allele, source, offset) to one; got %d"
        % len(deduped))


def test_drift_threshold_scales_with_source_length(caplog):
    """For long sources, the absolute drift threshold (3) is too
    tight — 5% of length is more permissive. A 5aa drift on a
    1000-aa source should NOT warn."""
    import logging
    # 200-aa source; 5% = 10aa threshold. Drift = 8aa → no warning.
    long_source = ("A" * 50) + "KLMNPV" + ("A" * 144)
    pred = _ep("KLMNPV", long_source, offset=42)  # declared 42, real 50 → drift 8
    with caplog.at_level(logging.WARNING):
        annotate_processing(
            [pred],
            predictor=StubPepsickle({long_source: [0.1] * 200}))
    # 8aa drift on 200-aa source: 5% threshold = 10, so no warning.
    assert not any(
        "re-located" in r.message and "positions" in r.message
        for r in caplog.records), (
        "8aa drift on 200-aa source should be below 5%% threshold and "
        "not emit a warning")
