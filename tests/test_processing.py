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
- Annotation is purely additive — never mutates ranking, and (post-2.23,
  see #272 Phase B) never mutates the flat record objects either.
  Per-(peptide, source) records land on
  :class:`vaxrank.processing_prediction.ProcessingPrediction` and reach
  report writers via the ``processing_predictions_by_key`` map.
- Continuous (not binary) per-position cleavage scores.
- Re-locates peptides in the source sequence when declared offset is
  off by one (off-by-one drift in offsets is common in upstream
  loaders).
- Composite ``ProcessingPrediction.processing_score`` =
  sqrt(c_term × (1 − max_internal)) (geometric mean of the two factors).
- Single pepsickle pass per unique source sequence (not per peptide).
"""

from mhctools.pred import Prediction

from vaxrank.candidate_epitope import CandidateEpitope, Peptide
from vaxrank.processing import (
    _component_probs,
    annotate_processing,
)


def _ep(peptide, source, offset, ic50=100.0, allele="HLA-A*02:01",
        predictor="stub"):
    """Build a minimal ``CandidateEpitope`` for tests — one mutant
    pMHC_affinity ``Prediction`` plus the (peptide, source, offset)
    layout the processing annotator joins on."""
    pred = Prediction(
        kind='pMHC_affinity',
        predictor_name=predictor,
        predictor_version='',
        allele=allele,
        peptide=peptide,
        value=ic50,
        score=0.0,
        percentile_rank=ic50 / 100.0,
    )
    return CandidateEpitope(
        mutant=Peptide(
            sequence=peptide,
            source=source,
            offset=offset,
            predictions=(pred,),
        ),
        overlaps_mutation=True,
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


# ---- Local guard around mhctools helpers ---------------------------------
#
# The component extraction (``c_term_prob`` / ``max_internal_prob``) and the
# composite ``c_term * (1 - max(internal))`` formula live in mhctools and are
# tested there. The only vaxrank-side concern is the out-of-range guard, so
# that's all we exercise locally.

def test_component_probs_out_of_range_returns_none_tuple():
    """Peptide span extending past the source returns ``(None, None)``
    so the caller can drop the row instead of raising IndexError from
    the mhctools slice helpers."""
    c_term, max_internal = _component_probs([0.1, 0.2, 0.3], start=2, length=5)
    assert c_term is None
    assert max_internal is None


# ---- Annotation (integration with flat records) ---------------------

def test_annotate_processing_attaches_continuous_scores():
    """Each (peptide, source) pair lands in the returned
    ProcessingPrediction map with c_term / max_internal /
    processing as floats. Post-2.23 (#272 Phase B) the
    flat record is *not* mutated — readers consume the map."""
    source = "AAAAKLMNPVAAAA"  # 14 aa
    # Peptide KLMNPV at offset 4, length 6
    probs = [0.0, 0.0, 0.0, 0.05,   # source positions 0-3
             0.10, 0.20, 0.05, 0.50, 0.30, 0.85,   # peptide positions 4-9
             0.0, 0.0, 0.0, 0.0]
    pred = _ep("KLMNPV", source, offset=4)
    n, by_key = annotate_processing(
        [pred], predictor=StubPepsickle({source: probs}))
    assert n == 1
    pp = by_key[(pred.sequence, source, 'pepsickle')]
    # C-term = probs at index 9 = 0.85 (clean release)
    assert abs(pp.c_term_cleavage_prob - 0.85) < 1e-9
    # max internal = max of probs[4..8] = max(0.10, 0.20, 0.05, 0.50, 0.30) = 0.50
    assert abs(pp.max_internal_cut_prob - 0.50) < 1e-9
    # processing = sqrt(0.85 * (1 - 0.50)) = sqrt(0.425) ≈ 0.6519
    # (geometric mean of c_term and (1 - max_internal); see
    # ``vaxrank.processing`` for the rationale).
    import math
    assert abs(pp.processing_score - math.sqrt(0.85 * 0.50)) < 1e-9


def test_annotate_processing_returns_processing_prediction_map():
    """Post-2.22 ``annotate_processing`` returns a tuple of
    ``(n_annotated, processing_predictions_by_key)`` — the map is
    the canonical record (decoupled from flat record).
    Each entry is a ``ProcessingPrediction`` keyed on
    ``(peptide, source, predictor_name)``."""
    from vaxrank.processing_prediction import ProcessingPrediction
    source = "AAAAKLMNPVAAAA"
    probs = [0.0] * 4 + [0.10, 0.20, 0.05, 0.50, 0.30, 0.85] + [0.0] * 4
    pred = _ep("KLMNPV", source, offset=4)
    n, by_key = annotate_processing(
        [pred], predictor=StubPepsickle({source: probs}))
    assert n == 1
    # Map carries one ProcessingPrediction keyed on
    # (peptide, source, 'pepsickle').
    assert len(by_key) == 1
    key = ('KLMNPV', source, 'pepsickle')
    assert key in by_key
    pp = by_key[key]
    assert isinstance(pp, ProcessingPrediction)
    assert pp.peptide_sequence == 'KLMNPV'
    assert pp.source_sequence == source
    assert pp.predictor_name == 'pepsickle'
    assert abs(pp.c_term_cleavage_prob - 0.85) < 1e-9
    assert abs(pp.max_internal_cut_prob - 0.50) < 1e-9
    # Composite is the geometric mean of c_term and (1 - max_internal).
    import math
    assert abs(pp.processing_score - math.sqrt(0.85 * 0.50)) < 1e-9


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
    """Epitopes with empty source_sequence are passed through
    untouched (annotation requires a source to slice probs against)
    and don't land in the ProcessingPrediction map."""
    pred = _ep("KLMNPV", source="", offset=0)
    n, by_key = annotate_processing(
        [pred], predictor=StubPepsickle({}))
    assert n == 0
    assert by_key == {}


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
    n, by_key = annotate_processing(
        [pred], predictor=StubPepsickle({source: probs}))
    assert n == 1
    pp = by_key[(pred.sequence, source, 'pepsickle')]
    # C-term should pick up probs[8] = 0.95 (re-located, not probs[1+5]=0.05)
    assert abs(pp.c_term_cleavage_prob - 0.95) < 1e-9


def test_annotate_processing_does_not_touch_ranking_score():
    """The annotation must not alter the underlying epitope's
    ranking-driving fields — pepsickle's job is purely additive."""
    source = "AAAAKLMNPVAAAA"
    pred = _ep("KLMNPV", source, offset=4, ic50=42.0)
    leaf = pred.mutant.best_affinity()
    pre_ic50 = leaf.value
    pre_rank = leaf.percentile_rank
    n, by_key = annotate_processing(
        [pred],
        predictor=StubPepsickle({source: [0.1] * len(source)}))
    # ProcessingPrediction landed in the map (post-2.23, the
    # canonical record).
    assert n == 1
    assert (pred.sequence, source, 'pepsickle') in by_key
    # Ranking-driving fields untouched (frozen CandidateEpitope, can't be mutated).
    leaf_after = pred.mutant.best_affinity()
    assert leaf_after.value == pre_ic50
    assert leaf_after.percentile_rank == pre_rank


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
    n, by_key = annotate_processing(
        [pred_ok, pred_fail], predictor=FlakyPredictor())
    # Only the OK one annotated; the failing source skipped, no crash.
    assert n == 1
    assert (pred_ok.sequence, "AAAAAAAAAA", 'pepsickle') in by_key
    assert (pred_fail.sequence, "FFAILFFFFFF", 'pepsickle') not in by_key


def test_annotate_processing_empty_input_returns_zero():
    n, _ = annotate_processing([], predictor=StubPepsickle({}))
    assert n == 0


# ---- Default-predictor construction -------------------------------------
#
# The default-path now goes straight through
# ``mhctools.Pepsickle(isolate_subprocess=True)`` (mhctools 3.13.3+).
# Vaxrank no longer ships its own subprocess launcher — those tests
# moved upstream alongside the implementation.

def test_load_default_predictor_returns_none_when_mhctools_missing(
        monkeypatch):
    """When mhctools isn't installed, ``_load_default_predictor`` logs
    a clear warning and returns None so the caller can degrade
    gracefully (no crash, no annotations)."""
    import builtins
    import sys
    from vaxrank.processing import _load_default_predictor

    real_import = builtins.__import__

    def _fail_on_mhctools(name, *args, **kwargs):
        if name == 'mhctools' or name.startswith('mhctools.'):
            raise ImportError("simulated: no module named 'mhctools'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, '__import__', _fail_on_mhctools)
    monkeypatch.delitem(sys.modules, 'mhctools', raising=False)
    assert _load_default_predictor() is None


# ---- Report integration --------------------------------------------------

def test_ascii_report_surfaces_each_predictor_in_per_epitope_table():
    """End-to-end render: when a VaccinePeptide carries two
    flat records for the same (peptide, allele) — one per
    predictor — both rows appear in the rendered ASCII table with
    a ``Predictor`` column distinguishing them. Pre-2.24 the
    VCF/BAM path collapsed via dict-overwrite so this test would
    have only seen one row."""
    import io
    from vaxrank.report import _make_report
    template_data = {
        'patient_info': {'Patient ID': 'TEST'},
        'package_versions': {},
        'reviewers': [], 'final_review': '',
        'input_json_file': None,
        'include_manufacturability': False,
        'include_wt_epitopes': False,
        'args': [], 'mrna_ranking': None,
        'variants': [{
            'num': 1,
            'short_description': 'chr1:100 A>T',
            'variant_data': {'Gene name': 'GENEA', 'Top score': 1.0},
            'effect_data': {'effect': 'missense'},
            'peptides': [{
                'header_display_data': {
                    'num': 1,
                    'aa_before_mutation': 'KLQ',
                    'aa_mutant': 'X',
                    'aa_after_mutation': 'GHS'},
                'peptide_data': {},
                'manufacturability_data': {},
                'wt_epitopes': [], 'ascii_wt_epitopes': None,
                'epitopes': [
                    {'Sequence': 'SIINFEKL', 'Predictor': 'mhcflurry',
                     'IC50': '50.00 nM', 'Allele': 'A*02:01'},
                    {'Sequence': 'SIINFEKL', 'Predictor': 'netmhcpan',
                     'IC50': '75.00 nM', 'Allele': 'A*02:01'},
                ],
                'ascii_epitopes': (
                    "Sequence Predictor IC50 Allele\n"
                    "SIINFEKL mhcflurry 50.00 A*02:01\n"
                    "SIINFEKL netmhcpan 75.00 A*02:01\n"),
            }],
        }],
    }
    f = io.StringIO()
    _make_report(template_data, f, 'templates/template.txt')
    rendered = f.getvalue()
    # Both predictors named in the per-epitope ASCII table; same
    # peptide+allele appears twice (one row per predictor).
    assert 'mhcflurry' in rendered
    assert 'netmhcpan' in rendered
    assert rendered.count('SIINFEKL') >= 2


def test_epitope_data_renders_one_row_per_predictor_for_same_pep_allele():
    """Issue #261: when LENS / pVACseq / multi-predictor VCF runs
    emit multiple Predictions per (peptide, allele) — one per
    predictor — the per-epitope report table renders one row
    per prediction with a ``Predictor`` column distinguishing them."""
    from collections import OrderedDict
    from vaxrank.report import TemplateDataCreator

    creator = TemplateDataCreator.__new__(TemplateDataCreator)
    creator.processing_predictions_by_key = {}

    def _pred(method, ic50):
        return Prediction(
            kind='pMHC_affinity', predictor_name=method,
            predictor_version='', allele='HLA-A*02:01',
            peptide='SIINFEKL', value=ic50, score=0.0,
            percentile_rank=0.5)
    mhcflurry, netmhcpan = _pred('mhcflurry', 50.0), _pred('netmhcpan', 75.0)
    epitope = CandidateEpitope(
        mutant=Peptide(
            sequence='SIINFEKL',
            source='SSIINFEKL', offset=1,
            predictions=(mhcflurry, netmhcpan)),
        comparators={'wt': Peptide(sequence='SIINFEKL')},
        overlaps_mutation=True, occurs_in_reference=False)

    rows = [creator._epitope_data(epitope, p) for p in (mhcflurry, netmhcpan)]
    assert len(rows) == 2
    # Each row preserves its own predictor name + IC50 — no collapsing.
    assert rows[0]['Predictor'] == 'mhcflurry'
    assert rows[1]['Predictor'] == 'netmhcpan'
    assert rows[0]['IC50'] == '50.00 nM'
    assert rows[1]['IC50'] == '75.00 nM'
    # Same allele + sequence (the join we want to surface, not collapse).
    assert rows[0]['Sequence'] == rows[1]['Sequence'] == 'SIINFEKL'
    assert rows[0]['Allele'] == rows[1]['Allele'] == 'A*02:01'
    # Result is an OrderedDict so the rendered HTML/PDF table
    # column order is stable across rows.
    assert all(isinstance(r, OrderedDict) for r in rows)
    assert list(rows[0].keys()) == list(rows[1].keys())


def test_epitope_data_surfaces_processing_columns_when_annotated():
    """``TemplateDataCreator._epitope_data`` adds three extra
    columns (Processing: C-term, Processing: max internal,
    Processing: combined) when ``include_processing=True`` — the
    values come from the ProcessingPrediction map by joining on
    ``(peptide, source, predictor_name)``. Default
    ``include_processing=False`` keeps the original 6-column shape
    so unannotated reports don't change."""
    from collections import OrderedDict

    source = "AAAAKLMNPVAAAA"
    pred = _ep("KLMNPV", source, offset=4)
    leaf = pred.mutant.best_affinity()

    from vaxrank.report import TemplateDataCreator
    creator = TemplateDataCreator.__new__(TemplateDataCreator)
    creator.processing_predictions_by_key = {}
    # Default: legacy 7-column shape (no Processing columns).
    pre = creator._epitope_data(pred, leaf)
    assert isinstance(pre, OrderedDict)
    assert 'Processing: C-term' not in pre
    assert 'Processing: combined' not in pre

    # After annotation, ``include_processing=True`` surfaces the
    # columns — read from the ProcessingPrediction map (#272 Phase B).
    n, by_key = annotate_processing(
        [pred],
        predictor=StubPepsickle({source: [0.1] * 9 + [0.85] + [0.0] * 4}))
    creator.processing_predictions_by_key = by_key
    post = creator._epitope_data(pred, leaf, include_processing=True)
    assert 'Processing: C-term' in post
    assert 'Processing: max internal' in post
    assert 'Processing: combined' in post
    # And the values are real (not '—' placeholders) since we have
    # a ProcessingPrediction for this peptide.
    assert post['Processing: C-term'] != '—'


# ---- entry_point integration --------------------------------------------

def test_manufacturability_default_on_for_peptide_only(monkeypatch):
    """Default manufacturability rendering depends on whether
    'peptide' is in the active --vaccine-type set: ON for runs that
    design peptides (synthesis-relevant metrics), OFF for mrna-only
    (those metrics don't apply to in-vivo translation)."""
    from vaxrank.cli.arg_parser import parse_vaxrank_args
    # Resolve mimics what main() does at startup. Test by parsing
    # args + invoking the same resolver block.
    cases = [
        (['peptide'], True),
        (['mrna'], False),
        (['peptide', 'mrna'], True),
        (['mrna', 'peptide'], True),
    ]
    for vaccine_type, expected in cases:
        args = parse_vaxrank_args([
            '--input-lens', '/dev/null',
            '--output-csv', '/dev/null',
            '--vaccine-type', *vaccine_type,
        ])
        # main() resolves args.manufacturability based on
        # _resolve_vaccine_types(args). Run that resolution here.
        from vaxrank.cli.entry_point import _resolve_vaccine_types
        if args.manufacturability is None:
            args.manufacturability = (
                'peptide' in _resolve_vaccine_types(args))
        assert args.manufacturability is expected, (
            "vaccine_type=%s: expected manufacturability=%s, got %s"
            % (vaccine_type, expected, args.manufacturability))


def test_manufacturability_explicit_flag_overrides_default():
    """User-supplied --include-manufacturability-in-report /
    --no-manufacturability-in-report wins over the modality-based
    default in either direction."""
    from vaxrank.cli.arg_parser import parse_vaxrank_args
    args = parse_vaxrank_args([
        '--input-lens', '/dev/null', '--output-csv', '/dev/null',
        '--vaccine-type', 'mrna',
        '--include-manufacturability-in-report',
    ])
    # explicit True survives the resolver
    assert args.manufacturability is True

    args = parse_vaxrank_args([
        '--input-lens', '/dev/null', '--output-csv', '/dev/null',
        '--vaccine-type', 'peptide',
        '--no-manufacturability-in-report',
    ])
    # explicit False survives the resolver
    assert args.manufacturability is False


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
    leaf_a = annotated.mutant.best_affinity()
    leaf_u = unannotated.mutant.best_affinity()
    n, by_key = annotate_processing(
        [annotated],
        predictor=StubPepsickle(
            {source: [0.1] * 9 + [0.85] + [0.0] * 4}))
    # Verify mixed state: only one of the two has a record in the map.
    assert (annotated.sequence, source, 'pepsickle') in by_key
    assert (unannotated.sequence, source, 'pepsickle') not in by_key
    creator.processing_predictions_by_key = by_key

    # When the caller turns include_processing on, BOTH rows have
    # the same key set — table renders cleanly.
    row_a = creator._epitope_data(annotated, leaf_a, include_processing=True)
    row_u = creator._epitope_data(unannotated, leaf_u, include_processing=True)
    assert list(row_a.keys()) == list(row_u.keys()), (
        "Annotated and unannotated rows should share identical "
        "column keys when include_processing=True; got "
        "annotated=%s vs unannotated=%s" % (
            list(row_a.keys()), list(row_u.keys())))
    # Unannotated row uses the '—' placeholder
    assert row_u['Processing: C-term'] == '—'
    assert row_u['Processing: max internal'] == '—'
    assert row_u['Processing: combined'] == '—'


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
    n, by_key = annotate_processing(
        [pred], predictor=StubPepsickle({source: probs}))
    pp = by_key[(pred.sequence, source, 'pepsickle')]
    # If re-location snapped to position 0 (first occurrence), c_term
    # would be probs[4] = 0.0; if it correctly snapped to position 8,
    # c_term = probs[12] = 0.95.
    assert pp.c_term_cleavage_prob == 0.95, (
        "Re-location should pick the closest occurrence to declared "
        "offset (8 → 12), not the first occurrence (0 → 4); got "
        "c_term=%s" % pp.c_term_cleavage_prob)


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


def test_warns_when_peptide_not_in_pep_context(caplog):
    """LENS rows occasionally carry a ``peptide`` that isn't a
    substring of its ``pep_context`` (peptide and pep_context built
    from different isoforms / annotation snapshots — see Pt02
    analysis, 8/2113 rows). Vaxrank skips those rows rather than
    fabricate an offset, and emits a single aggregate WARN at the
    end of ``annotate_processing`` with one example so the user can
    file an upstream LENS bug."""
    import logging
    # ``KLMXPV`` is one residue off from ``KLMNPV`` so it can't be
    # located in the source by substring or close-match search.
    source = "AAAAKLMNPVAAAA"
    pred = _ep("KLMXPV", source, offset=4)
    with caplog.at_level(logging.WARNING):
        n, _ = annotate_processing(
            [pred],
            predictor=StubPepsickle({source: [0.5] * len(source)}))
    assert n == 0
    skipped = [r for r in caplog.records
               if 'not a substring' in r.message]
    assert len(skipped) == 1, (
        "Expected exactly one aggregate skip-WARN; got %d" % len(skipped))
    assert "'KLMXPV'" in skipped[0].message
    assert source in skipped[0].message


def test_dedup_by_content_when_duplicate_objects(monkeypatch):
    """The CLI annotation dispatcher dedups by both id() AND a
    content key (peptide, allele, source, offset). A future loader
    that copies a flat record into a VP would produce two
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

    # Two distinct flat record objects with identical content.
    source = "AAAAKLMNPVAAAA"
    pred_a = _ep("KLMNPV", source, offset=4)
    pred_b = _ep("KLMNPV", source, offset=4)  # different object, same content
    assert id(pred_a) != id(pred_b)

    fragment = SimpleNamespace(amino_acids=source, gene_name='G')
    vp = SimpleNamespace(
        mutant_protein_fragment=fragment,
        mutant_epitopes=[pred_a])
    ranked = [(SimpleNamespace(), [vp])]
    lens_epitopes = [pred_b]

    captured = {}
    stub = StubPepsickle({source: [0.1] * 14})

    def _capture(pred_list, predictor=None, human_only=False, threshold=0.5):
        captured['list'] = list(pred_list)
        return real_annotate(pred_list, predictor=stub)

    monkeypatch.setattr(proc_mod, 'annotate_processing', _capture)
    _annotate_predictions_with_processing(ranked, lens_epitopes)

    assert 'list' in captured
    deduped = captured['list']
    assert len(deduped) == 1, (
        "Content-key dedup should collapse two distinct objects with "
        "identical (peptide, source, offset) to one; got %d"
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
