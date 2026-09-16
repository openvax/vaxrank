# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Tests for the small helpers that live in
:mod:`vaxrank.cli.entry_point` — the per-junction MHC resolver, the
grouped args summary printer, and the LENS provenance-marker
plumbing.

Note: importing :mod:`vaxrank.cli.entry_point` transitively triggers
isovar's ``logging.config.fileConfig`` which wipes pytest's caplog
handlers from the root logger. So these tests attach a local
``logging.Handler`` to the entry_point logger directly rather than
relying on the ``caplog`` fixture. The behavior we're pinning is
identical; only the capture mechanism differs.
"""

import logging
from types import SimpleNamespace
from contextlib import contextmanager

import pytest


@contextmanager
def _capture_logger(logger_name, level=logging.DEBUG):
    """Yield a list that fills with LogRecord objects emitted by
    ``logger_name`` (and its children, by propagation) at or above
    ``level``. Cleans up the handler on exit."""
    target = logging.getLogger(logger_name)
    records = []

    class _Capture(logging.Handler):
        def emit(self, r):
            records.append(r)

    handler = _Capture()
    handler.setLevel(level)
    prev_level = target.level
    target.setLevel(level)
    target.addHandler(handler)
    try:
        yield records
    finally:
        target.removeHandler(handler)
        target.setLevel(prev_level)


# -- resolve_mhc_for_linker_optimizer ----------------------------------


def _args_with_inferred_alleles(alleles=None):
    """Minimal Namespace shape that the LENS path actually produces.

    Real LENS-mode args don't have ``--mhc-predictor`` /
    ``--mhc-alleles`` (the external arg parser doesn't add them), so
    the mhctools helpers raise ``AttributeError`` when called against
    such a Namespace. The optimizer's job is to fall back to the
    inferred-alleles path.
    """
    return SimpleNamespace(
        _inferred_mhc_alleles_from_lens=list(alleles or []))


def test_resolve_mhc_defaults_to_mhcflurry_when_predictor_missing(monkeypatch):
    """LENS-path: args lacks ``--mhc-predictor`` but the loader
    stashed inferred alleles. Vaxrank should fall back to mhcflurry
    as a credible default rather than refusing to optimize.

    We don't actually instantiate mhcflurry here — loading it in
    the test process triggers the macOS libomp clash that issue
    #266 fixed for pepsickle (mhcflurry has the same constraint).
    Instead, we patch ``mhctools.MHCflurry`` to a sentinel and
    verify the function calls it with the inferred alleles + emits
    the "defaulting to mhcflurry" INFO line."""
    import mhctools
    from vaxrank.cli import entry_point as ep

    captured = {}

    class _SentinelPredictor:
        pass

    def _stub_ctor(alleles, **kw):
        captured['alleles'] = alleles
        captured['kwargs'] = kw
        return _SentinelPredictor()

    monkeypatch.setattr(mhctools, 'MHCflurry', _stub_ctor)
    args = _args_with_inferred_alleles(['HLA-A*02:01', 'HLA-B*07:02'])
    with _capture_logger('vaxrank.cli.entry_point') as records:
        predictor, alleles = ep.resolve_mhc_for_linker_optimizer(args)
    assert isinstance(predictor, _SentinelPredictor)
    assert alleles == ['HLA-A*02:01', 'HLA-B*07:02']
    assert captured['alleles'] == ['HLA-A*02:01', 'HLA-B*07:02']
    msgs = [r.getMessage() for r in records]
    assert any('inferred from the LENS / pVACseq report' in m for m in msgs)
    assert any('defaulting to mhcflurry' in m for m in msgs), \
        "Expected the 'defaulting to mhcflurry' INFO line; got %r" % msgs


def test_resolve_mhc_warns_when_mhcflurry_default_unavailable(monkeypatch):
    """When the mhcflurry default can't load (not installed, weights
    missing, …), the function returns ``(None, None)`` and surfaces
    the targeted "predictor missing, alleles available" warning so
    the operator knows what to fix.

    Patches ``mhctools.MHCflurry`` itself to raise — replacing the
    sys.modules entry would crater anything else that imports
    mhctools concurrently.
    """
    import mhctools
    from vaxrank.cli import entry_point as ep

    def _broken_ctor(*_a, **_kw):
        raise OSError("simulated: mhcflurry weights not available")
    monkeypatch.setattr(mhctools, 'MHCflurry', _broken_ctor)

    args = _args_with_inferred_alleles(['HLA-A*02:01'])
    with _capture_logger('vaxrank.cli.entry_point') as records:
        predictor, alleles = ep.resolve_mhc_for_linker_optimizer(args)
    assert predictor is None
    assert alleles is None
    msgs = [r.getMessage() for r in records]
    assert any(
        'alleles are available' in m and 'no MHC predictor' in m
        for m in msgs), \
        "Expected the 'predictor missing, alleles available' WARNING; got %r" % msgs


def test_resolve_mhc_no_inputs_at_all():
    """Pipeline path with neither flag set: both come back None and
    the warning explicitly mentions both flags."""
    from vaxrank.cli.entry_point import resolve_mhc_for_linker_optimizer
    args = SimpleNamespace()  # no inferred alleles either
    with _capture_logger('vaxrank.cli.entry_point') as records:
        predictor, alleles = resolve_mhc_for_linker_optimizer(args)
    assert predictor is None
    assert alleles is None
    msg = ' '.join(r.getMessage() for r in records)
    assert '--mhc-alleles' in msg
    assert '--mhc-predictor' in msg


def test_resolve_mhc_propagates_real_predictor_load_failure(monkeypatch):
    """A genuine model-load error inside mhctools (e.g. a missing
    weights file — surfaces as ``RuntimeError`` / ``OSError``) must
    propagate, not be swallowed into "optimizer disabled". Pin the
    narrowed exception scope so a future broadening of the catch
    falls red here."""
    from vaxrank.cli import entry_point as ep

    class _BadDay(RuntimeError):
        pass

    def _kaboom(_args):
        raise _BadDay("predictor went sideways")

    monkeypatch.setattr(ep, 'mhc_binding_predictor_from_args', _kaboom)
    args = SimpleNamespace(_inferred_mhc_alleles_from_lens=['HLA-A*02:01'])
    with pytest.raises(_BadDay):
        ep.resolve_mhc_for_linker_optimizer(args)


def test_resolve_target_alleles_uses_cli_then_external_fallback(monkeypatch):
    from vaxrank.cli import entry_point as ep

    monkeypatch.setattr(
        ep, 'mhc_alleles_from_args', lambda args: ['HLA-A*02:01'])
    args = SimpleNamespace(
        _inferred_mhc_alleles_from_lens=['HLA-B*07:02'])
    assert ep.resolve_target_alleles(args) == ['HLA-A*02:01']

    def unavailable(args):
        raise ValueError("no CLI alleles")

    monkeypatch.setattr(ep, 'mhc_alleles_from_args', unavailable)
    assert ep.resolve_target_alleles(args) == ['HLA-B*07:02']
    assert ep.resolve_target_alleles(SimpleNamespace()) == []


# -- log_args_summary ---------------------------------------------------


def _basic_args(**overrides):
    """Minimal args shape with a parser_defaults snapshot. The
    summary printer reads ``args._parser_defaults`` to decide what
    counts as default — so ``parser_defaults`` is captured *before*
    overrides are applied (mirroring how the real arg-parser
    snapshot works in :func:`vaxrank.cli.arg_parser.parse_vaxrank_args`).
    """
    base = dict(
        vcf=None, bam=None, input_lens=None, input_pvacseq=None,
        ensembl_release=None, genome=None, tumor_sample_name=None,
        input_json_file=None,
        mhc_predictor=None, mhc_alleles=None, mhc_alleles_file=None,
        mhc_epitope_lengths=None,
        vaccine_type=['peptide', 'mrna'],
        vaccine_peptide_length=25, padding_around_mutation=5,
        antigen_content=None, epitopes_per_antigen=None,
        peptide_mode='slp', peptide_linker='G4Sx3',
        mrna_signal_peptide='HLA_B', mrna_linker='G4Sx2',
        mrna_max_constructs=2, mrna_antigens_per_construct=5,
        output_dir='', output_csv='', output_xlsx_report='',
        output_ascii_report='', output_html_report='',
        output_pdf_report='', output_neoepitope_report='',
        output_json_file='', output_passing_variants_csv='',
        output_isovar_csv='', output_epitopes='',
        output_patient_id='', output_reviewed_by='',
        output_final_review='', pdf_backend='pdfkit',
        log_path='python.log',
        manufacturability=None, wt_epitopes=True,
        cosmic_vcf_filename='', max_mutations_in_report=None,
        config=None, config_set_overrides=None,
        config_expr_overrides=None,
        verbose=False,
        processing_aware_annotation=True,
        pepsickle_human_only=False, pepsickle_threshold=0.5,
    )
    parser_defaults = dict(base)  # snapshot before overrides
    base.update(overrides)
    args = SimpleNamespace(**base)
    args._parser_defaults = parser_defaults
    return args


def test_log_args_summary_smoke_no_user_overrides():
    """No CLI overrides → every value matches its default. Without
    ``--verbose`` we should still emit the header but nothing else
    that names a value."""
    from vaxrank.cli.entry_point import log_args_summary
    args = _basic_args()
    with _capture_logger('vaxrank.cli.entry_point', logging.INFO) as records:
        log_args_summary(args)
    msgs = [r.getMessage() for r in records
            if 'Vaxrank run configuration' in r.getMessage()]
    assert len(msgs) == 1, "Expected exactly one summary log call"
    msg = msgs[0]
    assert 'Vaxrank run configuration' in msg
    # No value lines should leak through when everything matches
    # defaults and verbose is off.
    assert 'mrna_linker' not in msg


def test_log_args_summary_surfaces_user_overrides():
    """A user-set value differs from parser default → it gets a
    line. Values that match defaults stay hidden."""
    from vaxrank.cli.entry_point import log_args_summary
    args = _basic_args(input_lens='/path/to/file.tsv', vaccine_type=['mrna'])
    with _capture_logger('vaxrank.cli.entry_point', logging.INFO) as records:
        log_args_summary(args)
    msg = [r.getMessage() for r in records
           if 'Vaxrank run configuration' in r.getMessage()][0]
    assert 'input_lens' in msg
    assert "/path/to/file.tsv" in msg
    assert 'vaccine_type' in msg
    assert "['mrna']" in msg
    # A value that didn't change stays hidden (peptide_mode default).
    assert 'peptide_mode' not in msg


def test_log_args_summary_shows_auto_inferred_block():
    """LENS-path inferred state lives on underscore-prefixed attrs
    on args. ``log_args_summary`` lifts them into a separate
    ``[auto-inferred]`` section."""
    from vaxrank.cli.entry_point import log_args_summary
    args = _basic_args()
    args._inferred_mhc_alleles_from_lens = ['HLA-A*02:01']
    with _capture_logger('vaxrank.cli.entry_point', logging.INFO) as records:
        log_args_summary(args)
    msg = [r.getMessage() for r in records
           if 'Vaxrank run configuration' in r.getMessage()][0]
    assert 'auto-inferred' in msg
    assert '_inferred_mhc_alleles_from_lens' in msg
    assert 'HLA-A*02:01' in msg


def test_log_args_summary_verbose_shows_defaults():
    """``--verbose`` flips every key into the visible set, including
    those equal to the parser default. The marker ``(default)``
    annotates them."""
    from vaxrank.cli.entry_point import log_args_summary
    args = _basic_args(verbose=True)
    with _capture_logger('vaxrank.cli.entry_point', logging.INFO) as records:
        log_args_summary(args)
    msg = [r.getMessage() for r in records
           if 'Vaxrank run configuration' in r.getMessage()][0]
    # Defaults visible with ``(default)`` annotation.
    assert 'peptide_mode' in msg
    assert '(default)' in msg


# -- populate_default_output_paths -------------------------------------


def test_auto_populate_pipeline_path_fills_csv_and_json():
    """Pipeline run (no LENS / pVACseq input) with just
    ``--output-dir`` set should auto-fill canonical CSV + JSON
    paths inside the directory."""
    from vaxrank.cli.entry_point import populate_default_output_paths
    args = SimpleNamespace(
        output_dir='/tmp/run',
        input_lens=None, input_pvacseq=None,
        output_csv='', output_json_file='')
    populate_default_output_paths(args)
    assert args.output_csv == '/tmp/run/ranked_vaccine_peptides.csv'
    assert args.output_json_file == '/tmp/run/ranked_vaccine_peptides.json'


def test_auto_populate_lens_path_fills_neoepitope_csv():
    """LENS path's natural rank report is per-(peptide, allele) —
    the auto-fill picks ``neoepitope_predictions.csv`` instead of
    ``ranked_vaccine_peptides.csv`` so the filename matches the
    content. JSON dump is skipped (the LENS path doesn't build the
    rich in-memory result that --output-json-file serializes)."""
    from vaxrank.cli.entry_point import populate_default_output_paths
    args = SimpleNamespace(
        output_dir='/tmp/lens-run',
        input_lens='/path/to/lens.tsv', input_pvacseq=None,
        output_csv='', output_json_file='')
    populate_default_output_paths(args)
    assert args.output_csv == '/tmp/lens-run/neoepitope_predictions.csv'
    assert args.output_json_file == ''


def test_auto_populate_explicit_paths_win():
    """When the user explicitly passes ``--output-csv`` /
    ``--output-json-file``, the auto-fill must not overwrite them."""
    from vaxrank.cli.entry_point import populate_default_output_paths
    args = SimpleNamespace(
        output_dir='/tmp/run',
        input_lens=None, input_pvacseq=None,
        output_csv='/elsewhere/custom.csv',
        output_json_file='/elsewhere/custom.json')
    populate_default_output_paths(args)
    assert args.output_csv == '/elsewhere/custom.csv'
    assert args.output_json_file == '/elsewhere/custom.json'


def test_auto_populate_no_output_dir_is_a_noop():
    """Without ``--output-dir`` the helper does nothing — the
    operator-passed flags determine output destinations as before."""
    from vaxrank.cli.entry_point import populate_default_output_paths
    args = SimpleNamespace(
        output_dir='', input_lens=None, input_pvacseq=None,
        output_csv='', output_json_file='')
    populate_default_output_paths(args)
    assert args.output_csv == ''
    assert args.output_json_file == ''


# -- LENS antigen-source breakdown ordering ---------------------------


def test_lens_funnel_summarizes_input_composition():
    """The consolidated ``LENS load funnel`` line replaces the old
    scattered per-stage count lines. It must report the input row
    count + antigen_source breakdown and the candidate-epitope count
    in one block (so the operator isn't left correlating drop counts
    across multiple lines)."""
    import os
    from vaxrank.epitope_io import read_lens_report
    from vaxrank.external_input import lens_ranking_result
    from .testing_helpers import data_path

    lens_path = data_path("epitope_fixtures/lens_example.tsv")
    if not os.path.exists(lens_path):
        pytest.skip("lens_example.tsv fixture not found")

    _loaded = read_lens_report(lens_path)
    predictions = list(_loaded.epitopes)
    with _capture_logger('vaxrank.external_input', logging.INFO) as records:
        lens_ranking_result(_loaded, predictions)

    funnel = [r.getMessage() for r in records
              if 'LENS load funnel' in r.getMessage()]
    assert len(funnel) == 1, (
        "Expected exactly one consolidated funnel line; got %d" % len(funnel))
    msg = funnel[0]
    assert 'rows in:' in msg
    assert 'candidate epitopes eligible for construct ranking' in msg
    # The old per-stage lines must be gone (no double-counting).
    assert not any('LENS report contains' in r.getMessage() for r in records)
    assert not any(
        'lack gene_name' in r.getMessage() for r in records)


def test_lens_antigen_source_breakdown_orders_snv_indel_first():
    """The breakdown's stable order: SNV / INDEL up-front, then
    other kinds by descending count, ``(missing)`` last. Pin the
    ordering so a future maintainer can't accidentally swap to
    pure-alphabetical."""
    import os
    import re
    from vaxrank.epitope_io import read_lens_report
    from vaxrank.external_input import lens_ranking_result
    from .testing_helpers import data_path

    lens_path = data_path("epitope_fixtures/lens_example.tsv")
    if not os.path.exists(lens_path):
        pytest.skip("lens_example.tsv fixture not found")

    _loaded = read_lens_report(lens_path)
    predictions = list(_loaded.epitopes)
    with _capture_logger('vaxrank.external_input', logging.INFO) as records:
        lens_ranking_result(_loaded, predictions)

    msgs = [
        r.getMessage() for r in records
        if 'LENS load funnel' in r.getMessage()
    ]
    if not msgs:
        pytest.skip("Test fixture didn't trigger the funnel line")
    line = msgs[0]
    m = re.search(r'rows in:\s*([^\n]+)', line)
    assert m, "Couldn't parse breakdown segment from: %r" % line
    segments = [s.strip() for s in m.group(1).split(',')]
    kinds_in_order = [s.split('=')[0].strip() for s in segments]
    snv_idx = kinds_in_order.index('SNV') if 'SNV' in kinds_in_order else -1
    indel_idx = (
        kinds_in_order.index('INDEL') if 'INDEL' in kinds_in_order else -1)
    missing_idx = (
        kinds_in_order.index('(missing)')
        if '(missing)' in kinds_in_order else -1)
    if snv_idx >= 0:
        for i, k in enumerate(kinds_in_order):
            if k in ('SNV', 'INDEL'):
                continue
            assert snv_idx < i, "SNV must precede %r" % k
    if missing_idx >= 0:
        assert missing_idx == len(kinds_in_order) - 1, \
            "'(missing)' must land last; got order %r" % kinds_in_order
    if snv_idx >= 0 and indel_idx >= 0:
        assert snv_idx < indel_idx


# -- output-dir overwrite confirmation --------------------------------


def test_confirm_overwrite_noop_for_missing_or_empty_dir(tmp_path):
    """Missing or empty --output-dir needs no confirmation (no exit)."""
    from vaxrank.cli.entry_point import confirm_output_dir_overwrite
    # missing
    confirm_output_dir_overwrite(
        SimpleNamespace(output_dir=str(tmp_path / "nope"), force_overwrite=False))
    # empty
    (tmp_path / "empty").mkdir()
    confirm_output_dir_overwrite(
        SimpleNamespace(output_dir=str(tmp_path / "empty"), force_overwrite=False))


def test_confirm_overwrite_force_proceeds(tmp_path, caplog):
    """--force-overwrite proceeds into a non-empty dir without prompting."""
    from vaxrank.cli.entry_point import confirm_output_dir_overwrite
    (tmp_path / "f.txt").write_text("x")
    with caplog.at_level(logging.INFO):
        confirm_output_dir_overwrite(
            SimpleNamespace(output_dir=str(tmp_path), force_overwrite=True))
    assert any('force-overwrite' in r.getMessage() for r in caplog.records)


def test_confirm_overwrite_non_interactive_warns_and_proceeds(
        tmp_path, caplog, monkeypatch):
    """Non-interactive (no TTY) runs warn and proceed rather than hang."""
    from vaxrank.cli import entry_point
    (tmp_path / "f.txt").write_text("x")
    monkeypatch.setattr(entry_point.sys.stdin, 'isatty', lambda: False)
    with caplog.at_level(logging.WARNING):
        entry_point.confirm_output_dir_overwrite(
            SimpleNamespace(output_dir=str(tmp_path), force_overwrite=False))
    assert any('already exists' in r.getMessage() for r in caplog.records)


def test_write_run_summary(tmp_path):
    """The top-level run_summary.txt records inputs, the MHC alleles in
    play (flagged inferred on the external path), antigen counts, and
    where each output landed."""
    from vaxrank.cli.entry_point import write_run_summary
    args = SimpleNamespace(
        output_dir=str(tmp_path),
        input_lens='patient.lens.tsv', input_pvacseq=None, vcf=None, bam=None,
        _inferred_mhc_alleles_from_lens=['HLA-A*02:01', 'HLA-B*07:02'],
        mhc_alleles=None, mhc_alleles_file=None,
        output_csv=str(tmp_path / 'neoepitope_predictions.csv'),
        output_ascii_report=str(tmp_path / 'vaccine_report.txt'),
        output_pdf_report=str(tmp_path / 'vaccine_report.pdf'),
        vaccine_type=['peptide', 'mrna'])
    patient = SimpleNamespace(
        num_somatic_variants=5, num_coding_effect_variants=5,
        num_variants_with_rna_support=4, num_variants_with_vaccine_peptides=5)

    write_run_summary(args, patient, source='external')

    text = (tmp_path / 'run_summary.txt').read_text()
    assert 'LENS report' in text and 'patient.lens.tsv' in text
    assert 'inferred from report' in text
    assert 'HLA-A*02:01' in text
    assert 'variants with antigens:' in text and '5' in text
    # Multi-type + reports requested → split layout, so each modality
    # subdir holds constructs + its own report.
    assert 'peptide/' in text and 'mrna/' in text
    assert 'constructs + vaccine_report' in text


def test_write_run_summary_noop_without_output_dir(tmp_path):
    """No --output-dir → no run summary written (ranking-only runs)."""
    from vaxrank.cli.entry_point import write_run_summary
    args = SimpleNamespace(output_dir='', input_lens=None, input_pvacseq=None)
    write_run_summary(args, None, source='external')
    assert not list(tmp_path.iterdir())


# -- output paths inside a not-yet-created directory -------------------


def test_ensure_parent_dir_creates_nested_path_and_tolerates_bare_name(tmp_path):
    """Every writer routes through this helper, so it has to handle a
    nested path, a bare filename (no parent to create) and an existing
    directory without complaining. Creation is logged exactly once, so a
    mistyped path shows up in the run log instead of quietly producing a
    new tree."""
    from vaxrank.epitope_io import ensure_parent_dir
    nested = tmp_path / "run" / "reports" / "out.csv"
    with _capture_logger("vaxrank.epitope_io", logging.INFO) as records:
        ensure_parent_dir(str(nested))
        assert nested.parent.is_dir()
        ensure_parent_dir(str(nested))  # idempotent
        ensure_parent_dir("out.csv")  # no parent component
    created = [r for r in records if "Created output directory" in r.getMessage()]
    assert len(created) == 1
    assert str(nested.parent) in created[0].getMessage()


def test_template_reports_create_their_parent_directory(tmp_path):
    """``--output-dir DIR`` auto-populates report paths inside DIR, and
    the writers run before any construct writer has created it. Writing
    a report into a missing directory used to raise FileNotFoundError
    after the whole pipeline had already run."""
    from vaxrank.patient_info import PatientInfo
    from vaxrank.report import (
        TemplateDataCreator, make_ascii_report, make_html_report,
    )

    template_data = TemplateDataCreator(
        ranked_variants_with_vaccine_peptides=[],
        patient_info=PatientInfo("TEST"),
        final_review="",
        reviewers="",
        args_for_report={"manufacturability": False, "wt_epitopes": False},
        input_json_file=None,
    ).compute_template_data()

    out_dir = tmp_path / "does-not-exist-yet"
    make_ascii_report(template_data, str(out_dir / "vaccine_report.txt"))
    make_html_report(template_data, str(out_dir / "vaccine_report.html"))
    assert (out_dir / "vaccine_report.txt").exists()
    assert (out_dir / "vaccine_report.html").exists()


def _stub_epitope():
    from mhctools.pred import Prediction
    from vaxrank.candidate_epitope import CandidateEpitope, Peptide
    prediction = Prediction(
        kind='pMHC_affinity', predictor_name='test', predictor_version='',
        allele='HLA-A*02:01', peptide='SIINFEKL', value=100.0, score=0.0,
        percentile_rank=0.5)
    return CandidateEpitope.from_peptide(
        Peptide(sequence='SIINFEKL', source_sequence='SIINFEKL', offset=0,
                predictions=(prediction,)),
        comparators={}, overlaps_mutation=True, occurs_in_reference=False)


def _stub_vaccine_peptide():
    """Smallest object ``make_csv_report`` / ``make_minimal_neoepitope_report``
    will actually emit a row for — both skip peptides with no target epitopes
    and return without writing anything when every peptide is skipped, so a
    writer test needs one that survives the filter."""
    from vaxrank.manufacturability import ManufacturabilityScores
    fragment = SimpleNamespace(
        gene_name="GENE1", amino_acids="SIINFEKLSIINFEKL",
        mutant_amino_acid_start_offset=3, mutant_amino_acid_end_offset=4,
        n_alt_reads=10)
    return SimpleNamespace(
        mutant_protein_fragment=fragment,
        contains_target_epitopes=lambda: True,
        target_epitopes=[_stub_epitope()],
        combined_score=1.0, expression_score=1.0, target_epitope_score=1.0,
        manufacturability_scores=ManufacturabilityScores(
            *([0.0] * len(ManufacturabilityScores._fields))))


def _stub_variant():
    return SimpleNamespace(
        contig="1", start=1, ref="A", alt="G",
        original_start=1, original_ref="A", original_alt="G",
        short_description="chr1:1 A>G")


def _write_ranked_csv(path):
    from vaxrank.report import make_csv_report
    make_csv_report(
        [(_stub_variant(), [_stub_vaccine_peptide()])], csv_report_path=path)


def _write_ranked_xlsx(path):
    from vaxrank.report import make_csv_report
    make_csv_report(
        [(_stub_variant(), [_stub_vaccine_peptide()])], excel_report_path=path)


def _write_minimal_neoepitope_xlsx(path):
    from vaxrank.report import make_minimal_neoepitope_report
    make_minimal_neoepitope_report(
        ranked_variants_with_vaccine_peptides=[
            (_stub_variant(), [_stub_vaccine_peptide()])],
        num_epitopes_per_peptide=None, excel_report_path=path)


def _write_predictions_tsv(path):
    from vaxrank.epitope_io import save_predictions
    save_predictions([], path)


@pytest.mark.parametrize("filename,writer", [
    ("ranked_vaccine_peptides.csv", _write_ranked_csv),
    ("ranked_vaccine_peptides.xlsx", _write_ranked_xlsx),
    ("neoepitopes.xlsx", _write_minimal_neoepitope_xlsx),
    ("epitopes.tsv", _write_predictions_tsv),
])
def test_tabular_writers_create_their_parent_directory(tmp_path, filename, writer):
    """Every ``--output-*`` writer has to create its own parent, not just
    the template reports: ``--output-csv sub/out.csv`` used to fail with
    pandas' "Cannot save file into a non-existent directory" only after the
    whole pipeline had run. One case per writer so removing any single
    ``ensure_parent_dir`` call fails a test."""
    target = tmp_path / "missing" / "deeper" / filename
    writer(str(target))
    assert target.exists()


def _pipeline_args(tmp_path, monkeypatch, **overrides):
    """Parsed args plus stubs for everything upstream of the writers, so
    ``ranked_vaccine_peptides_with_metadata_from_parsed_args`` runs the
    output block without a genome, a BAM or a predictor."""
    from vaxrank.cli import entry_point
    from vaxrank.cli.arg_parser import parse_vaxrank_args

    class _EmptyVariantCollection(list):
        sources = []

    variants = _EmptyVariantCollection()
    results = SimpleNamespace(
        variant_counts=lambda: {
            'num_total_variants': 0, 'num_coding_effect_variants': 0,
            'num_variants_with_rna_support': 0,
            'num_variants_with_vaccine_peptides': 0},
        variant_properties=lambda **kwargs: [],
        ranked_vaccine_peptides=[])
    monkeypatch.setattr(entry_point, "mhc_alleles_from_args", lambda args: [])
    monkeypatch.setattr(entry_point, "variant_collection_from_args", lambda args: [])
    monkeypatch.setattr(entry_point, "filter_unannotatable_variants", lambda v: variants)
    monkeypatch.setattr(
        entry_point, "extract_dna_vaf_by_variant", lambda v, **kwargs: {})
    monkeypatch.setattr(
        entry_point, "run_vaxrank_from_parsed_args", lambda args: results)
    monkeypatch.setattr(entry_point, "GenePathwayCheck", lambda: None)

    args = parse_vaxrank_args([
        "--vcf", "unused.vcf", "--bam", "unused.bam",
        "--mhc-predictor", "random", "--mhc-alleles", "HLA-A*02:01",
        "--output-csv", str(tmp_path / "ranked.csv"),
    ])
    for key, value in overrides.items():
        setattr(args, key, value)
    return entry_point, args


def test_json_dump_creates_its_parent_directory(tmp_path, monkeypatch):
    """The reported crash: ``--output-json-file`` is written straight from
    ``open(path, 'w')`` before any construct writer has run, so a path in a
    missing directory raised FileNotFoundError after the full pipeline."""
    entry_point, args = _pipeline_args(
        tmp_path, monkeypatch,
        output_json_file=str(tmp_path / "js" / "ranked_vaccine_peptides.json"),
        output_passing_variants_csv=str(tmp_path / "pv" / "passing.csv"))

    entry_point.ranked_vaccine_peptides_with_metadata_from_parsed_args(args)

    assert (tmp_path / "js" / "ranked_vaccine_peptides.json").exists()
    assert (tmp_path / "pv" / "passing.csv").exists()


def test_output_dir_bundle_is_written_into_a_missing_directory(
        tmp_path, monkeypatch):
    """``--output-dir DIR`` is the documented way to get the full bundle and
    ``populate_default_output_paths`` fills in paths inside DIR. Nothing
    creates DIR before the JSON dump, so the whole flow used to die on a
    directory the operator hadn't pre-created."""
    entry_point, args = _pipeline_args(tmp_path, monkeypatch)
    out_dir = tmp_path / "brand" / "new" / "dir"
    args.output_dir = str(out_dir)
    args.output_csv = ''
    args.input_lens = None
    args.input_pvacseq = None
    entry_point.populate_default_output_paths(args)
    assert args.output_json_file == str(out_dir / "ranked_vaccine_peptides.json")

    entry_point.ranked_vaccine_peptides_with_metadata_from_parsed_args(args)

    assert (out_dir / "ranked_vaccine_peptides.json").exists()


def test_isovar_csv_creates_its_parent_directory(tmp_path, monkeypatch):
    """``--output-isovar-csv`` is written from ``run_vaxrank_from_parsed_args``,
    the earliest of the output writers."""
    from unittest.mock import Mock
    from vaxrank.cli import entry_point
    from vaxrank.cli.arg_parser import parse_vaxrank_args

    monkeypatch.setattr(entry_point, "variant_collection_from_args", lambda args: [])
    monkeypatch.setattr(entry_point, "filter_unannotatable_variants", lambda v: v)
    monkeypatch.setattr(entry_point, "alignment_file_from_args", lambda args: None)
    monkeypatch.setattr(entry_point, "read_collector_from_args", lambda args: None)
    monkeypatch.setattr(entry_point, "predictors_from_args", lambda args: [Mock()])
    monkeypatch.setattr(entry_point, "run_isovar", lambda **kwargs: [])
    monkeypatch.setattr(entry_point, "run_vaxrank", lambda **kwargs: kwargs)

    target = tmp_path / "iso" / "nested" / "isovar.csv"
    args = parse_vaxrank_args([
        "--vcf", "unused.vcf", "--bam", "unused.bam",
        "--mhc-predictor", "random", "--mhc-alleles", "HLA-A*02:01",
        "--output-isovar-csv", str(target),
    ])
    entry_point.run_vaxrank_from_parsed_args(args)

    assert target.exists()
