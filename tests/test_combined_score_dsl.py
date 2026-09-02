# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Tests for ``vaxrank.combined_score_dsl``.

Pins:
- The grammar (allowed nodes / functions; rejected constructs).
- The fragment-aware default and the three legacy pre-3.1 mode formulas,
  all expressed as DSL strings.
- That a user-supplied expression is the only mechanism — the default
  formula travels the same code path.
"""

import math
import pytest


def _make_vp(
        n_alt_reads=10,
        n_alt_fragments=None,
        n_overlapping_reads=20,
        target_epitope_score=2.5,
        n_alt_reads_supporting_protein_sequence=8,
        combined_score_expr=None,
):
    """Minimal ``VaccinePeptide`` for scoring tests. Most knobs are
    irrelevant to the score arithmetic; the helper hides the
    boilerplate."""
    from vaxrank.vaccine_peptide import VaccinePeptide
    from vaxrank.mutant_protein_fragment import MutantProteinFragment
    from varcode import Variant

    frag = MutantProteinFragment(
        variant=Variant(contig='1', start=100, ref='A', alt='G'),
        gene_name='TEST',
        amino_acids='MASKLVPENTM',
        mutant_amino_acid_start_offset=2,
        mutant_amino_acid_end_offset=5,
        supporting_reference_transcripts=[],
        n_overlapping_reads=n_overlapping_reads,
        n_alt_reads=n_alt_reads,
        n_ref_reads=max(0, n_overlapping_reads - n_alt_reads),
        n_alt_reads_supporting_protein_sequence=
            n_alt_reads_supporting_protein_sequence,
        n_alt_fragments=n_alt_fragments,
    )
    vp = VaccinePeptide(
        mutant_protein_fragment=frag,
        epitopes=[],
        combined_score_expr=combined_score_expr,
    )
    # Fake target_epitope_score directly — the real one comes from
    # summed epitope predictions but we want to control the scalar.
    vp.target_epitope_score = target_epitope_score
    return vp


def test_dsl_legacy_formulas_round_trip_as_expressions():
    """Each pre-3.1 ``combined_score_mode`` value has an exact DSL
    expression equivalent. Pinning these numeric values guards
    against regressions in either the bindings (``expression_score``,
    ``n_alt_reads``) or the arithmetic semantics."""
    n_alt = 16
    epitope = 3.0
    # On a reads-only source, the historical reads expression and the
    # fragment-aware default see the same evidence count.
    vp_default = _make_vp(
        n_alt_reads=n_alt, target_epitope_score=epitope,
        combined_score_expr=None)
    vp_explicit = _make_vp(
        n_alt_reads=n_alt, target_epitope_score=epitope,
        combined_score_expr='sqrt(n_alt_reads) * target_epitope_score')
    assert vp_default.combined_score == pytest.approx(vp_explicit.combined_score)
    # Numeric value: sqrt(16) * 3 = 12
    assert vp_default.combined_score == pytest.approx(12.0)

    vp_reads = _make_vp(
        n_alt_reads=n_alt, target_epitope_score=epitope,
        combined_score_expr='n_alt_reads * target_epitope_score')
    assert vp_reads.combined_score == pytest.approx(48.0)

    vp_only = _make_vp(
        n_alt_reads=n_alt, target_epitope_score=epitope,
        combined_score_expr='target_epitope_score')
    assert vp_only.combined_score == pytest.approx(3.0)


def test_default_expression_uses_portable_rna_evidence():
    """The default weights independent RNA evidence, not raw reads."""
    from vaxrank.vaccine_config import DEFAULT_COMBINED_SCORE_EXPR
    assert DEFAULT_COMBINED_SCORE_EXPR == (
        "sqrt(n_rna_alt) * target_epitope_score")
    vp = _make_vp(n_alt_reads=9, target_epitope_score=2.0,
                  combined_score_expr=None)
    # __post_init__ resolved the default
    assert vp.combined_score_expr == DEFAULT_COMBINED_SCORE_EXPR
    # And the score is sqrt(9) * 2 = 6
    assert vp.combined_score == pytest.approx(6.0)


def test_default_reorders_duplicate_heavy_locus_behind_clean_locus():
    """Fragments prevent paired reads and duplicates receiving extra weight.

    Both candidates have the same epitope score. The duplicate-heavy locus
    wins under explicit read weighting (40 > 25), while the default correctly
    ranks the clean locus first because it has more independent fragments
    (16 > 9). This is the between-variant reordering that settled #394.
    """
    duplicate_heavy = _make_vp(
        n_alt_reads=40,
        n_alt_fragments=9,
        target_epitope_score=1.0,
    )
    clean = _make_vp(
        n_alt_reads=25,
        n_alt_fragments=16,
        target_epitope_score=1.0,
    )

    assert duplicate_heavy.expression_score == pytest.approx(3.0)
    assert clean.expression_score == pytest.approx(4.0)
    assert clean.combined_score > duplicate_heavy.combined_score

    read_weighted = "sqrt(n_alt_reads) * target_epitope_score"
    duplicate_heavy_reads = _make_vp(
        n_alt_reads=40,
        n_alt_fragments=9,
        target_epitope_score=1.0,
        combined_score_expr=read_weighted,
    )
    clean_reads = _make_vp(
        n_alt_reads=25,
        n_alt_fragments=16,
        target_epitope_score=1.0,
        combined_score_expr=read_weighted,
    )
    assert duplicate_heavy_reads.combined_score > clean_reads.combined_score


def test_dsl_supports_math_functions():
    """sqrt / log / log1p / exp / min / max / abs / pow are all
    available; arithmetic + unary minus + power compose freely."""
    vp = _make_vp(n_alt_reads=10, target_epitope_score=1.0)
    cases = [
        ('sqrt(n_alt_reads)', math.sqrt(10)),
        ('log(n_alt_reads)', math.log(10)),
        ('log1p(n_alt_reads)', math.log1p(10)),
        ('exp(0)', 1.0),
        ('min(n_alt_reads, n_overlapping_reads)', 10.0),
        ('max(n_alt_reads, n_overlapping_reads)', 20.0),
        ('abs(-n_alt_reads)', 10.0),
        ('pow(n_alt_reads, 0.5)', math.sqrt(10)),
        ('n_alt_reads ** 2', 100.0),
        ('n_alt_reads / n_overlapping_reads', 0.5),
        ('-n_alt_reads + 100', 90.0),
        ('log1p(n_alt_reads) * target_epitope_score + 0.1 * '
         'n_alt_reads_supporting_protein_sequence',
         math.log1p(10) * 1.0 + 0.1 * 8),
    ]
    for expr, expected in cases:
        from vaxrank.combined_score_dsl import (
            parse_combined_score_expr, evaluate_combined_score)
        tree = parse_combined_score_expr(expr)
        got = evaluate_combined_score(tree, vp)
        assert got == pytest.approx(expected), (
            "%s: expected %r, got %r" % (expr, expected, got))


def test_dsl_rejects_disallowed_constructs():
    """Anything outside the locked-down grammar fails at parse time
    so the operator gets the rejection at config-load, not at
    scoring time."""
    from vaxrank.combined_score_dsl import parse_combined_score_expr
    bad = [
        # Attribute access (could escape the namespace)
        'n_alt_reads.foo',
        # Subscript
        'n_alt_reads[0]',
        # Method call (function isn't a bare name)
        'n_alt_reads.bit_length()',
        # Unknown function
        'frobnicate(n_alt_reads)',
        # Comparison (nonsensical for a score)
        'n_alt_reads < 5',
        # Boolean op
        'n_alt_reads and target_epitope_score',
        # Statements
        'import os',
        '',
        '   ',
    ]
    for expr in bad:
        with pytest.raises(ValueError):
            parse_combined_score_expr(expr)


def test_dsl_unknown_identifier_raises_at_eval_time():
    """Unknown identifier survives parsing (we don't know the
    namespace at parse time) but fails cleanly at eval with a
    message that names the offending expression."""
    from vaxrank.combined_score_dsl import (
        parse_combined_score_expr, evaluate_combined_score)
    vp = _make_vp()
    tree = parse_combined_score_expr('not_a_real_field * 2')
    with pytest.raises(ValueError, match="not_a_real_field"):
        evaluate_combined_score(tree, vp)


def test_combined_score_bindings_exposes_documented_scalar_namespace():
    from vaxrank.combined_score_dsl import combined_score_bindings

    bindings = combined_score_bindings(_make_vp())

    assert bindings["target_epitope_score"] == pytest.approx(2.5)
    assert bindings["n_alt_reads"] == pytest.approx(10)
    assert bindings["expression_score"] == pytest.approx(math.sqrt(10))


def test_dsl_parsed_at_construction_time_not_scoring_time():
    """Bad ``combined_score_expr`` should raise when the
    VaccinePeptide is built (config-load surface), not later when
    ``combined_score`` is read. Pins that early-failure contract."""
    with pytest.raises(ValueError):
        _make_vp(combined_score_expr='import os')


def test_default_expr_is_pre_parsed_at_construction(monkeypatch):
    """When the user doesn't pass ``combined_score_expr``,
    ``__post_init__`` resolves to the default and pre-parses it. The
    AST is the SAME single artifact ``combined_score`` evaluates —
    no separate hardcoded Python branch for the default case."""
    vp = _make_vp(combined_score_expr=None)
    from vaxrank.vaccine_config import DEFAULT_COMBINED_SCORE_EXPR
    import vaxrank.vaccine_peptide as vaccine_peptide_module

    assert vp.combined_score_expr == DEFAULT_COMBINED_SCORE_EXPR

    def fail_on_reparse(_expr):
        raise AssertionError("unexpected reparse")

    monkeypatch.setattr(
        vaccine_peptide_module,
        "parse_combined_score_expr",
        fail_on_reparse,
    )
    assert vp.combined_score == pytest.approx(math.sqrt(10) * 2.5)


def test_dsl_parses_documented_default_yaml_examples():
    """The commented examples in ``vaxrank/config/default.yaml`` must
    parse and evaluate cleanly — drift between the documented
    snippets and the grammar is the easiest way for a config-knob
    doc to silently rot. Pin all four of the equivalences listed
    under ``vaccine_peptides.combined_score_expr`` in default.yaml,
    plus the documented ``log1p(...)`` example."""
    from vaxrank.combined_score_dsl import (
        parse_combined_score_expr, evaluate_combined_score)
    documented = [
        # The three legacy-mode equivalences.
        'target_epitope_score',
        'n_alt_reads * target_epitope_score',
        'expression_score * target_epitope_score',
        # The richer commented example.
        'log1p(n_alt_reads) * target_epitope_score',
    ]
    vp = _make_vp(
        n_alt_reads=10, target_epitope_score=2.0,
        n_alt_reads_supporting_protein_sequence=4)
    for expr in documented:
        tree = parse_combined_score_expr(expr)
        # Must produce a finite float — pinning that the bindings
        # used by the documented snippets are actually exposed by
        # the evaluator.
        result = evaluate_combined_score(tree, vp)
        assert math.isfinite(result), (
            "Documented example %r evaluated to non-finite: %r"
            % (expr, result))


def test_dsl_eval_error_truncates_bindings_in_user_message():
    """Evaluation-time errors should surface a *preview* of the
    bindings (the three headline scores), not the full dict — the
    full dict is logged at DEBUG. Keeps multi-VP failure stacks
    readable."""
    from vaxrank.combined_score_dsl import (
        parse_combined_score_expr, evaluate_combined_score)
    vp = _make_vp()
    tree = parse_combined_score_expr('log(-1)')  # math domain error
    with pytest.raises(ValueError) as excinfo:
        evaluate_combined_score(tree, vp)
    msg = str(excinfo.value)
    # Preview header is present; "more)" suffix covers the
    # truncated remainder.
    assert 'preview' in msg.lower()
    assert 'more)' in msg
    # And the message stays compact: under ~400 chars even though
    # there are 8 bindings of varying width.
    assert len(msg) < 400, "error message too long: %r" % msg
    # Opinionated preview: the three headline scores must be in
    # the message. Read counts and mutant-AA count must NOT —
    # those go to the DEBUG log.
    assert 'target_epitope_score' in msg
    assert 'self_epitope_score' in msg
    assert 'expression_score' in msg
    assert 'n_alt_reads' not in msg
    assert 'n_overlapping_reads' not in msg


def test_dsl_eval_error_warns_when_headline_binding_missing():
    """Soft contract: the binding extractor is supposed
    to populate every entry in ``HEADLINE_BINDINGS``. If a future
    refactor drops one, the user-facing preview shrinks silently —
    surface that through a warning so the maintainer notices."""
    import logging
    from vaxrank.combined_score_dsl import (
        HEADLINE_BINDINGS, combined_score_bindings, evaluate_combined_score,
        parse_combined_score_expr)
    import vaxrank.combined_score_dsl as dsl_mod

    # Monkeypatch the bindings extractor to drop one headline key.
    real_extractor = combined_score_bindings

    def _stripped(vp):
        b = real_extractor(vp)
        b.pop(HEADLINE_BINDINGS[0])  # drop target_epitope_score
        return b

    dsl_mod.combined_score_bindings = _stripped
    try:
        vp = _make_vp()
        tree = parse_combined_score_expr('log(-1)')
        # Trigger eval-time failure and capture the warning record.
        caplog_records = []

        class _Capture(logging.Handler):
            def emit(self, record):
                caplog_records.append(record)

        handler = _Capture(level=logging.WARNING)
        dsl_mod.logger.addHandler(handler)
        try:
            with pytest.raises(ValueError):
                evaluate_combined_score(tree, vp)
        finally:
            dsl_mod.logger.removeHandler(handler)
        warnings = [r for r in caplog_records if r.levelno == logging.WARNING]
        assert any(
            'missing expected headline binding' in r.getMessage()
            for r in warnings), (
                "Expected a warning about the missing headline "
                "binding; got: %r" % [r.getMessage() for r in warnings])
    finally:
        dsl_mod.combined_score_bindings = real_extractor


def _isovar_backed_peptide():
    """A stub vaccine peptide whose fragment came from isovar.

    Only the attributes ``combined_score_bindings`` reads are needed, and
    a real VaccinePeptide would drag in epitope selection irrelevant here.
    """
    import types

    from vaxrank.mutant_protein_fragment import MutantProteinFragment

    protein_sequence = types.SimpleNamespace(
        gene_name="BRAF", amino_acids="SIINFEKLAA",
        mutation_start_idx=2, mutation_end_idx=3,
        num_supporting_fragments=7, num_supporting_reads=13, transcripts=[])
    fragment = MutantProteinFragment.from_isovar_result(types.SimpleNamespace(
        variant=object(), top_protein_sequence=protein_sequence,
        num_total_fragments=20, num_alt_fragments=8, num_ref_fragments=12,
        num_total_reads=38, num_alt_reads=15, num_ref_reads=23))
    return types.SimpleNamespace(
        target_epitope_score=1.0, self_epitope_score=0.0,
        expression_score=1.0, mutant_protein_fragment=fragment)


def test_the_portable_name_reports_what_the_source_counted():
    """``n_rna_alt`` is the unit the run actually has; ``n_alt_reads`` is not.

    On the isovar path the two are different numbers for the same
    evidence — 8 fragments and 15 reads — and ``sqrt()`` weights whichever
    it is given identically. A threshold naming ``n_alt_reads`` therefore
    means a different bar depending on where the file came from, which is
    the portability problem (#385). ``n_rna_alt`` means one thing
    everywhere: what the source counted.
    """
    from vaxrank.combined_score_dsl import combined_score_bindings

    bindings = combined_score_bindings(_isovar_backed_peptide())

    assert bindings['n_alt_reads'] == 15.0     # reads, isovar's own
    assert bindings['n_rna_alt'] == 8.0        # fragments, what it counted


def test_a_fragment_with_no_rna_binds_zeros_rather_than_failing():
    """The varcode path consults no RNA, and must still score.

    n_rna_alt is None there rather than 0, so binding it unguarded would
    raise on exactly the runs that have no RNA — a config that works until
    someone runs without a BAM.
    """
    import types

    from vaxrank.combined_score_dsl import combined_score_bindings
    from vaxrank.mutant_protein_fragment import MutantProteinFragment

    fragment = MutantProteinFragment(
        variant=object(), gene_name="BRAF", amino_acids="SIINFEKL",
        mutant_amino_acid_start_offset=0, mutant_amino_acid_end_offset=1,
        supporting_reference_transcripts=[], n_overlapping_reads=0,
        n_alt_reads=0, n_ref_reads=0,
        n_alt_reads_supporting_protein_sequence=0)
    peptide = types.SimpleNamespace(
        target_epitope_score=1.0, self_epitope_score=0.0,
        expression_score=0.0, mutant_protein_fragment=fragment)

    bindings = combined_score_bindings(peptide)

    assert bindings['n_rna_alt'] == 0.0
    assert bindings['n_dna_alt'] == 0.0
    assert 'rna_evidence_subject' not in bindings
