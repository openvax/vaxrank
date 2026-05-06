# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Tests for ``vaxrank.peptide_context`` (replaces #282 / supersedes
the flat ``EpitopePrediction``).

Pins:
- Storage stays as a flat tuple of ``mhctools.Prediction``s, but
  every access path goes through methods that respect the
  (kind, predictor, version, allele) structure.
- ``best`` raises on multi-predictor ambiguity (cross-predictor
  ``max(score)`` is meaningless — each predictor's score scale is
  its own) and auto-resolves multi-version ambiguity to the most
  recent version (PEP 440).
- ``best`` accepts kind aliases (``'ba'`` / ``'el'`` / etc).
- ``predictors_for_kind`` / ``alleles_for`` / ``versions_for`` give
  callers the structured view they need to walk per-predictor.
- ``CandidateEpitope.comparators`` dict accommodates ``'wt'`` plus
  future ``nearest_self`` / ``nearest_vital_self`` / ``nearest_nonCTA``
  / ``nearest_oncovirus`` records (#254 / #257 / #258).
"""

import pytest

from mhctools.pred import Prediction

from vaxrank.peptide_context import (
    COMPARATOR_NEAREST_SELF,
    COMPARATOR_WT,
    CandidateEpitope,
    PeptideContext,
)


def _pred(kind, *, predictor='mhcflurry', version='', allele='HLA-A*02:01',
          score=0.7, value=100.0, percentile_rank=0.5):
    """Concise factory for an ``mhctools.Prediction``."""
    return Prediction(
        kind=kind, predictor_name=predictor, predictor_version=version,
        allele=allele, peptide='SIINFEKL',
        score=score, value=value, percentile_rank=percentile_rank,
    )


# ---- PeptideContext storage / structured views --------------------------


def test_peptide_context_predictions_by_kind_and_predictor_groups_correctly():
    """The flat storage tuple is the source of truth; the nested
    view is built on demand."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', predictor='mhcflurry',
                  allele='HLA-A*02:01'),
            _pred('pMHC_affinity', predictor='netmhcpan',
                  allele='HLA-A*02:01'),
            _pred('pMHC_presentation', predictor='mhcflurry',
                  allele='HLA-A*02:01'),
        ))
    nested = ctx.predictions_by_kind_and_predictor()
    assert set(nested.keys()) == {'pMHC_affinity', 'pMHC_presentation'}
    assert set(nested['pMHC_affinity'].keys()) == {'mhcflurry', 'netmhcpan'}
    assert 'HLA-A*02:01' in nested['pMHC_affinity']['mhcflurry']
    assert 'HLA-A*02:01' in nested['pMHC_presentation']['mhcflurry']


def test_predictors_for_kind_lists_emitting_predictors():
    """Drives multi-predictor disambiguation: callers ask "what
    predictors emitted ``kind`` for this peptide?" before deciding
    how to combine."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', predictor='mhcflurry'),
            _pred('pMHC_affinity', predictor='netmhcpan'),
            _pred('pMHC_presentation', predictor='mhcflurry'),
        ))
    assert ctx.predictors_for_kind('pMHC_affinity') == (
        'mhcflurry', 'netmhcpan')
    # Presentation only from mhcflurry — typical real data shape.
    assert ctx.predictors_for_kind('pMHC_presentation') == ('mhcflurry',)
    # No data → empty tuple.
    assert ctx.predictors_for_kind('immunogenicity') == ()


def test_alleles_for_kind_predictor_lists_alleles_only():
    """Alleles per (kind, predictor) — used by coverage code that
    walks per-allele evidence without picking a single 'best'."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', predictor='mhcflurry',
                  allele='HLA-A*02:01'),
            _pred('pMHC_affinity', predictor='mhcflurry',
                  allele='HLA-B*07:02'),
            _pred('pMHC_affinity', predictor='netmhcpan',
                  allele='HLA-A*02:01'),
        ))
    assert ctx.alleles_for('pMHC_affinity', 'mhcflurry') == (
        'HLA-A*02:01', 'HLA-B*07:02')
    assert ctx.alleles_for('pMHC_affinity', 'netmhcpan') == ('HLA-A*02:01',)
    assert ctx.alleles_for('pMHC_affinity', 'unknown_predictor') == ()


def test_kinds_returns_sorted_distinct_values():
    """Helps report code render only the kinds present without
    iterating predictions."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity'),
            _pred('pMHC_presentation'),
            _pred('pMHC_affinity', predictor='netmhcpan'),
        ))
    assert ctx.kinds() == ('pMHC_affinity', 'pMHC_presentation')


def test_kinds_empty_when_no_predictions():
    """A bare PeptideContext (e.g. a comparator that hasn't been
    binding-probed yet) has no kinds."""
    ctx = PeptideContext(peptide_sequence='SIINFEKL')
    assert ctx.kinds() == ()
    assert ctx.predictors_for_kind('pMHC_affinity') == ()


# ---- best: single-predictor unambiguous case ------------------


def test_best_single_predictor_picks_best_allele():
    """One predictor, multiple alleles → best-by-score across
    alleles. ``predictor=`` not required when there's only one."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', allele='HLA-A*02:01', score=0.4),
            _pred('pMHC_affinity', allele='HLA-B*07:02', score=0.9),
            _pred('pMHC_affinity', allele='HLA-C*03:04', score=0.3),
        ))
    best = ctx.best('pMHC_affinity')
    assert best is not None
    assert best.allele == 'HLA-B*07:02'
    assert best.score == 0.9


def test_best_returns_none_when_no_records():
    """Coverage / report code can call this safely; absence reads
    as ``None``."""
    ctx = PeptideContext(peptide_sequence='SIINFEKL')
    assert ctx.best('pMHC_affinity') is None
    assert ctx.best('immunogenicity') is None


# ---- best: multi-predictor disambiguation ---------------------


def test_best_multi_predictor_raises_without_predictor_arg():
    """The architectural fix from this design pass: cross-predictor
    ``max(score)`` is meaningless — score scales differ per
    predictor — so the API forces an explicit choice."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', predictor='mhcflurry', score=0.4),
            _pred('pMHC_affinity', predictor='netmhcpan', score=0.9),
        ))
    with pytest.raises(ValueError, match='multiple predictors'):
        ctx.best('pMHC_affinity')


def test_best_with_predictor_arg_picks_within_that_predictor():
    """Per-predictor disambiguation: caller picks the predictor,
    accessor returns the best allele for it."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', predictor='mhcflurry',
                  allele='HLA-A*02:01', score=0.4),
            _pred('pMHC_affinity', predictor='mhcflurry',
                  allele='HLA-B*07:02', score=0.7),
            _pred('pMHC_affinity', predictor='netmhcpan',
                  allele='HLA-A*02:01', score=0.9),
            _pred('pMHC_affinity', predictor='netmhcpan',
                  allele='HLA-B*07:02', score=0.5),
        ))
    # Best within mhcflurry: B*07 (0.7), not netmhcpan's A*02 (0.9).
    best_mhcf = ctx.best('pMHC_affinity', predictor='mhcflurry')
    assert best_mhcf.allele == 'HLA-B*07:02'
    assert best_mhcf.predictor_name == 'mhcflurry'
    # Best within netmhcpan: A*02.
    best_nmp = ctx.best('pMHC_affinity', predictor='netmhcpan')
    assert best_nmp.allele == 'HLA-A*02:01'
    assert best_nmp.predictor_name == 'netmhcpan'


def test_best_unknown_predictor_returns_none():
    """Asking for a predictor that didn't emit ``kind`` gives
    ``None`` (not an error). Lets callers iterate
    ``predictors_for_kind`` and call ``best`` without
    pre-checking."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', predictor='mhcflurry', score=0.4),
        ))
    assert ctx.best(
        'pMHC_affinity', predictor='netmhcpan') is None


# ---- Kind-named helpers ------------------------------------------------


def test_kind_named_helpers_forward_to_best():
    """``best_affinity`` / ``best_presentation`` / etc. are thin
    forwards. Verify each names the right ``kind`` string."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', score=0.7),
            _pred('pMHC_presentation', score=0.5),
            _pred('pMHC_stability', score=0.8),
            _pred('immunogenicity', score=0.3),
            _pred('proteasome_cleavage', score=0.6),
            _pred('antigen_processing', score=0.4),
        ))
    assert ctx.best_affinity().kind == 'pMHC_affinity'
    assert ctx.best_presentation().kind == 'pMHC_presentation'
    assert ctx.best_stability().kind == 'pMHC_stability'
    assert ctx.best_immunogenicity().kind == 'immunogenicity'
    assert ctx.best_cleavage().kind == 'proteasome_cleavage'
    assert ctx.best_antigen_processing().kind == 'antigen_processing'


def test_kind_named_helpers_raise_on_multi_predictor_ambiguity():
    """Same disambiguation contract as ``best``."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_presentation', predictor='mhcflurry'),
            _pred('pMHC_presentation', predictor='bigmhc'),
        ))
    with pytest.raises(ValueError):
        ctx.best_presentation()
    # Disambiguated form works.
    assert ctx.best_presentation(predictor='mhcflurry') is not None


# ---- Kind aliases ------------------------------------------------------


@pytest.mark.parametrize('alias,canonical', [
    ('ba', 'pMHC_affinity'),
    ('BA', 'pMHC_affinity'),
    ('affinity', 'pMHC_affinity'),
    ('binding', 'pMHC_affinity'),
    ('pMHC_affinity', 'pMHC_affinity'),
    ('PMHC_AFFINITY', 'pMHC_affinity'),
    ('el', 'pMHC_presentation'),
    ('EL', 'pMHC_presentation'),
    ('presentation', 'pMHC_presentation'),
    ('elution', 'pMHC_presentation'),
    ('stability', 'pMHC_stability'),
    ('cleavage', 'proteasome_cleavage'),
    ('proteasome', 'proteasome_cleavage'),
    ('processing', 'antigen_processing'),
    ('ap', 'antigen_processing'),
    ('tap', 'tap_transport'),
    ('erap', 'erap_trimming'),
])
def test_best_resolves_kind_aliases(alias, canonical):
    """mhcflurry / netmhcpan / pVACseq use BA / EL / etc. as
    shorthand for the canonical ``pMHC_*`` strings — accept both
    so call sites can stay terse without knowing which spelling
    vaxrank stores."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(_pred(canonical, score=0.7),))
    result = ctx.best(alias)
    assert result is not None
    assert result.kind == canonical


def test_predictors_for_kind_accepts_aliases():
    """Same alias contract on the structured-view accessors so
    callers don't have to switch spellings between methods."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_presentation', predictor='mhcflurry'),
            _pred('pMHC_presentation', predictor='bigmhc'),
        ))
    assert ctx.predictors_for_kind('el') == ('bigmhc', 'mhcflurry')
    assert ctx.predictors_for_kind('presentation') == ('bigmhc', 'mhcflurry')


def test_alleles_for_accepts_aliases():
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', allele='HLA-A*02:01'),
            _pred('pMHC_affinity', allele='HLA-B*07:02'),
        ))
    assert ctx.alleles_for('ba', 'mhcflurry') == (
        'HLA-A*02:01', 'HLA-B*07:02')


# ---- Version disambiguation --------------------------------------------


def test_best_multi_version_picks_most_recent():
    """When the same predictor ran at two versions (e.g. mhcflurry
    2.1.0 and 2.1.1 both on file), ``best`` defaults to the newest
    by PEP 440. Score scales are comparable across versions of the
    same predictor — only freshness differs — so this auto-resolves
    rather than raising."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', version='2.1.0', score=0.9),
            _pred('pMHC_affinity', version='2.1.1', score=0.4),
        ))
    result = ctx.best('pMHC_affinity')
    assert result is not None
    # Newer version wins even though its score is lower.
    assert result.predictor_version == '2.1.1'
    assert result.score == 0.4


def test_best_with_explicit_version_pins_to_that_version():
    """``version=`` overrides the most-recent default — useful
    when reproducing an older run or comparing version-over-version."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', version='2.1.0', score=0.9),
            _pred('pMHC_affinity', version='2.1.1', score=0.4),
        ))
    result = ctx.best('pMHC_affinity', version='2.1.0')
    assert result is not None
    assert result.predictor_version == '2.1.0'
    assert result.score == 0.9


def test_best_unknown_version_returns_none():
    """Asking for a version that isn't on file gives ``None``,
    not an error — same shape as unknown predictor."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', version='2.1.1', score=0.7),
        ))
    assert ctx.best('pMHC_affinity', version='99.0.0') is None


def test_best_picks_best_allele_within_most_recent_version():
    """Combination: multiple versions × multiple alleles. Pin to
    newest version, then take best-by-score across alleles."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', version='2.1.0',
                  allele='HLA-A*02:01', score=0.95),
            _pred('pMHC_affinity', version='2.1.1',
                  allele='HLA-A*02:01', score=0.4),
            _pred('pMHC_affinity', version='2.1.1',
                  allele='HLA-B*07:02', score=0.7),
        ))
    result = ctx.best('pMHC_affinity')
    assert result.predictor_version == '2.1.1'
    assert result.allele == 'HLA-B*07:02'  # best within 2.1.1


def test_best_per_predictor_with_version_disambiguation():
    """Two predictors, each with multiple versions. ``predictor=``
    picks the predictor; version then auto-resolves within it."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', predictor='mhcflurry',
                  version='2.1.0', score=0.9),
            _pred('pMHC_affinity', predictor='mhcflurry',
                  version='2.1.1', score=0.4),
            _pred('pMHC_affinity', predictor='netmhcpan',
                  version='4.0', score=0.5),
            _pred('pMHC_affinity', predictor='netmhcpan',
                  version='4.1', score=0.3),
        ))
    mhcf = ctx.best('pMHC_affinity', predictor='mhcflurry')
    assert mhcf.predictor_version == '2.1.1'
    nmp = ctx.best('pMHC_affinity', predictor='netmhcpan')
    assert nmp.predictor_version == '4.1'


def test_best_handles_empty_version_strings():
    """Pre-#261 records had empty ``predictor_version`` strings.
    Mixed with a real version, the real version wins (empty parses
    as InvalidVersion and falls below all PEP 440 versions)."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', version='', score=0.9),
            _pred('pMHC_affinity', version='2.1.1', score=0.4),
        ))
    result = ctx.best('pMHC_affinity')
    assert result.predictor_version == '2.1.1'


def test_best_all_invalid_version_strings_falls_back_to_lex_max():
    """If every ``predictor_version`` on file is unparseable
    (legacy / non-semver / pre-#261), there's no PEP 440 ordering
    to apply. Lexicographic ``max`` is the consistent fallback so
    behavior remains deterministic."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', version='legacy_a', score=0.4),
            _pred('pMHC_affinity', version='legacy_z', score=0.7),
        ))
    result = ctx.best('pMHC_affinity')
    assert result is not None
    assert result.predictor_version == 'legacy_z'


def test_versions_for_handles_invalid_version_strings():
    """Invalid versions sort before valid ones in the returned
    tuple (legacy first, semver after) so ``versions[-1]`` still
    yields the most recent valid version."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', version='legacy', score=0.5),
            _pred('pMHC_affinity', version='2.1.0', score=0.5),
            _pred('pMHC_affinity', version='2.1.1', score=0.5),
        ))
    versions = ctx.versions_for('pMHC_affinity', 'mhcflurry')
    assert versions == ('legacy', '2.1.0', '2.1.1')


def test_versions_for_returns_oldest_to_newest():
    """``versions_for`` exposes every version on file, sorted
    oldest → newest by PEP 440. ``versions[-1]`` is the latest."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', version='2.1.1'),
            _pred('pMHC_affinity', version='2.1.0'),
            _pred('pMHC_affinity', version='2.0.0'),
        ))
    versions = ctx.versions_for('pMHC_affinity', 'mhcflurry')
    assert versions == ('2.0.0', '2.1.0', '2.1.1')


def test_kind_named_helpers_forward_version_arg():
    """``best_affinity(version=…)`` forwards through to ``best``
    so the version override works on the convenience helpers too."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', version='2.1.0', score=0.9),
            _pred('pMHC_affinity', version='2.1.1', score=0.4),
        ))
    # No version → newest.
    assert ctx.best_affinity().predictor_version == '2.1.1'
    # Pinned version overrides.
    assert ctx.best_affinity(version='2.1.0').predictor_version == '2.1.0'


# ---- CandidateEpitope ---------------------------------------------------


def test_candidate_epitope_wt_shortcut():
    """``candidate.wt`` is a convenience for
    ``candidate.comparators['wt']``. Most callers reach for the WT
    pair without thinking about the dict."""
    mutant = PeptideContext(peptide_sequence='SIINFEKL')
    wt = PeptideContext(peptide_sequence='SIINFEKM')  # ref-allele variant
    candidate = CandidateEpitope(
        mutant=mutant,
        comparators={COMPARATOR_WT: wt},
        overlaps_mutation=True)
    assert candidate.wt is wt
    assert candidate.comparator('wt') is wt
    # Pass-through for the most common reads.
    assert candidate.peptide_sequence == 'SIINFEKL'
    assert candidate.length == 8


def test_candidate_epitope_no_wt_returns_none():
    """A WT-less candidate (e.g. a frameshift where there's no
    aligned reference position) has no ``'wt'`` entry — accessor
    is ``None`` rather than KeyError."""
    candidate = CandidateEpitope(
        mutant=PeptideContext(peptide_sequence='SIINFEKL'),
        comparators={})
    assert candidate.wt is None
    assert candidate.comparator('wt') is None


def test_candidate_epitope_holds_arbitrary_comparators():
    """Open-ended comparator names — future safety scorers (#254 /
    #257 / #258) populate ``nearest_self`` / ``nearest_vital_self``
    / ``nearest_nonCTA`` / ``nearest_oncovirus`` without changing
    the data shape."""
    mutant = PeptideContext(peptide_sequence='SIINFEKL')
    nearest_self = PeptideContext(
        peptide_sequence='SIINFEKM',
        source_name='HumanProteome:OVA',
        predictions=(
            _pred('pMHC_affinity', allele='HLA-A*02:01', value=80.0),
        ))
    candidate = CandidateEpitope(
        mutant=mutant,
        comparators={
            COMPARATOR_WT: PeptideContext(peptide_sequence='SIINFEKM'),
            COMPARATOR_NEAREST_SELF: nearest_self,
            'custom_comparator': PeptideContext(peptide_sequence='SIINFEEK'),
        })
    assert candidate.comparator('nearest_self') is nearest_self
    # The nearest-self comparator carries its own predictions —
    # which is what lets a safety scorer ask "does the patient's
    # MHC also bind the self peptide?"
    assert candidate.comparator('nearest_self').best_affinity().value == 80.0
    # Open-ended: arbitrary keys also work.
    assert candidate.comparator('custom_comparator') is not None
    assert candidate.comparator('not_present') is None


def test_candidate_epitope_freezes_mutation_context():
    """``overlaps_mutation`` / ``occurs_in_reference`` live at the
    candidate level (they're about the mutant-vs-source
    relationship, not the peptide itself)."""
    candidate = CandidateEpitope(
        mutant=PeptideContext(peptide_sequence='SIINFEKL'),
        comparators={},
        overlaps_mutation=True,
        occurs_in_reference=False)
    assert candidate.overlaps_mutation is True
    assert candidate.occurs_in_reference is False
    # Frozen dataclass — direct attribute mutation raises.
    with pytest.raises(Exception):
        candidate.overlaps_mutation = False  # type: ignore
