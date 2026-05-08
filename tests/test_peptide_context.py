# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0

"""Tests for ``vaxrank.peptide_context`` (replaces #282 / supersedes
the flat ``EpitopePrediction``).

Pins:
- Native nested storage: ``{kind: {predictor: {version: tuple}}}``.
  Constructor accepts flat ``Sequence[Prediction]`` for producer
  ergonomics; ``__post_init__`` groups.
- Allele is a *record property*, not a structural key — works
  cleanly for processing kinds (no allele) and for affinity /
  presentation (per-allele records).
- ``best`` / ``predictions_for`` raise on multi-predictor ambiguity
  (cross-predictor ranking is meaningless); auto-resolve multi-
  version to the most recent (PEP 440).
- Kind aliases defer to ``topiary.ranking._KIND_ALIASES``.
- ``CandidateEpitope.comparators`` dict accommodates ``'wt'`` plus
  future ``nearest_self`` / ``nearest_vital_self`` /
  ``nearest_nonCTA`` / ``nearest_oncovirus`` (#254 / #257 / #258).
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


# ---- Storage shape: native nested with auto-grouping --------------------


def test_constructor_auto_groups_flat_input():
    """Producers hand back flat ``list[Prediction]``; the
    constructor groups into the nested store. Pin so producer
    ergonomics stay intact even though storage is nested."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', predictor='mhcflurry',
                  version='2.1.1', allele='HLA-A*02:01'),
            _pred('pMHC_affinity', predictor='netmhcpan',
                  version='4.1', allele='HLA-A*02:01'),
            _pred('pMHC_presentation', predictor='mhcflurry',
                  version='2.1.1', allele='HLA-A*02:01'),
        ))
    assert set(ctx.predictions.keys()) == {
        'pMHC_affinity', 'pMHC_presentation'}
    assert set(ctx.predictions['pMHC_affinity'].keys()) == {
        'mhcflurry', 'netmhcpan'}
    assert set(ctx.predictions['pMHC_affinity']['mhcflurry'].keys()) == {
        '2.1.1'}
    # Leaf is a tuple of Predictions, not allele-keyed.
    leaf = ctx.predictions['pMHC_affinity']['mhcflurry']['2.1.1']
    assert isinstance(leaf, tuple)
    assert len(leaf) == 1
    assert leaf[0].allele == 'HLA-A*02:01'


def test_constructor_accepts_already_nested_input():
    """Pre-grouped nested dict passes through unchanged. Lets
    callers bypass auto-grouping when they've already done it."""
    nested = {
        'pMHC_affinity': {
            'mhcflurry': {
                '2.1.1': (_pred('pMHC_affinity', version='2.1.1'),)
            }
        }
    }
    ctx = PeptideContext(peptide_sequence='SIINFEKL', predictions=nested)
    assert ctx.predictions is nested


def test_storage_collects_per_allele_records_in_leaf_tuple():
    """Affinity emits one Prediction per allele under one
    (predictor, version). The leaf tuple holds them all — alleles
    aren't structural keys, so no collision."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', allele='HLA-A*02:01', score=0.4),
            _pred('pMHC_affinity', allele='HLA-B*07:02', score=0.9),
            _pred('pMHC_affinity', allele='HLA-C*03:04', score=0.3),
        ))
    leaf = ctx.predictions['pMHC_affinity']['mhcflurry']['']
    assert len(leaf) == 3
    assert {p.allele for p in leaf} == {
        'HLA-A*02:01', 'HLA-B*07:02', 'HLA-C*03:04'}


def test_storage_handles_no_allele_kinds():
    """Processing kinds (proteasome_cleavage, antigen_processing)
    emit Predictions with empty ``allele``. The nested store
    accommodates them naturally — leaf is a one-element tuple,
    no allele-keyed dict to collapse into."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('proteasome_cleavage', predictor='pepsickle',
                  allele='', score=0.6),
            _pred('antigen_processing', predictor='netchop',
                  allele='', score=0.8),
        ))
    cleavage_leaf = ctx.predictions['proteasome_cleavage']['pepsickle']['']
    assert len(cleavage_leaf) == 1
    assert cleavage_leaf[0].allele == ''


def test_predictions_flat_round_trips_through_constructor():
    """``predictions_flat()`` flattens the nested store back to a
    tuple. Feeding that tuple back through the constructor
    reproduces the same nested shape — round-trip pinned for
    serialization paths."""
    original = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', allele='HLA-A*02:01'),
            _pred('pMHC_presentation', allele='HLA-A*02:01'),
        ))
    flat = original.predictions_flat()
    rebuilt = PeptideContext(
        peptide_sequence=original.peptide_sequence,
        predictions=flat)
    assert rebuilt.predictions == original.predictions


# ---- Structured views --------------------------------------------------


def test_kinds_returns_sorted_distinct_values():
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity'),
            _pred('pMHC_presentation'),
            _pred('pMHC_affinity', predictor='netmhcpan'),
        ))
    assert ctx.kinds() == ('pMHC_affinity', 'pMHC_presentation')


def test_kinds_empty_when_no_predictions():
    ctx = PeptideContext(peptide_sequence='SIINFEKL')
    assert ctx.kinds() == ()
    assert ctx.predictors_for('pMHC_affinity') == ()


def test_predictors_for_lists_emitting_predictors():
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', predictor='mhcflurry'),
            _pred('pMHC_affinity', predictor='netmhcpan'),
            _pred('pMHC_presentation', predictor='mhcflurry'),
        ))
    assert ctx.predictors_for('pMHC_affinity') == ('mhcflurry', 'netmhcpan')
    assert ctx.predictors_for('pMHC_presentation') == ('mhcflurry',)
    assert ctx.predictors_for('immunogenicity') == ()


def test_alleles_for_lists_alleles_only():
    """Returns alleles among leaf records. Filters out empty
    strings (processing kinds)."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', allele='HLA-A*02:01'),
            _pred('pMHC_affinity', allele='HLA-B*07:02'),
            _pred('proteasome_cleavage', predictor='pepsickle', allele=''),
        ))
    assert ctx.alleles_for('pMHC_affinity') == (
        'HLA-A*02:01', 'HLA-B*07:02')
    # Processing kind: no alleles surfaced.
    assert ctx.alleles_for('proteasome_cleavage') == ()


# ---- Kind aliases (delegated to topiary) -------------------------------


@pytest.mark.parametrize('alias,canonical', [
    ('ba', 'pMHC_affinity'),
    ('BA', 'pMHC_affinity'),
    ('aff', 'pMHC_affinity'),
    ('ic50', 'pMHC_affinity'),
    ('affinity', 'pMHC_affinity'),
    ('pMHC_affinity', 'pMHC_affinity'),
    ('PMHC_AFFINITY', 'pMHC_affinity'),
    ('el', 'pMHC_presentation'),
    ('EL', 'pMHC_presentation'),
    ('presentation', 'pMHC_presentation'),
    ('stability', 'pMHC_stability'),
    ('processing', 'antigen_processing'),
])
def test_best_resolves_kind_aliases_via_topiary(alias, canonical):
    """vaxrank doesn't maintain its own alias table — defers to
    topiary's ``_KIND_ALIASES``. Pin the entries vaxrank cares
    about so a topiary change shows up here."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(_pred(canonical, score=0.7),))
    result = ctx.best(alias)
    assert result is not None
    assert result.kind == canonical


def test_predictors_for_accepts_aliases():
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_presentation', predictor='mhcflurry'),
            _pred('pMHC_presentation', predictor='bigmhc'),
        ))
    assert ctx.predictors_for('el') == ('bigmhc', 'mhcflurry')
    assert ctx.predictors_for('presentation') == ('bigmhc', 'mhcflurry')


# ---- _load_kind_aliases: defensive topiary import ----------------------


def test_load_kind_aliases_prefers_public_name(monkeypatch):
    """When topiary publishes a public ``KIND_ALIASES``, vaxrank
    prefers it over the legacy private ``_KIND_ALIASES``. Pins the
    upgrade path so when topiary promotes the name, vaxrank picks
    it up without code changes."""
    import topiary.ranking
    from vaxrank import peptide_context as pc

    sentinel = {'sentinel_alias': 'pMHC_affinity'}
    pc._load_kind_aliases.cache_clear()
    monkeypatch.setattr(
        topiary.ranking, 'KIND_ALIASES', sentinel, raising=False)
    try:
        assert pc._load_kind_aliases() is sentinel
    finally:
        pc._load_kind_aliases.cache_clear()


def test_load_kind_aliases_raises_when_neither_present(monkeypatch):
    """If a future topiary refactor removes both names, vaxrank
    raises ``ImportError`` with a clear diagnostic instead of
    letting an opaque ``AttributeError`` surface from inside
    PeptideContext construction. Proves the safety net actually
    fires + the message is well-formed."""
    import topiary.ranking
    from vaxrank import peptide_context as pc

    pc._load_kind_aliases.cache_clear()
    monkeypatch.delattr(topiary.ranking, 'KIND_ALIASES', raising=False)
    monkeypatch.delattr(topiary.ranking, '_KIND_ALIASES', raising=False)
    try:
        with pytest.raises(ImportError, match='KIND_ALIASES'):
            pc._load_kind_aliases()
    finally:
        pc._load_kind_aliases.cache_clear()


def test_load_kind_aliases_is_cached(monkeypatch):
    """``_load_kind_aliases`` is on the per-record hot path of
    ``best`` / ``predictions_for``; the result is module-constant
    so we cache it (``functools.lru_cache``) instead of redoing
    the ``getattr`` chain on every call."""
    import topiary.ranking
    from vaxrank import peptide_context as pc

    pc._load_kind_aliases.cache_clear()
    first = pc._load_kind_aliases()
    # Swap the underlying table; cached call should NOT see the swap.
    sentinel = {'sentinel_alias': 'pMHC_affinity'}
    monkeypatch.setattr(
        topiary.ranking, 'KIND_ALIASES', sentinel, raising=False)
    try:
        assert pc._load_kind_aliases() is first  # cached
        # After cache_clear, the new attr wins.
        pc._load_kind_aliases.cache_clear()
        assert pc._load_kind_aliases() is sentinel
    finally:
        pc._load_kind_aliases.cache_clear()


# ---- predictions_for: leaf access + disambiguation ---------------------


def test_predictions_for_returns_leaf_tuple():
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', allele='HLA-A*02:01', score=0.4),
            _pred('pMHC_affinity', allele='HLA-B*07:02', score=0.9),
        ))
    records = ctx.predictions_for('pMHC_affinity')
    assert len(records) == 2
    assert {p.allele for p in records} == {'HLA-A*02:01', 'HLA-B*07:02'}


def test_predictions_for_empty_when_no_records():
    ctx = PeptideContext(peptide_sequence='SIINFEKL')
    assert ctx.predictions_for('pMHC_affinity') == ()


def test_predictions_for_multi_predictor_raises():
    """Same disambiguation contract as ``best`` — predictor
    required when multiple emit ``kind``."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', predictor='mhcflurry', score=0.4),
            _pred('pMHC_affinity', predictor='netmhcpan', score=0.9),
        ))
    with pytest.raises(ValueError, match='multiple predictors'):
        ctx.predictions_for('pMHC_affinity')


def test_predictions_for_with_predictor_disambiguates():
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', predictor='mhcflurry',
                  allele='HLA-A*02:01', score=0.4),
            _pred('pMHC_affinity', predictor='netmhcpan',
                  allele='HLA-A*02:01', score=0.9),
        ))
    mhcf = ctx.predictions_for('pMHC_affinity', predictor='mhcflurry')
    assert len(mhcf) == 1
    assert mhcf[0].score == 0.4
    nmp = ctx.predictions_for('pMHC_affinity', predictor='netmhcpan')
    assert nmp[0].score == 0.9


def test_predictions_for_unknown_predictor_returns_empty():
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', predictor='mhcflurry', score=0.4),
        ))
    assert ctx.predictions_for(
        'pMHC_affinity', predictor='netmhcpan') == ()


def test_predictions_for_lets_caller_rank_on_arbitrary_axis():
    """The whole point of exposing the leaf tuple: callers rank
    on whatever axis they want. ``best`` does the common
    ``max(score)`` case; for "lowest IC50" / "lowest %-rank" /
    custom DSL, iterate ``predictions_for`` directly. This is
    *the* parametric-scoring affordance — no callable plumbed
    through ``best``."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', allele='HLA-A*02:01',
                  score=0.9, value=10.0, percentile_rank=2.5),
            _pred('pMHC_affinity', allele='HLA-B*07:02',
                  score=0.4, value=2.0, percentile_rank=0.3),
        ))
    records = ctx.predictions_for('pMHC_affinity')
    # Default best (by score): A*02 wins.
    assert ctx.best('pMHC_affinity').allele == 'HLA-A*02:01'
    # Caller-driven: rank by lowest IC50 (value).
    by_ic50 = min(records, key=lambda p: p.value)
    assert by_ic50.allele == 'HLA-B*07:02'
    # Caller-driven: rank by lowest percentile_rank.
    by_rank = min(records, key=lambda p: p.percentile_rank)
    assert by_rank.allele == 'HLA-B*07:02'


# ---- Version disambiguation --------------------------------------------


def test_predictions_for_multi_version_picks_most_recent():
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', version='2.1.0', score=0.9),
            _pred('pMHC_affinity', version='2.1.1', score=0.4),
        ))
    records = ctx.predictions_for('pMHC_affinity')
    assert len(records) == 1
    assert records[0].predictor_version == '2.1.1'


def test_predictions_for_with_explicit_version_pins():
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', version='2.1.0', score=0.9),
            _pred('pMHC_affinity', version='2.1.1', score=0.4),
        ))
    pinned = ctx.predictions_for('pMHC_affinity', version='2.1.0')
    assert len(pinned) == 1
    assert pinned[0].predictor_version == '2.1.0'


def test_predictions_for_unknown_version_returns_empty():
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(_pred('pMHC_affinity', version='2.1.1'),))
    assert ctx.predictions_for(
        'pMHC_affinity', version='99.0.0') == ()


def test_predictions_for_handles_empty_version_strings():
    """Pre-#261 records had empty ``predictor_version``. Mixed
    with a real version, the real version wins (empty parses as
    InvalidVersion and falls below all PEP 440 versions)."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', version='', score=0.9),
            _pred('pMHC_affinity', version='2.1.1', score=0.4),
        ))
    records = ctx.predictions_for('pMHC_affinity')
    assert records[0].predictor_version == '2.1.1'


def test_predictions_for_all_invalid_versions_falls_back_lex():
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', version='legacy_a', score=0.4),
            _pred('pMHC_affinity', version='legacy_z', score=0.7),
        ))
    records = ctx.predictions_for('pMHC_affinity')
    assert records[0].predictor_version == 'legacy_z'


def test_versions_for_returns_oldest_to_newest():
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', version='2.1.1'),
            _pred('pMHC_affinity', version='2.1.0'),
            _pred('pMHC_affinity', version='2.0.0'),
        ))
    assert ctx.versions_for('pMHC_affinity', 'mhcflurry') == (
        '2.0.0', '2.1.0', '2.1.1')


def test_versions_for_invalid_strings_sort_first():
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', version='legacy'),
            _pred('pMHC_affinity', version='2.1.1'),
            _pred('pMHC_affinity', version='2.1.0'),
        ))
    assert ctx.versions_for('pMHC_affinity', 'mhcflurry') == (
        'legacy', '2.1.0', '2.1.1')


# ---- best: thin wrapper over predictions_for ---------------------------


def test_best_single_predictor_picks_best_allele():
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
    ctx = PeptideContext(peptide_sequence='SIINFEKL')
    assert ctx.best('pMHC_affinity') is None


def test_best_multi_predictor_raises_without_predictor_arg():
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', predictor='mhcflurry', score=0.4),
            _pred('pMHC_affinity', predictor='netmhcpan', score=0.9),
        ))
    with pytest.raises(ValueError, match='multiple predictors'):
        ctx.best('pMHC_affinity')


def test_kind_named_helpers_forward_to_best():
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('pMHC_affinity', score=0.7),
            _pred('pMHC_presentation', score=0.5),
            _pred('pMHC_stability', score=0.8),
            _pred('immunogenicity', score=0.3),
            _pred('proteasome_cleavage', allele='', score=0.6),
            _pred('antigen_processing', allele='', score=0.4),
        ))
    assert ctx.best_affinity().kind == 'pMHC_affinity'
    assert ctx.best_presentation().kind == 'pMHC_presentation'
    assert ctx.best_stability().kind == 'pMHC_stability'
    assert ctx.best_immunogenicity().kind == 'immunogenicity'
    assert ctx.best_cleavage().kind == 'proteasome_cleavage'
    assert ctx.best_antigen_processing().kind == 'antigen_processing'


def test_best_for_processing_kind_works_without_allele():
    """Processing kinds emit ``allele=''``. ``best`` returns the
    one record without complaining about a missing allele key —
    proves allele-as-record-property handles this naturally."""
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        predictions=(
            _pred('proteasome_cleavage', predictor='pepsickle',
                  allele='', score=0.6),
        ))
    result = ctx.best_cleavage()
    assert result is not None
    assert result.score == 0.6
    assert result.allele == ''


# ---- Flanking residues + provenance -------------------------------------


def test_peptide_context_carries_flanks():
    ctx = PeptideContext(
        peptide_sequence='SIINFEKL',
        n_flank='AAAAA',
        c_flank='CCCCC',
        source_sequence='AAAAASIINFEKLCCCCC',
        source_name='OVA',
        offset=5)
    assert ctx.n_flank == 'AAAAA'
    assert ctx.c_flank == 'CCCCC'
    assert ctx.source_sequence == 'AAAAASIINFEKLCCCCC'
    assert ctx.source_name == 'OVA'
    assert ctx.offset == 5


def test_peptide_context_flanks_default_to_empty():
    ctx = PeptideContext(peptide_sequence='SIINFEKL')
    assert ctx.n_flank == ''
    assert ctx.c_flank == ''
    assert ctx.source_sequence == ''
    assert ctx.source_name == ''


# ---- CandidateEpitope ---------------------------------------------------


def test_candidate_epitope_wt_shortcut():
    mutant = PeptideContext(peptide_sequence='SIINFEKL')
    wt = PeptideContext(peptide_sequence='SIINFEKM')
    candidate = CandidateEpitope(
        mutant=mutant,
        comparators={COMPARATOR_WT: wt},
        overlaps_mutation=True)
    assert candidate.wt is wt
    assert candidate.comparator('wt') is wt
    assert candidate.peptide_sequence == 'SIINFEKL'
    assert candidate.length == 8


def test_candidate_epitope_no_wt_returns_none():
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
    assert candidate.comparator(
        'nearest_self').best_affinity().value == 80.0
    assert candidate.comparator('custom_comparator') is not None
    assert candidate.comparator('not_present') is None


def test_candidate_epitope_freezes_mutation_context():
    candidate = CandidateEpitope(
        mutant=PeptideContext(peptide_sequence='SIINFEKL'),
        comparators={},
        overlaps_mutation=True,
        occurs_in_reference=False)
    assert candidate.overlaps_mutation is True
    assert candidate.occurs_in_reference is False
    with pytest.raises(Exception):
        candidate.overlaps_mutation = False  # type: ignore
