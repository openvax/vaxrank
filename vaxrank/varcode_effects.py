from collections import OrderedDict

from varcode import MultiOutcomeEffect


OUTCOME_SELECTION_HIGHEST_PRIORITY = "highest_priority"
OUTCOME_SELECTION_MOST_LIKELY = "most_likely"
OUTCOME_SELECTION_MULTI_OUTCOME = "multi_outcome"

OUTCOME_SELECTIONS = {
    OUTCOME_SELECTION_HIGHEST_PRIORITY,
    OUTCOME_SELECTION_MOST_LIKELY,
    OUTCOME_SELECTION_MULTI_OUTCOME,
}


def is_multi_outcome_effect(effect):
    return isinstance(effect, MultiOutcomeEffect)


def select_varcode_effect_outcome(
        effect,
        outcome_selection=OUTCOME_SELECTION_HIGHEST_PRIORITY):
    """Collapse a varcode 5 multi-outcome effect to one concrete effect.

    ``most_likely_effect`` is producer-order-first; ``highest_priority_effect``
    is varcode's most protein-disruptive candidate. Vaxrank's peptide
    mechanics usually need a single concrete effect, so callers should choose
    the collapse that matches the question rather than relying on varcode's
    compatibility shims.
    """
    if outcome_selection not in OUTCOME_SELECTIONS:
        raise ValueError(
            "Unknown varcode outcome selection %r. Expected one of %s" % (
                outcome_selection,
                ", ".join(sorted(OUTCOME_SELECTIONS))))
    if effect is None:
        return None
    if outcome_selection == OUTCOME_SELECTION_MULTI_OUTCOME:
        return effect
    if not is_multi_outcome_effect(effect):
        return effect
    if outcome_selection == OUTCOME_SELECTION_MOST_LIKELY:
        return effect.most_likely_effect
    return effect.highest_priority_effect


def format_varcode_effect(effect):
    if effect is None:
        return "n/a"
    description = getattr(effect, "short_description", None) or "n/a"
    return "%s: %s" % (effect.__class__.__name__, description)


def format_varcode_candidate(candidate):
    text = format_varcode_effect(candidate.effect)
    source = getattr(candidate, "source", None)
    if source:
        text = "%s (%s)" % (text, source)
    return text


def summarize_varcode_effect_outcomes(effect):
    """Report rows for varcode 5 multi-outcome effects."""
    if not is_multi_outcome_effect(effect):
        return OrderedDict()
    return OrderedDict([
        ("Outcome set type", effect.__class__.__name__),
        ("Most likely effect", format_varcode_effect(effect.most_likely_effect)),
        ("Highest priority effect", format_varcode_effect(
            effect.highest_priority_effect)),
        ("Candidate effects", "; ".join(
            format_varcode_candidate(c)
            for c in effect.candidates)),
    ])
