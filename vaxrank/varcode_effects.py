from collections import OrderedDict
from varcode import get_effect_candidates


OUTCOME_SELECTION_HIGHEST_PRIORITY = "highest_priority"
OUTCOME_SELECTION_MOST_LIKELY = "most_likely"
OUTCOME_SELECTION_MULTI_OUTCOME = "multi_outcome"

OUTCOME_SELECTIONS = {
    OUTCOME_SELECTION_HIGHEST_PRIORITY,
    OUTCOME_SELECTION_MOST_LIKELY,
    OUTCOME_SELECTION_MULTI_OUTCOME,
}


def select_varcode_effect_outcome(
        effect,
        outcome_selection=OUTCOME_SELECTION_HIGHEST_PRIORITY):
    """Select a concrete comparison from an ordinary or wrapped effect.

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
    if not get_effect_candidates(effect):
        return effect
    name = ("most_likely" if outcome_selection == OUTCOME_SELECTION_MOST_LIKELY
            else "highest_priority")
    # Compatibility wrappers expose the effect; ordinary Varcode 11 results
    # expose its candidate. Respect the producer's selection in both cases.
    selected = getattr(effect, name + "_effect", None)
    if selected is not None:
        return selected
    candidate = getattr(effect, name + "_candidate", None)
    return candidate.effect if candidate is not None else effect


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
    """Report native alternatives without implying calibrated likelihoods."""
    candidates = get_effect_candidates(effect)
    if not candidates or (len(candidates) == 1 and candidates[0].effect is effect):
        return OrderedDict()
    return OrderedDict([
        ("Outcome set type", effect.__class__.__name__),
        ("Producer-preferred effect", format_varcode_effect(
            select_varcode_effect_outcome(effect, OUTCOME_SELECTION_MOST_LIKELY))),
        ("Highest priority effect", format_varcode_effect(
            select_varcode_effect_outcome(effect))),
        ("Outcome preference", "Producer ordering and consequence severity; not calibrated probabilities"),
        ("Candidate effects", "; ".join(
            format_varcode_candidate(c)
            for c in candidates)),
    ])
