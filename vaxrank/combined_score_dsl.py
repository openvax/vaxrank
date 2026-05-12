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

"""Vaccine-peptide-level combined-score DSL.

Mirrors the per-epitope DSL in :mod:`vaxrank.epitope_dsl` but at the
``VaccinePeptide`` scope. Lets users write expressions like::

    combined_score_expr: "sqrt(n_alt_reads) * target_epitope_score"
    combined_score_expr: "n_alt_reads * target_epitope_score"
    combined_score_expr: "log1p(n_alt_reads) * target_epitope_score
                          + 0.1 * n_alt_reads_supporting_protein_sequence"
    combined_score_expr: "target_epitope_score"

Bindings (read-only scalars from the active VaccinePeptide):

  ``target_epitope_score``
      Sum of per-epitope scores across the VP's mutant ligands
      (after EpitopeConfig filtering / scoring).

  ``self_epitope_score``
      Same metric for the WT-aligned ligands (when present).

  ``expression_score``
      ``sqrt(n_alt_reads)`` — the canonical "expression strength"
      surrogate; available so the legacy
      ``sqrt_reads_times_epitope`` mode is exactly
      ``expression_score * target_epitope_score``.

  ``n_alt_reads``
      RNA reads supporting the variant allele.

  ``n_ref_reads``
      RNA reads supporting the reference allele at the locus.

  ``n_overlapping_reads``
      Total reads spanning the locus (ref + alt + other).

  ``n_alt_reads_supporting_protein_sequence``
      Subset of ``n_alt_reads`` that fully span the assembled
      mutant protein context — the "structurally supports the
      ligand" count.

  ``n_mutant_amino_acids``
      Number of mutated residues in the antigen window.

Functions:

  ``sqrt``, ``log``, ``log1p``, ``exp``, ``min``, ``max``, ``abs``,
  ``pow``.

Operators:

  ``+ - * / // % **`` and unary ``-`` / ``+``.

Anything outside this set raises ``ValueError`` at parse time —
bare attribute access (``something.foo``), subscripts, function
attribute access (``math.sqrt``), comparisons, and statements are
all rejected. The whole expression must be a single
``ast.Expression``.

Same precedence relationship as the per-epitope DSL: when
``vaccine_peptides.combined_score_expr`` is set, it supersedes
``combined_score_mode``. The three legacy modes are equivalent to
these expressions:

  ``epitope_only``               ≡ ``target_epitope_score``
  ``reads_times_epitope``        ≡ ``n_alt_reads * target_epitope_score``
  ``sqrt_reads_times_epitope``   ≡ ``expression_score * target_epitope_score``
                                 ≡ ``sqrt(n_alt_reads) * target_epitope_score``
"""

import ast
import logging
import math


logger = logging.getLogger(__name__)


# When evaluation fails, these bindings get shown verbatim in the
# user-facing error message — they're the per-VP score axes most
# useful for diagnosing "why did my expression evaluate weirdly"
# (the read counts and mutant-AA count are useful too, but more
# numerous and less specific; they go in the DEBUG log). Order is
# the order they appear in the preview.
_HEADLINE_BINDINGS = (
    'target_epitope_score',
    'self_epitope_score',
    'expression_score',
)


# Functions exposed to the expression namespace. All real-valued and
# safe to call on user input — no I/O, no side effects.
_FUNCTIONS = {
    'sqrt': math.sqrt,
    'log': math.log,
    'log1p': math.log1p,
    'exp': math.exp,
    'min': min,
    'max': max,
    'abs': abs,
    'pow': pow,
}

# AST node types the parser permits. Every other node type raises at
# parse time so the user gets the rejection at config-load, not at
# scoring time.
#
# Trust boundary: ``combined_score_expr`` reaches us from a YAML
# config the operator wrote — same channel as ``epitopes.score_expr``
# and the rest of the vaxrank knobs. We are NOT defending against a
# hostile user trying to wedge the process; the threat model is
# "operator typo / pasted-in expression that happens to be expensive
# to evaluate." Power and ``pow()`` are both allowed because they're
# legitimately useful for score expressions; an operator who writes
# ``2 ** (10 ** 9)`` will eat their own CPU and notice. If/when the
# DSL ever takes input from an untrusted channel (e.g. a web form),
# revisit: drop ``ast.Pow``, drop ``pow``, and add a per-node count
# / arithmetic-magnitude cap. Until then the simplicity wins.
_ALLOWED_NODES = (
    ast.Expression,
    ast.BinOp, ast.UnaryOp, ast.Constant, ast.Name, ast.Call,
    ast.Add, ast.Sub, ast.Mult, ast.Div, ast.FloorDiv,
    ast.Mod, ast.Pow, ast.USub, ast.UAdd,
    ast.Load,
)


def _validate_ast(tree, expr):
    """Walk ``tree`` and reject anything outside ``_ALLOWED_NODES``.

    Function calls additionally must target a name in
    ``_FUNCTIONS`` (no method calls, no attribute access on bound
    names like ``math.sqrt``).
    """
    for node in ast.walk(tree):
        if not isinstance(node, _ALLOWED_NODES):
            raise ValueError(
                "Disallowed expression element %s in combined_score_expr "
                "%r. Allowed: identifiers, numeric literals, +-*/// %% **, "
                "unary +/-, and calls to %s." % (
                    type(node).__name__, expr,
                    ', '.join(sorted(_FUNCTIONS))))
        if isinstance(node, ast.Call):
            if not isinstance(node.func, ast.Name):
                raise ValueError(
                    "combined_score_expr %r contains a non-name function "
                    "call (e.g. attribute access or method call). Only "
                    "bare-name calls to %s are allowed." % (
                        expr, ', '.join(sorted(_FUNCTIONS))))
            if node.func.id not in _FUNCTIONS:
                raise ValueError(
                    "Unknown function %r in combined_score_expr %r. "
                    "Allowed: %s." % (
                        node.func.id, expr, ', '.join(sorted(_FUNCTIONS))))


def parse_combined_score_expr(expr):
    """Parse + validate ``expr`` once. Returns the AST root for
    later evaluation. Raises ``ValueError`` on syntax errors or
    disallowed constructs."""
    if not isinstance(expr, str) or not expr.strip():
        raise ValueError("combined_score_expr must be a non-empty string")
    try:
        tree = ast.parse(expr, mode='eval')
    except SyntaxError as e:
        raise ValueError(
            "combined_score_expr %r failed to parse: %s" % (expr, e))
    _validate_ast(tree, expr)
    return tree


def _bindings_from_vaccine_peptide(vp):
    """Read-only namespace of the scalars exposed to the DSL.

    Pulls from both the ``VaccinePeptide`` and its
    ``MutantProteinFragment`` so the user writes flat names
    (``n_alt_reads``) without needing to know which struct holds
    each value.
    """
    frag = vp.mutant_protein_fragment
    return {
        'target_epitope_score': float(vp.target_epitope_score),
        'self_epitope_score': float(getattr(
            vp, 'self_epitope_score', 0.0) or 0.0),
        'expression_score': float(vp.expression_score),
        'n_alt_reads': float(frag.n_alt_reads or 0),
        'n_ref_reads': float(frag.n_ref_reads or 0),
        'n_overlapping_reads': float(frag.n_overlapping_reads or 0),
        'n_alt_reads_supporting_protein_sequence': float(
            frag.n_alt_reads_supporting_protein_sequence or 0),
        'n_mutant_amino_acids': float(
            frag.mutant_amino_acid_end_offset
            - frag.mutant_amino_acid_start_offset),
    }


def evaluate_combined_score(expr_or_tree, vaccine_peptide):
    """Score ``vaccine_peptide`` via ``expr_or_tree``.

    ``expr_or_tree`` accepts either a raw string or a pre-parsed
    AST returned by :func:`parse_combined_score_expr`. Pre-parsing
    saves work when scoring many VPs with the same expression.

    Returns a float. Re-raises any evaluation-time error
    (``ZeroDivisionError``, ``ValueError`` from ``log(<=0)``, etc.)
    after wrapping in ``ValueError`` with the expression text for
    context — these tend to come from edge-case inputs (empty
    fragments, zero reads) and the operator needs to know which
    expression triggered them.
    """
    if isinstance(expr_or_tree, str):
        tree = parse_combined_score_expr(expr_or_tree)
        expr_text = expr_or_tree
    else:
        tree = expr_or_tree
        # Re-render for error messages; stays cheap because this
        # only runs on failure.
        expr_text = ast.unparse(tree)
    bindings = _bindings_from_vaccine_peptide(vaccine_peptide)
    namespace = dict(_FUNCTIONS)
    namespace.update(bindings)
    try:
        # Compile the AST and eval against the locked-down
        # namespace. ``__builtins__`` is set to {} so even if a
        # validator bug let an attribute access through, the
        # surface area is empty.
        return float(eval(
            compile(tree, '<combined_score_expr>', 'eval'),
            {'__builtins__': {}},
            namespace,
        ))
    except Exception as e:
        # Opinionated preview: surface the per-VP *score* axes
        # (mutant / wildtype / expression) since those are what
        # users actually look at when a score expression goes
        # sideways. Read counts and mutant-AA count live in the
        # DEBUG dump; per-VP exceptions can fan out across the
        # whole ranked list and a tight preview keeps the stack
        # readable.
        preview = {
            k: bindings[k] for k in _HEADLINE_BINDINGS if k in bindings
        }
        missing = [k for k in _HEADLINE_BINDINGS if k not in bindings]
        if missing:
            # Soft contract: ``_bindings_from_vaccine_peptide`` is
            # supposed to populate every headline binding. If a
            # future refactor renames or drops one, the preview will
            # silently shrink — surface that through a warning so
            # the maintainer notices instead of users wondering why
            # their error message lost a field.
            logger.warning(
                "combined_score_expr error preview is missing "
                "expected headline binding(s) %s; the bindings "
                "contract from _bindings_from_vaccine_peptide may "
                "have changed.", missing)
        other_count = len(bindings) - len(preview)
        suffix = (
            ", ... (%d more)" % other_count if other_count > 0 else "")
        logger.debug(
            "combined_score_expr %r evaluation failure; full "
            "bindings: %s", expr_text, bindings)
        raise ValueError(
            "combined_score_expr %r failed during evaluation: %s. "
            "Bindings preview: %s%s (set vaxrank logging to DEBUG "
            "for the full bindings dict)." % (
                expr_text, e, preview, suffix)) from e
