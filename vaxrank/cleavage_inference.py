# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Non-mutating full-context processing through existing mhctools backends."""

import hashlib
from importlib import metadata, util
import math
from numbers import Real
from pathlib import Path

import mhctools
from mhctools.cleavage import CleavageInput, CleavageModel

from .cleavage_profile import (
    CleavageProfile, profile_from_cleavage_result, profile_from_residue_scores,
)


def _inputs(peptides):
    peptides = tuple(peptides)
    if any(not isinstance(p, CleavageInput) for p in peptides):
        raise ValueError("Processing audit requires explicit mhctools CleavageInput records")
    return peptides


def _pepsickle_model_identity(human_only):
    """Fingerprint the weights actually found on this interpreter's import path."""
    identity_error = ""
    try:
        spec = util.find_spec("pepsickle")
    except (ImportError, ValueError) as error:
        spec, identity_error = None, str(error)
    digest = ""
    version = "unknown"
    if spec is not None:
        try:
            version = metadata.version("pepsickle")
        except metadata.PackageNotFoundError:
            pass
        if spec.origin:
            path = Path(spec.origin).parent / "trained_model_dict.pickle"
            try:
                if path.is_file():
                    hasher = hashlib.sha256()
                    with path.open("rb") as handle:
                        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                            hasher.update(chunk)
                    digest = hasher.hexdigest()
            except OSError as error:
                identity_error = str(error)
    model = CleavageModel(
        name="pepsickle-epitope-human" if human_only else "pepsickle-epitope-mammal",
        version=version, enzyme="proteasome", uniprot="", species="",
        compartments=("proteasomal",), evidence="quantitative_model",
        references=("https://doi.org/10.1093/bioinformatics/btab628",
                    "https://github.com/pdxgx/pepsickle"),
        assay="in-vivo epitope-trained sequence model",
        limitations="Proteasome-type agnostic. Missing terminal context is padded; "
                    "not proof of cleavage, presentation, serum stability or clinical safety. "
                    "Human-only training is experimental and uses fewer observations.",
        score_name="cleavage_probability", score_units="dimensionless",
        scored_endpoint="site_cleavage")
    return model, digest, identity_error


def audit_pepsickle_inputs(peptides, *, human_only=False, threshold=.5):
    """One isolated mhctools batch per set of distinct, unmodified sequences.

    Return one immutable profile per input occurrence, preserving original source
    offsets and chemistry. Missing installation/assets and invalid outputs stay
    explicitly unassessed/failed. No ranking scores are read or changed.
    """
    peptides = _inputs(peptides)
    if (type(human_only) is not bool or isinstance(threshold, bool)
            or not isinstance(threshold, Real) or not math.isfinite(threshold)
            or not 0 <= threshold <= 1):
        raise ValueError("Invalid Pepsickle model settings")
    if not peptides:
        return ()
    model, digest, identity_error = _pepsickle_model_identity(human_only)
    provenance = dict(
        backend_version=mhctools.__version__, model_asset_sha256=digest,
        settings=(("human_only", str(human_only)), ("threshold", str(threshold)),
                  ("model_type", "epitope"), ("proteasome_type", "agnostic"),
                  ("upstream_context", "8"), ("downstream_context", "8"),
                  ("isolate_subprocess", "True")))
    eligible = tuple(dict.fromkeys(p.sequence for p in peptides
                                  if p.n_term == p.c_term == "free"))
    outputs, batch_error = {}, None
    if eligible and digest:
        try:
            predictor = mhctools.Pepsickle(
                human_only=human_only, threshold=threshold, isolate_subprocess=True)
            outputs = predictor.cleavage_probs_many(eligible)
            if not isinstance(outputs, dict) or set(outputs) != set(eligible):
                raise ValueError("Processing batch did not return exactly the requested contexts")
        except Exception as error:
            batch_error = str(error) or type(error).__name__
    profiles = []
    for peptide in peptides:
        if peptide.n_term != "free" or peptide.c_term != "free":
            profile = CleavageProfile(
                peptide, model, "unassessed", reason_codes=("unsupported_terminal_chemistry",),
                **provenance)
        elif not digest:
            profile = CleavageProfile(
                peptide, model, "unassessed", reason_codes=("pepsickle_model_assets_unavailable",),
                error_message=identity_error,
                **provenance)
        elif batch_error is not None:
            profile = CleavageProfile(
                peptide, model, "failed", reason_codes=("processing_prediction_failed",),
                error_message=batch_error, **provenance)
        else:
            try:
                profile = profile_from_residue_scores(
                    peptide, model, outputs[peptide.sequence],
                    upstream_context=8, downstream_context=8, **provenance)
            except (TypeError, ValueError) as error:
                profile = CleavageProfile(
                    peptide, model, "failed", reason_codes=("invalid_processing_output",),
                    error_message=str(error), **provenance)
        profiles.append(profile)
    return tuple(profiles)


def audit_peptidase_inputs(peptides, predictor):
    """Retain each enzyme model's native bond evidence without filling gaps.

    A recognition rule's non-match is not resistance; a native kinetic score
    is not a cleavage probability or whole-peptide half-life.
    """
    peptides = _inputs(peptides)
    model = getattr(predictor, "model", None)
    if not isinstance(model, CleavageModel):
        raise ValueError("Peptidase predictor requires canonical mhctools model metadata")
    by_input = {}
    for peptide in dict.fromkeys(peptides):
        try:
            result = predictor.predict(peptide)
            if result.peptide != peptide or result.model != model:
                raise ValueError("Peptidase result changed the requested input or model identity")
            by_input[peptide] = profile_from_cleavage_result(
                result, backend_version=mhctools.__version__)
        except Exception as error:
            by_input[peptide] = CleavageProfile(
                peptide, model, "failed", reason_codes=("peptidase_prediction_failed",),
                error_message=str(error) or type(error).__name__,
                backend_version=mhctools.__version__)
    return tuple(by_input[peptide] for peptide in peptides)
