# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""Unfiltered MHC and site evidence for complete, source-aware contexts."""

from dataclasses import dataclass, field, replace
import mhctools
from mhctools.cleavage import CleavageInput, CleavageModel
from serializable import DataclassSerializable
from topiary import TopiaryPredictor

from .cleavage_inference import audit_pepsickle_inputs, audit_peptidase_inputs
from .cleavage_profile import CleavageProfile
from .allele_validation import parse_prediction_allele, validate_prediction_allele, validate_prediction_frame
from .safety_assessment import SafetyPrediction, prediction_occurrences_from_frame
from .sequence_context import SequenceContext


@dataclass(frozen=True)
class MHCRequest(DataclassSerializable):
    """One explicitly requested model/kind/HLA/length, not observed coverage."""

    kind: str
    predictor_name: str
    predictor_version: str
    allele: str
    peptide_length: int

    def __post_init__(self):
        if any(not isinstance(v, str) or not v for v in (
                self.kind, self.predictor_name, self.predictor_version, self.allele)):
            raise ValueError("MHC request needs complete model/kind/allele identity")
        if type(self.peptide_length) is not int or self.peptide_length < 1:
            raise ValueError("MHC request peptide length must be a positive integer")
        validate_prediction_allele(self.kind, self.allele, "MHC request")
        object.__setattr__(self, "allele", parse_prediction_allele(self.allele))

    @property
    def identity(self):
        return (self.kind, self.predictor_name, self.predictor_version, self.allele)


@dataclass(frozen=True)
class ContextLigand(DataclassSerializable):
    """Every model observation for one exact ligand occurrence."""

    peptide: str
    start: int
    predictions: tuple[SafetyPrediction, ...]

    def __post_init__(self):
        object.__setattr__(self, "predictions", tuple(self.predictions))
        if (not isinstance(self.peptide, str) or not self.peptide
                or type(self.start) is not int or self.start < 0
                or not self.predictions
                or any(not isinstance(p, SafetyPrediction) for p in self.predictions)):
            raise ValueError("Context ligand requires exact coordinates and typed predictions")
        if len({p.identity for p in self.predictions}) != len(self.predictions):
            raise ValueError("Duplicate prediction identity for one context ligand")

    @property
    def end(self):
        return self.start + len(self.peptide)


@dataclass(frozen=True)
class SequenceContextAudit(DataclassSerializable):
    """Exact inference graph, with missing combinations derived from requests.

    This is an evidence inventory, not ranking, historical selection agreement,
    a self-safety verdict, or proof of cleavage/presentation. Individual-source
    self assessments remain available via ``audit_construct_sequence``.
    """

    context: SequenceContext
    profiles: tuple[CleavageProfile, ...] = field(default_factory=tuple)
    mhc_status: str = "unassessed"
    mhc_requests: tuple[MHCRequest, ...] = field(default_factory=tuple)
    ligands: tuple[ContextLigand, ...] = field(default_factory=tuple)
    reason_codes: tuple[str, ...] = ("mhc_predictor_not_requested",)
    error_message: str = ""
    mhc_backend_version: str = ""

    def __post_init__(self):
        for name in ("profiles", "mhc_requests", "ligands", "reason_codes"):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        self.validate()

    def validate(self):
        if not isinstance(self.context, SequenceContext):
            raise ValueError("Context audit requires a typed sequence context")
        self.context.validate()
        if (any(not isinstance(p, CleavageProfile) for p in self.profiles)
                or any(not isinstance(r, MHCRequest) for r in self.mhc_requests)
                or any(not isinstance(ligand, ContextLigand) for ligand in self.ligands)):
            raise ValueError("Context evidence requires typed profiles, requests and ligands")
        if len(set(self.mhc_requests)) != len(self.mhc_requests):
            raise ValueError("Duplicate MHC request")
        if len({(p.model, p.settings) for p in self.profiles}) != len(self.profiles):
            raise ValueError("Duplicate cleavage model/settings in one context")
        for profile in self.profiles:
            sequence_only = dict(profile.settings).get("input_scope") == "sequence_context_only"
            if sequence_only and (self.context.n_term != "unknown" or self.context.c_term != "unknown"
                                  or not self.context.is_sequence_context):
                raise ValueError("Invalid sequence-only processing chemistry assumption")
            if profile.peptide != self.context.cleavage_input(sequence_only=sequence_only):
                raise ValueError("Cleavage evidence belongs to a different context")
            if self.context.prediction_unassessed_reason and profile.sites:
                raise ValueError("Unsupported context cannot claim sequence-model observations")
        if self.mhc_status not in {"unassessed", "failed", "no_predictions", "predictions_returned"}:
            raise ValueError("Unknown context MHC status")
        if bool(self.ligands) != (self.mhc_status == "predictions_returned"):
            raise ValueError("Context MHC status disagrees with ligand inventory")
        if self.mhc_status != "predictions_returned" and not self.reason_codes:
            raise ValueError("Missing MHC evidence requires an explanation")
        if self.mhc_status == "failed" and not self.error_message:
            raise ValueError("Failed MHC inference requires an error")
        if self.mhc_status in {"predictions_returned", "no_predictions"} and not self.mhc_requests:
            raise ValueError("MHC inference must retain explicit requested coverage")
        if self.context.prediction_unassessed_reason and self.ligands:
            raise ValueError("Unsupported context cannot claim sequence-model MHC observations")
        if any(not isinstance(r, str) or not r for r in self.reason_codes):
            raise ValueError("Context reasons must be nonempty text")
        if len({(ligand.peptide, ligand.start) for ligand in self.ligands}) != len(self.ligands):
            raise ValueError("Duplicate context ligand occurrence")
        for ligand in self.ligands:
            if self.context.sequence[ligand.start:ligand.end] != ligand.peptide:
                raise ValueError("MHC ligand does not match its exact context")

    @property
    def coverage(self):
        """Each requested combination, including every missing occurrence offset."""
        observed = {}
        for ligand in self.ligands:
            for prediction in ligand.predictions:
                observed.setdefault((prediction.identity, len(ligand.peptide)), set()).add(ligand.start)
        return tuple((request, tuple(sorted(observed.get(
            (request.identity, request.peptide_length), set()))),
            tuple(start for start in range(max(0, len(self.context.sequence) - request.peptide_length + 1))
                  if start not in observed.get((request.identity, request.peptide_length), set())))
            for request in self.mhc_requests)

    @property
    def unrequested_predictions(self):
        requested = {(r.identity, r.peptide_length) for r in self.mhc_requests}
        return tuple((ligand, prediction) for ligand in self.ligands for prediction in ligand.predictions
                     if (prediction.identity, len(ligand.peptide)) not in requested)

    def overlays(self, profile):
        """Target masks and all returned ligands, never a best-only subset."""
        if profile not in self.profiles:
            raise ValueError("Overlay profile is not part of this context audit")
        targets = tuple((target, profile.interval_evidence(target.start, target.end))
                        for target in self.context.targets)
        ligands = tuple((ligand, tuple(t for t in self.context.targets if t.overlaps(ligand.start, ligand.end)),
                         profile.interval_evidence(ligand.start, ligand.end)) for ligand in self.ligands)
        return {"targets": targets, "ligands": ligands}


def _mhc_inventories(contexts, predictor, requests):
    if predictor is None:
        return tuple(dict(mhc_requests=requests) for _ in contexts)
    if not requests:
        raise ValueError("MHC prediction requires explicit model/HLA/length requests")
    from topiary import __version__ as topiary_version
    backend = f"topiary={topiary_version};mhctools={mhctools.__version__}"
    # Batch unique sequences, then reattach each complete occurrence's source
    # graph. Sharing raw inference never equates different manufacturing records.
    sequences = tuple(dict.fromkeys(c.sequence for c in contexts if not c.prediction_unassessed_reason))
    names = {sequence: f"context_{i}" for i, sequence in enumerate(sequences)}
    frames, batch_error = {}, None
    if sequences:
        try:
            model = predictor if isinstance(predictor, TopiaryPredictor) else TopiaryPredictor(models=[predictor])
            frame = model.predict_from_named_sequences({names[s]: s for s in sequences})
            if frame is None or not frame.columns.is_unique:
                raise ValueError("MHC output is missing or has duplicate columns")
            validate_prediction_frame(frame, "Context MHC output")
            if not frame.empty:
                if not set(frame["source_sequence_name"]) <= set(names.values()):
                    raise ValueError("MHC output includes an unrequested source")
                frames = {s: frame[frame["source_sequence_name"] == names[s]] for s in sequences}
            else:
                frames = {s: frame for s in sequences}
        except Exception as error:
            batch_error = str(error) or type(error).__name__
    results = []
    for context in contexts:
        common = dict(mhc_requests=requests, mhc_backend_version=backend)
        reason = context.prediction_unassessed_reason
        if reason:
            result = dict(mhc_status="unassessed", reason_codes=(reason,))
        elif batch_error:
            result = dict(mhc_status="failed", reason_codes=("mhc_prediction_failed",), error_message=batch_error)
        else:
            try:
                groups, _ = prediction_occurrences_from_frame(
                    frames[context.sequence], context.sequence, expected_source_name=names[context.sequence])
                ligands = tuple(ContextLigand(peptide, start, tuple(group["predictions"]))
                                for (peptide, start), group in sorted(groups.items(), key=lambda item: item[0][1]))
                result = dict(mhc_status="predictions_returned" if ligands else "no_predictions", ligands=ligands,
                              reason_codes=() if ligands else ("no_predictions_emitted",))
            except (TypeError, ValueError, RuntimeError) as error:
                result = dict(mhc_status="failed", reason_codes=("invalid_mhc_output",), error_message=str(error))
        results.append(dict(common, **result))
    return tuple(results)


def audit_sequence_contexts(contexts, *, mhc_predictor=None, mhc_requests=(),
                            pepsickle=False, human_only=False, threshold=.5,
                            peptidase_predictors=()):
    """Inventory whole-context evidence without changing selection or ranking.

    Pepsickle native/product inputs are explicitly sequence-context-only. Enzyme
    rules cannot interpret reconstructed-window edges as molecular termini.
    Final peptide chemistry defaults unknown; free termini require a recorded
    conditional or documented basis. Other chemistry is not silently stripped.
    """
    contexts, requests, enzymes = tuple(contexts), tuple(mhc_requests), tuple(peptidase_predictors)
    if any(not isinstance(c, SequenceContext) for c in contexts):
        raise ValueError("Audit requires typed sequence contexts")
    if any(not isinstance(r, MHCRequest) for r in requests) or len(set(requests)) != len(requests):
        raise ValueError("Audit requires distinct typed MHC requests")
    for context in contexts:
        context.validate()
    for enzyme in enzymes:
        if not isinstance(getattr(enzyme, "model", None), CleavageModel):
            raise ValueError("Enzyme requires canonical mhctools model metadata")
    if not contexts:
        return ()
    profiles = [[] for _ in contexts]
    if pepsickle:
        # No clinical chemistry claim is made for native/translated context.
        sequence_only = tuple(c.is_sequence_context and not c.prediction_unassessed_reason
                              and c.n_term == c.c_term == "unknown" for c in contexts)
        inputs = tuple(c.cleavage_input(sequence_only=only) for c, only in zip(contexts, sequence_only))
        eligible = [i for i, c in enumerate(contexts) if not c.prediction_unassessed_reason]
        returned = audit_pepsickle_inputs([inputs[i] for i in eligible], human_only=human_only, threshold=threshold)
        for i, profile in zip(eligible, returned):
            if sequence_only[i]:
                profile = replace(profile, settings=profile.settings + (("input_scope", "sequence_context_only"),),
                                  reason_codes=profile.reason_codes + ("molecular_termini_not_established",))
            profiles[i].append(profile)
        # An unsupported construct must still show the requested model as missing.
        if len(eligible) != len(contexts):
            from .cleavage_inference import _pepsickle_model_identity
            model, digest, error = _pepsickle_model_identity(human_only)
            for i, c in enumerate(contexts):
                if i not in eligible:
                    profiles[i].append(CleavageProfile(inputs[i], model, "unassessed",
                        reason_codes=(c.prediction_unassessed_reason,), model_asset_sha256=digest,
                        error_message=error, backend_version=mhctools.__version__))
    for enzyme in enzymes:
        by_input, pending = {}, []
        for i, context in enumerate(contexts):
            peptide = context.cleavage_input()
            reason = context.prediction_unassessed_reason
            if reason is None and context.kind == "native_context":
                reason = "native_window_molecular_termini_unestablished"
            if reason is None and context.kind == "final_construct" and context.construct.modality != "peptide":
                reason = "encoded_segment_molecular_termini_unestablished"
            if reason:
                profiles[i].append(CleavageProfile(peptide, enzyme.model, "unassessed",
                    reason_codes=(reason,), backend_version=mhctools.__version__))
            else:
                # Source metadata cannot change the native enzyme prediction.
                # Score each chemistry/sequence once, then restore occurrences.
                key = CleavageInput(peptide.sequence, n_term=peptide.n_term, c_term=peptide.c_term)
                by_input.setdefault(key, None)
                pending.append((i, peptide, key))
        if by_input:
            by_input = dict(zip(by_input, audit_peptidase_inputs(by_input, enzyme)))
            for i, peptide, key in pending:
                profiles[i].append(replace(by_input[key], peptide=peptide))
    inventories = _mhc_inventories(contexts, mhc_predictor, requests)
    return tuple(SequenceContextAudit(context, profiles=tuple(evidence), **inventory)
                 for context, evidence, inventory in zip(contexts, profiles, inventories))
