# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

"""Exact native, manufactured and complete translated contexts for auditing."""

from copy import deepcopy
from dataclasses import dataclass
import hashlib
from typing import Optional

from mhctools.cleavage import CleavageInput
from serializable import DataclassSerializable

from .construct_sequence import ConstructEvidence, ConstructSequence
from .mrna import RNAConstruct, validate_translated_construct
from .vaccine_antigen import VaccineAntigen


@dataclass(frozen=True)
class ContextTarget(DataclassSerializable):
    """A target-mask interval, not an independently documented MHC ligand."""

    source_name: str
    start: int
    end: int

    def __post_init__(self):
        if (not isinstance(self.source_name, str) or not self.source_name
                or type(self.start) is not int or type(self.end) is not int
                or not 0 <= self.start <= self.end):
            raise ValueError("Invalid context target interval")

    def overlaps(self, start, end):
        if self.start == self.end:
            return start < self.start < end
        return start < self.end and self.start < end


@dataclass(frozen=True)
class SequenceContext(DataclassSerializable):
    """One exact input plus its original source graph; no synthetic antigen.

    Native contexts can be cropped reconstruction windows. Complete translation
    is a CDS fact, not evidence of a mature protein's termini or localization.
    Unknown chemistry is the default. A conditional free-terminal analysis needs
    an explicit ``chemistry_basis`` and does not revise manufacture provenance.
    """

    name: str
    kind: str
    sequence: str
    native_antigen: Optional[VaccineAntigen] = None
    construct: Optional[ConstructSequence] = None
    product: Optional[RNAConstruct] = None
    evidence: Optional[ConstructEvidence] = None
    n_term: str = "unknown"
    c_term: str = "unknown"
    chemistry_basis: str = ""

    def __post_init__(self):
        # RNAConstruct is a public mutable assembly record. Snapshot its whole
        # dataclass, not a hand-maintained subset of fields; revalidate on export.
        if self.product is not None:
            object.__setattr__(self, "product", deepcopy(self.product))
        self.validate()

    def validate(self):
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("Sequence context name is required")
        if not isinstance(self.chemistry_basis, str):
            raise ValueError("Context chemistry basis must be text")
        CleavageInput(self.sequence, n_term=self.n_term, c_term=self.c_term)
        if (self.n_term != "unknown" or self.c_term != "unknown") and not self.chemistry_basis:
            raise ValueError("Terminal chemistry requires an explicit evidence/assumption basis")
        if self.evidence is not None and not isinstance(self.evidence, ConstructEvidence):
            raise ValueError("Context evidence must be attributed")
        if self.kind == "native_context":
            if (not isinstance(self.native_antigen, VaccineAntigen)
                    or self.construct is not None or self.product is not None
                    or self.sequence != self.native_antigen.amino_acids):
                raise ValueError("Native context must preserve the complete supplied antigen")
            if self.n_term != "unknown" or self.c_term != "unknown":
                raise ValueError("Native reconstruction windows do not establish molecular termini")
        elif self.kind == "final_construct":
            if (not isinstance(self.construct, ConstructSequence)
                    or self.native_antigen is not None or self.product is not None
                    or self.sequence != self.construct.sequence):
                raise ValueError("Final context must preserve its complete construct")
        elif self.kind == "translated_product":
            if (not isinstance(self.product, RNAConstruct) or self.evidence is None
                    or self.native_antigen is not None or self.construct is not None
                    or self.sequence != self.product.cds_aa):
                raise ValueError("Translated context requires its complete attributed RNA product")
            validate_translated_construct(self.product)
        else:
            raise ValueError("Unknown sequence context kind")

    @property
    def source_id(self):
        # Same bare sequence may share inference, but never source attribution.
        from .native_serialization import to_native_json
        return "sequence-context:" + hashlib.sha256(to_native_json(self).encode()).hexdigest()

    @property
    def targets(self):
        if self.kind == "native_context":
            masks = ((self.name, 0, self.native_antigen.targetable_mask),)
        elif self.kind == "final_construct":
            masks = ((self.construct.name, 0, self.construct.targetable_mask),)
        else:
            masks = tuple((p.construct.name, p.offset, p.construct.targetable_mask)
                          for p in self.product.construct_placements)
        return tuple(ContextTarget(name, shift + interval.start, shift + interval.end)
                     for name, shift, mask in masks for interval in mask.intervals)

    @property
    def prediction_unassessed_reason(self):
        if self.kind == "final_construct":
            return self.construct.sequence_prediction_unassessed_reason
        if self.product is not None and any(
                p.construct.sequence_prediction_unassessed_reason
                for p in self.product.construct_placements):
            return "unresolved_construct_placement"
        return None

    @property
    def target_mapping_status(self):
        if self.prediction_unassessed_reason == "unresolved_native_mapping":
            return "unresolved"
        if self.product is not None:
            if not self.product.construct_placements:
                return "unassessed_no_placements"
            if any(p.construct.mapping_status != "resolved" for p in self.product.construct_placements):
                return "unresolved"
        return "mapped"

    @property
    def is_sequence_context(self):
        return self.kind != "final_construct" or self.construct.modality != "peptide"

    @property
    def source_antigens(self):
        """All occurrence-specific sources; never choose one gene for a product."""
        if self.native_antigen is not None:
            return (self.native_antigen,)
        if self.construct is not None:
            return (self.construct.native_antigen,) if self.construct.native_antigen else ()
        return tuple(p.construct.native_antigen for p in self.product.construct_placements
                     if p.construct.native_antigen is not None)

    def cleavage_input(self, *, sequence_only=False):
        # Pepsickle accepts only bare sequence. Free here is an explicit model
        # input assumption for a context window, NOT exposed peptidase termini.
        n_term, c_term = ("free", "free") if sequence_only else (self.n_term, self.c_term)
        return CleavageInput(self.sequence, n_term=n_term, c_term=c_term, source_id=self.source_id)


def native_sequence_context(antigen, *, name=None):
    return SequenceContext(name or antigen.source_identifier or antigen.gene_name or "native",
                           "native_context", antigen.amino_acids, native_antigen=antigen)


def final_sequence_context(construct, **chemistry):
    return SequenceContext(construct.name, "final_construct", construct.sequence,
                           construct=construct, **chemistry)


def translated_sequence_context(product, evidence, **chemistry):
    return SequenceContext(product.name, "translated_product", product.cds_aa,
                           product=product, evidence=evidence, **chemistry)
