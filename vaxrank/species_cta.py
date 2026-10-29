"""Experimental species-scoped CTA expression and frozen OncoRef evidence."""

from dataclasses import dataclass, replace
import hashlib
import json
import re

from serializable import DataclassSerializable

from .cta_admission import CTAAdmissionPolicy, CTAOverrideEvidence, PatientTumorExpressionEvidence
from .cta_reference_adapter import feature_id, reference_namespace, require_common_cta_data, validate_host_mapping
from .cta_restriction import (
    CTARestrictionAssessment, CTARestrictionEvidence, CTARestrictionPolicy,
    assess_restriction, host_occurrences, validate_restriction_evidence,
)
from .amino_acids import STANDARD_AMINO_ACIDS
from .cta_expression import (
    CTAExpressionDecision, CTAExpressionInput, CTAExpressionResult,
    CTATargetExclusionPolicy, _expression_rows,
)
from .vaccine_antigen import (
    AminoAcidInterval, TargetableMask, TumorSpecificityAttestation,
    TumorSpecificityEvidence, VaccineAntigen,
)


@dataclass(frozen=True)
class SpeciesCTAExpressionInput(CTAExpressionInput):
    """Reference-native gene or exact versioned RNA IDs, never gene symbols."""

    namespace: str = 'ncbi'

    def __post_init__(self):
        super().__post_init__()
        if self.namespace not in ('ncbi', 'ensembl'):
            raise ValueError('Unsupported CTA expression namespace')

    @property
    def id_namespace(self):
        if self.namespace == 'ensembl':
            return f'ensembl_{self.measurement_level}_id'
        return 'ncbi_gene_id' if self.measurement_level == 'gene' else 'refseq_transcript_id'


@dataclass(frozen=True)
class SpeciesCTAExpressionDecision(CTAExpressionDecision):
    """A selected reference protein, including sequences unsuitable for prediction."""

    protein_id: str = ''
    protein_sequence_sha256: str = ''
    unsupported_amino_acids: str = ''

    def __post_init__(self):
        super().__post_init__()
        if bool(self.protein_id) != bool(self.protein_sequence_sha256):
            raise ValueError('Species CTA protein selection requires its sequence checksum')
        if self.protein_id and not re.fullmatch('[0-9a-f]{64}', self.protein_sequence_sha256):
            raise ValueError('Species CTA protein selection requires a SHA-256')
        if (bool(self.unsupported_amino_acids) != (self.status == 'unsupported_protein_sequence')
                or (self.unsupported_amino_acids and (not self.protein_id
                    or any(a in STANDARD_AMINO_ACIDS for a in self.unsupported_amino_acids)))):
            raise ValueError('Species CTA unresolved sequence status disagrees with its residues')


@dataclass(frozen=True)
class SpeciesCTAAdmission(DataclassSerializable):
    """A reviewed research target with independent patient-expression evidence."""

    antigen: VaccineAntigen
    policy: CTAAdmissionPolicy
    tumor_expression: PatientTumorExpressionEvidence
    candidate_ids: tuple[str, ...]
    snapshot_id: str
    override_evidence: CTAOverrideEvidence | None = None
    restriction: CTARestrictionAssessment | None = None

    def __post_init__(self):
        expression, antigen = self.tumor_expression, self.antigen
        if (antigen.kind != 'CTA' or expression.gene_id != antigen.gene_id
                or expression.unit != self.policy.expression_unit):
            raise ValueError('Species CTA antigen and expression contracts disagree')
        if not self.candidate_ids or not self.snapshot_id:
            raise ValueError('Species CTA admission requires mapped candidate evidence')
        passes = expression.value >= self.policy.min_tumor_expression
        reviewed = self.override_evidence is not None
        if self.restriction is not None:
            if reviewed or not isinstance(self.restriction, CTARestrictionAssessment):
                raise ValueError('Common CTA restriction cannot be replaced by an override')
            expected = 'admitted' if passes and self.restriction.passed else 'held_out'
        else:
            expected = 'overridden' if passes and reviewed else 'held_out'
        if antigen.tumor_specificity.status != expected:
            raise ValueError('Species CTA selection requires expression and explicit evidence review')
        if (not antigen.tumor_specificity.requires_review or not antigen.tumor_specificity.patient_specific
                or antigen.tumor_specificity.override_reason != (self.override_evidence.rationale
                    if expected == 'overridden' else '')):
            raise ValueError('Species CTA review attestation disagrees with its evidence')
        expected_exclusions = (antigen.gene_id,) if expected in ('admitted', 'overridden') else ()
        if antigen.self_reference_excluded_gene_ids != expected_exclusions:
            raise ValueError('Species CTA self exclusions must contain only the reviewed target')
        if expression.measurement_level == 'transcript' and antigen.transcript_ids != (expression.transcript_id,):
            raise ValueError('Species CTA must select its measured transcript')


@dataclass(frozen=True)
class SpeciesCTAExpressionResult(CTAExpressionResult):
    """Complete detached evidence; restoring it never queries current OncoRef."""

    evidence_snapshot: dict | None = None
    host_reference_sha256: str = ''

    def __post_init__(self):
        from oncoref.species_antigens import CTAEvidenceSnapshot
        from orfdata.reference import species_taxon

        if self.annotation_reference is not None:
            raise ValueError('Species CTA evidence cannot claim the historical human annotation contract')
        snapshot = CTAEvidenceSnapshot.from_dict(self.evidence_snapshot)
        reference, candidates = _snapshot_reference(snapshot)
        namespace = reference_namespace(reference)
        if self.input_contract.namespace != namespace:
            raise ValueError('CTA expression namespace differs from its reference')
        if (reference['assembly_accession'] != self.reference_assembly
                or reference['annotation_source'] != self.annotation_name
                or reference['annotation_release'] != self.annotation_version):
            raise ValueError('Species CTA result and evidence annotation disagree')
        if not re.fullmatch('[0-9a-f]{64}', self.host_reference_sha256):
            raise ValueError('Species CTA result requires its host reference identity')
        by_gene = _unique_candidates(candidates)
        for decision in self.decisions:
            feature_id(decision.input_identifier, self.input_contract.measurement_level, namespace)
            selected = by_gene.get(decision.gene_id, ())
            pairs = {(s['transcript_id'], s['protein_id']) for c in selected
                     for t in c['mapping']['targets'] if t['id'] == decision.gene_id
                     for s in t['sequences']}
            if (isinstance(decision, SpeciesCTAExpressionDecision) and decision.protein_id
                    and (decision.transcript_id, decision.protein_id) not in pairs):
                raise ValueError('Species CTA selected protein disagrees with its frozen mapping')
            assessment = decision.assessment
            if assessment is None:
                if decision.status == 'admitted':
                    raise ValueError('Admitted species CTA requires its assessment')
                continue
            if not isinstance(assessment, SpeciesCTAAdmission):
                raise ValueError('Species CTA result contains a different admission contract')
            antigen = assessment.antigen
            expression = assessment.tumor_expression
            if (assessment.snapshot_id != snapshot.snapshot_id
                    or assessment.candidate_ids != tuple(c['id'] for c in selected)
                    or assessment.policy != self.admission_policy
                    or expression.value != decision.value or expression.input_sha256 != self.input_sha256
                    or expression.input_identifier != decision.input_identifier
                    or expression.sample_id != self.input_contract.sample_id
                    or expression.measurement_level != self.input_contract.measurement_level
                    or expression.evidence_source != self.input_contract.evidence_source
                    or expression.evidence_version != self.input_contract.evidence_version
                    or expression.assay != self.input_contract.assay
                    or decision.gene_identity != dict(source_gene_id=decision.gene_id,
                        canonical_gene_id=decision.gene_id, taxon_id=snapshot.taxon_id,
                        status='species_reference', reference_id=reference['reference_id'])
                    or antigen.gene_id != decision.gene_id
                    or species_taxon(antigen.species) != snapshot.taxon_id
                    or len(antigen.transcript_ids) != 1 or len(antigen.protein_ids) != 1
                    or (decision.transcript_id, antigen.protein_ids[0]) not in pairs
                    or (isinstance(decision, SpeciesCTAExpressionDecision) and (
                        decision.protein_id != antigen.protein_ids[0]
                        or decision.protein_sequence_sha256 != hashlib.sha256(antigen.amino_acids.encode()).hexdigest()))
                    or antigen.transcript_ids != (decision.transcript_id,)
                    or decision.status != ('admitted' if antigen.tumor_specificity.admits_construct else 'held_out')):
                raise ValueError('Species CTA decision disagrees with frozen mapping or expression evidence')
            metadata = dict(antigen.source_metadata)
            if (metadata.get('host_reference_sha256') != self.host_reference_sha256
                    or metadata.get('protein_sequence_sha256') != hashlib.sha256(antigen.amino_acids.encode()).hexdigest()):
                raise ValueError('Species CTA sequence differs from its saved reference evidence')
            if (not antigen.tumor_specificity.evidence_records
                    or antigen.tumor_specificity.evidence_records[0] != _reference_evidence(
                        snapshot, selected, decision.gene_id, reference)):
                raise ValueError('Species CTA antigen lost or changed its candidate observations')
            restriction = assessment.restriction
            if restriction is not None:
                evidence = restriction.evidence
                if evidence is not None and (
                        evidence.snapshot_id != snapshot.snapshot_id
                        or evidence.reference_id != reference['reference_id']
                        or evidence.host_reference_sha256 != self.host_reference_sha256
                        or evidence.taxon_id != snapshot.taxon_id):
                    raise ValueError('Saved CTA restriction reference identity differs')
                if (any(o.sequence_sha256 != decision.protein_sequence_sha256 for o in restriction.occurrences)
                        or not any((o.gene_id, o.transcript_id, o.entry_id) ==
                            (antigen.gene_id, decision.transcript_id, decision.protein_id)
                            for o in restriction.occurrences)
                        or restriction != assess_restriction(restriction.policy, evidence,
                            restriction.occurrences, snapshot)
                        or _restriction_record(restriction, antigen.gene_id) not in
                            antigen.tumor_specificity.evidence_records
                        or metadata.get('cta_restriction_sha256') != _restriction_digest(restriction)):
                    raise ValueError('Saved CTA restriction assessment disagrees with its frozen evidence')


def _snapshot_reference(snapshot):
    if snapshot.coverage['status'] == 'unavailable':
        raise ValueError(f'OncoRef CTA evidence is unavailable for species {snapshot.taxon_id}')
    references = snapshot.as_dict()['references']
    if len(references) != 1:
        raise ValueError('Species CTA input requires one explicit mapped reference')
    reference = references[0]['identity']
    reference_namespace(reference)
    return reference, snapshot.candidates(reference=references[0]['id'])


def _unique_candidates(candidates):
    eligible = {c['gene_ids'][0] for c in candidates
                if c['mapping']['status'] == 'unique' and c['mapping']['complete']}
    # Eligibility requires a resolved mapping; retain family observations too.
    return {gene: [c for c in candidates if gene in c['gene_ids']] for gene in eligible}


def _reference_evidence(snapshot, candidates, gene, reference):
    """Keep relevant raw facts on the antigen as well as the full run snapshot."""
    identifiers = {candidate['id'] for candidate in candidates}
    observations = [row for row in snapshot.evidence() if identifiers.intersection(row['candidate_ids'])]
    def encoded(value):
        return json.dumps(value, sort_keys=True, separators=(',', ':'))
    return TumorSpecificityEvidence(
        evidence_kind='species_cta_candidate', evidence_source='OncoRef',
        evidence_version=snapshot.as_dict()['source_version'], subject_id=gene,
        details=(('snapshot_id', snapshot.snapshot_id), ('reference_id', reference['reference_id']),
                 ('candidate_ids', ','.join(c['id'] for c in candidates)),
                 ('coverage', encoded(snapshot.coverage)),
                 ('candidate_records', encoded(candidates)), ('observations', encoded(observations))))


def _restriction_record(assessment, gene):
    return TumorSpecificityEvidence(
        evidence_kind='cta_rna_restriction', evidence_source=assessment.policy.policy_id,
        evidence_version=assessment.policy.policy_version, subject_id=gene,
        passed=assessment.passed, details=(('assessment', _restriction_json(assessment)),))


def _restriction_json(assessment):
    from .native_serialization import to_native_json
    return json.dumps(json.loads(to_native_json(assessment)), sort_keys=True, separators=(',', ':'))


def _restriction_digest(assessment):
    return hashlib.sha256(_restriction_json(assessment).encode()).hexdigest()


def cta_overrides(options):
    """Parse target-specific evidence; no run-wide switch admits all candidates."""
    records = options.get('overrides', {})
    if not isinstance(records, dict):
        raise ValueError('CTA overrides must map gene IDs to versioned evidence')
    result = {}
    for gene, record in records.items():
        if not isinstance(gene, str) or not gene or not isinstance(record, dict):
            raise ValueError('CTA overrides require string gene IDs and evidence records')
        if any(not isinstance(value, str) or (key != 'evidence_id' and not value.strip())
               for key, value in record.items()):
            raise ValueError('CTA override evidence fields must be nonempty strings')
        try:
            result[gene] = CTAOverrideEvidence(**record)
        except TypeError as error:
            raise ValueError('Invalid CTA override evidence fields') from error
    return result


def admit_species_cta_expression(path, *, input_contract, admission_policy, host_reference,
                                 exclusion_policy=None, overrides=None, evidence_snapshot=None,
                                 restriction_policy=None, restriction_evidence=None):
    """Use mapped host sequences and retain all OncoRef evidence, including limits."""
    from pathlib import Path
    import pandas as pd
    from oncoref.species_antigens import cta_snapshot

    if restriction_policy is not None:
        require_common_cta_data()
    reference = host_reference.proteome.reference
    snapshot = cta_snapshot(reference.taxon_id, path=evidence_snapshot)
    mapped_reference, candidates = _snapshot_reference(snapshot)
    namespace = validate_host_mapping(host_reference, mapped_reference)
    if restriction_policy is not None:
        if not isinstance(restriction_policy, CTARestrictionPolicy) or overrides:
            raise ValueError('Common CTA restriction requires a typed policy and no overrides')
        if restriction_evidence is not None and not isinstance(restriction_evidence, CTARestrictionEvidence):
            restriction_evidence = CTARestrictionEvidence.read(restriction_evidence)
        validate_restriction_evidence(restriction_evidence, host_reference, snapshot, mapped_reference)
    elif restriction_evidence is not None:
        raise ValueError('CTA restriction evidence requires its explicit policy')
    if input_contract.expression_unit != admission_policy.expression_unit:
        raise ValueError('Expression input units do not match the CTA admission policy')
    input_contract = SpeciesCTAExpressionInput(**{**input_contract.to_dict(), 'namespace': namespace})
    exclusions = exclusion_policy or CTATargetExclusionPolicy((), ())
    overrides = overrides or {}
    fingerprint = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    rows = _expression_rows(path, input_contract, normalize=lambda value, level: feature_id(value, level, namespace))
    by_gene = _unique_candidates(candidates)
    ambiguous = {gid for c in candidates if c['mapping']['status'] != 'unique'
                 or not c['mapping']['complete'] for gid in c['gene_ids']}
    proteins = {}
    for entry in host_reference.proteome:
        if entry.kind == 'annotated':
            for occurrence in entry.occurrences:
                proteins.setdefault(occurrence.gene_id, []).append((entry, occurrence))
    decisions, used_overrides = [], set()
    for identifier, feature, abundance in rows:
        value = None if pd.isna(abundance) else float(abundance)
        measured_transcript = feature if input_contract.measurement_level == 'transcript' else ''
        genes = {t['id'] for c in candidates for t in c['mapping']['targets']
                 if any(s['transcript_id'] == measured_transcript for s in t['sequences'])}
        gene = feature if not measured_transcript else next(iter(genes)) if len(genes) == 1 else ''
        decision = SpeciesCTAExpressionDecision(identifier, value, gene_id=gene, transcript_id=measured_transcript)
        selected = by_gene.get(gene, ())
        if not selected:
            status = 'ambiguous_mapping' if gene in ambiguous or len(genes) > 1 else 'not_cta'
            decisions.append(replace(decision, status=status))
            continue
        target = next(t for c in selected for t in c['mapping']['targets'] if t['id'] == gene)
        decision = replace(decision, gene_name=target['symbol'], gene_identity={
            'source_gene_id': gene, 'canonical_gene_id': gene, 'taxon_id': reference.taxon_id,
            'status': 'species_reference', 'reference_id': mapped_reference['reference_id']})
        if gene in overrides:
            used_overrides.add(gene)
        if exclusions.excludes(decision.gene_name):
            decisions.append(replace(decision, status='excluded_target'))
            continue
        if value is None:
            decisions.append(replace(decision, status='unknown_expression'))
            continue
        pairs = {(s['transcript_id'], s['protein_id']) for c in selected
                 for t in c['mapping']['targets'] if t['id'] == gene for s in t['sequences']}
        available = [(e, o) for e, o in proteins.get(gene, ())
                     if (o.transcript_id, e.accession) in pairs
                     and (not measured_transcript or o.transcript_id == measured_transcript)]
        if not available:
            decisions.append(replace(decision, status='no_mapped_protein_sequence'))
            continue
        # Gene RNA selects a reference protein, never an inferred expressed isoform.
        entry, occurrence = min(available, key=lambda pair: (-len(pair[0].sequence), pair[1].transcript_id, pair[0].accession))
        unsupported = ''.join(sorted(set(entry.sequence) - STANDARD_AMINO_ACIDS))
        decision = replace(decision, transcript_id=occurrence.transcript_id, protein_id=entry.accession,
                           protein_sequence_sha256=hashlib.sha256(entry.sequence.encode()).hexdigest())
        if unsupported:
            decisions.append(replace(decision, status='unsupported_protein_sequence',
                                     unsupported_amino_acids=unsupported))
            continue
        expression = PatientTumorExpressionEvidence(
            gene, input_contract.sample_id, value, input_contract.expression_unit,
            input_contract.evidence_source, input_contract.evidence_version, input_contract.assay,
            input_contract.measurement_level, measured_transcript, identifier, fingerprint)
        override = overrides.get(gene)
        passes = value >= admission_policy.min_tumor_expression
        restriction = (assess_restriction(restriction_policy, restriction_evidence,
            host_occurrences(host_reference, entry.sequence), snapshot) if restriction_policy is not None else None)
        admitted = passes and (restriction.passed if restriction is not None else override is not None)
        records = [_reference_evidence(snapshot, selected, gene, mapped_reference),
            TumorSpecificityEvidence(
                evidence_kind='patient_tumor_expression', evidence_source=expression.evidence_source,
                evidence_version=expression.evidence_version, subject_id=gene, sample_id=expression.sample_id,
                patient_specific=True, passed=passes, numeric_value=value, unit=expression.unit,
                threshold=admission_policy.min_tumor_expression,
                details=(('measurement_level', input_contract.measurement_level), ('input_sha256', fingerprint)))]
        if override:
            records.append(TumorSpecificityEvidence(
                evidence_kind='cta_admission_override', evidence_source=override.evidence_source,
                evidence_version=override.evidence_version, subject_id=gene, passed=True,
                details=(('evidence_id', override.evidence_id), ('rationale', override.rationale))))
        if restriction is not None:
            records.append(_restriction_record(restriction, gene))
        antigen = VaccineAntigen(
            kind='CTA', amino_acids=entry.sequence,
            targetable_mask=TargetableMask((AminoAcidInterval(0, len(entry.sequence)),)),
            tumor_specificity=TumorSpecificityAttestation(
                status=('admitted' if restriction is not None else 'overridden') if admitted else 'held_out',
                evidence_kind='species_cta_and_patient_tumor_expression', evidence_source='OncoRef plus patient expression',
                evidence_version=snapshot.as_dict()['source_version'], patient_specific=True,
                rationale_code=('cta_restriction_policy_pass' if restriction is not None else
                    'cta_explicit_evidence_override') if admitted else 'cta_expression_below_threshold' if not passes else
                    'cta_restriction_policy_held' if restriction is not None else 'cta_species_evidence_review_required',
                requires_review=True, override_reason=override.rationale if admitted and override else '',
                evidence_records=tuple(records)),
            self_reference_excluded_gene_ids=(gene,) if admitted else (),
            gene_id=gene, gene_name=decision.gene_name, transcript_ids=(occurrence.transcript_id,),
            protein_ids=(entry.accession,), species=reference.species,
            source_identifier=f'{input_contract.sample_id}:{identifier}',
            source_metadata=(('host_reference_sha256', host_reference.identity),
                             ('protein_sequence_sha256', hashlib.sha256(entry.sequence.encode()).hexdigest()),
                             ('sequence_selection', 'measured_transcript' if measured_transcript else 'annotation_longest_protein'),
                             ('expression_id_namespace', input_contract.id_namespace),
                             ('reference_assembly', mapped_reference['assembly_accession']),
                             ('annotation_name', mapped_reference['annotation_source']),
                             ('annotation_version', mapped_reference['annotation_release'])) +
                            ((('cta_restriction_sha256', _restriction_digest(restriction)),)
                             if restriction is not None else ()))
        assessment = SpeciesCTAAdmission(antigen, admission_policy, expression,
                                        tuple(c['id'] for c in selected), snapshot.snapshot_id, override, restriction)
        decisions.append(replace(decision, status='admitted' if admitted else 'held_out',
                                 transcript_id=occurrence.transcript_id, assessment=assessment))
    if set(overrides) != used_overrides:
        raise ValueError('CTA overrides do not match resolved expression targets: ' + ', '.join(sorted(set(overrides) - used_overrides)))
    if fingerprint != hashlib.sha256(Path(path).read_bytes()).hexdigest():
        raise ValueError('Expression input changed during admission')
    return SpeciesCTAExpressionResult(
        input_contract, admission_policy, exclusions, fingerprint,
        mapped_reference['assembly_accession'], mapped_reference['annotation_source'],
        mapped_reference['annotation_release'], tuple(decisions),
        evidence_snapshot=snapshot.as_dict(), host_reference_sha256=host_reference.identity)
