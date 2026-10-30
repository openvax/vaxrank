"""Explicit, frozen RNA restriction evidence for experimental CTA selection.

Thresholds and panel membership are supplied by the caller. This module does not
calibrate them or infer protein abundance, malignant-cell expression or immunity.
"""

from dataclasses import dataclass, fields
import hashlib
import json
import math
from pathlib import Path
import re

from serializable import DataclassSerializable


def _text(value, name):
    if not isinstance(value, str) or not value.strip() or value != value.strip():
        raise ValueError(f'CTA restriction requires {name}')


def _sha(value, name):
    if not isinstance(value, str) or not re.fullmatch('[0-9a-f]{64}', value):
        raise ValueError(f'CTA restriction requires a SHA-256 for {name}')


def _number(value, name):
    if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
        raise ValueError(f'CTA restriction requires a finite nonnegative {name}')


def _tuple(record, name, cls=None):
    values = getattr(record, name)
    if not isinstance(values, (tuple, list)):
        raise ValueError(f'CTA restriction requires a sequence for {name}')
    values = tuple(values)
    for value in values:
        if cls is None:
            _text(value, name)
        elif not isinstance(value, cls):
            raise ValueError(f'Invalid CTA restriction {name} record')
    if len(set(values)) != len(values):
        raise ValueError(f'Duplicate CTA restriction {name}')
    object.__setattr__(record, name, values)


@dataclass(frozen=True)
class CTARestrictionPolicy(DataclassSerializable):
    """A caller's versioned panel and bounds; no taxon-dependent thresholds."""

    policy_id: str
    policy_version: str
    rna_unit: str
    rna_scale: str
    max_somatic_expression: float
    min_reproductive_expression: float
    required_somatic_tissues: tuple[str, ...]
    required_reproductive_tissues: tuple[str, ...]
    min_independent_donors: int
    min_independent_samples: int
    rationale: str
    schema_version: int = 1

    def __post_init__(self):
        if type(self.schema_version) is not int or self.schema_version != 1:
            raise ValueError('Unsupported CTA restriction policy schema')
        for name in ('policy_id', 'policy_version', 'rna_unit', 'rna_scale', 'rationale'):
            _text(getattr(self, name), name)
        for name in ('max_somatic_expression', 'min_reproductive_expression'):
            _number(getattr(self, name), name)
        if self.min_reproductive_expression <= self.max_somatic_expression:
            raise ValueError('Reproductive RNA lower bound must exceed the somatic upper bound')
        for name in ('min_independent_donors', 'min_independent_samples'):
            if type(getattr(self, name)) is not int or getattr(self, name) < 1:
                raise ValueError(f'CTA restriction requires a positive {name}')
        for name in ('required_somatic_tissues', 'required_reproductive_tissues'):
            _tuple(self, name)
            if not getattr(self, name):
                raise ValueError(f'CTA restriction requires an explicit {name}')
        if set(self.required_somatic_tissues) & set(self.required_reproductive_tissues):
            raise ValueError('CTA somatic and reproductive panels overlap')


@dataclass(frozen=True)
class CTARestrictionOccurrence(DataclassSerializable):
    """One full-protein occurrence, including noncanonical host ORFs."""

    gene_id: str
    transcript_id: str | None
    entry_id: str
    sequence_sha256: str
    kind: str = 'annotated'

    def __post_init__(self):
        for name in ('gene_id', 'entry_id'):
            _text(getattr(self, name), name)
        if self.transcript_id is not None:
            _text(self.transcript_id, 'transcript_id')
        _sha(self.sequence_sha256, 'occurrence sequence')
        if self.kind not in ('annotated', 'extra_orf'):
            raise ValueError('CTA restriction occurrence must belong to the host')


@dataclass(frozen=True)
class CTARestrictionObservation(DataclassSerializable):
    """Native RNA bounds and explicit biological identities, never donor counts.

Sample/donor IDs are globally reconciled within the supplied bundle. Repeated
assays use the same IDs. Pools retain all donor IDs but cannot establish an
independent negative observation for each donor.
"""

    evidence_id: str
    occurrences: tuple[CTARestrictionOccurrence, ...]
    tissue: str
    tissue_role: str
    feature_level: str
    rna_unit: str
    rna_scale: str
    lower_bound: float | None
    upper_bound: float | None
    allocation: str
    allocation_provenance: str
    sample_id: str
    donor_ids: tuple[str, ...]
    pooled: bool
    source_id: str
    source_version: str
    source_sha256: str
    source_locator: str
    lineage_ids: tuple[str, ...]
    exclusion_reason: str = ''

    def __post_init__(self):
        for name in ('evidence_id', 'tissue', 'rna_unit', 'rna_scale', 'allocation_provenance',
                     'sample_id', 'source_id', 'source_version', 'source_locator'):
            _text(getattr(self, name), name)
        _tuple(self, 'occurrences', CTARestrictionOccurrence)
        _tuple(self, 'donor_ids')
        _tuple(self, 'lineage_ids')
        if not self.occurrences or not self.lineage_ids:
            raise ValueError('CTA restriction observations require targets and study lineage')
        if self.tissue_role not in ('somatic', 'reproductive'):
            raise ValueError('CTA restriction RNA must declare a normal tissue role')
        if self.feature_level not in ('gene', 'transcript'):
            raise ValueError('CTA restriction requires RNA feature resolution')
        if self.allocation not in ('independent', 'group', 'unknown'):
            raise ValueError('Unknown CTA RNA allocation')
        if self.allocation == 'independent' and len(self.occurrences) != 1:
            raise ValueError('Independent CTA RNA bounds require one occurrence')
        if type(self.pooled) is not bool or (not self.pooled and len(self.donor_ids) > 1):
            raise ValueError('CTA donor identities disagree with pooling')
        if not isinstance(self.exclusion_reason, str):
            raise ValueError('CTA exclusion reason must be text')
        _sha(self.source_sha256, 'RNA source')
        for name in ('lower_bound', 'upper_bound'):
            if getattr(self, name) is not None:
                _number(getattr(self, name), name)
        if (self.lower_bound is not None and self.upper_bound is not None
                and self.lower_bound > self.upper_bound):
            raise ValueError('CTA RNA lower bound exceeds its upper bound')


@dataclass(frozen=True)
class CTARestrictionEvidence(DataclassSerializable):
    """Reference-bound supplemental evidence, retained alongside OncoRef facts."""

    taxon_id: int
    reference_id: str
    host_reference_sha256: str
    snapshot_id: str
    observations: tuple[CTARestrictionObservation, ...]
    schema_version: int = 1

    def __post_init__(self):
        if type(self.schema_version) is not int or self.schema_version != 1:
            raise ValueError('Unsupported CTA restriction evidence schema')
        if type(self.taxon_id) is not int or self.taxon_id < 1:
            raise ValueError('CTA restriction requires a taxon')
        for name, prefix in (('reference_id', 'openvax:reference:v1:sha256:'),
                             ('snapshot_id', 'openvax:cta-evidence:v2:sha256:')):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.startswith(prefix):
                raise ValueError(f'CTA restriction requires a versioned {name}')
            _sha(value[len(prefix):], name)
        _sha(self.host_reference_sha256, 'host reference')
        _tuple(self, 'observations', CTARestrictionObservation)
        if len({r.evidence_id for r in self.observations}) != len(self.observations):
            raise ValueError('Duplicate CTA restriction evidence IDs')
        samples = {}
        for row in self.observations:
            identity = (tuple(sorted(row.donor_ids)), row.pooled, row.tissue, row.tissue_role)
            if row.sample_id in samples and samples[row.sample_id] != identity:
                raise ValueError('Conflicting CTA sample/donor identities')
            samples[row.sample_id] = identity

    @classmethod
    def read(cls, path):
        """Read an ordinary JSON bundle with a closed field vocabulary."""
        def construct(kind, record):
            if not isinstance(record, dict) or set(record) - {f.name for f in fields(kind)}:
                raise ValueError(f'Unknown {kind.__name__} fields')
            try:
                return kind(**record)
            except TypeError as error:
                raise ValueError(f'Invalid {kind.__name__} record') from error
        data = json.loads(Path(path).read_text())
        if not isinstance(data, dict) or not isinstance(data.get('observations'), list):
            raise ValueError('CTA restriction evidence requires observations')
        for row in data['observations']:
            if not isinstance(row, dict) or not isinstance(row.get('occurrences'), list):
                raise ValueError('CTA restriction observation requires occurrences')
            row['occurrences'] = tuple(construct(CTARestrictionOccurrence, x) for x in row['occurrences'])
        data['observations'] = tuple(construct(CTARestrictionObservation, x) for x in data['observations'])
        return construct(cls, data)


@dataclass(frozen=True)
class CTARestrictionFinding(DataclassSerializable):
    code: str
    subject: str
    evidence_id: str = ''
    tissue: str = ''

    def __post_init__(self):
        for name in ('code', 'subject'):
            _text(getattr(self, name), name)
        if not isinstance(self.evidence_id, str) or not isinstance(self.tissue, str):
            raise ValueError('CTA restriction finding identities must be text')


@dataclass(frozen=True)
class CTARestrictionAssessment(DataclassSerializable):
    policy: CTARestrictionPolicy
    evidence: CTARestrictionEvidence | None
    occurrences: tuple[CTARestrictionOccurrence, ...]
    findings: tuple[CTARestrictionFinding, ...]
    native_evidence_json: str

    def __post_init__(self):
        if not isinstance(self.policy, CTARestrictionPolicy):
            raise ValueError('CTA restriction assessment requires its frozen policy')
        if self.evidence is not None and not isinstance(self.evidence, CTARestrictionEvidence):
            raise ValueError('CTA restriction assessment requires typed evidence')
        _tuple(self, 'occurrences', CTARestrictionOccurrence)
        _tuple(self, 'findings', CTARestrictionFinding)
        if not self.occurrences:
            raise ValueError('CTA restriction assessment requires host occurrences')
        native = json.loads(self.native_evidence_json)
        if (not isinstance(native, dict) or set(native) != {
                'snapshot_id', 'taxon_id', 'source_version', 'references', 'candidates',
                'evidence', 'facts', 'contracts', 'sources'}):
            raise ValueError('CTA restriction assessment requires its detached native facts')
        if self.evidence is not None and (self.evidence.snapshot_id != native['snapshot_id']
                                         or self.evidence.taxon_id != native['taxon_id']):
            raise ValueError('CTA restriction native and supplemental identities disagree')
        if self.findings != _restriction_findings(self.policy, self.evidence, self.occurrences, native):
            raise ValueError('CTA restriction findings disagree with their frozen evidence')

    @property
    def passed(self):
        return not self.findings


def host_occurrences(host_reference, sequence=None):
    """Enumerate full sequence equivalents; partial sharing is audited downstream."""
    return tuple(dict.fromkeys(CTARestrictionOccurrence(
        o.gene_id, o.transcript_id or None, e.accession, hashlib.sha256(e.sequence.encode()).hexdigest(), e.kind)
        for e in host_reference.proteome if e.kind in ('annotated', 'extra_orf')
        and (sequence is None or e.sequence == sequence) for o in e.occurrences))


def validate_restriction_evidence(evidence, host_reference, snapshot, reference):
    """Reject contradictory inputs before sequence selection or predictors."""
    if evidence is None:
        return
    if (evidence.taxon_id != snapshot.taxon_id
            or evidence.reference_id != reference['reference_id']
            or evidence.host_reference_sha256 != host_reference.identity
            or evidence.snapshot_id != snapshot.snapshot_id):
        raise ValueError('CTA restriction evidence differs from its host/reference/snapshot')
    available = set(host_occurrences(host_reference))
    if any(o not in available for r in evidence.observations for o in r.occurrences):
        raise ValueError('CTA restriction occurrence or sequence differs from its host inventory')


def assess_restriction(policy, evidence, occurrences, snapshot):
    """Evaluate all occurrences without summing pooled or repeated measurements."""
    genes = {o.gene_id for o in occurrences}
    candidates = [c for c in snapshot.candidates() if genes.intersection(c['gene_ids'])]
    candidate_ids = {c['id'] for c in candidates}
    record = snapshot.as_dict()
    parents = [e for e in record['evidence'] if candidate_ids.intersection(e['candidate_ids'])]
    parent_ids = {e['id'] for e in parents}
    contract_ids = {e.get('contract_id') for e in parents}
    source_ids = {e.get('source_id') for e in parents}
    native = dict(snapshot_id=snapshot.snapshot_id, taxon_id=snapshot.taxon_id,
        source_version=record.get('source_version', 'not reported'),
        references=record.get('references', []), candidates=candidates, evidence=parents,
        facts=[f for f in record.get('facts', ()) if f['evidence_id'] in parent_ids],
        contracts=[c for c in record.get('contracts', ()) if c['contract_id'] in contract_ids],
        sources=[s for s in record.get('sources', ()) if s['id'] in source_ids])
    findings = _restriction_findings(policy, evidence, occurrences, native)
    return CTARestrictionAssessment(policy, evidence, tuple(occurrences), findings,
                                    json.dumps(native, sort_keys=True, separators=(',', ':')))


def _restriction_findings(policy, evidence, occurrences, native):
    findings = []
    def add(code, subject, row=None, tissue=''):
        finding = CTARestrictionFinding(code, subject, row.evidence_id if row else '', tissue)
        if finding not in findings:
            findings.append(finding)

    if evidence is None:
        add('missing_restriction_evidence', native['snapshot_id'])
    for occurrence in occurrences:
        subject = f'{occurrence.gene_id}:{occurrence.transcript_id}:{occurrence.entry_id}'
        if occurrence.transcript_id is None:
            add('unresolved_host_transcript', subject)
        rows = [r for r in evidence.observations if occurrence in r.occurrences] if evidence else []
        for row in rows:
            compatible = (row.rna_unit, row.rna_scale) == (policy.rna_unit, policy.rna_scale)
            if ((row.tissue in policy.required_somatic_tissues and row.tissue_role != 'somatic')
                    or (row.tissue in policy.required_reproductive_tissues and row.tissue_role != 'reproductive')):
                raise ValueError('CTA observation tissue role conflicts with the declared panel')
            # Exclusion cannot erase observed somatic counterevidence. Group RNA
            # and incompatible scales remain unresolved, rather than numeric zero.
            if row.tissue_role == 'somatic':
                if not compatible or row.allocation != 'independent' or row.feature_level != 'transcript':
                    add('unresolved_somatic_rna', subject, row, row.tissue)
                elif row.lower_bound is not None and row.lower_bound > policy.max_somatic_expression:
                    add('somatic_rna_above_bound', subject, row, row.tissue)
                elif row.upper_bound is None or row.upper_bound > policy.max_somatic_expression:
                    add('unresolved_somatic_rna', subject, row, row.tissue)
        for role, tissues in (('somatic', policy.required_somatic_tissues),
                              ('reproductive', policy.required_reproductive_tissues)):
            for tissue in tissues:
                donors, samples = set(), set()
                for row in rows:
                    if (row.tissue, row.tissue_role) != (tissue, role):
                        continue
                    if (row.exclusion_reason or row.pooled or len(row.donor_ids) != 1
                            or row.feature_level != 'transcript' or row.allocation != 'independent'
                            or (row.rna_unit, row.rna_scale) != (policy.rna_unit, policy.rna_scale)):
                        continue
                    passes = (row.upper_bound is not None and row.upper_bound <= policy.max_somatic_expression
                              if role == 'somatic' else row.lower_bound is not None
                              and row.lower_bound >= policy.min_reproductive_expression)
                    if passes:
                        donors.update(row.donor_ids)
                        samples.add(row.sample_id)
                if len(donors) < policy.min_independent_donors:
                    add('insufficient_independent_donors', subject, tissue=tissue)
                if len(samples) < policy.min_independent_samples:
                    add('insufficient_independent_samples', subject, tissue=tissue)

    # The supplemental panel cannot replace contrary/uncertain native facts from
    # any sequence-equivalent candidate locus. Historical summary fields retain
    # their original contract; they are not independent normal observations.
    contracts = {c['contract_id']: c for c in native['contracts']}
    parents = {p['id']: p for p in native['evidence']}
    for fact in native['facts']:
        if fact['tissue_role'] != 'somatic' or fact['observation']['status'] == 'not_reported':
            continue
        # Unquantified negatives do not supply coverage. A native quantitative
        # RNA fact can be compared only with the same reference and measurement
        # contract. It still cannot supply independent transcript/donor coverage:
        # the supplemental bundle must retain those allocations and identities.
        if fact['observation']['status'] == 'no_amplification':
            continue
        contract = contracts.get(parents[fact['evidence_id']].get('contract_id'), {})
        measurement = fact.get('measurement', {})
        value = fact['observation'].get('value')
        comparable = (evidence is not None and fact.get('statistic') == 'native_observation'
            and contract.get('context') == 'normal' and contract.get('reference_status') == 'verified'
            and contract.get('taxon_id') == evidence.taxon_id
            and (contract.get('reference') or {}).get('reference_id') == evidence.reference_id
            and measurement.get('feature_level') == 'transcript'
            and (measurement.get('unit'), measurement.get('scale')) == (policy.rna_unit, policy.rna_scale)
            and measurement.get('is_proxy') is False and measurement.get('transformations') == []
            and fact['observation']['status'] in ('reported', 'below_limit')
            and type(value) in (int, float) and math.isfinite(value) and value >= 0)
        if comparable and value <= policy.max_somatic_expression:
            continue
        code = ('somatic_fact_above_bound' if comparable and fact['observation']['status'] == 'reported'
                else 'somatic_fact_requires_review')
        finding = CTARestrictionFinding(code, fact['id'],
                                        fact['evidence_id'], fact['tissue'] or '')
        if finding not in findings:
            findings.append(finding)
    return tuple(findings)
