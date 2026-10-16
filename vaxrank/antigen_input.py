"""Experimental sequence/evidence inputs shared by all tumor antigen kinds."""

from dataclasses import asdict, dataclass, field, fields, replace
import hashlib
import json
from pathlib import Path

import msgspec
from topiary import TopiaryResult

from .allele_validation import parse_prediction_allele, validate_prediction_frame
from .cells import text as cell_text
from .epitope_config import EpitopeConfig
from .epitope_dataset import ANTIGEN_COLUMN, EpitopeDataset
from .input_scope import (
    InputProvenance, InputScope, SHARED_FIELDS, combine_scopes,
    scope_from_mapping, scope_dict, provenance_dict, validate_input_scopes, validate_host_taxon,
)
from .native_serialization import _json_loads, to_native_json
from .reference_proteome import self_reference_matches
from .vaccine_antigen import (
    AminoAcidInterval, TargetableMask, TumorSpecificityAttestation,
    TumorSpecificityEvidence, VaccineAntigen,
)


ANTIGEN_DATASET_SCHEMA = 'vaxrank.antigen_dataset.v1'


def _typed_record(cls, value):
    if not isinstance(value, dict):
        raise ValueError(cls.__name__ + ' must be an object')
    unknown = set(value) - {f.name for f in fields(cls)}
    if unknown:
        raise ValueError('Unknown ' + cls.__name__ + ' fields: ' + ', '.join(sorted(unknown)))
    return msgspec.convert(value, type=cls, strict=True)


def _antigen_from_dict(value):
    value = dict(value)
    value['targetable_mask'] = _typed_record(TargetableMask, {
        **value['targetable_mask'],
        'intervals': tuple(_typed_record(AminoAcidInterval, i) for i in value['targetable_mask']['intervals']),
    })
    attestation = dict(value['tumor_specificity'])
    attestation['evidence_records'] = tuple(
        _typed_record(TumorSpecificityEvidence, r) for r in attestation.get('evidence_records', ()))
    value['tumor_specificity'] = _typed_record(TumorSpecificityAttestation, attestation)
    return _typed_record(VaccineAntigen, value)


@dataclass
class AntigenDataset:
    """Named antigen sequences and admission evidence before prediction.

    ``scope`` declares the patient and host reference, independently of the
    source species on each antigen. ``provenance`` retains input history;
    it does not assign measurements from one source to another antigen.
    Admission decisions are supplied evidence, not inferred by this container.
    """

    antigens: dict[str, VaccineAntigen]
    scope: InputScope = field(default_factory=InputScope)
    provenance: tuple[InputProvenance, ...] = ()

    def __post_init__(self):
        self.antigens = dict(self.antigens)
        if any(not isinstance(k, str) or not k.strip() for k in self.antigens):
            raise ValueError('Antigen record identifiers must be nonempty strings')
        if any(not isinstance(a, VaccineAntigen) for a in self.antigens.values()):
            raise TypeError('Antigen records must contain VaccineAntigen objects')
        self.scope = scope_from_mapping(scope_dict(self.scope), 'Antigen dataset')
        self.provenance = tuple(self.provenance)
        for source in self.provenance:
            combine_scopes(self.scope, source.scope, source.path, SHARED_FIELDS)

    def to_dict(self):
        """Return versioned JSON data without Python class metadata."""
        return dict(schema=ANTIGEN_DATASET_SCHEMA, scope=scope_dict(self.scope),
                    antigens={key: asdict(a) for key, a in self.antigens.items()},
                    provenance=[provenance_dict(p) for p in self.provenance])

    def save(self, path):
        """Write the same portable evidence format accepted by the CLI."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2, sort_keys=True, allow_nan=False) + '\n')

    @classmethod
    def load(cls, path):
        """Load a versioned antigen file without imports, prediction or annotation."""
        document = _json_loads(Path(path).read_text())
        if not isinstance(document, dict) or document.get('schema') != ANTIGEN_DATASET_SCHEMA:
            raise ValueError('Expected antigen schema: ' + ANTIGEN_DATASET_SCHEMA)
        if set(document) - {'schema', 'scope', 'antigens', 'provenance'}:
            raise ValueError('Unknown antigen dataset fields')
        if not isinstance(document.get('antigens'), dict):
            raise ValueError('Antigens must be an object keyed by record identifier')
        try:
            provenance = tuple(_typed_record(InputProvenance, {
                **p, 'scope': scope_from_mapping(p['scope'], 'Antigen provenance'),
                'observed_mhc_alleles': tuple(p.get('observed_mhc_alleles', ())),
            }) for p in document.get('provenance', ()))
            return cls({key: _antigen_from_dict(a) for key, a in document['antigens'].items()},
                       scope_from_mapping(document.get('scope', {}), 'Antigen dataset'), provenance)
        except (TypeError, KeyError) as error:
            raise ValueError('Malformed antigen dataset: ' + str(error)) from error

    @classmethod
    def from_epitope_dataset(cls, dataset):
        """Export resolved source antigens, retaining distinct evidence and input history.

        This exports sequences and admission evidence, not prediction values,
        wild-type comparators or resolved design policies. Keep the native
        EpitopeDataset for replay. Unresolved table rows are not invented as
        antigen sequences.
        """
        antigens = {}
        catalogs = dataset.selection.get('antigen_inputs', {})
        for input_id, records in catalogs.items():
            for name, antigen in records.items():
                key = name if len(catalogs) == 1 else input_id + ':' + name
                antigens[key] = _antigen_from_dict(antigen)
        captured = {to_native_json(a) for a in antigens.values()}
        for antigen in dataset.antigens.values():
            payload = to_native_json(antigen)
            if payload in captured:
                continue
            identity = hashlib.sha256(payload.encode()).hexdigest()
            antigens[identity] = antigen
        scope = InputScope()
        for source in dataset.provenance:
            scope = combine_scopes(scope, source.scope, source.path, SHARED_FIELDS)
        return cls(antigens, scope, dataset.provenance)

    def validate_for_prediction(self, genome, *, host_reference=None):
        """Check the common prediction contract before acquiring models."""
        if not self.antigens:
            raise ValueError('Antigen input contains no resolved sequences')
        scope = self.scope
        if host_reference is not None:
            scope = host_reference.validate(scope, genome=genome)
        elif scope.host_reference_sha256 is not None:
            raise ValueError('Antigen input requires its declared host_reference snapshot')
        missing = [name for name in ('patient_id', 'reference_assembly', 'mhc_alleles')
                   if getattr(scope, name) is None]
        if missing:
            raise ValueError('Antigen prediction requires declared ' + ', '.join(missing))
        if genome is None and host_reference is None:
            raise ValueError('New antigen predictions require a host annotation for exact-self assessment')
        validate_input_scopes([InputProvenance('validation', 'antigens', '', '', scope)],
                              genome=genome, host_reference=host_reference)
        if genome is None:
            return scope
        release = getattr(genome, 'release', None)
        annotation = ('ensembl:' + str(release) if release is not None else
                      str(genome.annotation_name) + ':' + str(genome.annotation_version))
        return combine_scopes(scope, InputScope(annotation=annotation), 'Selected host annotation')

    def prediction_provenance(self, *, genome=None, path='', manifest_path=None,
                              manifest_declarations=None, host_reference=None):
        """Resolve identity and scope before acquiring prediction models."""
        scope = self.validate_for_prediction(genome, host_reference=host_reference)
        resolved = dict(self.to_dict(), scope=scope_dict(scope))
        digest = hashlib.sha256(json.dumps(resolved, sort_keys=True, allow_nan=False).encode()).hexdigest()
        return InputProvenance(
            input_id='antigens:' + digest, source_format='antigens', path=str(path),
            content_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest() if path else digest,
            scope=scope,
            manifest_path=str(Path(manifest_path).resolve()) if manifest_path else None,
            manifest_declarations=dict(manifest_declarations or {}),
            report_declarations={'input_history': [provenance_dict(p) for p in self.provenance]})

    def predict(self, predictor, *, genome=None, config=None, path='', host_reference=None):
        """Predict native occurrences through Topiary and return shared evidence."""
        if host_reference is not None:
            from .host_reference import HostReference
            host_reference = HostReference.load(host_reference)
        provenance = self.prediction_provenance(genome=genome, path=path, host_reference=host_reference)
        return predict_antigens(self.antigens, predictor, genome=genome,
                                provenance=provenance, config=config, host_reference=host_reference)


def predict_antigens(antigens, predictor, *, genome, provenance, config=None,
                     annotate=None, self_matcher=None, host_reference=None):
    """Shared native-protein prediction for source-specific admission adapters.

    ``annotate`` may add source evidence columns to the prediction frame before
    shared adaptation. It must not change sequence coordinates or measurements.
    Distinct source names and repeated native occurrences remain distinct.
    """
    if host_reference is not None:
        from .host_reference import HostReference
        host_reference = HostReference.load(host_reference)
        provenance = replace(provenance, scope=host_reference.validate(provenance.scope, genome=genome))
        def self_matcher(peptides, antigen, genome):
            return host_reference.matches(peptides, antigen)
    elif provenance.scope.host_reference_sha256 is not None:
        raise ValueError('Prediction requires its declared host_reference snapshot')
    validate_input_scopes([provenance], genome=genome, host_reference=host_reference)
    self_matcher = self_matcher or self_reference_matches
    alleles = provenance.scope.mhc_alleles
    if not alleles:
        raise ValueError('Antigen prediction requires a declared patient genotype')
    frame = predictor.predict_from_named_sequences({name: a.amino_acids for name, a in antigens.items()})
    if frame.empty:
        raise ValueError('Antigen predictor returned no peptide measurements')
    frame = frame.reset_index(drop=True)
    validate_prediction_frame(frame, 'Native antigen prediction')
    reference = host_reference.proteome.reference if host_reference is not None else genome
    validate_host_taxon(provenance.scope, reference, 'Native antigen predictor output',
                        [a for a in frame.allele if cell_text(a)])
    if any(cell_text(a) and parse_prediction_allele(cell_text(a)) not in alleles for a in frame.allele):
        raise ValueError('Antigen predictor returned an allele outside the declared patient genotype')
    if not set(frame.source_sequence_name) <= set(antigens):
        raise ValueError('Antigen predictor returned an unknown source sequence')
    for column, values in {
        'source_sequence': [antigens[n].amino_acids for n in frame.source_sequence_name],
        ANTIGEN_COLUMN: [to_native_json(antigens[n]) for n in frame.source_sequence_name],
        'gene_id': [antigens[n].gene_id for n in frame.source_sequence_name],
        'gene_name': [antigens[n].display_gene_name for n in frame.source_sequence_name],
        'source_species': [antigens[n].species for n in frame.source_sequence_name],
    }.items():
        frame[column] = values
    frame['sample_name'] = provenance.scope.sample_id
    if annotate is not None:
        annotate(frame)
    adapted = EpitopeDataset.from_topiary(TopiaryResult(frame), label=provenance.input_id)

    def occurrence_id(name, peptide, offset):
        return provenance.input_id + ':' + hashlib.sha256(
            json.dumps([name, peptide, int(offset)]).encode()).hexdigest()

    frame['prediction_id'] = [occurrence_id(r.source_sequence_name, r.peptide, r.peptide_offset)
                              for r in frame.itertuples()]
    candidates = tuple(replace(e, prediction_id=occurrence_id(e.source_name, e.sequence, e.offset))
                       for e in adapted.epitopes)
    dataset = EpitopeDataset.from_predictions(candidates, result=TopiaryResult(frame))
    dataset.antigens = {e.prediction_id: antigens[e.source_name] for e in candidates}
    # Preserve sequences for which the predictor emitted no observations too.
    observed = {e.source_name for e in candidates}
    dataset.antigens.update({provenance.input_id + ':unmeasured:' + name: a
                            for name, a in antigens.items() if name not in observed})
    matches = {name: self_matcher(frame.loc[frame.source_sequence_name.eq(name), 'peptide'], antigen, genome)
               for name, antigen in antigens.items() if name in observed}
    host_matches = {name: self_matcher(
        frame.loc[frame.source_sequence_name.eq(name), 'peptide'],
        replace(antigen, self_reference_excluded_gene_ids=()), genome)
        if antigen.self_reference_excluded_gene_ids else matches[name]
        for name, antigen in antigens.items() if name in observed}
    dataset.epitopes = tuple(replace(
        e, input_provenance=provenance, patient_alleles=alleles,
        occurs_in_reference=host_matches[e.source_name][e.sequence].occurs,
        occurs_in_non_CTA_reference=matches[e.source_name][e.sequence].occurs,
        self_reference_match=matches[e.source_name][e.sequence]) for e in dataset.epitopes)
    dataset.result.df['occurs_in_reference'] = [host_matches[r.source_sequence_name][r.peptide].occurs
                                               for r in frame.itertuples()]
    dataset.result.df['occurs_in_non_CTA_reference'] = [matches[r.source_sequence_name][r.peptide].occurs
                                                       for r in frame.itertuples()]
    dataset.provenance = (provenance,)
    dataset.host_reference = host_reference
    dataset.config = config or EpitopeConfig()
    dataset.selection['antigen_inputs'] = {provenance.input_id: {name: asdict(a) for name, a in antigens.items()}}
    return dataset


def validate_antigen_args(args):
    from .external_rescoring import external_input_specs
    specs = external_input_specs(args)
    if not any(fmt == 'antigens' for fmt, _, _ in specs):
        return
    conflicts = ('vcf', 'bam')
    if any(getattr(args, name, None) for name in conflicts):
        raise ValueError('Combine new antigen sequences with saved VCF/BAM evidence '
                         'using epitopes inputs in --input-manifest')
    if not getattr(args, 'mhc_predictor', None):
        raise ValueError('--input-antigens requires --mhc-predictor for new predictions')
    if getattr(args, 'external_predictions', 'input') != 'input' or getattr(args, 'external_peptide_only', False):
        raise ValueError('--input-antigens predicts native sequences; external prediction modes do not apply')
    if getattr(args, 'prediction_cache', None):
        raise ValueError('Antigen context prediction caches are not yet supported')


def antigen_report(inputs, predictor, *, genome, provenance, config=None, host_reference=None):
    """Adapt validated sequences into the common report representation."""
    from .external_report import ExternalReport
    dataset = predict_antigens(inputs.antigens, predictor, genome=genome,
                               provenance=provenance, config=config, host_reference=host_reference)
    return ExternalReport('antigens', provenance.path, report_df=dataset.report_frame(),
                          epitopes=dataset.epitopes, dataset=dataset,
                          input_provenance=provenance)
