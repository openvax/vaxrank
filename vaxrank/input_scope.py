"""Patient/reference declarations for external reports, separate from evidence.

The ordinary VCF/BAM entry point does not use this module. Manifest values are
user declarations, not facts inferred from filenames or prediction coverage.
"""

from dataclasses import asdict, dataclass, field, fields, replace
from functools import lru_cache
import json
from numbers import Integral, Real
from pathlib import Path
import re

import pandas as pd
from mhcgnomes import ParseError
from serializable import DataclassSerializable
import yaml

from . import cells
from .allele_validation import parse_prediction_allele


MANIFEST_SCHEMA = "vaxrank.input_manifest.v1"
PROVENANCE_SCHEMA = "vaxrank.input_provenance.v1"


@dataclass(frozen=True)
class InputScope(DataclassSerializable):
    patient_id: str | None = None
    reference_assembly: str | None = None
    annotation: str | None = None
    mhc_alleles: tuple[str, ...] | None = None
    sample_id: str | None = None
    timepoint: str | None = None
    library_id: str | None = None
    host_taxon_id: int | None = None
    mhc_species: tuple[str, ...] | None = None

    def __post_init__(self):
        if self.host_taxon_id is not None and (
                type(self.host_taxon_id) is not int or self.host_taxon_id < 1):
            raise ValueError("host_taxon_id must be a positive integer NCBI taxonomy ID")
        if self.mhc_species is not None:
            object.__setattr__(self, 'mhc_species', normalize_mhc_species(self.mhc_species))

    @classmethod
    def from_dict(cls, values):
        return scope_from_mapping(values, 'Input scope')


SCOPE_FIELDS = tuple(f.name for f in fields(InputScope))
DECLARATION_FIELDS = (*SCOPE_FIELDS, "host_species")
SPECIES_FIELDS = ("host_taxon_id", "mhc_species")
SPECIES_DECLARATIONS = (*SPECIES_FIELDS, "host_species")
SHARED_FIELDS = ("patient_id", "reference_assembly", "annotation", "mhc_alleles", *SPECIES_FIELDS)
REQUIRED_POOLED_FIELDS = ("patient_id", "reference_assembly", "mhc_alleles")


@dataclass(frozen=True)
class InputProvenance(DataclassSerializable):
    input_id: str
    source_format: str
    path: str
    content_sha256: str
    scope: InputScope
    manifest_path: str | None = None
    manifest_declarations: dict = field(default_factory=dict)
    report_declarations: dict = field(default_factory=dict)
    observed_mhc_alleles: tuple[str, ...] = ()


def scope_dict(scope):
    """Keep historical identities unchanged when no host taxon was declared."""
    result = asdict(scope)
    for name in SPECIES_FIELDS:
        if result[name] is None:
            del result[name]
    return result


def provenance_dict(provenance):
    result = asdict(provenance)
    result['scope'] = scope_dict(provenance.scope)
    return result


def normalize_alleles(values, location):
    if not isinstance(values, (list, tuple)) or not values:
        raise ValueError(f"{location}: mhc_alleles must be a nonempty list, or null if unknown")
    try:
        return tuple(sorted({parse_prediction_allele(a) for a in values}))
    except (ParseError, TypeError, ValueError) as error:
        raise ValueError(f"{location}: invalid mhc_alleles: {error}") from error


def normalize_host_species(value):
    """Resolve a name or taxonomy ID to the canonical host taxon."""
    from orfdata.reference import species_taxon
    if isinstance(value, str):
        value = value.strip()
        value = int(value) if value.isdecimal() else species_taxon(value)
    if type(value) is not int or value < 1:
        raise ValueError('host_species requires a known species name or positive NCBI taxonomy ID')
    return value


def normalize_mhc_species(values):
    """Resolve an explicit repertoire, including transgenic/chimeric animals."""
    from mhcgnomes import Species
    from orfdata.reference import canonical_species
    if isinstance(values, (str, int)):
        values = [values]
    if not isinstance(values, (list, tuple)) or not values:
        raise ValueError('mhc_species requires a species name or taxonomy ID, or a nonempty list')
    result = set()
    for value in values:
        if isinstance(value, str):
            value = value.strip()
            if value.isdecimal():
                value = int(value)
        species = (_species_for_taxon(value) if type(value) is int else
                   Species.get(canonical_species(value).replace('_', ' ').capitalize())
                   if isinstance(value, str) else None)
        if species is None or species.is_group:
            raise ValueError(f'mhc_species requires a known concrete species, got {value!r}')
        result.add(species.name)
    return tuple(sorted(result))


def species_scope_from_args(args):
    return scope_from_mapping(
        {name: getattr(args, name, None) for name in SPECIES_DECLARATIONS}, 'CLI species')


def validate_prediction_species(args, references, alleles, location):
    """Check live candidate/junction queries before acquiring a model."""
    scope = combine_scopes(scope_from_mapping(getattr(args, '_input_species_scope', {}), location),
                           species_scope_from_args(args), location)
    references = [reference for reference in references if reference is not None]
    if not references and scope.reference_assembly is not None:
        references = [scope.reference_assembly]
    for reference in references or (None,):
        validate_host_taxon(scope, resolve_input_genome(reference), location, alleles)
    return scope


def with_species_declarations(declarations, scope, location):
    """Apply CLI assertions without replacing manifest or adapter options."""
    existing = scope_from_mapping(
        {k: v for k, v in declarations.items() if k in SPECIES_DECLARATIONS}, location)
    combined = combine_scopes(existing, scope, location, SPECIES_FIELDS)
    return {**declarations, **{k: list(v) if isinstance(v, tuple) else v
                              for k, v in scope_dict(combined).items()
                              if k in SPECIES_FIELDS}}


def scope_from_mapping(values, location):
    """Validate declaration types without guessing missing identifiers."""
    values = dict(values)
    host_species = values.pop('host_species', None)
    if host_species is not None:
        taxon = normalize_host_species(host_species)
        declared = values.get('host_taxon_id')
        if declared is not None and declared != taxon:
            raise ValueError(f'{location}: conflicting host_species and host_taxon_id')
        if declared is None:
            values['host_taxon_id'] = taxon
    normalized = {}
    for name, value in values.items():
        if name not in SCOPE_FIELDS:
            raise ValueError(f"{location}: unknown scope field {name!r}")
        if value is None:
            normalized[name] = None
        elif name == "mhc_alleles":
            normalized[name] = normalize_alleles(value, location)
        elif name == "mhc_species":
            normalized[name] = normalize_mhc_species(value)
        elif name == "host_taxon_id":
            if type(value) is not int or value < 1:
                raise ValueError(f"{location}: host_taxon_id must be a positive integer NCBI taxonomy ID")
            normalized[name] = value
        elif not isinstance(value, str) or not value.strip():
            raise ValueError(f"{location}: {name} must be a nonempty string, or null if unknown")
        else:
            normalized[name] = value.strip()
    return InputScope(**normalized)


def combine_scopes(left, right, location, names=SCOPE_FIELDS):
    """Fill unknowns, rejecting contradictory stated declarations."""
    values = asdict(left)
    for name in names:
        a, b = getattr(left, name), getattr(right, name)
        if a is not None and b is not None and a != b:
            raise ValueError(f"{location}: conflicting {name}: {a!r} versus {b!r}")
        if b is not None:
            values[name] = b
    return InputScope(**values)


class _ManifestLoader(yaml.SafeLoader):
    pass


def _unique_mapping(loader, node):
    pairs = loader.construct_pairs(node, deep=True)
    result = {}
    for key, value in pairs:
        if not isinstance(key, str):
            raise ValueError("Input manifest keys must be strings")
        if key in result:
            raise ValueError(f"Duplicate input manifest key: {key}")
        result[key] = value
    return result


_ManifestLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _unique_mapping)


def read_input_manifest(path, *, include_direct=False):
    """Return paths and declarations; resolve paths beside the manifest."""
    path = Path(path).resolve()
    with path.open() as stream:
        document = yaml.load(stream, Loader=_ManifestLoader)
    if not isinstance(document, dict) or document.get("schema") != MANIFEST_SCHEMA:
        raise ValueError(f"{path}: expected schema: {MANIFEST_SCHEMA}")
    unknown = set(document) - {"schema", "inputs", "direct", *DECLARATION_FIELDS}
    if unknown:
        raise ValueError(f"{path}: unknown manifest fields: {sorted(unknown)}")
    shared = {k: v for k, v in document.items() if k in DECLARATION_FIELDS}
    shared_scope = scope_from_mapping(shared, str(path))
    inputs = document.get("inputs")
    if not isinstance(inputs, list) or not inputs:
        raise ValueError(f"{path}: inputs must be a nonempty list")
    specs = []
    for index, item in enumerate(inputs, start=1):
        location = f"{path}, input {index}"
        if not isinstance(item, dict):
            raise ValueError(f"{location}: expected an input object")
        unknown = set(item) - {"format", "path", "exacto", "cta_expression", *DECLARATION_FIELDS}
        if unknown:
            raise ValueError(f"{location}: unknown input fields: {sorted(unknown)}")
        from .external_rescoring import INPUT_FORMATS
        if item.get("format") not in INPUT_FORMATS:
            raise ValueError(f"{location}: format must be one of {', '.join(INPUT_FORMATS)}")
        if not isinstance(item.get("path"), str) or not item["path"].strip():
            raise ValueError(f"{location}: path must be a nonempty string")
        local = {k: v for k, v in item.items() if k in DECLARATION_FIELDS}
        local_scope = scope_from_mapping(local, location)
        # Shared patient/reference/genotype assertions cannot be overridden.
        # Sample/timepoint/library defaults can differ between inputs.
        combine_scopes(shared_scope, local_scope, location, SHARED_FIELDS)
        declarations = {**shared, **{k: v for k, v in local.items()
                                   if v is not None or k not in shared}}
        if item['format'] == 'cta_expression' or 'cta_expression' in item:
            from .cta_input import cta_expression_contract
            if item['format'] != 'cta_expression':
                raise ValueError(f'{location}: cta_expression options require a CTA expression input')
            options = item.get('cta_expression')
            cta_expression_contract(options, scope_from_mapping(declarations, location), location)
            options = dict(options)
            if 'hitlist_evidence_bundle' in options:
                options['hitlist_evidence_bundle'] = str((path.parent / options['hitlist_evidence_bundle']).resolve())
            declarations['cta_expression'] = options
        if 'exacto' in item:
            options = item['exacto']
            if item['format'] != 'exacto' or not isinstance(options, dict):
                raise ValueError(f"{location}: exacto options require an Exacto input object")
            allowed = {'schema', 'primary_structures', 'transcript_read_support', 'read_set_id', 'tag'}
            if set(options) - allowed:
                raise ValueError(f"{location}: unknown Exacto options: {sorted(set(options) - allowed)}")
            if any(not isinstance(v, str) or not v.strip() for v in options.values()):
                raise ValueError(f"{location}: Exacto options must be nonempty strings")
            options = dict(options)
            for name in ('primary_structures', 'transcript_read_support'):
                if name in options:
                    options[name] = str((path.parent / options[name]).resolve())
            declarations['exacto'] = options
        specs.append((item["format"], str((path.parent / item["path"]).resolve()),
                      declarations))
    direct = document.get("direct", {})
    if not isinstance(direct, dict):
        raise ValueError(f"{path}: direct must be a scope object")
    direct_scope = scope_from_mapping(direct, f"{path}, direct input")
    combine_scopes(shared_scope, direct_scope, str(path), SHARED_FIELDS)
    if include_direct:
        declarations = {**shared, **{k: v for k, v in direct.items()
                                    if v is not None or k not in shared}}
        return specs, declarations
    return specs


def normalize_scope_columns(frame, location):
    """Restore nullable integer taxon IDs at the table ingestion boundary.

    Pandas promotes an integer column with blanks to floats. Normalize that
    representation here, leaving manifest/API validation strictly typed.
    """
    for name in (*SPECIES_DECLARATIONS, *("input_" + n for n in SPECIES_DECLARATIONS)):
        if name not in frame:
            continue
        values = []
        for value in frame[name]:
            if cells.missing(value):
                values.append(None)
                continue
            if isinstance(value, str):
                try:
                    parsed = json.loads(value)
                except ValueError as error:
                    if name.endswith('host_taxon_id'):
                        raise ValueError(f'{location}: invalid {name}: {value!r}') from error
                else:
                    if name.endswith('host_taxon_id') or isinstance(parsed, (Real, bool)):
                        value = parsed
            if not name.endswith('host_taxon_id') and isinstance(value, (str, list, tuple)):
                values.append(value)
                continue
            if isinstance(value, Integral) and not isinstance(value, bool):
                value = int(value)
            elif (isinstance(value, Real) and not isinstance(value, bool)
                  and value > 0 and value < 2**53 and int(value) == value):
                value = int(value)
            if type(value) is not int or value < 1:
                raise ValueError(f'{location}: {name} must be a positive integer NCBI taxonomy ID')
            values.append(value)
        frame[name] = pd.Series(values, index=frame.index, dtype=object)
    return frame


def report_declarations(path, source_format, *, frame=None):
    """Read explicit scope columns before a format normalizer drops them.

    These exact column names are also accepted in annotated producer tables.
    Per-row HLA prediction columns are deliberately not genotype declarations.
    LENS ERV origin identifiers independently state an assembly, not a patient.
    """
    # Reading the path preserves pandas' support for compressed TSVs.
    if frame is None:
        native = source_format == 'epitopes'
        skiprows = 0
        if native:
            with open(path) as stream:
                for line in stream:
                    if not line.startswith('#'):
                        break
                    skiprows += 1
        frame = pd.read_csv(path, sep=',' if native and not str(path).endswith('.tsv') else '\t',
                           skiprows=skiprows, dtype=str,
                           keep_default_na=False,
                           usecols=lambda name: name in (*DECLARATION_FIELDS, "origin_descriptor"))
    normalize_scope_columns(frame, str(path))
    result = {}
    resolved = InputScope()
    for name in DECLARATION_FIELDS:
        if name not in frame:
            continue
        values = []
        for value in frame[name]:
            if not cells.missing(value) and value not in values:
                values.append(value)
        for value in values:
            parsed = value
            if name in {"mhc_alleles", "mhc_species"} and isinstance(value, str):
                if name == 'mhc_alleles' or value.lstrip().startswith(('[', '"')):
                    parsed = json.loads(value)
            claim = scope_from_mapping({name: parsed}, str(path))
            resolved = combine_scopes(resolved, claim, str(path))
        if values:
            result[name] = values
    if source_format == "lens" and "origin_descriptor" in frame:
        markers = {}
        for value in frame["origin_descriptor"].unique():
            match = re.match(r"^Hsap(37|38)\.", value)
            if match:
                assembly = "GRCh" + match.group(1)
                # The prefix is the assembly declaration; retaining every
                # genomic coordinate here would duplicate the entire source
                # table on every candidate's provenance record.
                markers[match.group(0)] = assembly
                resolved = combine_scopes(
                    resolved, InputScope(reference_assembly=assembly), str(path))
        if markers:
            result["origin_descriptor"] = markers
    return resolved, result


def resolve_input_genome(genome):
    """Resolve a configured name/release with Varcode's existing reference API."""
    if isinstance(genome, (str, int)):
        from varcode.reference import infer_genome
        return infer_genome(genome)[0]
    return genome


@lru_cache(maxsize=128)
def _species_for_taxon(taxon_id):
    """Resolve taxonomy without requiring an annotation for that species."""
    from mhcgnomes import Species
    from orfdata.reference import species_name_for_taxon
    name = species_name_for_taxon(taxon_id)
    return Species.get(name.replace('_', ' ').capitalize()) if name is not None else None


def validate_host_taxon(scope, reference, location, observed_alleles=()):
    """Validate an explicit host independently of antigen source species.

    ORFdata owns taxon aliases; mhcgnomes owns allele species and group
    compatibility (for example DLA's Canis group). An unknown allele taxon
    cannot establish host identity. New prediction inputs also validate the
    actual annotation, while saved evidence can be replayed without one.
    """
    from mhcgnomes import Species, parse
    from orfdata import ReferenceContext
    from orfdata.reference import canonical_species

    host_species = _species_for_taxon(scope.host_taxon_id) if scope.host_taxon_id else None
    if reference is not None:
        if scope.host_taxon_id is not None:
            try:
                context = ReferenceContext.from_genome(reference)
            except (TypeError, ValueError) as error:
                raise ValueError(f"{location}: cannot validate host_taxon_id against annotation: {error}") from error
            if context.taxon_id != scope.host_taxon_id:
                raise ValueError(f"{location}: conflicting host_taxon_id: "
                                 f"{scope.host_taxon_id} versus annotation {context.taxon_id}")
            host_species = Species.get(canonical_species(context.species).replace('_', ' ').capitalize())
        else:
            name = getattr(reference, 'species', None)
            if name is None and getattr(reference, 'reference_name', None):
                from pyensembl import find_species_by_reference
                try:
                    name = find_species_by_reference(reference.reference_name)
                except ValueError:
                    pass
            name = getattr(name, 'latin_name', name)
            host_species = Species.get(canonical_species(name).replace('_', ' ').capitalize()) if name else None
    if scope.host_taxon_id is not None and host_species is None:
        raise ValueError(f'{location}: cannot resolve host_taxon_id {scope.host_taxon_id} to an MHC species')
    expected = ([Species.get(name) for name in scope.mhc_species]
                if scope.mhc_species is not None else [host_species] if host_species else [])
    for allele in set(scope.mhc_alleles or ()) | set(observed_alleles):
        species = parse(parse_prediction_allele(allele), infer_class2_pairing=False).species
        if expected and not any(item.compatible_with(species) for item in expected):
            label = (f'mhc_species {scope.mhc_species}' if scope.mhc_species is not None else
                     f'host_taxon_id {scope.host_taxon_id}' if scope.host_taxon_id else
                     f'host species {host_species.name}')
            raise ValueError(f'{location}: allele {allele} conflicts with {label}. '
                             'For a chimeric/transgenic animal, declare its intended mhc_species.')


def validate_input_scopes(provenance, *, genome=None, prediction_alleles=(),
                          retained_genomes=None):
    """Validate the whole batch before DSL scoring, output or live prediction."""
    shared = InputScope()
    genome = resolve_input_genome(genome)
    retained_genomes = retained_genomes or {}
    for source in provenance:
        shared = combine_scopes(shared, source.scope, source.path, SHARED_FIELDS)
    for source in provenance:
        # A shared host assertion applies to every source, including older
        # evidence which did not itself declare a host taxon.
        host_scope = replace(source.scope, **{name: getattr(shared, name) for name in SPECIES_FIELDS})
        if len(provenance) > 1:
            missing = [name for name in REQUIRED_POOLED_FIELDS
                       if getattr(source.scope, name) is None]
            if missing:
                raise ValueError(
                    f"{source.path}: combining reports requires declared "
                    f"{', '.join(missing)}. Use --input-manifest to declare the "
                    "shared patient, reference_assembly and mhc_alleles; "
                    "--output-patient-id is only an output label.")
        genotype = source.scope.mhc_alleles
        if genotype is not None:
            outside = set(source.observed_mhc_alleles) - set(genotype)
            if outside:
                raise ValueError(f"{source.path}: reported alleles outside declared "
                                 f"mhc_alleles: {sorted(outside)}")
        references = list(retained_genomes.get(source.input_id, ()))
        if genome is not None:
            references.append(genome)
        observed_alleles = (*source.observed_mhc_alleles, *prediction_alleles)
        if not references and source.scope.reference_assembly is not None:
            # Variant adapters historically skip unresolvable rows. Reject a
            # bad declared assembly here instead of silently dropping them.
            try:
                declared_reference = resolve_input_genome(source.scope.reference_assembly)
            except (TypeError, ValueError) as error:
                raise ValueError(f"{source.path}: unsupported reference_assembly "
                                 f"{source.scope.reference_assembly!r}: {error}") from error
            # Assembly metadata identifies species, not the saved annotation release.
            validate_host_taxon(
                host_scope, declared_reference, source.path,
                observed_alleles)
        elif not references:
            validate_host_taxon(host_scope, None, source.path, observed_alleles)
        for reference in references:
            validate_host_taxon(host_scope, reference, source.path, observed_alleles)
            assembly = getattr(reference, "reference_name", None)
            if source.scope.reference_assembly is not None and assembly is None:
                raise ValueError(f"{source.path}: configured genome has no reference_name "
                                 "to validate against reference_assembly")
            release = getattr(reference, "release", None)
            configured = InputScope(
                reference_assembly=assembly,
                annotation=f"ensembl:{release}" if release is not None else None)
            combine_scopes(source.scope, configured,
                           f"{source.path}, configured or retained genome",
                           ("reference_assembly", "annotation"))
    if prediction_alleles and shared.mhc_alleles is not None:
        requested = normalize_alleles(prediction_alleles, "Fresh prediction")
        outside = set(requested) - set(shared.mhc_alleles)
        if outside:
            raise ValueError(f"Fresh prediction alleles outside declared mhc_alleles: {sorted(outside)}")
    return shared


def provenance_columns(provenance):
    """Plain table columns plus a complete, reloadable provenance record."""
    scope = scope_dict(provenance.scope)
    scope["mhc_alleles"] = (json.dumps(scope["mhc_alleles"])
                            if scope["mhc_alleles"] is not None else None)
    if scope.get('mhc_species') is not None:
        scope['mhc_species'] = json.dumps(scope['mhc_species'])
    return {
        **{"input_" + name: value for name, value in scope.items()},
        "input_observed_mhc_alleles": json.dumps(provenance.observed_mhc_alleles),
        "input_provenance_json": json.dumps(provenance_dict(provenance), sort_keys=True),
    }
