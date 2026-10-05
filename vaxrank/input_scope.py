"""Patient/reference declarations for external reports, separate from evidence.

The ordinary VCF/BAM entry point does not use this module. Manifest values are
user declarations, not facts inferred from filenames or prediction coverage.
"""

from dataclasses import asdict, dataclass, field, fields
import json
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


SCOPE_FIELDS = tuple(f.name for f in fields(InputScope))
SHARED_FIELDS = ("patient_id", "reference_assembly", "annotation", "mhc_alleles")
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


def normalize_alleles(values, location):
    if not isinstance(values, (list, tuple)) or not values:
        raise ValueError(f"{location}: mhc_alleles must be a nonempty list, or null if unknown")
    try:
        return tuple(sorted({parse_prediction_allele(a) for a in values}))
    except (ParseError, TypeError, ValueError) as error:
        raise ValueError(f"{location}: invalid mhc_alleles: {error}") from error


def scope_from_mapping(values, location):
    """Validate declaration types without guessing missing identifiers."""
    normalized = {}
    for name, value in values.items():
        if name not in SCOPE_FIELDS:
            raise ValueError(f"{location}: unknown scope field {name!r}")
        if value is None:
            normalized[name] = None
        elif name == "mhc_alleles":
            normalized[name] = normalize_alleles(value, location)
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
    unknown = set(document) - {"schema", "inputs", "direct", *SCOPE_FIELDS}
    if unknown:
        raise ValueError(f"{path}: unknown manifest fields: {sorted(unknown)}")
    shared = {k: v for k, v in document.items() if k in SCOPE_FIELDS}
    shared_scope = scope_from_mapping(shared, str(path))
    inputs = document.get("inputs")
    if not isinstance(inputs, list) or not inputs:
        raise ValueError(f"{path}: inputs must be a nonempty list")
    specs = []
    for index, item in enumerate(inputs, start=1):
        location = f"{path}, input {index}"
        if not isinstance(item, dict):
            raise ValueError(f"{location}: expected an input object")
        unknown = set(item) - {"format", "path", "exacto", *SCOPE_FIELDS}
        if unknown:
            raise ValueError(f"{location}: unknown input fields: {sorted(unknown)}")
        from .external_rescoring import READERS
        if item.get("format") not in READERS:
            raise ValueError(f"{location}: format must be one of {', '.join(READERS)}")
        if not isinstance(item.get("path"), str) or not item["path"].strip():
            raise ValueError(f"{location}: path must be a nonempty string")
        local = {k: v for k, v in item.items() if k in SCOPE_FIELDS}
        local_scope = scope_from_mapping(local, location)
        # Shared patient/reference/genotype assertions cannot be overridden.
        # Sample/timepoint/library defaults can differ between inputs.
        combine_scopes(shared_scope, local_scope, location, SHARED_FIELDS)
        declarations = {**shared, **{k: v for k, v in local.items()
                                   if v is not None or k not in shared}}
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


def report_declarations(path, source_format, *, frame=None):
    """Read explicit scope columns before a format normalizer drops them.

    These exact column names are also accepted in annotated producer tables.
    Per-row HLA prediction columns are deliberately not genotype declarations.
    LENS ERV origin identifiers independently state an assembly, not a patient.
    """
    # Reading the path preserves pandas' support for compressed TSVs.
    if frame is None:
        frame = pd.read_csv(path, sep="\t", dtype=str,
                           keep_default_na=False,
                           usecols=lambda name: name in (*SCOPE_FIELDS, "origin_descriptor"))
    result = {}
    resolved = InputScope()
    for name in SCOPE_FIELDS:
        if name not in frame:
            continue
        values = []
        for value in frame[name]:
            if not cells.missing(value) and value not in values:
                values.append(value)
        for value in values:
            parsed = (json.loads(value) if name == "mhc_alleles" and isinstance(value, str)
                      else value)
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


def validate_input_scopes(provenance, *, genome=None, prediction_alleles=()):
    """Validate the whole batch before DSL scoring, output or live prediction."""
    shared = InputScope()
    genome = resolve_input_genome(genome)
    for source in provenance:
        shared = combine_scopes(shared, source.scope, source.path, SHARED_FIELDS)
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
        if genome is None and source.scope.reference_assembly is not None:
            # Variant adapters historically skip unresolvable rows. Reject a
            # bad declared assembly here instead of silently dropping them.
            try:
                resolve_input_genome(source.scope.reference_assembly)
            except (TypeError, ValueError) as error:
                raise ValueError(f"{source.path}: unsupported reference_assembly "
                                 f"{source.scope.reference_assembly!r}: {error}") from error
        elif genome is not None:
            assembly = getattr(genome, "reference_name", None)
            if source.scope.reference_assembly is not None and assembly is None:
                raise ValueError(f"{source.path}: configured genome has no reference_name "
                                 "to validate against reference_assembly")
            release = getattr(genome, "release", None)
            configured = InputScope(
                reference_assembly=assembly,
                annotation=f"ensembl:{release}" if release is not None else None)
            combine_scopes(source.scope, configured,
                           f"{source.path}, configured genome", ("reference_assembly", "annotation"))
    if prediction_alleles and shared.mhc_alleles is not None:
        requested = normalize_alleles(prediction_alleles, "Fresh prediction")
        outside = set(requested) - set(shared.mhc_alleles)
        if outside:
            raise ValueError(f"Fresh prediction alleles outside declared mhc_alleles: {sorted(outside)}")
    return shared


def provenance_columns(provenance):
    """Plain table columns plus a complete, reloadable provenance record."""
    scope = asdict(provenance.scope)
    scope["mhc_alleles"] = (json.dumps(scope["mhc_alleles"])
                            if scope["mhc_alleles"] is not None else None)
    return {
        **{"input_" + name: value for name, value in scope.items()},
        "input_observed_mhc_alleles": json.dumps(provenance.observed_mhc_alleles),
        "input_provenance_json": json.dumps(asdict(provenance), sort_keys=True),
    }
