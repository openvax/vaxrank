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


def read_input_manifest(path):
    """Return paths and declarations; resolve paths beside the manifest."""
    path = Path(path).resolve()
    with path.open() as stream:
        document = yaml.load(stream, Loader=_ManifestLoader)
    if not isinstance(document, dict) or document.get("schema") != MANIFEST_SCHEMA:
        raise ValueError(f"{path}: expected schema: {MANIFEST_SCHEMA}")
    unknown = set(document) - {"schema", "inputs", *SCOPE_FIELDS}
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
        unknown = set(item) - {"format", "path", *SCOPE_FIELDS}
        if unknown:
            raise ValueError(f"{location}: unknown input fields: {sorted(unknown)}")
        if item.get("format") not in ("lens", "pvacseq"):
            raise ValueError(f"{location}: format must be lens or pvacseq")
        if not isinstance(item.get("path"), str) or not item["path"].strip():
            raise ValueError(f"{location}: path must be a nonempty string")
        local = {k: v for k, v in item.items() if k in SCOPE_FIELDS}
        local_scope = scope_from_mapping(local, location)
        # Shared patient/reference/genotype assertions cannot be overridden.
        # Sample/timepoint/library defaults can differ between inputs.
        combine_scopes(shared_scope, local_scope, location, SHARED_FIELDS)
        declarations = {**shared, **{k: v for k, v in local.items()
                                   if v is not None or k not in shared}}
        specs.append((item["format"], str((path.parent / item["path"]).resolve()),
                      declarations))
    return specs


def report_declarations(path, source_format):
    """Read explicit scope columns before a format normalizer drops them.

    These exact column names are also accepted in annotated producer tables.
    Per-row HLA prediction columns are deliberately not genotype declarations.
    LENS ERV origin identifiers independently state an assembly, not a patient.
    """
    # Reading the path preserves pandas' support for compressed TSVs.
    frame = pd.read_csv(path, sep="\t", dtype=str,
                        keep_default_na=False,
                        usecols=lambda name: name in (*SCOPE_FIELDS, "origin_descriptor"))
    result = {}
    resolved = InputScope()
    for name in SCOPE_FIELDS:
        if name not in frame:
            continue
        values = [v for v in frame[name].unique() if not cells.missing(v)]
        for value in values:
            parsed = json.loads(value) if name == "mhc_alleles" else value
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
