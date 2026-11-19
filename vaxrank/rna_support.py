"""Source-independent RNA measurements, with explicit protein-window scope."""

from dataclasses import dataclass, replace
from functools import lru_cache
import hashlib
import json

from isovar import validate_rna_feature_support
from isovar.rna_feature_support import SCHEMA as ISOVAR_RNA_SCHEMA
from serializable import DataclassSerializable

from .native_serialization import _json_loads


MEASUREMENT_REF_FIELD = "__vaxrank_rna_measurement__"
PROFILE_REF_FIELD = "__vaxrank_sampling_profile__"
UNION_REF_FIELD = "__vaxrank_shared_rna_union__"
SHARED_RNA_UNION_SCHEMA = "vaxrank.shared_sequence_rna_union.v1"
_HEX64 = frozenset("0123456789abcdef")

RNA_BINDINGS = frozenset(
    "rna_" + name
    for unit in ("reads", "fragments")
    for name in (unit, "library_" + unit, "effective_length_" + unit,
                 unit + "_per_million", unit + "_per_effective_kb",
                 unit + "_per_effective_kb_per_million"))


def _canonical_json(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def is_sha256_digest(value) -> bool:
    return isinstance(value, str) and len(value) == 64 and set(value) <= _HEX64


@lru_cache(maxsize=256)
def _validate_canonical_measurement(canonical: str) -> None:
    """Validate one canonical Isovar measurement JSON payload once."""
    validate_rna_feature_support(_json_loads(canonical))


def _contains_reserved_marker(value) -> bool:
    if isinstance(value, list):
        return any(_contains_reserved_marker(item) for item in value)
    if isinstance(value, dict):
        if (
            MEASUREMENT_REF_FIELD in value
            or PROFILE_REF_FIELD in value
            or UNION_REF_FIELD in value
        ):
            return True
        return any(_contains_reserved_marker(item) for item in value.values())
    return False


@dataclass(frozen=True)
class RNAFeatureEvidence(DataclassSerializable):
    """One Isovar feature measurement and its declared applicability.

    ``event`` attaches local event evidence to this source sequence and retains
    it through slicing; it does not assert observation of downstream coding
    bases. ``window`` describes exactly ``protein_interval`` (zero-based,
    half-open amino acids). Its score binds only for that exact window.
    Multiple libraries/features remain separate; no pooling or best-sample
    policy is inferred. Normalized values are computed by the producer before
    any explicitly chosen longitudinal comparison.
    """

    measurement: dict
    scope: str
    protein_interval: tuple[int, int] | None = None

    def __post_init__(self):
        if not isinstance(self.measurement, dict):
            raise ValueError("RNA feature support measurement must be an object")
        canonical = _canonical_json(self.measurement)
        _validate_canonical_measurement(canonical)
        object.__setattr__(self, "measurement", json.loads(canonical))
        if self.scope not in ("event", "window"):
            raise ValueError("RNA evidence scope must be event or window")
        if self.scope == "event":
            if self.protein_interval is not None:
                raise ValueError("Event RNA support does not assert a protein-window span")
        else:
            interval = tuple(self.protein_interval or ())
            if (len(interval) != 2 or any(type(x) is not int for x in interval)
                    or not 0 <= interval[0] < interval[1]):
                raise ValueError("Window RNA support requires an amino-acid interval")
            object.__setattr__(self, "protein_interval", interval)

    def __hash__(self):
        return hash((self.measurement["measurement_id"], self.scope, self.protein_interval))


class RNASupportArchive:
    """Content-addressed in-file deduplication for RNA measurements and profiles."""

    def __init__(
        self,
        profiles=None,
        measurements=None,
        resolved_canonical=None,
        unions=None,
        resolved_unions=None,
    ):
        self.profiles = {} if profiles is None else profiles
        self.measurements = {} if measurements is None else measurements
        self._resolved_canonical = {} if resolved_canonical is None else resolved_canonical
        self.unions = {} if unions is None else unions
        self._resolved_unions = {} if resolved_unions is None else resolved_unions

    def encode_tree(self, value):
        if isinstance(value, list):
            return [self.encode_tree(item) for item in value]
        if isinstance(value, tuple):
            return tuple(self.encode_tree(item) for item in value)
        if not isinstance(value, dict):
            return value
        if (
            MEASUREMENT_REF_FIELD in value
            or PROFILE_REF_FIELD in value
            or UNION_REF_FIELD in value
        ):
            raise ValueError("Reserved native RNA reference marker in unencoded payload")
        if value.get("schema") == ISOVAR_RNA_SCHEMA:
            if _contains_reserved_marker(value):
                raise ValueError("Reserved native RNA reference marker inside RNA measurement")
            canonical = _canonical_json(value)
            _validate_canonical_measurement(canonical)
            digest = _sha256(canonical)
            if digest not in self.measurements:
                stored = json.loads(canonical)
                profile = stored.get("sampling_profile")
                if isinstance(profile, dict):
                    profile_canon = _canonical_json(profile)
                    profile_digest = _sha256(profile_canon)
                    self.profiles[profile_digest] = json.loads(profile_canon)
                    stored["sampling_profile"] = {PROFILE_REF_FIELD: profile_digest}
                self.measurements[digest] = stored
            return {MEASUREMENT_REF_FIELD: digest}
        if value.get("schema") == SHARED_RNA_UNION_SCHEMA:
            if _contains_reserved_marker(value):
                raise ValueError("Reserved native RNA reference marker inside shared RNA union")
            canonical = _canonical_json(value)
            digest = _sha256(canonical)
            if digest not in self.unions:
                self.unions[digest] = json.loads(canonical)
            return {UNION_REF_FIELD: digest}
        return {key: self.encode_tree(item) for key, item in value.items()}

    def encode_json(self, serialized):
        if serialized is None or (
            ISOVAR_RNA_SCHEMA not in serialized
            and SHARED_RNA_UNION_SCHEMA not in serialized
            and MEASUREMENT_REF_FIELD not in serialized
            and PROFILE_REF_FIELD not in serialized
            and UNION_REF_FIELD not in serialized
        ):
            return serialized
        return json.dumps(self.encode_tree(_json_loads(serialized)), allow_nan=False)

    def write_to_payload(self, payload: dict) -> None:
        if self.profiles:
            payload["rna_sampling_profiles"] = self.profiles
        if self.measurements:
            payload["rna_measurements"] = self.measurements
        if self.unions:
            payload["rna_shared_unions"] = self.unions

    @classmethod
    def from_payload(cls, payload: dict) -> "RNASupportArchive":
        raw_profiles = payload.pop("rna_sampling_profiles", None)
        raw_measurements = payload.pop("rna_measurements", None)
        raw_unions = payload.pop("rna_shared_unions", None)
        profiles = {}
        if raw_profiles is not None:
            if not isinstance(raw_profiles, dict):
                raise ValueError("Native RNA sampling profiles must be an object")
            for digest, profile in raw_profiles.items():
                if (not is_sha256_digest(digest) or not isinstance(profile, dict)
                        or _contains_reserved_marker(profile)):
                    raise ValueError("Invalid native RNA sampling profile entry")
                canonical = _canonical_json(profile)
                if _sha256(canonical) != digest:
                    raise ValueError("Native RNA sampling profile checksum mismatch")
                profiles[digest] = json.loads(canonical)
        measurements = {}
        resolved_canonical = {}
        if raw_measurements is not None:
            if not isinstance(raw_measurements, dict):
                raise ValueError("Native RNA measurements must be an object")
            for digest, entry in raw_measurements.items():
                if not is_sha256_digest(digest) or not isinstance(entry, dict):
                    raise ValueError("Invalid native RNA measurement entry")
                resolved = dict(entry)
                profile_ref = resolved.get("sampling_profile")
                if isinstance(profile_ref, dict) and PROFILE_REF_FIELD in profile_ref:
                    if set(profile_ref) != {PROFILE_REF_FIELD}:
                        raise ValueError("Malformed native RNA sampling profile reference")
                    profile_digest = profile_ref[PROFILE_REF_FIELD]
                    if not isinstance(profile_digest, str) or profile_digest not in profiles:
                        raise ValueError("Unknown native RNA sampling profile reference")
                    resolved["sampling_profile"] = json.loads(_canonical_json(profiles[profile_digest]))
                if _contains_reserved_marker(resolved):
                    raise ValueError("Unexpected native RNA reference marker in measurement")
                canonical = _canonical_json(resolved)
                if _sha256(canonical) != digest:
                    raise ValueError("Native RNA measurement checksum mismatch")
                _validate_canonical_measurement(canonical)
                measurements[digest] = entry
                resolved_canonical[digest] = canonical
        unions = {}
        resolved_unions = {}
        if raw_unions is not None:
            if not isinstance(raw_unions, dict):
                raise ValueError("Native shared RNA unions must be an object")
            for digest, entry in raw_unions.items():
                if (not is_sha256_digest(digest) or not isinstance(entry, dict)
                        or entry.get("schema") != SHARED_RNA_UNION_SCHEMA
                        or _contains_reserved_marker(entry)):
                    raise ValueError("Invalid native shared RNA union entry")
                canonical = _canonical_json(entry)
                if _sha256(canonical) != digest:
                    raise ValueError("Native shared RNA union checksum mismatch")
                unions[digest] = json.loads(canonical)
                resolved_unions[digest] = canonical
        return cls(
            profiles=profiles,
            measurements=measurements,
            resolved_canonical=resolved_canonical,
            unions=unions,
            resolved_unions=resolved_unions,
        )

    def decode_tree(self, value):
        if isinstance(value, list):
            return [self.decode_tree(item) for item in value]
        if isinstance(value, tuple):
            return tuple(self.decode_tree(item) for item in value)
        if not isinstance(value, dict):
            return value
        if PROFILE_REF_FIELD in value:
            raise ValueError("Unexpected native RNA sampling profile reference")
        if MEASUREMENT_REF_FIELD in value:
            if set(value) != {MEASUREMENT_REF_FIELD}:
                raise ValueError("Malformed native RNA measurement reference")
            digest = value[MEASUREMENT_REF_FIELD]
            if not isinstance(digest, str) or digest not in self._resolved_canonical:
                raise ValueError("Unknown native RNA measurement reference")
            return json.loads(self._resolved_canonical[digest])
        if UNION_REF_FIELD in value:
            if set(value) != {UNION_REF_FIELD}:
                raise ValueError("Malformed native shared RNA union reference")
            digest = value[UNION_REF_FIELD]
            if not isinstance(digest, str) or digest not in self._resolved_unions:
                raise ValueError("Unknown native shared RNA union reference")
            return json.loads(self._resolved_unions[digest])
        return {key: self.decode_tree(item) for key, item in value.items()}

    def decode_json(self, serialized):
        if serialized is None or (
            MEASUREMENT_REF_FIELD not in serialized
            and PROFILE_REF_FIELD not in serialized
            and UNION_REF_FIELD not in serialized
        ):
            return serialized
        return json.dumps(self.decode_tree(_json_loads(serialized)), allow_nan=False)


def validate_rna_support(records, protein_length):
    records = tuple(records)
    identities = []
    for record in records:
        if not isinstance(record, RNAFeatureEvidence):
            raise TypeError("RNA support requires RNAFeatureEvidence records")
        if not isinstance(record.measurement, dict):
            raise ValueError("RNA feature support measurement must be an object")
        _validate_canonical_measurement(_canonical_json(record.measurement))
        if record.protein_interval is not None and record.protein_interval[1] > protein_length:
            raise ValueError("RNA window evidence extends beyond its protein sequence")
        identities.append((record.measurement["measurement_id"], record.scope, record.protein_interval))
    if len(identities) != len(set(identities)):
        raise ValueError("Duplicate RNA feature measurement attachment")
    return records


def slice_rna_support(records, start, end):
    """Retain local event evidence and rebase contained, exactly named windows."""
    result = []
    for record in records:
        if record.scope == "event":
            result.append(record)
        elif start <= record.protein_interval[0] and record.protein_interval[1] <= end:
            result.append(replace(record, protein_interval=tuple(x - start for x in record.protein_interval)))
    return tuple(result)


def applicable_rna_support(antigen):
    return tuple(record for record in getattr(antigen, "rna_support", ())
                 if record.scope == "event" or record.protein_interval == (0, len(antigen.amino_acids)))


def rna_support_bindings(antigen):
    """Missing or multiply scoped evidence remains unavailable, never zero."""
    records = applicable_rna_support(antigen)
    values = records[0].measurement["values"] if len(records) == 1 else {}
    return {name: values.get(name.removeprefix("rna_")) for name in RNA_BINDINGS}
