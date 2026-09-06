# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Allowlisted JSON serialization for native epitope object graphs."""

import json

from serializable import to_json


# Resolve ``serializable`` type metadata through a static schema so native
# files have a deterministic, constrained object graph. Adding a new non-
# primitive CandidateEpitope field fails loudly until its value type is
# reviewed and added here; the field can never silently disappear.

_CLASS_FIELD = "__class__"
_VALUE_FIELD = "__value__"
_SERIALIZED_KEYS_FIELD = "__serialized_keys__"
_SERIALIZED_KEY_PREFIX = _SERIALIZED_KEYS_FIELD + "element_"


def _native_serializable_classes():
    """Return the native class registry without creating an import cycle."""
    from mhctools.pred import Prediction

    from .allele_evidence import AlleleAttribution
    from .candidate_epitope import CandidateEpitope, Peptide
    from .vaccine_antigen import SelfReferenceMatch, SelfReferenceSource

    return {
        ("builtins", "set"): set,
        ("builtins", "tuple"): tuple,
        ("mhctools.pred", "Prediction"): Prediction,
        ("vaxrank.allele_evidence", "AlleleAttribution"): AlleleAttribution,
        ("vaxrank.candidate_epitope", "CandidateEpitope"): CandidateEpitope,
        ("vaxrank.candidate_epitope", "Peptide"): Peptide,
        ("vaxrank.vaccine_antigen", "SelfReferenceMatch"): SelfReferenceMatch,
        ("vaxrank.vaccine_antigen", "SelfReferenceSource"): SelfReferenceSource,
    }


def _json_loads(payload):
    """Parse JSON while rejecting ambiguous duplicate object keys."""
    def unique_object(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("Duplicate key in native serialized JSON: %r" % key)
            result[key] = value
        return result

    return json.loads(payload, object_pairs_hook=unique_object)


def _class_from_metadata(metadata, registry):
    """Resolve exact class metadata through the static native registry."""
    if not isinstance(metadata, dict):
        raise ValueError("Native serialized class metadata must be an object")
    if set(metadata) != {"__module__", "__name__"}:
        raise ValueError("Malformed native serialized class metadata")
    class_key = (metadata["__module__"], metadata["__name__"])
    try:
        return registry[class_key]
    except KeyError:
        raise ValueError(
            "Native serialized class %r is not permitted"
            % "%s.%s" % class_key
        ) from None


def _decode_native_repr(value, registry):
    """Decode one representation without importing input-selected objects."""
    if value is None or isinstance(value, (bool, float, int, str)):
        return value
    if isinstance(value, list):
        return [_decode_native_repr(item, registry) for item in value]
    if not isinstance(value, dict):
        raise TypeError("Invalid native serialized value %r" % (value,))

    # Native epitope graphs contain no standalone class or function values.
    if "__module__" in value or "__name__" in value:
        raise ValueError("Standalone class/function metadata is not permitted")

    class_metadata = value.get(_CLASS_FIELD)
    class_object = (
        _class_from_metadata(class_metadata, registry)
        if class_metadata is not None else None
    )

    serialized_keys = value.get(_SERIALIZED_KEYS_FIELD, [])
    if not isinstance(serialized_keys, list):
        raise ValueError("Native serialized dictionary keys must be a list")
    decoded_keys = []
    for serialized_key in serialized_keys:
        if not isinstance(serialized_key, str):
            raise ValueError("A native serialized dictionary key must be JSON")
        decoded_keys.append(_decode_native_repr(
            _json_loads(serialized_key), registry))

    state = {}
    used_key_indexes = set()
    for key, item in value.items():
        if key in {_CLASS_FIELD, _SERIALIZED_KEYS_FIELD}:
            continue
        if key.startswith(_SERIALIZED_KEY_PREFIX):
            index_text = key[len(_SERIALIZED_KEY_PREFIX):]
            try:
                index = int(index_text)
                decoded_key = decoded_keys[index]
            except (IndexError, ValueError):
                raise ValueError(
                    "Invalid native serialized dictionary key reference %r" % key
                ) from None
            if index < 0 or index in used_key_indexes:
                raise ValueError(
                    "Invalid native serialized dictionary key reference %r" % key
                )
            key = decoded_key
            used_key_indexes.add(index)
        decoded_item = _decode_native_repr(item, registry)
        if key in state:
            raise ValueError("Duplicate decoded native dictionary key %r" % (key,))
        state[key] = decoded_item
    if len(used_key_indexes) != len(decoded_keys):
        raise ValueError("Unused native serialized dictionary key metadata")

    if class_object is None:
        return state
    if class_object in (set, tuple):
        if set(state) != {_VALUE_FIELD}:
            raise ValueError("Malformed native serialized collection")
        return class_object(state[_VALUE_FIELD])
    if _VALUE_FIELD in state:
        raise ValueError("Unexpected value wrapper for native serialized class")
    if hasattr(class_object, "from_dict"):
        return class_object.from_dict(state)
    return class_object(**state)


def to_native_json(value):
    """Encode one native object after validating its complete class graph."""
    payload = to_json(value)
    restored = _decode_native_repr(
        _json_loads(payload), _native_serializable_classes())
    if type(restored) is not type(value):
        raise ValueError(
            "Native payload decoded to %s, expected %s"
            % (type(restored).__name__, type(value).__name__)
        )
    return payload


def from_native_json(payload, expected_type):
    """Decode one allowlisted native object and enforce its root type."""
    value = _decode_native_repr(
        _json_loads(payload), _native_serializable_classes())
    if type(value) is not expected_type:
        raise ValueError(
            "Native payload decoded to %s, expected %s"
            % (type(value).__name__, expected_type.__name__)
        )
    return value
