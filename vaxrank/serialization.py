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

"""JSON/dict serialization for @dataclass-based vaxrank objects.

This is the replacement path for inheriting from ``serializable.Serializable``.
``@dataclass`` already gives us ``__init__``, ``__eq__``, and ``__repr__``
(and ``__hash__`` when frozen); the ``Serializable`` base class in the
``serializable`` package also injected those, which meant classes couldn't
cleanly be dataclasses. This mixin keeps the wire-level JSON format
identical to what ``Serializable`` produced — the underlying ``serializable``
helpers pick up ``to_dict()``, ``from_dict``, ``__class__``, and
``__module__``, so migrated classes remain interoperable with any persisted
JSON (and with classes that still inherit from ``Serializable``).
"""

from dataclasses import fields
from typing import Any, ClassVar

import serializable as _serializable


class DataclassSerializable:
    """Mixin granting @dataclass classes the ``to_dict`` / ``from_dict`` /
    ``to_json`` / ``from_json`` API used elsewhere in vaxrank, without
    inheriting from the old ``Serializable`` base.

    Subclasses may set ``_SERIALIZABLE_KEYWORD_ALIASES`` to migrate old
    field names (mirrors the same hook on ``Serializable``). Map an old
    name to the new name, or to ``None`` to drop it.
    """

    _SERIALIZABLE_KEYWORD_ALIASES: ClassVar[dict[str, str | None]] = {}

    def to_dict(self) -> dict[str, Any]:
        return {f.name: getattr(self, f.name) for f in fields(self)}

    @classmethod
    def from_dict(cls, d: dict[str, Any]):
        d = dict(d)
        for old, new in cls._SERIALIZABLE_KEYWORD_ALIASES.items():
            if old in d:
                value = d.pop(old)
                if new is not None and new not in d:
                    d[new] = value
        return cls(**d)

    def to_json(self) -> str:
        return _serializable.to_json(self)

    @classmethod
    def from_json(cls, json_string: str):
        return _serializable.from_json(json_string)
