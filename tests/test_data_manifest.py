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

"""The checked-in fixtures under ``tests/data`` must match their digests.

Osteosarc supplies the Sid reads and variant identities, but it publishes the
human dataset only and has no epitope-table product, so the B16-F10 mouse data
and the LENS/pVACseq tables stay in this repository. Nothing verified them: a
silent edit, a truncated checkout or a stray file went unnoticed, and two tests
used to skip themselves when a fixture was missing, so an absent fixture read as
a pass.

``tests/data/manifest.json`` pins every one of those files by sha256. Regenerate
it deliberately, and only when the fixture change is the point of the commit:

    python tests/test_data_manifest.py
"""

import hashlib
import json
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"
MANIFEST_PATH = DATA_DIR / "manifest.json"


def _load_manifest():
    return json.loads(MANIFEST_PATH.read_text())


def _digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _present_files():
    """Every fixture file actually on disk, as manifest-relative posix paths.

    Dotfiles (``.DS_Store``) and ``__pycache__`` are editor/interpreter debris
    rather than fixtures, so they are not expected in the manifest.
    """
    found = set()
    for path in DATA_DIR.rglob("*"):
        if not path.is_file():
            continue
        relative = path.relative_to(DATA_DIR)
        if relative.as_posix() == "manifest.json":
            continue
        if any(part.startswith(".") or part == "__pycache__"
               for part in relative.parts):
            continue
        found.add(relative.as_posix())
    return found


def test_manifest_lists_every_fixture_and_nothing_else():
    """The manifest and the directory must agree exactly.

    Listing a file that is gone hides a deletion; carrying a file nobody
    listed is how unpinned test data crept in to begin with.
    """
    listed = set(_load_manifest()["files"])
    present = _present_files()

    missing = sorted(listed - present)
    assert not missing, (
        "manifest.json lists %d file(s) that are not on disk: %s"
        % (len(missing), ", ".join(missing)))

    unlisted = sorted(present - listed)
    assert not unlisted, (
        "%d fixture file(s) under tests/data are not pinned in manifest.json: "
        "%s. Add them (see this module's docstring) or delete them."
        % (len(unlisted), ", ".join(unlisted)))


def test_every_fixture_matches_its_recorded_digest():
    """Content, not just presence. A fixture edited by accident — or a
    partial checkout of a binary BAM — fails here rather than surfacing as
    an unrelated assertion in whichever test happens to read it."""
    changed = []
    for relative, expected in sorted(_load_manifest()["files"].items()):
        actual = _digest(DATA_DIR / relative)
        if actual != expected:
            changed.append((relative, expected, actual))

    assert not changed, "\n".join(
        "%s: expected sha256 %s, found %s" % entry for entry in changed)


def test_manifest_explains_why_each_group_is_not_osteosarc_sourced():
    """Each top-level fixture group carries a note. These files are the
    documented exceptions to "osteosarc is the test-data source", so the
    reason travels with the digests instead of living only in a commit
    message."""
    manifest = _load_manifest()
    notes = manifest["notes"]
    groups = {relative.split("/")[0] for relative in manifest["files"]}
    undocumented = sorted(groups - set(notes))
    assert not undocumented, (
        "fixture group(s) with no note in manifest.json: %s"
        % ", ".join(undocumented))
    for group, note in notes.items():
        assert note.strip(), "empty note for fixture group %s" % group


def _regenerate():
    """Rewrite the digests in place, keeping the existing notes.

    Deliberately does not invent a note for a new fixture group — the
    coverage test fails until a maintainer writes one.
    """
    manifest = _load_manifest()
    manifest["files"] = {
        relative: _digest(DATA_DIR / relative)
        for relative in sorted(_present_files())
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n")
    print("pinned %d fixture(s) in %s" % (len(manifest["files"]), MANIFEST_PATH))


if __name__ == "__main__":
    _regenerate()
