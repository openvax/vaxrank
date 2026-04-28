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

"""Guard the invariant that vaxrank never imports pysam directly.

Issue #220 dropped the standalone ``pysam>=0.23.0`` line from
requirements.txt: vaxrank reads BAMs only through isovar, which owns
its own pysam pin. If a future change adds ``import pysam`` to a
vaxrank module, the package will still work in dev (pysam is
transitively installed) but requirements.txt no longer documents the
dependency. This test fails before that drift can ship.
"""

import ast
import pathlib


VAXRANK_ROOT = pathlib.Path(__file__).resolve().parent.parent / "vaxrank"


def _python_files(root):
    for path in root.rglob("*.py"):
        if "__pycache__" in path.parts:
            continue
        yield path


def _imports_pysam(source):
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == "pysam" or alias.name.startswith("pysam."):
                    return True
        elif isinstance(node, ast.ImportFrom):
            if node.module == "pysam" or (
                    node.module and node.module.startswith("pysam.")):
                return True
    return False


def test_no_direct_pysam_import():
    offenders = []
    for path in _python_files(VAXRANK_ROOT):
        if _imports_pysam(path.read_text()):
            offenders.append(str(path.relative_to(VAXRANK_ROOT.parent)))
    assert not offenders, (
        "vaxrank must not import pysam directly (see #220 — isovar owns the "
        "pysam dep). Offending files:\n  " + "\n  ".join(offenders))
