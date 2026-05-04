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

"""Pepsickle invocation that runs in an isolated subprocess.

Imported and executed via ``python -m vaxrank._pepsickle_subprocess``
from :mod:`vaxrank.processing`. The subprocess is a *fresh*
interpreter — only ``mhctools`` (which imports torch + pepsickle) and
``json`` get loaded. The parent process keeps its pandas / numpy /
pyarrow imports without colliding with torch's bundled libomp.
That's the entire reason this lives in its own module: see issue
openvax/vaxrank#266 for the macOS libomp clash diagnosis and
openvax/mhctools#200 for the upstream cleanup that will let us drop
this wrapper.

The body delegates *all* per-position scoring to ``mhctools.Pepsickle``
— this module owns no scoring logic, only the I/O contract:

  stdin (JSON):
    {
      "sequences": [str, ...],
      "human_only": bool,
      "threshold": float
    }

  stdout (JSON):
    {sequence: [float, ...], ...}      # per-position cleavage probs

  stderr:
    Diagnostic messages on import / instantiate / inference failures.

  exit codes:
    0  — success (stdout is parseable JSON)
    2  — import_error (mhctools / pepsickle / torch not installed)
    3  — instantiate_error (Pepsickle constructor failed)
"""

import json
import sys


def main(stdin=None, stdout=None, stderr=None):
    """Read JSON from ``stdin``, run mhctools.Pepsickle, write JSON to
    ``stdout``. Defaults to ``sys.stdin`` / ``sys.stdout`` / ``sys.stderr``.

    Exposed as a function so unit tests can drive it without spawning
    a subprocess (pass StringIO / BytesIO buffers and call directly).
    """
    if stdin is None:
        stdin = sys.stdin
    if stdout is None:
        stdout = sys.stdout
    if stderr is None:
        stderr = sys.stderr
    data = json.loads(stdin.read())
    try:
        from mhctools import Pepsickle
    except Exception as e:
        stderr.write("import_error: %s\n" % e)
        return 2
    try:
        p = Pepsickle(
            human_only=data.get("human_only", False),
            threshold=data.get("threshold", 0.5))
    except Exception as e:
        stderr.write("instantiate_error: %s\n" % e)
        return 3
    out = {}
    for seq in data.get("sequences", []):
        try:
            out[seq] = [float(x) for x in p.cleavage_probs(seq)]
        except Exception as e:
            # Per-source inference failure: drop this sequence from
            # the result, don't fail the whole run. Caller treats
            # missing keys as "no signal."
            stderr.write("inference_error %r: %s\n" % (seq[:40], e))
    stdout.write(json.dumps(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
