"""Inspect actual distributions and open installed Sid resources offline."""

import gzip
import io
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
import tarfile
import zipfile


ROOT = Path(__file__).resolve().parents[1]


def test_distributions_include_only_the_selected_sid_bundle(tmp_path):
    result = subprocess.run([
        sys.executable, "setup.py", "--quiet", "sdist", "--dist-dir", str(tmp_path),
        "bdist_wheel", "--dist-dir", str(tmp_path)], cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    sdist_path, = tmp_path.glob("vaxrank-*.tar.gz")
    wheel_path, = tmp_path.glob("vaxrank-*.whl")
    with tarfile.open(sdist_path, "r:gz") as source:
        sdist = {str(PurePosixPath(*PurePosixPath(m.name).parts[1:])): source.extractfile(m).read()
                 for m in source if m.isfile()}
    with zipfile.ZipFile(wheel_path) as source:
        wheel = {name: source.read(name) for name in source.namelist()}
    assert not any(name.startswith("tests/") for name in sdist)
    assert {"LICENSE", "MANIFEST.in", "setup.py", "requirements.txt", "TEST_DATA.md",
            "vaxrank/sid_test_data.py", "examples/osteosarc_test_data/build.py",
            "examples/osteosarc_test_data/pin_catalog.py",
            "examples/osteosarc_test_data/recipe/selection.json.gz"} <= set(sdist)
    bundle = "vaxrank/data/sid-test-data.zip"
    assert sdist[bundle] == wheel[bundle] == (ROOT / bundle).read_bytes()
    for payload in (sdist, wheel):
        for name, data in payload.items():
            if name.endswith((".bam", ".bai", ".cram", ".fastq", ".fastq.gz", ".sam", ".sam.gz")):
                # The generator's original SAM headers contain no read records.
                assert name.startswith("examples/osteosarc_test_data/recipe/headers/")
                assert all(line.startswith(b"@") for line in gzip.decompress(data).splitlines())
        with zipfile.ZipFile(io.BytesIO(payload[bundle])) as source:
            assert "provenance.json" in source.namelist()
            assert len([n for n in source.namelist() if n.endswith(".bam")]) == 54
    installed = tmp_path / "installed"
    installed.mkdir()
    with zipfile.ZipFile(wheel_path) as source:
        source.extractall(installed)
    # Import from the extracted wheel outside the checkout with every network
    # connection forbidden and no persistent data/reference caches available.
    env = dict(os.environ, PYTHONPATH=str(installed) + os.pathsep + os.environ.get("PYTHONPATH", ""),
               OPENVAX_DATA_CACHE=str(tmp_path / "empty-cache"), OSTEOSARC_CACHE=str(tmp_path / "empty-cache"))
    script = '''
import json, socket
from pathlib import Path
socket.socket.connect = lambda *a, **k: (_ for _ in ()).throw(AssertionError("network"))
import vaxrank
from vaxrank.sid_test_data import sid_test_data, sid_reads
assert Path(vaxrank.__file__).is_relative_to(Path.cwd() / "installed")
root = sid_test_data()
m = json.loads((root / "osteosarc/shared-v1/manifest.json").read_text())
assert len(m["cases"]) == 49
for case in m["cases"]:
    with sid_reads("osteosarc/shared-v1/" + case["bam"]).open() as bam:
        assert bam.check_index() and sum(1 for _ in bam) == case["selected_record_count"]
'''
    result = subprocess.run([sys.executable, "-c", script], cwd=tmp_path, env=env,
                            text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert not (tmp_path / "empty-cache").exists()
