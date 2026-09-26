"""Inspect actual distributions and build installed Sid test data offline."""

import gzip
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
import tarfile
import zipfile

import osteosarc


ROOT = Path(__file__).resolve().parents[1]


def test_distributions_include_the_sid_recipe_and_no_reads(tmp_path):
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
            "vaxrank/sid_test_data.py", "examples/osteosarc_test_data/pin_catalog.py"} <= set(sdist)
    recipe = ROOT / "vaxrank/data/sid-recipe"
    shipped = {p.relative_to(ROOT).as_posix(): p.read_bytes() for p in recipe.rglob("*") if p.is_file()}
    assert shipped
    for name, data in shipped.items():
        assert sdist[name] == wheel[name] == data, name
    for payload in (sdist, wheel):
        for name, data in payload.items():
            if name.endswith((".bam", ".bai", ".cram", ".fastq", ".fastq.gz", ".sam", ".sam.gz", ".zip")):
                # The recipe's fixture SAM headers contain no read records.
                assert name.startswith("vaxrank/data/sid-recipe/headers/"), name
                assert all(line.startswith(b"@") for line in gzip.decompress(data).splitlines())
    installed = tmp_path / "installed"
    installed.mkdir()
    with zipfile.ZipFile(wheel_path) as source:
        source.extractall(installed)
    # Import from the extracted wheel outside the checkout with every network
    # connection forbidden: only the cached openvax-v1 reads are available.
    reads = Path(osteosarc.fetch_bundle("openvax-v1"))
    cache = reads.parents[2]
    assert reads.is_relative_to(cache / "osteosarc" / "bundles")
    env = dict(os.environ, PYTHONPATH=str(installed) + os.pathsep + os.environ.get("PYTHONPATH", ""),
               OPENVAX_DATA_CACHE=str(tmp_path / "empty-cache"), OSTEOSARC_CACHE=str(cache))
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
