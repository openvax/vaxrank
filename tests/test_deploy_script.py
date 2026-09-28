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

import os
import shutil
import subprocess
from pathlib import Path

import pytest


def _write_fake_python(path):
    path.write_text(
        "#!/bin/sh\n"
        "printf '%s\\n' \"$@\" > \"$PYTHON_LOG\"\n",
        encoding="utf-8",
    )
    path.chmod(0o755)


def _deploy_resume_harness(tmp_path, *, include_sdist):
    repo_root = Path(__file__).resolve().parents[1]
    for filename in ("deploy.sh", "lint.sh", "test.sh"):
        shutil.copy2(repo_root / filename, tmp_path / filename)

    fake_bin = tmp_path / "fake-bin"
    fake_bin.mkdir()
    python_log = tmp_path / "python-log"
    fake_python = fake_bin / "python"
    fake_python.write_text(
        "#!/bin/sh\n"
        "if [ \"$1\" = \"-\" ]; then\n"
        "  echo 3.13.9\n"
        "  exit 0\n"
        "fi\n"
        "printf '%s\\n' \"$*\" >> \"$PYTHON_LOG\"\n",
        encoding="utf-8",
    )
    fake_python.chmod(0o755)

    fake_git = fake_bin / "git"
    fake_git.write_text(
        "#!/bin/sh\n"
        "case \"$1 $2\" in\n"
        "  'status --porcelain') exit 0 ;;\n"
        "  'rev-parse HEAD') echo release-commit; exit 0 ;;\n"
        "  'rev-parse @{upstream}') echo release-commit; exit 0 ;;\n"
        "  'rev-parse -q') exit 0 ;;\n"
        "  'rev-list -n') echo release-commit; exit 0 ;;\n"
        "esac\n"
        "exit 0\n",
        encoding="utf-8",
    )
    fake_git.chmod(0o755)

    dist = tmp_path / "dist"
    dist.mkdir()
    (dist / "vaxrank-3.13.9-py3-none-any.whl").write_bytes(b"wheel")
    if include_sdist:
        (dist / "vaxrank-3.13.9.tar.gz").write_bytes(b"sdist")

    env = os.environ.copy()
    env.update(
        PATH=f"{fake_bin}{os.pathsep}{env['PATH']}",
        PYTHON=str(fake_python),
        PYTHON_LOG=str(python_log),
        VAXRANK_DEPLOY_BRANCH="main",
    )
    result = subprocess.run(
        ["bash", "deploy.sh", "--dry-run"],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    return result, python_log


def test_deploy_rejects_non_release_branch():
    repo_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["VAXRANK_DEPLOY_BRANCH"] = "yaml-config"
    result = subprocess.run(
        ["bash", "deploy.sh", "--dry-run"],
        cwd=repo_root,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 1
    assert "only allowed from main or master" in result.stderr


@pytest.mark.parametrize(
    ("script", "extra_args", "expected_args"),
    [
        (
            "lint.sh",
            [],
            ["-m", "ruff", "check", "vaxrank", "tests", "release_upload.py",
             "check_release_version.py", "examples/osteosarc_read_corpus"],
        ),
        ("test.sh", ["-q"], ["-m", "pytest", "tests", "-q"]),
    ],
)
def test_quality_scripts_use_the_configured_python(
        tmp_path, script, extra_args, expected_args):
    repo_root = Path(__file__).resolve().parents[1]
    fake_python = tmp_path / "release-python"
    log_path = tmp_path / "python-args"
    _write_fake_python(fake_python)
    env = os.environ.copy()
    env.update(PYTHON=str(fake_python), PYTHON_LOG=str(log_path))

    subprocess.run(
        ["bash", script, *extra_args],
        cwd=repo_root,
        env=env,
        check=True,
    )

    assert log_path.read_text(encoding="utf-8").splitlines() == expected_args


def test_quality_scripts_prefer_the_active_environment_over_local_venv(
        tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    active_venv = tmp_path / "active-venv"
    fake_python = active_venv / "bin" / "python"
    fake_python.parent.mkdir(parents=True)
    log_path = tmp_path / "python-args"
    _write_fake_python(fake_python)
    env = os.environ.copy()
    env.pop("PYTHON", None)
    env.update(VIRTUAL_ENV=str(active_venv), PYTHON_LOG=str(log_path))

    subprocess.run(
        ["bash", "lint.sh"],
        cwd=repo_root,
        env=env,
        check=True,
    )

    assert log_path.read_text(encoding="utf-8").splitlines()[:3] == [
        "-m", "ruff", "check",
    ]


def test_deploy_resume_reuses_the_original_artifacts(tmp_path):
    result, python_log = _deploy_resume_harness(
        tmp_path,
        include_sdist=True,
    )

    assert result.returncode == 0, result.stderr
    assert "Resuming release from v3.13.9" in result.stdout
    assert "-m build" not in python_log.read_text(encoding="utf-8")


def test_deploy_resume_rejects_a_missing_original_artifact(tmp_path):
    result, _ = _deploy_resume_harness(tmp_path, include_sdist=False)

    assert result.returncode == 1
    assert "Cannot safely resume v3.13.9" in result.stderr
    assert "vaxrank-3.13.9.tar.gz" in result.stderr


@pytest.mark.parametrize("status", [0, 1])
@pytest.mark.parametrize("explicit_root", [False, True])
def test_test_script_owns_only_its_automatic_temp_root(tmp_path, status, explicit_root):
    fake_python = tmp_path / "python"
    fake_python.write_text(
        '#!/bin/sh\n'
        'printf "%s" "$PYTEST_DEBUG_TEMPROOT" > "$PYTHON_LOG"\n'
        'touch "$PYTEST_DEBUG_TEMPROOT/evidence"\n'
        'exit "$TEST_STATUS"\n')
    fake_python.chmod(0o755)
    env = os.environ.copy()
    env.pop("PYTEST_DEBUG_TEMPROOT", None)
    env.update(PYTHON=str(fake_python), PYTHON_LOG=str(tmp_path / "root"),
               TMPDIR=str(tmp_path), TEST_STATUS=str(status))
    if explicit_root:
        supplied = tmp_path / "supplied"
        supplied.mkdir()
        env["PYTEST_DEBUG_TEMPROOT"] = str(supplied)
    result = subprocess.run(
        ["bash", "test.sh"], cwd=Path(__file__).resolve().parents[1],
        env=env, capture_output=True, text=True)
    root = Path((tmp_path / "root").read_text())
    assert result.returncode == status
    assert (root / "evidence").exists() == (explicit_root or status != 0)
    if explicit_root:
        assert root == supplied
    else:
        assert root.parent == tmp_path
        assert root.name.startswith("vaxrank-pytest.")
