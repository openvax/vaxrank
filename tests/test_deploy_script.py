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
import sys
from pathlib import Path

import pytest


def _write_fake_python(path):
    path.write_text(
        "#!/bin/sh\n"
        "printf '%s\\n' \"$@\" > \"$PYTHON_LOG\"\n",
        encoding="utf-8",
    )
    path.chmod(0o755)


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
            ["-m", "ruff", "check", "vaxrank", "tests", "release_deploy.py", "release_upload.py",
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


@pytest.mark.parametrize("arguments,valid", [
    (["--dry-run"], True),
    (["--dry-run", "1.0.0"], True),
    (["v1.0.0", "--dry-run"], True),
    (["--dry-run", "1.0.1"], False),
    (["--dry-run", "v1.0.1"], False),
    (["1.0.1"], False),
    (["v1.0.1"], False),
])
def test_deploy_version_arguments_never_mutate_git(tmp_path, arguments, valid):
    """Exercise real version reads and local/remote refs without package services."""
    repo_root = Path(__file__).resolve().parents[1]
    repo = tmp_path / "checkout"
    repo.mkdir()
    remote = tmp_path / "origin.git"

    def git(*args):
        return subprocess.run(
            ["git", *args], cwd=repo, check=True, capture_output=True,
            text=True).stdout

    git("init", "-b", "main")
    git("init", "--bare", str(remote))
    git("config", "user.name", "Deploy Test")
    git("config", "user.email", "deploy-test@example.invalid")
    for filename in ("deploy.sh", "lint.sh", "test.sh"):
        shutil.copy2(repo_root / filename, repo / filename)
    (repo / "vaxrank").mkdir()
    version_path = repo / "vaxrank/version.py"
    version_path.write_text('__version__ = "1.0.0"\n', encoding="utf-8")
    git("add", ".")
    git("commit", "--no-gpg-sign", "-m", "Initial release")
    git("remote", "add", "origin", str(remote))
    git("push", "-u", "origin", "main")
    before = (version_path.read_bytes(), git("show-ref"), git("ls-remote", "origin"))

    fake_python = tmp_path / "release-python"
    python_log = tmp_path / "python-log"
    fake_python.write_text(
        "#!" + sys.executable + "\n"
        "import os, sys\n"
        "if sys.argv[1:] == ['-']:\n"
        "    os.execv(sys.executable, [sys.executable, '-'])\n"
        "with open(os.environ['PYTHON_LOG'], 'a') as log:\n"
        "    log.write(' '.join(sys.argv[1:]) + '\\n')\n"
        "if sys.argv[1:2] == ['release_deploy.py']:\n"
        "    sys.exit(0)\n"
        "if sys.argv[1:3] not in (['-m', 'pip'], ['-m', 'build'], "
        "['-m', 'ruff'], ['-m', 'pytest']):\n"
        "    sys.exit('Unexpected package publication command')\n",
        encoding="utf-8")
    fake_python.chmod(0o755)
    env = dict(os.environ, PYTHON=str(fake_python), PYTHON_LOG=str(python_log))
    env.pop("VAXRANK_DEPLOY_BRANCH", None)
    result = subprocess.run(
        ["bash", "deploy.sh", *arguments], cwd=repo, env=env,
        capture_output=True, text=True, check=False)

    assert (version_path.read_bytes(), git("show-ref"), git("ls-remote", "origin")) == before
    assert git("status", "--porcelain") == ""
    if valid:
        assert result.returncode == 0, result.stderr
        commands = python_log.read_text().splitlines()
        assert "release_deploy.py --project vaxrank --version 1.0.0 --dry-run" in commands
        assert any(command.startswith("-m ruff check") for command in commands)
        assert "-m pytest tests" in commands
    else:
        assert result.returncode == 1, result.stderr
        assert "does not match vaxrank/version.py (1.0.0)" in result.stderr
        assert "Bump the version in a PR" in result.stderr
        assert not python_log.exists()


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
