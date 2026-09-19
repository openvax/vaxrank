"""Exercise the release check against committed files and real Git history."""

import os
from pathlib import Path
import subprocess
import sys

import pytest
import yaml


ROOT = Path(__file__).resolve().parents[1]
CHECKER = ROOT / "check_release_version.py"


def git(repo, *args):
    return subprocess.run(
        ["git", "-C", str(repo), *args], check=True, capture_output=True,
        text=True).stdout.strip()


def commit_version(repo, version):
    (repo / "vaxrank/version.py").write_text(
        '__version__ = "%s"\n' % version, encoding="utf-8")
    git(repo, "add", ".")
    git(repo, "commit", "--no-gpg-sign", "--allow-empty", "-m", "Set " + version)
    return git(repo, "rev-parse", "HEAD")


@pytest.fixture
def repo(tmp_path):
    git(tmp_path, "init", "-b", "main")
    git(tmp_path, "config", "user.email", "version-test@example.invalid")
    git(tmp_path, "config", "user.name", "Version Check Test")
    (tmp_path / "vaxrank").mkdir()
    commit_version(tmp_path, "3.19.6")
    git(tmp_path, "checkout", "-b", "topic")
    return tmp_path


def check(repo, base="main", head="HEAD"):
    return subprocess.run(
        [sys.executable, str(CHECKER), "--base", base, "--head", head],
        cwd=repo, capture_output=True, text=True, check=False)


@pytest.mark.parametrize("base,head", [
    ("3.19.6", "3.19.7"),
    ("3.9.9", "3.10.0"),
    ("3.20.0rc1", "3.20.0rc2"),
    ("3.20.0rc2", "3.20.0"),
])
def test_unreleased_version_increase_passes(repo, base, head):
    base_commit = commit_version(repo, base)
    commit_version(repo, head)
    git(repo, "tag", "v" + base, base_commit)
    git(repo, "tag", "validation-only")
    result = check(repo, base=base_commit)
    assert result.returncode == 0, result.stderr
    assert "has no existing tag" in result.stdout


@pytest.mark.parametrize("version", ["3.19.6", "3.19.6.0", "3.19.5", "3.19.6rc1"])
def test_missing_or_stale_bump_fails(repo, version):
    commit_version(repo, version)
    result = check(repo)
    assert result.returncode == 1
    assert "must be greater than target branch version 3.19.6" in result.stderr


@pytest.mark.parametrize("annotated", [False, True])
@pytest.mark.parametrize("tag", ["v3.19.7", "v3.19.7.0"])
def test_existing_release_tag_fails_even_on_another_commit(repo, annotated, tag):
    args = ["tag", "--no-sign"]
    if annotated:
        args += ["-a", "-m", "Already released"]
    git(repo, *args, tag, "main")
    commit_version(repo, "3.19.7")
    result = check(repo)
    assert result.returncode == 1
    assert "already tagged as " + tag in result.stderr


@pytest.mark.parametrize("source,message", [
    ("other = '3.19.7'\n", "must assign __version__ exactly once"),
    ("__version__ = 3197\n", "must assign __version__ exactly once"),
    ("__version__ = '3.19.' + '7'\n", "must assign __version__ exactly once"),
    ("__version__ = '3.19.7'\n__version__ = '3.19.8'\n", "exactly once"),
    ("__version__ = 'not-a-version'\n", "Invalid release version"),
    ("__version__ = (\n", "Invalid version file"),
    ("  __version__ = '3.19.7'\n", "Invalid version file"),
])
def test_invalid_version_assignment_fails(repo, source, message):
    (repo / "vaxrank/version.py").write_text(source, encoding="utf-8")
    git(repo, "add", ".")
    git(repo, "commit", "--no-gpg-sign", "-m", "Invalid version")
    result = check(repo)
    assert result.returncode == 1
    assert message in result.stderr
    assert "Traceback" not in result.stderr


def test_version_file_is_never_executed(repo):
    (repo / "vaxrank/version.py").write_text(
        "from pathlib import Path\n"
        "Path('executed').touch()\n"
        "__version__: str = '3.19.7'\n", encoding="utf-8")
    git(repo, "add", ".")
    git(repo, "commit", "--no-gpg-sign", "-m", "Version with side effects")
    result = check(repo)
    assert result.returncode == 0, result.stderr
    assert not (repo / "executed").exists()


def test_uses_committed_version_not_working_tree(repo):
    (repo / "vaxrank/version.py").write_text(
        "__version__ = '3.19.7'\n", encoding="utf-8")
    result = check(repo)
    assert result.returncode == 1
    assert "PR version 3.19.6" in result.stderr


def test_missing_version_file_fails(repo):
    git(repo, "rm", "vaxrank/version.py")
    git(repo, "commit", "--no-gpg-sign", "-m", "Remove version")
    result = check(repo)
    assert result.returncode == 1
    assert "Git inspection failed" in result.stderr
    assert "vaxrank/version.py" in result.stderr


def test_unknown_base_fails(repo):
    commit_version(repo, "3.19.7")
    result = check(repo, base="does-not-exist")
    assert result.returncode == 1
    assert "Git inspection failed" in result.stderr


def test_clean_duplicate_version_merge_is_rejected(repo):
    old_base = git(repo, "rev-parse", "main")
    (repo / "first-change").write_text("first PR\n", encoding="utf-8")
    first = commit_version(repo, "3.19.7")
    assert check(repo).returncode == 0
    git(repo, "checkout", "-b", "second", old_base)
    (repo / "second-change").write_text("second PR\n", encoding="utf-8")
    second = commit_version(repo, "3.19.7")
    assert check(repo).returncode == 0
    git(repo, "checkout", "main")
    git(repo, "merge", "--ff-only", first)
    result = check(repo, head=second)
    assert result.returncode == 1
    assert "PR version 3.19.7 must be greater than target branch version 3.19.7" in result.stderr
    # This is the bug: Git itself accepts both PRs without a conflict.
    git(repo, "merge", "--no-gpg-sign", "--no-edit", second)
    assert (repo / "first-change").exists()
    assert (repo / "second-change").exists()


@pytest.mark.parametrize("update", ["base", "tag"])
def test_workflow_fetches_current_target_and_tags(repo, tmp_path_factory, update):
    """Run the actual workflow step against a local remote advanced after checkout."""
    (repo / "check_release_version.py").write_bytes(CHECKER.read_bytes())
    head = commit_version(repo, "3.19.7")
    checkout = tmp_path_factory.mktemp("version-check") / "checkout"
    git(repo, "clone", "--no-local", str(repo), str(checkout))
    assert check(checkout, base="origin/main", head=head).returncode == 0
    if update == "base":
        git(repo, "checkout", "main")
        commit_version(repo, "3.19.7")
    else:
        git(repo, "tag", "v3.19.7", "main")

    workflow = yaml.safe_load(
        (ROOT / ".github/workflows/release-version.yml").read_text(encoding="utf-8"))
    command = workflow["jobs"]["check"]["steps"][-1]["run"]
    result = subprocess.run(
        ["bash", "-eo", "pipefail", "-c", command], cwd=checkout,
        env=dict(os.environ, BASE_REF="main", HEAD_SHA=head,
                 PATH=str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"]),
        capture_output=True, text=True, check=False)
    assert result.returncode == 1
    assert ("must be greater" if update == "base" else "already tagged") in result.stderr
