"""Release ownership and recovery against independent clones of a real remote."""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import sys
from threading import Barrier, Event

import pytest

from release_deploy import Release, deploy
import release_deploy
from release_upload import ReleaseUploadError, file_sha256


def git(repo, *args):
    return subprocess.run(["git", "-C", str(repo), *args], check=True,
                          capture_output=True, text=True).stdout.strip()


@pytest.fixture
def releases(tmp_path):
    remote = tmp_path / "origin.git"
    git(tmp_path, "init", "--bare", str(remote))
    first = tmp_path / "first"
    git(tmp_path, "clone", str(remote), str(first))
    git(first, "checkout", "-b", "main")
    git(first, "config", "user.name", "Release test")
    git(first, "config", "user.email", "release@example.invalid")
    (first / ".gitignore").write_text("dist/\n")
    git(first, "add", ".")
    git(first, "commit", "--no-gpg-sign", "-m", "Release source")
    git(first, "push", "-u", "origin", "main")
    second = tmp_path / "second"
    git(tmp_path, "clone", "--branch", "main", str(remote), str(second))
    git(second, "config", "user.name", "Release test")
    git(second, "config", "user.email", "release@example.invalid")
    return [Release(path, "vaxrank", "1.0.0") for path in (first, second)]


@pytest.fixture
def package_io():
    published, builds, uploads = {}, [], []

    def build(output):
        builds.append(output)
        output.mkdir(exist_ok=True, parents=True)
        for name in ("vaxrank-1.0.0-py3-none-any.whl", "vaxrank-1.0.0.tar.gz"):
            # Different builds really have different bytes, like archive timestamps.
            (output / name).write_text("build %s: %s" % (len(builds), name))

    def upload(path):
        uploads.append(path)
        assert path.name not in published
        published[path.name] = file_sha256(path)

    return dict(fetch_artifacts=lambda: dict(published), build_artifacts=build,
                check_artifacts=lambda paths: None, upload_file=upload), published, builds, uploads


def test_concurrent_clones_build_and_publish_exactly_once(releases, package_io, monkeypatch, capsys):
    io, published, builds, uploads = package_io
    start, collision = Barrier(2), Event()
    original_git = Release.git

    def observed_git(self, *args, **kwargs):
        result = original_git(self, *args, **kwargs)
        if "--force-with-lease=%s:" % self.lock_ref in args and result.returncode:
            collision.set()
        return result

    monkeypatch.setattr(Release, "git", observed_git)
    original_build = io["build_artifacts"]

    def build(output):
        assert collision.wait(10), "The other clone must contend before this build finishes"
        original_build(output)

    def run(release):
        start.wait(timeout=10)
        deploy(release, **dict(io, build_artifacts=build), poll_seconds=0.02)

    with ThreadPoolExecutor(max_workers=2) as executor:
        list(executor.map(run, releases))
    assert len(builds) == 1 and len(uploads) == 2
    for release in releases:
        assert release.remote_ref(release.lock_ref) is None
        assert release.completed(release.preflight(), io["fetch_artifacts"])
    assert set(published) == releases[0].filenames
    assert "Already published and verified" in capsys.readouterr().out


def test_partial_upload_resumes_only_original_bytes(releases, package_io):
    io, published, builds, uploads = package_io
    release = releases[0]
    original_upload = io["upload_file"]

    def interrupted(path):
        if published:
            raise RuntimeError("interrupted upload")
        original_upload(path)

    with pytest.raises(RuntimeError, match="interrupted upload"):
        deploy(release, **dict(io, upload_file=interrupted))
    assert len(published) == 1
    assert release.local_tag() and release.remote_ref(release.tag_ref) is None
    assert release.remote_ref(release.lock_ref) is None
    with pytest.raises(ReleaseUploadError, match="Verify the published artifacts and tag"):
        deploy(releases[1], **io)
    assert len(builds) == 1 and len(uploads) == 1
    deploy(release, **io)
    assert len(builds) == 1 and len(uploads) == 2
    assert len(published) == 2


@pytest.mark.parametrize("damage", ["missing", "changed"])
@pytest.mark.parametrize("dry_run", [False, True])
def test_resume_rejects_missing_or_changed_originals(releases, package_io, damage, dry_run):
    io, published, builds, uploads = package_io
    release = releases[0]

    def interrupted(path):
        raise RuntimeError("interrupted")

    with pytest.raises(RuntimeError):
        deploy(release, **dict(io, upload_file=interrupted))
    artifact = next((release.repo / "dist").iterdir())
    if damage == "missing":
        artifact.unlink()
    else:
        artifact.write_text("rebuilt with different bytes")
    with pytest.raises(ReleaseUploadError, match="Cannot safely resume"):
        deploy(release, **io, dry_run=dry_run)
    assert len(builds) == 1 and not published and not uploads
    assert release.remote_ref(release.lock_ref) is None


def test_completed_release_rejects_pypi_hash_mismatch(releases, package_io):
    io, published, builds, uploads = package_io
    deploy(releases[0], **io)
    published[next(iter(published))] = "0" * 64
    with pytest.raises(ReleaseUploadError, match="SHA-256"):
        deploy(releases[1], **io)
    assert len(builds) == 1 and len(uploads) == 2


def test_preexisting_publication_without_tag_is_not_rebuilt(releases, package_io):
    io, published, builds, uploads = package_io
    published["vaxrank-1.0.0.tar.gz"] = "0" * 64
    with pytest.raises(ReleaseUploadError, match="Verify the published artifacts and tag"):
        deploy(releases[0], **io)
    assert not builds and not uploads
    assert releases[0].remote_ref(releases[0].lock_ref) is None


def test_legacy_release_tag_requires_manual_verification(releases, package_io):
    io, _, builds, uploads = package_io
    release = releases[0]
    git(release.repo, "tag", release.tag)
    git(release.repo, "push", "origin", release.tag_ref)
    with pytest.raises(ReleaseUploadError, match="no recorded artifact hashes"):
        deploy(releases[1], **io)
    assert not builds and not uploads


def test_remote_branch_is_rechecked_after_quality_gates(releases, package_io):
    io, _, builds, uploads = package_io
    git(releases[0].repo, "commit", "--no-gpg-sign", "--allow-empty", "-m", "New main")
    git(releases[0].repo, "push", "origin", "main")
    with pytest.raises(ReleaseUploadError, match="Remote release branch changed"):
        deploy(releases[1], **io)
    assert not builds and not uploads


def test_source_change_during_build_prevents_upload(releases, package_io):
    io, _, _, uploads = package_io
    release = releases[0]
    original_build = io["build_artifacts"]

    def build(output):
        original_build(output)
        (release.repo / "unexpected-source-change").write_text("changed")

    with pytest.raises(ReleaseUploadError, match="Working tree not clean"):
        deploy(release, **dict(io, build_artifacts=build))
    assert not uploads and release.local_tag() is None
    assert release.remote_ref(release.lock_ref) is None


def test_stale_claim_is_never_stolen(releases, package_io):
    io, _, builds, uploads = package_io
    first, second = releases
    owner = first.acquire(first.preflight(), io["fetch_artifacts"], 0, 0)
    try:
        with pytest.raises(ReleaseUploadError, match="only then remove the stale claim"):
            deploy(second, **io, wait_seconds=0)
        assert second.remote_ref(second.lock_ref) == owner
        assert not builds and not uploads
    finally:
        first.release_owner(owner)


@pytest.mark.parametrize("timing", ["before_check", "during_push"])
def test_cleanup_cannot_delete_another_owner(releases, package_io, monkeypatch, timing):
    io, _, _, _ = package_io
    release = releases[0]
    owner = release.acquire(release.preflight(), io["fetch_artifacts"], 0, 0)
    # Simulate an operator replacing a stale claim while the old process resumes.
    original_git = release.git

    def replace_claim():
        git(release.repo, "push", "--force", "origin", "HEAD:" + release.lock_ref)

    def raced_push(*args, **kwargs):
        if ":" + release.lock_ref in args:
            replace_claim()
        return original_git(*args, **kwargs)

    if timing == "before_check":
        replace_claim()
    else:
        monkeypatch.setattr(release, "git", raced_push)
    with pytest.raises(ReleaseUploadError, match="ownership changed|Could not release ownership"):
        release.release_owner(owner)
    assert release.remote_ref(release.lock_ref) == git(release.repo, "rev-parse", "HEAD")


def test_dry_run_keeps_refs_and_existing_dist_untouched(releases, package_io):
    io, _, builds, uploads = package_io
    release = releases[0]
    original = release.repo / "dist" / "original-build"
    original.parent.mkdir()
    original.write_bytes(b"keep")
    before = (git(release.repo, "show-ref"), git(release.repo, "ls-remote", "origin"))
    deploy(release, **io, dry_run=True)
    assert (git(release.repo, "show-ref"), git(release.repo, "ls-remote", "origin")) == before
    assert original.read_bytes() == b"keep"
    assert len(builds) == 1 and not builds[0].exists() and not uploads


def test_lost_lock_push_acknowledgement_is_reconciled(releases, package_io, monkeypatch):
    io, _, builds, uploads = package_io
    release = releases[0]
    original_git = release.git

    def lost_response(*args, **kwargs):
        result = original_git(*args, **kwargs)
        if "--force-with-lease=%s:" % release.lock_ref in args:
            result.returncode = 1
            result.stderr = "connection lost after push"
        return result

    monkeypatch.setattr(release, "git", lost_response)
    deploy(release, **io)
    assert len(builds) == 1 and len(uploads) == 2
    assert release.remote_ref(release.lock_ref) is None


def test_remote_claim_rejection_fails_before_build(releases, package_io):
    io, _, builds, uploads = package_io
    release = releases[0]
    remote = release.repo.parent / "origin.git"
    hook = remote / "hooks" / "pre-receive"
    hook.write_text('#!/bin/sh\necho "claim denied" >&2\nexit 1\n')
    hook.chmod(0o755)
    with pytest.raises(ReleaseUploadError, match="claim denied"):
        deploy(release, **io)
    assert not builds and not uploads and release.remote_ref(release.lock_ref) is None


@pytest.mark.parametrize("stage", ["build_artifacts", "check_artifacts"])
def test_build_failure_releases_ownership_without_creating_tag(releases, package_io, stage):
    io, _, _, uploads = package_io

    def failure(*args):
        raise RuntimeError("artifact failure")

    with pytest.raises(RuntimeError, match="artifact failure"):
        deploy(releases[0], **dict(io, **{stage: failure}))
    assert not uploads
    assert releases[0].local_tag() is None
    assert releases[0].remote_ref(releases[0].lock_ref) is None


def test_release_tag_for_another_commit_is_not_accepted(releases, package_io):
    io, _, builds, uploads = package_io
    release = releases[0]
    deploy(release, **io)
    git(release.repo, "commit", "--no-gpg-sign", "--allow-empty", "-m", "Different source")
    git(release.repo, "push", "origin", "main")
    with pytest.raises(ReleaseUploadError, match="does not point to HEAD"):
        deploy(release, **io)
    assert len(builds) == 1 and len(uploads) == 2


def test_completed_tag_with_missing_public_artifact_fails(releases, package_io):
    io, published, builds, uploads = package_io
    deploy(releases[0], **io)
    published.pop(next(iter(published)))
    releases[1].verify_timeout_seconds = 0
    with pytest.raises(ReleaseUploadError, match="Published artifact set"):
        deploy(releases[1], **io)
    assert len(builds) == 1 and len(uploads) == 2


@pytest.mark.parametrize("message", ["not JSON", "{}"])
def test_invalid_release_manifest_is_not_trusted(releases, package_io, message):
    io, _, builds, uploads = package_io
    first, second = releases
    git(first.repo, "tag", "--no-sign", "-a", "-m", message, first.tag)
    git(first.repo, "push", "origin", first.tag_ref)
    with pytest.raises(ReleaseUploadError, match="Invalid artifact manifest"):
        deploy(second, **io)
    assert not builds and not uploads


def test_cli_dry_run_drives_build_and_check_without_publication(releases, monkeypatch, capsys):
    release = releases[0]
    monkeypatch.chdir(release.repo)
    monkeypatch.setattr(release_deploy, "pypi_release_artifacts", lambda *args: {})
    original_run = subprocess.run
    commands = []

    def run(command, *args, **kwargs):
        if command[:3] == [sys.executable, "-m", "build"]:
            commands.append(command)
            assert command[3] == "--outdir"
            assert kwargs["cwd"] != release.repo
            assert (kwargs["cwd"] / ".gitignore").read_text() == "dist/\n"
            assert not (kwargs["cwd"] / ".git").exists()
            for name in release.filenames:
                (Path(command[4]) / name).write_bytes(b"artifact")
            return subprocess.CompletedProcess(command, 0)
        if command[:3] == [sys.executable, "-m", "twine"]:
            commands.append(command)
            assert command[3] == "check"
            assert {Path(path).name for path in command[4:]} == release.filenames
            assert all(Path(path).read_bytes() == b"artifact" for path in command[4:])
            return subprocess.CompletedProcess(command, 0)
        return original_run(command, *args, **kwargs)

    monkeypatch.setattr(subprocess, "run", run)
    before = (git(release.repo, "show-ref"), git(release.repo, "ls-remote", "origin"))
    assert release_deploy.main(["--version", "1.0.0", "--dry-run"]) == 0
    assert "Dry run: ready to release v1.0.0" in capsys.readouterr().out
    assert len(commands) == 2
    assert (git(release.repo, "show-ref"), git(release.repo, "ls-remote", "origin")) == before
