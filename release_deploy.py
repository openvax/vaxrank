"""Coordinate release builds across clones using an owned remote Git ref."""

import argparse
from dataclasses import dataclass
import json
import os
from pathlib import Path
import re
from secrets import token_hex
import shutil
import subprocess
import sys
import tempfile
import time

from release_upload import (
    DEFAULT_UPLOAD_TIMEOUT_SECONDS,
    DEFAULT_VERIFY_TIMEOUT_SECONDS,
    ReleaseUploadError,
    expected_release_filenames,
    file_sha256,
    publish_release,
    pypi_release_artifacts,
    upload_distribution,
    wait_for_complete_release,
)


@dataclass
class Release:
    repo: Path
    project: str
    version: str
    verify_timeout_seconds: float = DEFAULT_VERIFY_TIMEOUT_SECONDS

    @property
    def tag(self):
        return "v" + self.version

    @property
    def tag_ref(self):
        return "refs/tags/" + self.tag

    @property
    def lock_ref(self):
        return "refs/tags/release-locks/" + self.tag

    @property
    def filenames(self):
        return expected_release_filenames(self.project, self.version)

    def git(self, *args, input=None, check=True):
        result = subprocess.run(
            ["git", "-C", str(self.repo), *args], input=input,
            capture_output=True, text=True, check=False)
        if check and result.returncode:
            raise ReleaseUploadError(result.stderr.strip() or result.stdout.strip())
        return result

    def remote_ref(self, ref):
        output = self.git("ls-remote", "--refs", "origin", ref).stdout.splitlines()
        return next((line.split()[0] for line in output if line.split()[1] == ref), None)

    def preflight(self):
        branch = self.git("branch", "--show-current").stdout.strip()
        if branch not in ("main", "master"):
            raise ReleaseUploadError("Deploys are only allowed from main or master.")
        if self.git("status", "--porcelain").stdout:
            raise ReleaseUploadError("Working tree not clean after release checks.")
        head = self.git("rev-parse", "HEAD").stdout.strip()
        if self.remote_ref("refs/heads/" + branch) != head:
            raise ReleaseUploadError("Remote release branch changed; pull main and rerun ./deploy.sh.")
        return head

    def local_tag(self):
        result = self.git("rev-parse", "--verify", "--quiet", self.tag_ref, check=False)
        if result.returncode not in (0, 1):
            raise ReleaseUploadError(result.stderr.strip())
        return result.stdout.strip() or None

    def manifest(self, oid, head):
        if self.git("rev-parse", oid + "^{commit}").stdout.strip() != head:
            raise ReleaseUploadError("Release tag %s does not point to HEAD %s." % (self.tag, head))
        if self.git("cat-file", "-t", oid).stdout.strip() != "tag":
            raise ReleaseUploadError(
                "Release %s has no recorded artifact hashes. Verify the published artifacts "
                "and tag manually; do not rebuild or overwrite this release." % self.tag)
        try:
            data = json.loads(self.git("cat-file", "tag", oid).stdout.split("\n\n", 1)[1])
            hashes = data["artifacts"]
            valid = (
                data["schema"] == 1 and data["project"] == self.project
                and data["version"] == self.version and data["commit"] == head
                and set(hashes) == self.filenames
                and all(isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value)
                        for value in hashes.values()))
        except (KeyError, IndexError, TypeError, ValueError, AttributeError):
            valid = False
        if not valid:
            raise ReleaseUploadError("Invalid artifact manifest in release tag %s." % self.tag)
        return hashes

    def completed(self, head, fetch_artifacts):
        oid = self.remote_ref(self.tag_ref)
        if oid is None:
            return False
        local = self.local_tag()
        if local and local != oid:
            raise ReleaseUploadError("Local and remote release tags differ; inspect %s." % self.tag)
        # Fetch objects only; never replace a local tag or unrelated remote refs.
        self.git("fetch", "--no-tags", "origin", self.tag_ref)
        hashes = self.manifest(oid, head)
        published = wait_for_complete_release(
            hashes, fetch_artifacts, timeout_seconds=self.verify_timeout_seconds)
        if published != hashes:
            raise ReleaseUploadError("Published artifact set does not match release tag %s." % self.tag)
        print("Already published and verified %s: remote tag and PyPI SHA-256 hashes match; "
              "no build or upload needed." % self.tag, flush=True)
        return True

    def acquire(self, head, fetch_artifacts, wait_seconds, poll_seconds):
        ident = self.git("var", "GIT_COMMITTER_IDENT").stdout.strip()
        tag = "object %s\ntype commit\ntag release-locks/%s\ntagger %s\n\nOwner %s\n" % (
            head, self.tag, ident, token_hex(16))
        owner = self.git("mktag", input=tag).stdout.strip()
        deadline = time.monotonic() + wait_seconds
        while True:
            if self.completed(head, fetch_artifacts):
                return None
            result = self.git(
                "-c", "push.followTags=false", "push", "--porcelain",
                "--force-with-lease=%s:" % self.lock_ref,
                "origin", "%s:%s" % (owner, self.lock_ref), check=False)
            observed = self.remote_ref(self.lock_ref)
            # Reconcile a lost acknowledgement without starting a second owner.
            if observed == owner:
                return owner
            if observed is None:
                raise ReleaseUploadError("Could not acquire release ownership: " + result.stderr.strip())
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise ReleaseUploadError(
                    "Another deploy owns %s at %s. Retry after it finishes. If the owner "
                    "has stopped, recover its original artifacts; only then remove the stale "
                    "claim with: git push --force-with-lease=%s:%s origin :%s" % (
                        self.lock_ref, observed, self.lock_ref, observed, self.lock_ref))
            print("Waiting for release owner %s (%s)." % (observed, self.tag), flush=True)
            time.sleep(min(poll_seconds, remaining))

    def check_owner(self, owner):
        if self.remote_ref(self.lock_ref) != owner:
            raise ReleaseUploadError("Release ownership changed; refusing to build or publish.")

    def release_owner(self, owner):
        observed = self.remote_ref(self.lock_ref)
        if observed is None:
            return
        self.check_owner(owner)
        result = self.git(
            "-c", "push.followTags=false", "push", "--porcelain",
            "--force-with-lease=%s:%s" % (self.lock_ref, owner),
            "origin", ":" + self.lock_ref, check=False)
        # A successful deletion can be followed by a transient SSH read failure.
        for attempt in range(3):
            try:
                observed = self.remote_ref(self.lock_ref)
            except ReleaseUploadError:
                if attempt == 2:
                    raise
                time.sleep(1)
                continue
            if observed is not None:
                raise ReleaseUploadError("Could not release ownership: " + result.stderr.strip())
            return

    def record_artifacts(self, head, paths):
        hashes = {path.name: file_sha256(path) for path in paths}
        message = json.dumps(dict(schema=1, project=self.project, version=self.version,
                                  commit=head, artifacts=hashes), sort_keys=True)
        self.git("tag", "--no-sign", "-a", "-F", "-", self.tag, head, input=message)
        return hashes


def verify_local_artifacts(paths, hashes):
    for path in paths:
        if not path.is_file():
            raise ReleaseUploadError("Cannot safely resume: original artifact is missing: %s" % path)
        if file_sha256(path) != hashes[path.name]:
            raise ReleaseUploadError("Cannot safely resume: original artifact SHA-256 changed: %s" % path)


def deploy(release, *, fetch_artifacts, build_artifacts, check_artifacts, upload_file,
           dry_run=False, wait_seconds=120, poll_seconds=2):
    """Build and publish once; another clone verifies the completed release."""
    head = release.preflight()
    if release.completed(head, fetch_artifacts):
        return
    if dry_run:
        local = release.local_tag()
        if local:
            paths = [release.repo / "dist" / name for name in sorted(release.filenames)]
            verify_local_artifacts(paths, release.manifest(local, head))
            check_artifacts(paths)
            print("Resuming release from %s with its original artifacts." % release.tag)
        else:
            if fetch_artifacts():
                raise ReleaseUploadError("PyPI already contains %s; verify its artifacts and tag "
                                         "before building." % release.tag)
            # Never replace a live deploy's original dist/ files, even in this clone.
            with tempfile.TemporaryDirectory(prefix="vaxrank-dry-run-") as directory:
                output = Path(directory)
                build_artifacts(output)
                check_artifacts([output / name for name in sorted(release.filenames)])
        print("Dry run: ready to release %s" % release.tag)
        return

    owner = release.acquire(head, fetch_artifacts, wait_seconds, poll_seconds)
    if owner is None:
        return
    try:
        if release.preflight() != head:
            raise ReleaseUploadError("HEAD changed while waiting; rerun ./deploy.sh.")
        # Recheck after ownership: an earlier owner may have just finished.
        if release.completed(head, fetch_artifacts):
            return
        published = fetch_artifacts()
        local = release.local_tag()
        if published and not local:
            raise ReleaseUploadError(
                "PyPI already contains %s without a verifiable release tag. Verify the "
                "published artifacts and tag; do not rebuild this version." % release.tag)
        output = release.repo / "dist"
        paths = [output / name for name in sorted(release.filenames)]
        release.check_owner(owner)
        if local:
            hashes = release.manifest(local, head)
            verify_local_artifacts(paths, hashes)
            print("Resuming release from %s with its original artifacts." % release.tag)
        else:
            if output.exists():
                shutil.rmtree(output)
            build_artifacts(output)
        check_artifacts(paths)
        if release.preflight() != head:
            raise ReleaseUploadError("HEAD changed during build; rerun ./deploy.sh.")
        if not local:
            hashes = release.record_artifacts(head, paths)
        verify_local_artifacts(paths, hashes)
        release.check_owner(owner)
        # The upload helper retains its exact-artifact guard and partial-upload recovery.
        publish_release(paths, expected_filenames=release.filenames,
                        fetch_release_artifacts=fetch_artifacts, upload_file=upload_file,
                        verify_timeout_seconds=release.verify_timeout_seconds)
        verify_local_artifacts(paths, hashes)
        release.check_owner(owner)
        pushed = release.git("-c", "push.followTags=false", "push", "origin",
                             release.tag_ref, check=False)
        if release.remote_ref(release.tag_ref) != release.local_tag():
            raise ReleaseUploadError("Release tag push failed: " + pushed.stderr.strip())
        release.completed(head, fetch_artifacts)
    finally:
        release.release_owner(owner)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", default="vaxrank")
    parser.add_argument("--version", required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    release = Release(Path.cwd(), args.project, args.version,
                      float(os.environ.get("VAXRANK_VERIFY_TIMEOUT_SECONDS", DEFAULT_VERIFY_TIMEOUT_SECONDS)))
    upload_timeout = float(os.environ.get("VAXRANK_UPLOAD_TIMEOUT_SECONDS", DEFAULT_UPLOAD_TIMEOUT_SECONDS))

    def build_artifacts(output):
        command = [sys.executable, "-m", "build", "--outdir", str(output)]
        if args.dry_run:
            # Build backends also write build/ and egg-info beside the sources.
            # Isolate those writes from a real deploy in the same checkout.
            with tempfile.TemporaryDirectory(prefix="vaxrank-dry-run-source-") as directory:
                root = Path(directory)
                source = root / "source"
                source.mkdir()
                archive = root / "source.tar"
                release.git("archive", "--format=tar", "--output=" + str(archive), "HEAD")
                subprocess.run(["tar", "-xf", str(archive), "-C", str(source)], check=True)
                subprocess.run(command, cwd=source, check=True)
        else:
            subprocess.run(command, cwd=release.repo, check=True)

    def check_artifacts(paths):
        subprocess.run([sys.executable, "-m", "twine", "check", *map(str, paths)], check=True)

    try:
        deploy(release, dry_run=args.dry_run,
               fetch_artifacts=lambda: pypi_release_artifacts(args.project, args.version),
               build_artifacts=build_artifacts, check_artifacts=check_artifacts,
               upload_file=lambda path: upload_distribution(
                   path, python_executable=sys.executable, timeout_seconds=upload_timeout))
    except (ReleaseUploadError, subprocess.CalledProcessError) as error:
        print("Release failed: %s" % error, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
