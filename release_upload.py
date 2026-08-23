"""Bounded, idempotent PyPI uploads for vaxrank releases."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Callable, Iterable
from urllib.error import HTTPError
from urllib.parse import quote
from urllib.request import urlopen


PYPI_JSON_BASE_URL = "https://pypi.org/pypi"
DEFAULT_UPLOAD_TIMEOUT_SECONDS = 60.0
DEFAULT_VERIFY_TIMEOUT_SECONDS = 30.0
DEFAULT_VERIFY_POLL_SECONDS = 1.0


class ReleaseUploadError(RuntimeError):
    """Raised when the complete release artifact set cannot be verified."""


def expected_release_filenames(project: str, version: str) -> frozenset[str]:
    """Return the wheel and source-distribution names built by vaxrank."""
    wheel_project = project.replace("-", "_")
    return frozenset({
        f"{wheel_project}-{version}-py3-none-any.whl",
        f"{project}-{version}.tar.gz",
    })


def pypi_release_filenames(
    project: str,
    version: str,
    *,
    json_base_url: str = PYPI_JSON_BASE_URL,
    request_timeout_seconds: float = 10.0,
) -> frozenset[str]:
    """Return filenames currently published for one PyPI release."""
    url = "%s/%s/%s/json" % (
        json_base_url.rstrip("/"),
        quote(project, safe=""),
        quote(version, safe=""),
    )
    try:
        with urlopen(url, timeout=request_timeout_seconds) as response:
            payload = json.load(response)
    except HTTPError as error:
        if error.code == 404:
            return frozenset()
        raise
    return frozenset(item["filename"] for item in payload.get("urls", ()))


def upload_distribution(
    distribution_path: str | Path,
    *,
    python_executable: str = sys.executable,
    timeout_seconds: float = DEFAULT_UPLOAD_TIMEOUT_SECONDS,
    repository_url: str | None = None,
) -> None:
    """Upload one distribution with a finite timeout and quiet progress."""
    command = [
        python_executable,
        "-m",
        "twine",
        "upload",
        "--disable-progress-bar",
    ]
    if repository_url:
        command.extend(["--repository-url", repository_url])
    command.append(str(distribution_path))
    subprocess.run(command, check=True, timeout=timeout_seconds)


def wait_for_release_file(
    filename: str,
    fetch_release_filenames: Callable[[], Iterable[str]],
    *,
    timeout_seconds: float = DEFAULT_VERIFY_TIMEOUT_SECONDS,
    poll_seconds: float = DEFAULT_VERIFY_POLL_SECONDS,
) -> bool:
    """Poll release metadata until ``filename`` is visible or time expires."""
    deadline = time.monotonic() + timeout_seconds
    while True:
        if filename in set(fetch_release_filenames()):
            return True
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            return False
        time.sleep(min(poll_seconds, remaining))


def publish_release(
    distribution_paths: Iterable[str | Path],
    *,
    expected_filenames: Iterable[str],
    fetch_release_filenames: Callable[[], Iterable[str]],
    upload_file: Callable[[Path], None],
    verify_timeout_seconds: float = DEFAULT_VERIFY_TIMEOUT_SECONDS,
    verify_poll_seconds: float = DEFAULT_VERIFY_POLL_SECONDS,
) -> frozenset[str]:
    """Upload missing artifacts and reconcile every response with PyPI state.

    A timeout or nonzero Twine exit is ambiguous: the server may have accepted
    the immutable file before the client lost its response. This operation
    therefore verifies server state after every attempt and treats the file as
    published only when its exact filename appears in release metadata.
    """
    paths_by_name = {
        Path(distribution_path).name: Path(distribution_path)
        for distribution_path in distribution_paths
    }
    expected = frozenset(expected_filenames)
    if frozenset(paths_by_name) != expected:
        raise ReleaseUploadError(
            "Distribution files do not match the expected release set: "
            f"expected={sorted(expected)}, observed={sorted(paths_by_name)}"
        )

    published = frozenset(fetch_release_filenames())
    for filename in sorted(expected):
        if filename in published:
            print(f"Already published: {filename}")
            continue

        upload_error = None
        try:
            upload_file(paths_by_name[filename])
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
            upload_error = error

        if not wait_for_release_file(
            filename,
            fetch_release_filenames,
            timeout_seconds=verify_timeout_seconds,
            poll_seconds=verify_poll_seconds,
        ):
            message = f"PyPI did not publish {filename} after the upload attempt"
            if upload_error is not None:
                raise ReleaseUploadError(message) from upload_error
            raise ReleaseUploadError(message)
        if upload_error is not None:
            print(f"Reconciled ambiguous upload from PyPI state: {filename}")
        else:
            print(f"Published: {filename}")
        published = frozenset(fetch_release_filenames())

    published = frozenset(fetch_release_filenames())
    missing = expected - published
    if missing:
        raise ReleaseUploadError(
            f"PyPI release is missing expected files: {sorted(missing)}"
        )
    return published


def main(argv: list[str] | None = None) -> int:
    """Upload and verify one complete vaxrank release."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True)
    parser.add_argument("--version", required=True)
    parser.add_argument("--json-base-url", default=PYPI_JSON_BASE_URL)
    parser.add_argument("--repository-url")
    parser.add_argument(
        "--upload-timeout-seconds",
        type=float,
        default=float(os.environ.get(
            "VAXRANK_UPLOAD_TIMEOUT_SECONDS",
            DEFAULT_UPLOAD_TIMEOUT_SECONDS,
        )),
    )
    parser.add_argument(
        "--verify-timeout-seconds",
        type=float,
        default=float(os.environ.get(
            "VAXRANK_VERIFY_TIMEOUT_SECONDS",
            DEFAULT_VERIFY_TIMEOUT_SECONDS,
        )),
    )
    parser.add_argument("distributions", nargs="+")
    args = parser.parse_args(argv)

    expected = expected_release_filenames(args.project, args.version)

    def fetch_release_filenames():
        return pypi_release_filenames(
            args.project,
            args.version,
            json_base_url=args.json_base_url,
        )

    def upload_file(path):
        upload_distribution(
            path,
            python_executable=sys.executable,
            timeout_seconds=args.upload_timeout_seconds,
            repository_url=args.repository_url,
        )

    published = publish_release(
        args.distributions,
        expected_filenames=expected,
        fetch_release_filenames=fetch_release_filenames,
        upload_file=upload_file,
        verify_timeout_seconds=args.verify_timeout_seconds,
    )
    print(
        "Verified PyPI release files: %s"
        % ", ".join(sorted(expected & published))
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
