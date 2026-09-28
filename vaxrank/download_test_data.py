"""Export of the Sid retrieval test cases.

The default exports the 49 retrieval cases from the Sid test data, which is
built from the shared openvax-v1 reads (downloaded once into the osteosarc
cache; ``--offline`` never downloads them). Explicit custom manifests retain the
generic verified-download API for existing callers.

Run it as ``vaxrank-test-data`` or ``python -m vaxrank.download_test_data``.
Both print one JSON object describing what was exported or inspected, so a
caller can read a per-file status instead of parsing prose. From Python,
``inspect_dataset`` reports those statuses without raising, ``verify_dataset``
raises on the first problem, and ``download_test_data`` exports.
"""

import argparse
import ctypes
import errno
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
from urllib.parse import urlsplit

from datacache import Cache, FileValidationError, get_data_dir, inspect_files
from osteosarc import OsteosarcError


DEFAULT_MANIFEST = None
CACHE_NAMESPACE = "openvax"
CACHE_ENVIRONMENT = "OPENVAX_DATA_CACHE"
# Osteosarc resolves its cache as OSTEOSARC_CACHE, then OPENVAX_DATA_CACHE, then
# the platform directory. The Sid half of this command goes through osteosarc, so
# honouring the same order here keeps one invocation from reading bundled reads
# out of one root and custom assets out of another.
OSTEOSARC_ENVIRONMENT = "OSTEOSARC_CACHE"


def cache_root_for(cache_root=None):
    """The shared OpenVax cache root, resolved as osteosarc resolves it."""
    return Path(
        cache_root
        or os.environ.get(OSTEOSARC_ENVIRONMENT)
        or os.environ.get(CACHE_ENVIRONMENT)
        or get_data_dir(CACHE_NAMESPACE))


def _sid_test_data(offline):
    """The Sid test files; offline, only if openvax-v1 is already cached."""
    import osteosarc
    from .sid_test_data import READS, sid_test_data
    if offline:
        osteosarc.fetch_bundle(READS, offline=True)
    return sid_test_data()


def load_manifest(path=DEFAULT_MANIFEST, *, offline=False):
    """Validate an explicit manifest, or the Sid test data's retrieval manifest."""
    if path is None:
        path = _sid_test_data(offline) / "osteosarc/shared-v1/manifest.json"
    manifest = json.loads(Path(path).read_text())
    if manifest.get("schema_version") != 1 or not manifest.get("dataset") or not manifest.get("data_version"):
        raise ValueError("Unsupported test-data manifest")
    assets = manifest.get("assets")
    if not isinstance(assets, list) or not assets:
        raise ValueError("Manifest must contain assets")
    names = set()
    for asset in assets:
        name = asset["filename"]
        if not isinstance(name, str) or not name or name in (".", "..", "manifest.json") or "/" in name or "\\" in name or ":" in name:
            raise ValueError("Asset filename must be a plain, safe basename")
        if name in names:
            raise ValueError("Duplicate asset filename: " + name)
        names.add(name)
        if not re.fullmatch("[0-9a-f]{64}", asset["sha256"]):
            raise ValueError("Asset needs a lowercase SHA-256 digest")
        if type(asset["size_bytes"]) is not int or asset["size_bytes"] < 0:
            raise ValueError("Asset needs a nonnegative byte size")
        if "bundle_path" in asset:
            if asset["bundle_path"] != "osteosarc/shared-v1/" + name:
                raise ValueError("Invalid bundled asset path")
        elif urlsplit(asset["url"]).scheme not in ("https", "http", "file"):
            raise ValueError("Unsupported asset URL scheme")
    cases = manifest.get("cases", [])
    if len({c["case_id"] for c in cases}) != len(cases):
        raise ValueError("Duplicate test case identity")
    if any(c["bam"] not in names or c["bam"] + ".bai" not in names for c in cases):
        raise ValueError("Test case BAM/index absent from manifest")
    return manifest


def cache_filename(asset):
    """Content identity shared across consumers, URLs and dataset revisions."""
    # Keep the original suffix so datacache never infers decompression from a
    # suffix-less destination for a compressed upstream object.
    return "%s%s" % (asset["sha256"], "".join(Path(asset["filename"]).suffixes))


def _inspect_dataset(output, manifest):
    """``(report, datacache inspection or None)``; hashes each asset once."""

    def report(status, detail, files=None):
        return dict(path=str(output), status=status, detail=detail,
                    verified=status == "available", files=files or {})

    if output.is_symlink() or not output.is_dir():
        return report("not_a_directory", "Dataset must be a real directory"), None
    expected_names = {a["filename"] for a in manifest["assets"]} | {"manifest.json"}
    present = {p.name for p in output.iterdir()}
    if present != expected_names:
        missing = sorted(expected_names - present)
        unexpected = sorted(present - expected_names)
        return report("unexpected_contents", "Dataset has missing or unexpected files"
                      + ("; missing: " + ", ".join(missing) if missing else "")
                      + ("; unexpected: " + ", ".join(unexpected) if unexpected else "")), None
    if (output / "manifest.json").is_symlink():
        return report("symlink", "Dataset manifest must not be a symlink"), None
    if json.loads((output / "manifest.json").read_text()) != manifest:
        return report("manifest_mismatch", "Dataset manifest does not match requested dataset"), None
    symlinked = sorted(a["filename"] for a in manifest["assets"]
                       if (output / a["filename"]).is_symlink())
    if symlinked:
        return report("symlink", "Dataset assets must not be symlinks: "
                      + ", ".join(symlinked)), None

    inspection = inspect_files(str(output), {
        a["filename"]: dict(expected_sha256=a["sha256"], expected_size=a["size_bytes"])
        for a in manifest["assets"]})
    files = {name: dict(status=inspected.status, verified=inspected.verified)
             for name, inspected in inspection.files.items()}
    unverified = sorted(name for name, inspected in inspection.files.items()
                        if not inspected.verified)
    if unverified:
        return report(inspection.status, "Assets failed verification: "
                      + ", ".join("%s (%s)" % (name, files[name]["status"])
                                  for name in unverified), files), inspection
    return report("available", "Every asset matched its recorded digest and size",
                  files), inspection


def inspect_dataset(output, manifest=None):
    """Report a materialized dataset's condition offline, without raising.

    Returns ``{"path", "status", "verified", "files", "detail"}``. ``status`` is
    ``"available"`` only when every asset's digest and size were checked and
    matched; otherwise it names what is wrong, and ``detail`` says it in prose.
    ``files`` maps each asset to datacache's per-file ``status``/``verified``.

    Structural problems are checked here rather than delegated: datacache's
    ``inspect_files`` follows symlinks and only looks at the inventory it is
    given, so on its own it would call a symlinked asset available and would
    never notice an extra file.
    """
    manifest = load_manifest() if manifest is None else manifest
    return _inspect_dataset(Path(output), manifest)[0]


def verify_dataset(output, manifest=None):
    """Validate a materialized dataset offline, without repair or cache writes.

    Raises on the first problem, so callers that only care whether the dataset
    is usable keep one failure path. ``inspect_dataset`` reports the same
    findings without raising.
    """
    manifest = load_manifest() if manifest is None else manifest
    report, inspection = _inspect_dataset(Path(output), manifest)
    if report["status"] == "available":
        return Path(output)
    if inspection is not None:
        # Re-raise datacache's own error so callers keep distinguishing a
        # corrupt asset (FileValidationError) from an absent one.
        for inspected in inspection.files.values():
            if inspected.error is not None:
                raise inspected.error
    raise ValueError("%s: %s" % (report["detail"], output))


def _publish_no_replace(staged, output):
    """Atomically publish a directory, refusing even an empty destination.

    Python's rename replaces empty directories on POSIX. Use the native
    exclusive operation on Linux/macOS, and fail closed if it is unavailable.
    """
    libc = ctypes.CDLL(None, use_errno=True)
    if sys.platform == "linux":
        rename = getattr(libc, "renameat2", None)
        argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
        # AT_FDCWD = -100; RENAME_NOREPLACE = 1 (Linux renameat2(2)).
        args = (-100, os.fsencode(staged), -100, os.fsencode(output), 1)
    elif sys.platform == "darwin":
        rename = getattr(libc, "renamex_np", None)
        argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint]
        # RENAME_EXCL = 0x4 (Darwin sys/stdio.h).
        args = (os.fsencode(staged), os.fsencode(output), 0x4)
    else:
        rename = None
    if rename is None:
        raise OSError(errno.ENOTSUP, "Atomic no-replace directory publication is unavailable", str(output))
    rename.argtypes = argtypes
    rename.restype = ctypes.c_int
    if rename(*args) != 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error), str(output))


def download_test_data(output=None, *, manifest_path=DEFAULT_MANIFEST,
                       cache_root=None, offline=False, repair_cache=False,
                       show_progress=False, timeout=60, max_retries=2):
    """Export the Sid retrieval cases, or fetch an explicitly supplied custom manifest.

    Existing output directories are only validated, never overwritten. Failed
    downloads leave verified cache objects reusable but do not publish output.
    ``repair_cache`` explicitly permits replacing invalid cached objects.
    ``show_progress`` draws datacache's download progress; ``timeout`` and
    ``max_retries`` bound each request.
    """
    if offline and repair_cache:
        raise ValueError("Offline mode cannot repair the download cache")
    manifest = load_manifest(manifest_path, offline=offline)
    if output is not None:
        output = Path(output).absolute()
        if output.exists() or output.is_symlink():
            return verify_dataset(output, manifest)
    if all("bundle_path" in asset for asset in manifest["assets"]):
        bundled = _sid_test_data(offline) / "osteosarc/shared-v1"
        verify_dataset(bundled, manifest)
        paths = {asset["filename"]: bundled / asset["filename"] for asset in manifest["assets"]}
    else:
        root = cache_root_for(cache_root)
        cache = Cache(CACHE_NAMESPACE, cache_root=root / "objects" / "sha256")
        paths = {}
        for asset in manifest["assets"]:
            name = cache_filename(asset)
            expected = dict(expected_sha256=asset["sha256"], expected_size=asset["size_bytes"])
            # Inspection never downloads, creates directories or repairs files,
            # so it is the read-only probe both branches below need.
            cached = cache.inspect(filename=name, **expected)
            if offline:
                if cached.error is not None:
                    raise cached.error
                path = cached.path
            else:
                path = cache.fetch(asset["url"], filename=name, timeout=timeout,
                                   force=repair_cache and cached.status == "corrupt",
                                   show_progress=show_progress,
                                   max_retries=max_retries, **expected)
            paths[asset["filename"]] = Path(path)
    if output is None:
        return paths
    output.parent.mkdir(parents=True, exist_ok=True)
    # A one-time output export, not a general mutable bundle registry. Siblings
    # permit atomic exclusive rename; no destination tree is deleted or replaced.
    with tempfile.TemporaryDirectory(prefix=".osteosarc-", dir=output.parent) as temporary:
        staged = Path(temporary) / "dataset"
        staged.mkdir()
        for name, path in paths.items():
            shutil.copyfile(path, staged / name)
        (staged / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        verify_dataset(staged, manifest)
        try:
            _publish_no_replace(staged, output)
        except FileExistsError:
            # A concurrent equivalent publisher is fine; a foreign/modified
            # output is rejected without deleting or repairing it.
            return verify_dataset(output, manifest)
    return verify_dataset(output, manifest)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST,
                        help="Custom test-data manifest; omit for the Sid retrieval cases")
    parser.add_argument("--output", type=Path, help="New directory to export, or an existing identical dataset to verify")
    parser.add_argument("--cache-root", type=Path,
                        help="Shared OpenVax cache root (else OSTEOSARC_CACHE, OPENVAX_DATA_CACHE, the platform cache)")
    parser.add_argument("--offline", action="store_true", help="Use verified local objects only; never access the network")
    parser.add_argument("--repair-cache", action="store_true", help="Explicitly refetch invalid cached objects")
    parser.add_argument("--verify-only", action="store_true", help="Validate --output without fetching or touching the cache")
    parser.add_argument("--progress", action="store_true", help="Show download progress for custom manifests")
    parser.add_argument("--timeout", type=float, default=60, help="Seconds allowed per request (default: 60)")
    parser.add_argument("--max-retries", type=int, default=2, help="Retries per failed request (default: 2)")
    args = parser.parse_args(argv)
    if args.verify_only and (args.output is None or args.repair_cache):
        parser.error("--verify-only requires --output and cannot repair the cache")

    report = dict(dataset=None, data_version=None, assets=None,
                  action="verify" if args.verify_only else "export",
                  status="unavailable", output=None, cached_paths=None,
                  files={}, error=None)

    def emit(code):
        """One JSON object on stdout whatever happened; exit non-zero on failure.

        Success returns normally so importing callers of ``main`` do not have to
        catch ``SystemExit`` for the good case.
        """
        print(json.dumps(report, indent=2, sort_keys=True))
        if code:
            parser.exit(code)

    try:
        manifest = load_manifest(args.manifest, offline=args.offline or args.verify_only)
        report.update(dataset=manifest["dataset"], data_version=manifest["data_version"],
                      assets=len(manifest["assets"]))
        if args.verify_only:
            inspected = inspect_dataset(args.output, manifest)
            report.update(status=inspected["status"], output=inspected["path"],
                          files=inspected["files"],
                          error=None if inspected["verified"] else inspected["detail"])
            return emit(0 if inspected["verified"] else 1)
        result = download_test_data(
            args.output, manifest_path=args.manifest, cache_root=args.cache_root,
            offline=args.offline, repair_cache=args.repair_cache,
            show_progress=args.progress, timeout=args.timeout,
            max_retries=args.max_retries)
    except (OSError, ValueError, FileValidationError, OsteosarcError) as error:
        report["error"] = "Test data unavailable: %s" % error
        return emit(1)
    if args.output is None:
        report.update(status="cached", cached_paths={k: str(v) for k, v in result.items()})
    else:
        inspected = inspect_dataset(result, manifest)
        report.update(status=inspected["status"], output=inspected["path"],
                      files=inspected["files"])
    return emit(0)


if __name__ == "__main__":
    main()
