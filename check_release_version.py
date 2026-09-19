"""Reject missing, stale, or already-tagged PR versions without importing code."""

import argparse
import ast
import subprocess
import sys

from packaging.version import InvalidVersion, Version


def _git(*args):
    result = subprocess.run(
        ["git", *args], capture_output=True, text=True, check=False)
    if result.returncode:
        raise ValueError("Git inspection failed: %s" % result.stderr.strip())
    return result.stdout.strip()


def _read_version(revision):
    commit = _git("rev-parse", "--verify", "--end-of-options", revision + "^{commit}")
    source = _git("show", commit + ":vaxrank/version.py")
    try:
        module = ast.parse(source)
    except SyntaxError as error:
        raise ValueError("Invalid version file at %s: %s" % (revision, error)) from error
    values = []
    for statement in module.body:
        if isinstance(statement, ast.Assign):
            targets = statement.targets
        elif isinstance(statement, ast.AnnAssign):
            targets = [statement.target]
        else:
            continue
        if any(isinstance(target, ast.Name) and target.id == "__version__"
               for target in targets):
            values.append(statement.value)
    if (len(values) != 1 or not isinstance(values[0], ast.Constant)
            or not isinstance(values[0].value, str)):
        raise ValueError(
            "%s:vaxrank/version.py must assign __version__ exactly once "
            "to a string literal" % revision)
    try:
        return Version(values[0].value)
    except InvalidVersion as error:
        raise ValueError("Invalid release version at %s: %s" % (revision, error)) from error


def check_release_version(base, head="HEAD"):
    """Compare committed versions and local tags; callers must fetch them first."""
    base_version = _read_version(base)
    head_version = _read_version(head)
    if head_version <= base_version:
        raise ValueError(
            "PR version %s must be greater than target branch version %s. "
            "Update vaxrank/version.py before merging." % (head_version, base_version))
    for tag in _git("tag", "--list", "v*").splitlines():
        try:
            tagged_version = Version(tag[1:])
        except InvalidVersion:
            continue
        if tagged_version == head_version:
            raise ValueError(
                "PR version %s is already tagged as %s. Choose a new version "
                "in vaxrank/version.py." % (head_version, tag))
    return "Release version %s is newer than %s and has no existing tag." % (
        head_version, base_version)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, help="Fetched target branch revision")
    parser.add_argument("--head", default="HEAD", help="Exact PR head revision")
    args = parser.parse_args(argv)
    try:
        message = check_release_version(args.base, args.head)
    except ValueError as error:
        print("Release version check failed: %s" % error, file=sys.stderr)
        return 1
    print(message)
    return 0


if __name__ == "__main__":
    sys.exit(main())
