#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: ./deploy.sh [--dry-run] [version]" >&2
}

DRY_RUN=0
VERSION_INPUT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      if [[ -z "${VERSION_INPUT}" ]]; then
        VERSION_INPUT="$1"
        shift
      else
        echo "Unexpected argument: $1" >&2
        usage
        exit 1
      fi
      ;;
  esac
done

PYTHON="${PYTHON:-python3}"
if [[ -x ".venv/bin/python" ]]; then
  PYTHON=".venv/bin/python"
fi

CURRENT_VERSION="$(
  "$PYTHON" - <<'PY'
import re
from pathlib import Path

text = Path("vaxrank/version.py").read_text()
match = re.search(r'__version__ = "([^"]+)"', text)
if not match:
    raise SystemExit("Failed to read vaxrank/version.py; unexpected format.")
print(match.group(1))
PY
)"

VERSION="${CURRENT_VERSION}"
if [[ -n "${VERSION_INPUT}" ]]; then
  VERSION="${VERSION_INPUT#v}"
fi

TAG="v${VERSION}"

# Deploys should only run from the release branch.
CURRENT_BRANCH="${VAXRANK_DEPLOY_BRANCH:-$(git rev-parse --abbrev-ref HEAD)}"
if [[ "${CURRENT_BRANCH}" != "main" && "${CURRENT_BRANCH}" != "master" ]]; then
  echo "Deploys are only allowed from main or master (current: ${CURRENT_BRANCH})." >&2
  exit 1
fi

# Ensure a clean tree before deploy.
if [[ -n "$(git status --porcelain)" ]]; then
  echo "Working tree not clean; commit or stash changes before deploy." >&2
  exit 1
fi

# Ensure the tag doesn't already exist.
if git rev-parse -q --verify "refs/tags/${TAG}" >/dev/null; then
  echo "Tag ${TAG} already exists; aborting." >&2
  exit 1
fi

./lint.sh
./test.sh

# Bump version, commit, and push if version was provided.
if [[ -n "${VERSION_INPUT}" && "${VERSION}" != "${CURRENT_VERSION}" ]]; then
  VERSION="${VERSION}" "$PYTHON" - <<'PY'
import os
import re
from pathlib import Path

version = os.environ["VERSION"]
path = Path("vaxrank/version.py")
text = path.read_text()
new_text, n = re.subn(
    r'__version__ = "([^"]+)"',
    f'__version__ = "{version}"',
    text,
    count=1,
)
if n != 1:
    raise SystemExit("Failed to update vaxrank/version.py; unexpected format.")
path.write_text(new_text)
PY
  git add vaxrank/version.py
  git commit -m "Bump version to ${VERSION}"
  git push
fi

"$PYTHON" -m pip install --upgrade build twine
rm -rf dist
"$PYTHON" -m build
git --version

if [[ "${DRY_RUN}" -eq 1 ]]; then
  if [[ -n "${VERSION_INPUT}" && "${VERSION}" != "${CURRENT_VERSION}" ]]; then
    echo "Dry run: would bump version to ${VERSION}"
  fi
  echo "Dry run: ready to release ${TAG}"
  exit 0
fi

"$PYTHON" -m twine check dist/*
"$PYTHON" -m twine upload dist/*
git tag "${TAG}"
git push --tags
