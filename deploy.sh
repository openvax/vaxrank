#!/usr/bin/env bash
set -euo pipefail

# macOS SIP strips DYLD_* env vars when exec'ing /bin/bash (how
# `#!/usr/bin/env bash` resolves). WeasyPrint — used by tests/ and by
# --pdf-backend=weasyprint — loads Pango/Cairo via ctypes.find_library,
# which on Apple Silicon can't locate Homebrew's /opt/homebrew/lib unless
# dyld's fallback search path is told about it. Re-export the var inside
# the script so ./test.sh → pytest → weasyprint see it.
if [[ "$(uname)" == "Darwin" && -d /opt/homebrew/lib \
      && ":${DYLD_FALLBACK_LIBRARY_PATH:-}:" != *":/opt/homebrew/lib:"* ]]; then
  export DYLD_FALLBACK_LIBRARY_PATH="/opt/homebrew/lib${DYLD_FALLBACK_LIBRARY_PATH:+:$DYLD_FALLBACK_LIBRARY_PATH}"
fi

usage() {
  echo "Usage: ./deploy.sh [--dry-run] [version]" >&2
  echo "The optional version must match vaxrank/version.py; bump versions in a PR." >&2
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

if [[ -n "${PYTHON:-}" ]]; then
  :
elif [[ -n "${VIRTUAL_ENV:-}" && -x "${VIRTUAL_ENV}/bin/python" ]]; then
  PYTHON="${VIRTUAL_ENV}/bin/python"
elif [[ -x ".venv/bin/python" ]]; then
  PYTHON=".venv/bin/python"
else
  PYTHON="python3"
fi
export PYTHON

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
if [[ -n "${VERSION_INPUT}" && "${VERSION_INPUT#v}" != "${CURRENT_VERSION}" ]]; then
  echo "Requested version ${VERSION_INPUT} does not match vaxrank/version.py (${CURRENT_VERSION})." >&2
  echo "Bump the version in a PR, merge it, then run ./deploy.sh without a version argument." >&2
  exit 1
fi

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

./lint.sh
./test.sh

HEAD_COMMIT="$(git rev-parse HEAD)"
if ! UPSTREAM_COMMIT="$(git rev-parse '@{upstream}')"; then
  echo "Release branch has no upstream; pull and verify main before deploying." >&2
  exit 1
fi
if [[ "${HEAD_COMMIT}" != "${UPSTREAM_COMMIT}" ]]; then
  echo "Release branch is not synchronized with its upstream." >&2
  echo "HEAD=${HEAD_COMMIT} upstream=${UPSTREAM_COMMIT}" >&2
  exit 1
fi

# The coordinator owns build/upload/tag across independent clones.
"$PYTHON" -m pip install --upgrade build twine
DEPLOY_ARGS=(--project vaxrank --version "${VERSION}")
if [[ "${DRY_RUN}" -eq 1 ]]; then
  DEPLOY_ARGS+=(--dry-run)
fi
"$PYTHON" release_deploy.py "${DEPLOY_ARGS[@]}"
