#!/usr/bin/env bash
set -eo pipefail

# Apple Silicon Homebrew + macOS SIP combine to strip dyld's search
# path when bash runs this script, so WeasyPrint can't find Pango.
# Re-export here so direct `./test.sh` runs work the same as
# deploy-driven ones.
if [[ "$(uname)" == "Darwin" && -d /opt/homebrew/lib \
      && ":${DYLD_FALLBACK_LIBRARY_PATH:-}:" != *":/opt/homebrew/lib:"* ]]; then
  export DYLD_FALLBACK_LIBRARY_PATH="/opt/homebrew/lib${DYLD_FALLBACK_LIBRARY_PATH:+:$DYLD_FALLBACK_LIBRARY_PATH}"
fi

if [[ -z "${PYTHON:-}" ]]; then
  if [[ -n "${VIRTUAL_ENV:-}" && -x "${VIRTUAL_ENV}/bin/python" ]]; then
    PYTHON="${VIRTUAL_ENV}/bin/python"
  elif [[ -x ".venv/bin/python" ]]; then
    PYTHON=".venv/bin/python"
  else
    PYTHON="python3"
  fi
fi

# Isolate this run from pytest cleanup in sibling repositories. A caller's
# explicit root is preserved, including when nested tests invoke this script.
own_temproot=0
if [[ -z "${PYTEST_DEBUG_TEMPROOT:-}" ]]; then
  PYTEST_DEBUG_TEMPROOT="$(mktemp -d "${TMPDIR:-/tmp}/vaxrank-pytest.XXXXXX")"
  export PYTEST_DEBUG_TEMPROOT
  own_temproot=1
fi

status=0
"${PYTHON}" -m pytest tests "$@" || status=$?
if (( own_temproot )); then
  if (( status == 0 )); then
    rm -rf -- "${PYTEST_DEBUG_TEMPROOT}"
  else
    echo "Kept pytest temp root: ${PYTEST_DEBUG_TEMPROOT}" >&2
  fi
fi
exit "$status"
