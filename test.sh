#!/usr/bin/env bash
set -euo pipefail

# See deploy.sh — Apple Silicon Homebrew + macOS SIP combine to strip dyld's
# search path when bash runs this script, so WeasyPrint can't find pango.
# Re-export here so direct `./test.sh` runs work the same as deploy-driven ones.
if [[ "$(uname)" == "Darwin" && -d /opt/homebrew/lib ]]; then
  export DYLD_FALLBACK_LIBRARY_PATH="/opt/homebrew/lib:${DYLD_FALLBACK_LIBRARY_PATH:-}"
fi

pytest --cov=vaxrank/ --cov-report=term-missing tests
