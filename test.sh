#!/usr/bin/env bash
set -euo pipefail

# See deploy.sh — Apple Silicon Homebrew + macOS SIP combine to strip dyld's
# search path when bash runs this script, so WeasyPrint can't find pango.
# Re-export here so direct `./test.sh` runs work the same as deploy-driven ones.
if [[ "$(uname)" == "Darwin" && -d /opt/homebrew/lib \
      && ":${DYLD_FALLBACK_LIBRARY_PATH:-}:" != *":/opt/homebrew/lib:"* ]]; then
  export DYLD_FALLBACK_LIBRARY_PATH="/opt/homebrew/lib${DYLD_FALLBACK_LIBRARY_PATH:+:$DYLD_FALLBACK_LIBRARY_PATH}"
fi

# Use pytest-xdist with 2 workers when available. Each worker loads
# its own pyensembl / mhctools / topiary state (~3 GB resident), so
# naive `-n auto` on an 8-core box can balloon to ~24 GB. Two
# workers is a sweet spot: meaningful parallelism without the memory
# blowup. xdist is not a required dependency, so fall back to serial
# pytest when it isn't installed (e.g., CI).
XDIST_FLAG=""
if python -c "import xdist" 2>/dev/null; then
  XDIST_FLAG="-n 2"
fi
pytest ${XDIST_FLAG} --cov=vaxrank/ --cov-report=term-missing tests
