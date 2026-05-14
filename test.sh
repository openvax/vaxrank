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

exec pytest tests "$@"
