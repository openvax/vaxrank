set -e

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
PYTHON="${PYTHON:-python3}"

if command -v uv >/dev/null 2>&1; then
  uv pip install -e "$ROOT_DIR"
else
  "$PYTHON" -m pip install -e "$ROOT_DIR"
fi
# Prefer local sources when running dev/test commands from this shell.
export PYTHONPATH="${ROOT_DIR}${PYTHONPATH:+:$PYTHONPATH}"
