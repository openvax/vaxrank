set -e

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"

uv pip install -e "$ROOT_DIR"
# Prefer local sources when running dev/test commands from this shell.
export PYTHONPATH="${ROOT_DIR}${PYTHONPATH:+:$PYTHONPATH}"
