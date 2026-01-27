#!/bin/bash
set -euo pipefail

PYTHON="${PYTHON:-python3}"
if [[ -x ".venv/bin/python" ]]; then
  PYTHON=".venv/bin/python"
fi

./lint.sh
./test.sh

"$PYTHON" -m pip install --upgrade build twine
rm -rf dist
"$PYTHON" -m build
git --version
"$PYTHON" -m twine upload dist/*
git tag "$("$PYTHON" vaxrank/version.py)"
git push --tags
