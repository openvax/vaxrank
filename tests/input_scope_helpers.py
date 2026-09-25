"""Explicit scope assertions for synthetic external-input fixtures."""

import json
from pathlib import Path


FIXTURE_SCOPE = {
    "patient_id": "fixture-patient",
    "reference_assembly": "GRCh38",
    "mhc_alleles": ["HLA-A*02:01", "HLA-B*07:02"],
}


def write_input_manifest(path, inputs, **scope):
    """Write an actual manifest for CLI and loader integration tests."""
    document = {
        "schema": "vaxrank.input_manifest.v1",
        **FIXTURE_SCOPE, **scope,
        "inputs": [{"format": fmt, "path": str(Path(filename).resolve())}
                   for fmt, filename in inputs],
    }
    Path(path).write_text(json.dumps(document))
    return str(path)
