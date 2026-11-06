# Local Hitlist evidence across candidate inputs

`--hitlist-local-evidence` opts into exact public MS evidence for candidate
inputs: antigen/Isovar sequences, direct variant predictions, Topiary/native
evidence, and LENS/pVACseq tables. It applies before DSL filtering and scoring,
including fresh and additive prediction modes. CTA expression inputs continue
to query automatically. An explicit CTA evidence bundle keeps precedence.

```bash
vaxrank --input-topiary predictions.tsv --hitlist-local-evidence \
  --output-epitopes annotated.tsv
```

Lookup alone adds annotations. To require measured monoallelic support for the
candidate allele, explicitly select that facet in the existing DSL:

```bash
vaxrank --input-topiary predictions.tsv --hitlist-local-evidence \
  --config-text 'epitopes.filter_expr=hitlist_ms_monoallelic > 0' \
  --output-epitopes selected.tsv
```

Vaxrank reads Hitlist's current local index (`HITLIST_DATA_DIR` can select it).
It never downloads data or builds an index. Missing, stale, or unverifiable data
remain unknown; a completed query with no positive MS observations yields zero.
See [evidence semantics and restriction facets](design/expression-cta-inputs.md#hitlist-evidence)
for the distinction between any-MHC observations, measured restriction, predicted
restriction, and donor typing. Local MS evidence does not establish patient
expression, CTA specificity, or tissue safety. Exact matching does not transfer
support from a longer observed ligand to its nested peptides.

Native export retains per-source query scopes, index metadata, returned rows,
and unavailable states in `selection["source_evidence"]` under
`hitlist_local_queries`. Direct runs also retain these queries at the selection
root; automatic CTA captures retain their existing `hitlist_local_evidence` key.
Each query keeps its own peptide/species scope. Original occurrence IDs, antigen
source records, prediction values, and patient allele membership are preserved.
Explicit filtering limits scored alleles without removing the original evidence.

Reloading `--input-epitopes annotated.tsv` uses saved annotations and raw evidence
without requiring the local index. Supplying the flag again only queries rows
that lack a captured annotation; unavailable results are frozen too. To query a
new index for existing candidates, start from the original unannotated input.
Legacy `--input-json-file` stores selected vaccines rather than complete candidate
evidence and rejects this option.

The shared Python scoring API exposes the same opt-in:

```python
from vaxrank.external_rescoring import prepare_reports

reports, frame, candidates = prepare_reports(
    [("topiary", "predictions.tsv")], hitlist_local_evidence=True)
dataset = frame.attrs["epitope_dataset"]
dataset.save("annotated.tsv")
```

For low-level direct prediction, pass a `LocalHitlistLookup` instance as the
`annotate` callback to `predict_epitopes` or `run_vaxrank`. Its `queries` property
contains the captures for the caller to persist. The CLI does this automatically.
