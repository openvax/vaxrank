# Experimental CTA expression in the shared manifest

Add `format: cta_expression` to the existing manifest, alongside new antigen
sequences and saved evidence. Each input supplies its expression contract,
minimum expression and optional target exclusions/Hitlist bundle. Patient,
sample, host reference and genotype use the common scope declarations.

Separate CTA admission from prediction. The CLI alias and manifest use the
same adapter, existing OncoRef admission policy and shared Topiary predictor.
Validate every input before acquiring models. Preserve gene/transcript units,
annotation identity, held-out and unresolved decisions, and captured evidence.
An entirely held-out CTA input can remain as evidence beside usable sources.

Retain the common EpitopeDataset and peptide/mRNA outputs with offline replay.
A single CTA source keeps existing admission/product reports. Multiple inputs
retain per-source admission records; an old product audit is not an audit of
a new mixed design. Write separate admission JSON files under output-dir when
there are multiple CTA sources. The explicit single-file admission flag keeps
its existing one-source meaning.

Validate alias/manifest equivalence, mixed sources and samples, conflicts before
models, explicit contracts and relative bundle paths, all-held-out sources,
source-specific evidence and exact replay. All new interfaces remain experimental.
