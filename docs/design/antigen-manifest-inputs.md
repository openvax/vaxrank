# Experimental antigen inputs in a shared manifest

Allow `antigens` inputs beside saved `epitopes`, Topiary and imported prediction
tables in the existing manifest. A sequence file keeps its admission decisions;
only its missing patient/reference/genotype declarations may be filled by the
manifest. Relative paths resolve beside the manifest.

Load and validate every source before acquiring prediction models. Predict only
new sequence inputs through the shared Topiary adapter, then rank all evidence
through the existing EpitopeDataset pipeline. Preserve original prediction
values, source identities, held-out records and per-input admission evidence.
The CLI and Python preparation API use the same loader.

New sequence inputs require a host annotation and predictor. Saved-only replay
keeps its evidence-driven requirements. Fresh/additive rescoring of mixed raw
and saved inputs is outside this increment. Direct VCF/BAM and CTA expression
remain their existing adapters; their native outputs can join this manifest.

Verify mixed new/saved sources in peptide and mRNA designs, offline replay,
manifest/file genotype and reference conflicts before models, duplicate input
rejection, and single-file/manifest equivalence. Keep new help and docs concise
and explicitly experimental.
