# Pinned original RNA for downstream context and selection tests

Copied without modification from Isovar 1.8.0, commit
`7bee9bef690a0184a41a0d8397879915dd86fe37`:
https://github.com/openvax/isovar/tree/7bee9bef690a0184a41a0d8397879915dd86fe37/tests/data/osteosarc

`manifest.json` pins original-source and uncompressed fixture checksums,
selection rules, read names, and source provenance. `selection.json` records
the GRCh38 alleles and the two distinct RNA samples/timepoints. The 482
bulk STAR T0 and 242 nanopore T1 alignments retain their original bases,
qualities, CIGARs and tags. This deliberately small subset is not a VAF or
independent-molecule estimator; do not pool these samples as replicates.
The source is public CC0 data (see the unmodified `source_registry.yaml`):
https://registry.opendata.aws/sid-osteosarc/.

`protein_reference/` contains six original Ensembl 87 transcript models,
cDNAs and reference proteins, plus their source URLs and SHA-256 hashes.
This is a partial GRCh38 annotation, given its own reference/cache identity
in tests; it must not populate the full-genome GRCh38 annotation cache.
Transcript versions remain in the original FASTA headers. No reference
bases are added to RNA-derived protein candidates.

The context regressions use the real Vaxrank/Isovar factory and read
collector (including their assembly/mate-merging defaults), followed by
fragment conversion and mutation-overlapping window generation. They are
not claims of historical MHC ranking or clinical vaccine agreement.
The independent reference-protein edits checked here are documented on:

- https://osteosarc.com/variant/DYNC1H1-chr14-101980529/ (p.Val314Ile)
- https://osteosarc.com/variant/EXOC4-chr7-133274996/ (p.Ser34Ile)
- https://osteosarc.com/variant/H1_2-chr6-26055824/ (p.Ala197_Lys201del)
- https://osteosarc.com/variant/GTF3C5-chr9-133057893/ (p.Glu503_Glu506del)

MAP2 has no exact alternate deletion support in these two RNA fixtures;
that remains a no-RNA-peptide case, not evidence of absent mutant expression
or grounds for silently enabling DNA fallback. PIP5K1A's frameshift and the
historically selected vaccine sequences require separate selection-level
expectations under Vaxrank #414.

The small fixtures intentionally do not reproduce full-region read counts:
ONT DYNC1H1 has 16 alternate names, but the 20-aa candidate has 11 compatible
names and all full 25mer contexts have at most 9 (below 85% of 11). Bulk H1-2
has three alternate read objects / two names; after mate merging its longest
two-object-supported context is 24 aa. These produce no eligible default
25-aa vaccine peptide. Explicitly permitting a shorter peptide is distinct
from weakening either RNA support threshold (Vaxrank #419).

Tests check fixture integrity and build their small reference/BAM indices in
temporary directories. They perform no network requests or model downloads.
