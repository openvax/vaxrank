# Native Isovar effective genetic codes

Consume Isovar's `variant_rna_orfs.v2` candidate code through its public
`variant_orf_translation_table` reader. Use the returned code for protein
validation, frame grouping and candidate provenance. Keep source reference
metadata intact. Historical v1 imports retain their original nullable-code
decoding and frame identity; do not infer new codes for saved results.

Consolidate `_translation` onto `translate_cdna` for all selected tables.
Require Isovar 1.89.0, which supplies the corrected reconstruction and public
reader. Validate explicit and context-derived mitochondrial candidates,
nuclear controls, invalid/missing codes, source conflicts, and saved replay.
Coordinate this release with openvax/isovar#700; do not claim NuMT exclusion
from a transcript code or mapped RNA alone.

CI reference contract: use the already provisioned human Ensembl release 93
and stable ND5, ND6 and TP53 transcript IDs. The initial release-87 tests
passed against a local cache but failed on clean CI because that undeclared
reference was absent. Keep all three native reconstruction/replay controls;
verify the corrected fixtures locally and across the complete CI matrix.
