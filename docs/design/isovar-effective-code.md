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
