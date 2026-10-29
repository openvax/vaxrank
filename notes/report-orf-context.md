# Coding context alongside antigen event kind (#717)

Keep the event kind (`mutation`) independent of its coding reference. Carry
the producer's reference kind, ORF type and ID, catalog, reference identity,
and transcript in the existing Isovar candidate metadata. A recorded
`extra_orf` is a noncanonical ORF; an `annotated` reference is an annotated
CDS. Missing historical context stays unknown, and an unfamiliar kind is
shown as recorded rather than guessed from an identifier or input family.

Use one small context view for candidate cards/ranking, candidate CSV/XLSX,
source outcome rows, and actual product mappings. Preserve all distinct
contexts in products shared by multiple candidates. Do not change event
identity, product names/sequences, ranking, admission, or native replay
assembly fingerprints just to improve a report label.

Historical native inputs may recover context from an unambiguous candidate
model ID and its linked source catalog. If a merged record does not uniquely
pair candidate and source evidence, retain the unknown instead of selecting
the first reference. New product source mappings retain their exact context;
old saved products can display context by matching the existing complete
source mapping to ranked candidates without rewriting their manifests.

The producer contract is Isovar's variant_rna_orfs model reference, built
from ORFdata's reference kind/ORF type metadata. This is a provenance display
change, not a claim of sample-specific translation or tumor specificity.

Validate annotated CDS, noncanonical dORF, absent historical metadata,
ambiguous merged provenance, and identical products with mixed contexts.
Check native round-trip, HTML/PDF/text and CSV/XLSX semantics, product parity,
lint/full tests, review, CI and release. Evidence grouping remains #719.
