# Isovar support compatibility

Vaxrank requires Isovar 1.42.1 through 1.44.x. That release corrects
flanking insertions/deletions in the anchored RNA candidate matcher. The fusion
adapter still consumes and validates `isovar.fusion_rna.v3`; compatibility is
tested with both the minimum 1.42.1 release and 1.44.2. Isovar 1.45 requires
Osteosarc 0.15, while Vaxrank and Topiary currently share the 0.14 window;
the coordinated update is tracked in [#571](https://github.com/openvax/vaxrank/issues/571).

## Counts and candidate scores

Vaxrank's direct mutation fragments retain the following distinct measurements:

| Fields | Meaning |
| --- | --- |
| `n_alt_reads`, `n_alt_fragments` | All alternate-supporting reads/fragments at the focal locus |
| `n_alt_reads_supporting_protein_sequence`, `n_alt_fragments_supporting_protein_sequence` | Reads/fragments used to assemble the translated cDNA |
| `n_rna_supporting_protein_sequence` | Assembly fragment count when available, otherwise assembly read count |

Partial assembly reads need not span a full vaccine peptide. Fragment identity
is scoped by SAM read group, rather than a QNAME alone. These counts are not
independent molecule counts, TPM or cancer-cell fractions.

Isovar's separately enabled `partial_read_support` export reconsiders original
reads against a specified RNA/protein hypothesis catalog. One scored fragment
contributes one unit split across compatible hypotheses. The resulting
`fractional_fragments` is a catalog-relative evidence score, not a raw count.
Shorter local hypotheses remain alternatives and can change that allocation.

The corrected `anchored_edit_compatibility.v2` matcher charges internal flanking
indels once and treats unobserved distal overhangs as free. It conservatively
reports coverage across tied alignments. Its tolerance is an explicit
compatibility rule, not a calibrated error probability. See the producer's
[partial-read support contract](https://github.com/openvax/isovar/blob/v1.42.1/docs/partial-read-support.md)
and [matching specification](https://github.com/openvax/isovar/blob/v1.42.1/docs/flanking-indel-matching-spec.md).

## Ranking and provenance

This adoption changes no frozen `openvax-v1` formula or selection golden.
Vaxrank continues to consume Isovar's assembly-support counts and ordering;
it does not enable fractional allocation or substitute those scores into RNA
counts. Fusion admission still requires its own tumor-specificity evidence.

New direct fragments store the actual runtime Isovar version as
`sequence_source_version`, carried into mixed-input evidence and native replay.
An older fragment without this field retains an unknown version. Historical
fixture producer versions remain immutable and are never rewritten to the
current installation.
