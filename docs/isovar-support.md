# Isovar support compatibility

Vaxrank 3.47 requires Isovar 1.61.0 or later, below 2.0, with ORFdata
0.39.1 or later and the compatible Osteosarc 0.15 series. Native
[Isovar inputs](isovar-inputs.md) cover variant/ncORFs, SV ORFs, fusion,
viral coding, targeted isoforms and ordinary protein-hypothesis exports.
That guide documents exact supported schemas, reconstruction/admission controls,
source-specific evidence boundaries and replay. The historical fusion API
continues to accept `isovar.fusion_rna.v3`.

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

The ordinary VCF/BAM workflow keeps its frozen `openvax-v1` formula and
selection policy. It continues to consume assembly-support counts and ordering;
it does not enable fractional allocation or substitute those scores into RNA
counts. New Isovar source inputs retain every candidate before shared ranking;
their admission requires an explicit tumor-specificity attestation. Reports
distinguish that assertion from measured sample evidence and catalog attribution.

New direct fragments store the actual runtime Isovar version as
`sequence_source_version`, carried into mixed-input evidence and native replay.
An older fragment without this field retains an unknown version. Historical
fixture producer versions remain immutable and are never rewritten to the
current installation.

## Evidence in reports

Full candidate details (`--report-candidate-details`) distinguish **construct
admission** from biological evidence and workflow review flags. An admitted
research assumption does not establish tumor specificity, and an unset admission
review flag does not record scientific sign-off. The report shows the supplied
evidence type, rationale, override, and effective partial/ambiguous admission
controls. Each retained attestation lists its subject, source, sample and units.

Isovar candidates show coding completeness and their recorded RNA evidence:
reconstruction reads/fragments, coding-interval depth, ORF-occurrence support or
direct junction support, as applicable. These are distinct subjects; none is
automatically exact support for a chosen peptide. A catalog's ORF attribution
and literature evidence appear separately from sample evidence. Missing fields
remain **Not recorded**, including ambiguity flags absent from older exports.
Multiple source contexts remain separate; a missing candidate-to-source pairing
is reported explicitly instead of assigning one origin's evidence to another.

HTML/PDF and text use the same evidence sections. Native prediction replay reads
the linked source catalog already stored in the dataset; it needs no new
sidecar, reconstruction or prediction. API callers can pass that catalog as
`TemplateDataCreator(source_evidence=dataset.selection['source_evidence'], ...)`.
