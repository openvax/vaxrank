# Isovar support compatibility

Vaxrank requires Isovar 1.90.9 or later, below 2.0, with ORFdata 0.52.0 or
later and Osteosarc 0.15.6 or later in the 0.15 series (the pins in
`requirements.txt`). Native
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
their admission requires an explicit tumor-specificity attestation or a somatic
variant declaration for a reconstructed changed sequence. Reports
distinguish that assertion from measured sample evidence and catalog attribution.

New direct fragments store the actual runtime Isovar version as
`sequence_source_version`, carried into mixed-input evidence and native replay.
An older fragment without this field retains an unknown version. Historical
fixture producer versions remain immutable and are never rewritten to the
current installation.

## Shared-sequence RNA set unions

By default (`vaccine_peptides.shared_sequence_rna_policy: source`), each source
gene or transcript is scored from its own `RNAFeatureEvidence` record, keeping
`openvax-v1` and historical replay unchanged. When multiple origins encode the
exact same candidate vaccine amino-acid window (for example paralogs or
overlapping isoforms), adding per-origin counts can double-count shared reads
and selecting one origin's count omits distinct reads supporting other origins.

Setting `--shared-sequence-rna-policy within_library_union` (or
`vaccine_peptides.shared_sequence_rna_policy: within_library_union` in YAML)
computes `isovar.union_rna_support` across applicable `RNAFeatureEvidence`
records sharing that exact candidate sequence within each
`evidence_scope = (sample_id, source_id)`:

- `rna_shared_reads`, `rna_shared_fragments`, `rna_shared_reads_per_million`,
  and `rna_shared_fragments_per_million` bind in `combined_score_expr` when a
  single library scope is active and Contributing measurements agree on
  `fragment_mode` and whole-library denominators (`library_reads`,
  `library_fragments`).
- Individually spanning reads and paired-mate-only fragments remain distinct;
  distinct libraries and timepoints are never pooled.
- Group effective-length quantities (`rna_shared_effective_length_*`,
  `rna_shared_*_per_effective_kb*`) remain unavailable
  (`shared_sequence_effective_length_undefined`) because a read-set union across
  alternative transcript paths has no shared effective sampling length.
- Candidate and shared-product reports record the union counts, naive sum across
  origins, per-origin shared/exclusive counts, `evidence_set_id`, and
  attribution limitations, and native exports deduplicate shared union records
  in-file.

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
sidecar, reconstruction or prediction. Shared RNA feature measurements,
whole-library sampling profiles and multi-epitope antigen/fragment records are
stored once per unique SHA-256 content identity inside native exports and
validated on load before restoring per-occurrence attachments. API callers can
pass the source catalog as
`TemplateDataCreator(source_evidence=dataset.selection['source_evidence'], ...)`.
