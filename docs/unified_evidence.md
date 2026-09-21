# Combining ORF hypotheses and their evidence

The combined view should link biological hypotheses to observations, rather
than flattening every caller's result into one row per variant. This is the
design contract for the generalized Topiary/Vaxrank workflow; the external
rescoring entry point currently implements the peptide/context portion only.

## Identities and relationships

| Entity | Identity and information retained |
| --- | --- |
| Biological event | Patient/sample, reference identity, normalized variant or structural event; uncertainty and linked events remain explicit. One event can have several transcript/ORF hypotheses. |
| ORF hypothesis | Transcript path, coding nucleotide sequence, frame, initiation/termination and completeness when known, linked events and source references. Missing boundaries are unknown, not inferred from a short peptide. |
| Protein sequence | Exact translated amino-acid sequence and its extent (complete protein or partial window). Several ORFs can encode the same protein sequence, including synonymous nucleotide differences. |
| Peptide occurrence | Peptide, location within a specific protein/window and available flanks. The same peptide can occur at several positions and come from several genes or ORFs. |
| Candidate pMHC | Patient/sample, peptide and normalized MHC allele, with links to all source occurrences. Context-dependent predictions remain attached to their occurrence. |
| Observation | Input artifact/row, tool/model/version, sample/library, measured entity, quantity, unit, value and missingness. RNA abundance and predictions are separate observations. |

Exact full protein matches may share a protein-sequence node while retaining
distinct ORF hypotheses and every source observation. Different proteins for
the same event remain alternatives; overlapping short windows only establish
agreement in the overlap. A matching peptide never establishes ORF identity.
Do not merge contradictory novelty or tumor-specificity annotations.

For example, LENS and Exacto might nominate one protein sequence with separate
RNA estimates, while a pVACseq table supplies a matching peptide-HLA prediction.
The peptide links to that protein only when its occurrence is established; a
sequence match alone does not select its gene of origin. A second Exacto ORF
with a different frame remains a second hypothesis. Shared peptides may reuse
compatible binding predictions, while distinct flanks require distinct
processing predictions.

## Combining evidence and choosing vaccines

Keep RNA quantities at their stated level: gene, transcript, ORF, variant,
reads, fragments or molecules. Retain assay, sample, library and quantifier
provenance. Do not promote gene TPM to ORF abundance, add two pipelines' TPMs,
or count rediscovery of the same reads as independent support. Combine counts
only when their underlying evidence units and overlap are known; otherwise
retain separate estimates and select one with an explicit ranking policy.

Original predictions act as sparse cached evidence. A missing pMHC prediction
on an Exacto ORF is unknown, not zero binding. Generating peptides from that
ORF expands the candidate universe and must be requested separately. Fresh
prediction adds a new observation carrying the actual model/version and input
context. Only compatible peptide, allele/genotype and context queries can
share prediction work. Ranking selects which observations to use through the
Topiary DSL; mixed raw scores are not automatically calibrated.

Expose both hypothesis-level comparisons and peptide-HLA coverage across
hypotheses. Selecting a vaccine window must name its actual supporting
hypothesis and tumor-specificity evidence. Multiple callers agreeing on it do
not multiply its score or coverage. Ambiguous ORFs remain visible even when a
selection policy chooses one construct. Full ORF-aware coverage/selection is
upstream integration work, not implemented by the current external entry point.

## What the current PR provides

LENS/pVACseq source observations retain separate prediction identities, values
and context. The report exposes exact peptide/context sequence hashes and the
reported sequence extent. These hashes are sequence-only grouping aids, not
sample-scoped biological identities or a full ORF registry. All contexts remain
in the report; the existing vaccine path selects one source-derived window per
biological source. No RNA counts or scores are summed across input reports.

The common-model mode scores reported peptides on the explicitly requested
HLA set; it can therefore add peptide-HLA combinations, but does not scan new
peptides. The historical mode preserves source-local model selection and
scoring. Native exports retain the original candidate predictions and context;
the CSV report retains source-specific evidence alongside fresh values.

The reusable table workflow is tracked in
[Topiary #366](https://github.com/openvax/topiary/issues/366), native Exacto
ingestion in [#365](https://github.com/openvax/topiary/issues/365), and consumer
adoption in [Vaxrank #497](https://github.com/openvax/vaxrank/issues/497).
These should own the general identity and evidence model rather than adding
independent ORF reconciliation rules to each consumer.
The detailed reconciliation contract is
[Topiary #370](https://github.com/openvax/topiary/issues/370); the translation
and RNA evidence export is [Isovar #324](https://github.com/openvax/isovar/issues/324).

## Source format references

[Exacto](https://github.com/pirl-unc/exacto) describes full-transcript
translation with amino-acid variant annotations. The
[pVACseq output documentation](https://pvactools.readthedocs.io/en/latest/pvacseq/output_files.html)
distinguishes peptide tables, aggregated selections, and the separate
mutant/wild-type FASTA artifacts. Importers must preserve these differences in
extent instead of claiming every table supplies a complete ORF.
