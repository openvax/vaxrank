# Evidence and provenance

Vaxrank retains source observations while selecting a representative for a
vaccine design. The experimental shared model distinguishes:

| Entity | Identity |
| --- | --- |
| Event | Patient/sample, reference and variant or structural event |
| ORF hypothesis | Transcript path, coding sequence, frame and completeness |
| Protein | Exact translated sequence and its known extent |
| Peptide occurrence | Sequence, position in a protein/window and flanks |
| Candidate pMHC | Peptide and normalized MHC allele, linked to occurrences |
| Observation | Source file/row, model/version or assay, sample, value and unit |

Several ORFs may encode one protein, and several proteins may contain the same
peptide. Sequence matches retain every source mapping. Context-dependent
predictions remain attached to their occurrence.

## Combining inputs

RNA evidence keeps its measured level: gene, transcript, ORF, variant, read,
fragment or molecule. Counts from separate reports are retained separately.
Ranking chooses among observations with a [scoring policy](selection-policies.md).
It does not sum callers' expression estimates or calibrate different models.

Fresh prediction evaluates reported peptides on the requested alleles. Additive
prediction attaches named features to existing candidates. Neither mode scans
unreported peptide windows. See [prediction modes](input-workflows.md#generate-common-predictions-and-re-rank).

## Inspect and replay

`EpitopeDataset.evidence_views()` exposes event, ORF, occurrence and RNA
relationships supplied through Topiary's reconciliation API. Native exports
retain source observations, scoring policies, antigen references and saved
products. Replayed designs use that evidence without rerunning upstream tools.

Construction still selects one source-derived window per grouped external
source. Full ORF-aware selection remains integration work. The
[input guide](input-workflows.md#replay-a-saved-vaccine-design) documents current
replay behavior and the [native reference guide](native-references.md) explains
which files to move with an export.
