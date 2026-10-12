# Core delivery for all candidate SLP windows

Consume saved, scored EpitopeDataset evidence and its explicit antigen contexts.
Enumerate all target-containing windows in a requested length range before the
usual score-fraction/construct-count selection. Preserve each scored target
occurrence, allele, score expression and source provenance. Use Topiary's existing
representative policy for combined observations; never rescore binding evidence.

Pass candidates to mhctools 3.53's dose-normalized core exposure adapter. A profile
supplies explicit compartment rates and fragment cleavage assessments, with
per-candidate overrides. Missing context or kinetics remains visible. Export the
editable request, complete JSON result, per-core CSV and readable report. Keep
this an opt-in analysis command; do not change existing design selection.

The target interval is the scored ligand when no minimal binding core was
supplied. Neither serum motif scores nor assay half-life predictors become IM
rate constants. SLP terminal chemistry is explicit; unsupported chemistry stays
unassessed. Actual mRNA constructs and transcript doses use mhctools' mRNA input
schema; do not invent a transcript from a protein window.

Validate enumeration before pruning, multiple cores and alleles, repeated
occurrences, retained source metadata, missing context, conditional analytic
arrival, saved-dataset/CLI composition and no output overwrites. Run lint, full
tests, CI, then merge and deploy; verify the published wheel.
