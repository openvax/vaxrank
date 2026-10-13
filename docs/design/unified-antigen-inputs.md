# Experimental unified antigen inputs

Vaxrank should use one sequence/evidence model, prediction path and design
representation for mutation, fusion, CTA, ERV, splice and viral antigens.
Source adapters resolve evidence; vaccine type controls downstream assembly.

This increment adds `AntigenDataset`: named `VaccineAntigen` records with a
declared patient, host reference and genotype. Its versioned JSON is accepted
by `--input-antigens`. A shared prediction API returns `EpitopeDataset`, so
selection, native replay, reports and peptide/mRNA outputs use existing code.
CTA expression prediction also uses this API. Admission decisions remain
explicit evidence supplied by the adapter or caller, never inferred from a
source label or sequence alone. Held-out records remain in the saved evidence.

`--output-antigens` exports resolved antigen records and their input provenance
in the same JSON format. Prediction-free input and saved prediction evidence
remain separate: replay uses `--input-epitopes`, while `--input-antigens` requests
new predictions against an explicitly selected host annotation.

Validation covers all six kinds, repeated sequences with distinct sources,
held-out records, patient/reference/genotype conflicts, shared CTA prediction,
and peptide/mRNA design followed by offline native replay. Existing mutation
and CTA policies remain unchanged. All new interfaces are experimental.

Later increments can put expression and other upstream adapters behind one
manifest, preserve unresolved sequence hypotheses in the common inventory,
add source allocation policies, and consolidate final-product audit outputs.
