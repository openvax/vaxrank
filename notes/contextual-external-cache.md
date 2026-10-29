# Contextual external prediction caches (#624)

Use Topiary's released CachedPredictor through the existing public
predict_peptide_occurrences path. Remove the obsolete peptide-only guard and
explicitly select the patient's genotype from potentially broader cache data.
Do not infer historical model versions from the installed fallback.

Keep fresh/additive inference opt-in and preserve primary/comparator occurrence
identity and the distinction between unknown and known-empty flanks. Topiary
owns coverage checks, contextual lookup, and fallback model/version validation;
Vaxrank must not recreate that machinery. Protein-scan and explicit peptide-only
semantics are unchanged.

Exercise the real CLI with a partial contextual cache and extra cached alleles,
context-specific misses and comparator hits, mismatched fallback versions,
unknown versus terminal contexts, and prediction-free native replay. Check
stored historical evidence and duplicate source occurrences survive unchanged.

Native frames must carry the exact source window and local peptide offset.
Topiary derives current inference flanks from those fields; Vaxrank preserves
recorded flank fields as historical invocation metadata. This keeps fresh,
native and additive paths on one occurrence-context rule without a sidecar.
