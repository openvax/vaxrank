# One resolved host for species-dependent defaults

Resolve the host from the configured/live reference, its retained native scope,
or the mapped host snapshot, falling back to an explicit host declaration.
Reject contradictory evidence. Store only the existing canonical host taxon;
do not add another organism field or infer the host from antigen source species
or MHC alleles (viral targets and transgenic animals are valid).

Use that host for new mRNA codon defaults and CTA routing/exclusions. Explicit
codon targets remain supported because expression targets can intentionally
differ. Unknown hosts must not silently select human codons or human CTA
behavior. Retained assemblies keep their recorded effective codon configuration
on replay. Frozen named configurations retain their explicit historical choices.

Pepsickle's mammalian model remains the default: the upstream human-only model
is experimental, trained on less data, and is not automatically preferable for
a human host. Validate an explicit human-only request against a known host.
Keep per-antigen source species and explicit MHC repertoires independent.

Implementation: consolidate host resolution in input_scope; use it at direct,
shared-input and construct boundaries; remove the human codon fallback in the
CLI/default configuration and require context for automatic API assembly;
resolve CTA host before defaults or admission; clarify help and user docs.

Verification: regression cases for human/mouse references, host-only inputs,
conflicts, unknown hosts, viral sources, HLA-transgenic mice, explicit codon
choices, native replay and CTA rejection. Run lint, the full suite, a CLI smoke,
GitHub CI, then release from clean main and verify PyPI. Finish existing cleanup
PR #779 first, preserving unrelated work and repository-local environments.

Primary contracts: DNA Chisel's CodonOptimize uses the target species codon
table and preserves translation with EnforceTranslation; Pepsickle's README
labels human-only experimental. References:
https://edinburgh-genome-foundry.github.io/DnaChisel/ref/builtin_specifications.html
https://github.com/pdxgx/pepsickle
