# Experimental species CTA expression

Use OncoRef's frozen species evidence with ORFdata host sequences in the existing
cta_expression manifest/API. Human defaults and archives retain their contract.
Nonhuman inputs require the exact mapped annotation, a protein-source checksum,
versioned transcript/protein pairs, and explicit expression units/provenance.

Retain all candidate facts, ambiguous mappings, normal-tissue counterevidence and
coverage gaps in native outputs. A complete unique mapping permits consideration,
not selection. Require target-specific CTAOverrideEvidence plus the expression
threshold for nonhuman design. Exclude only that reviewed target gene from its
own self comparison, preserving every other gene and all extra ORFs. Unsupported
species and wrong references fail before predictor acquisition.

Test the real released canine snapshot against synthetic host sequences, checked
annotation identities, mixed manifests, wrong species, missing/low expression,
ambiguous mappings, self-sharing, native export and offline replay. Run lint,
full tests, strict docs and a real CLI smoke before PR merge/deploy.
