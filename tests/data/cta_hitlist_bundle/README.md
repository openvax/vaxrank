# Synthetic portable Hitlist evidence fixture

This is software test data, **not clinical or published MS evidence**. Hitlist
1.65.1's public `write_cta_evidence_bundle` produced the bundle, and its unmocked
`verify_evidence_bundle` verified it. The fixture includes a tiny synthetic
IEDB-shaped input (fake PMID 99999999): two source copies of one MS observation
and one fluorescence observation. Only the MS observation counts as presentation;
all three contributors and the rejected observation remain auditable.

The synthetic PRAME/non-CTA mappings and atlas-shaped benign-tissue donor rows
exercise shared-source and tissue-risk provenance. They make no biological claim
about SLYNTVATL, PRAME or its presence in any patient. CTA membership/identifier
resolution was mocked during export; the writer, observation builder, hashes,
relationships, positive-MS selection and verifier were real. The original paths
in provenance deliberately do not exist on CI: verification is portable.

`input-expression.tsv` is the exact hashed expression input. The manifest records
OncoRef 1.8.207 and Ensembl 93. Regenerate the bundle when changing the producer or
the tested reference version rather than relabeling captured provenance.
