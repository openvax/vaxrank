# Scoped SV placement support (#707)

Consume the additive `occurrences[].rna_support` field released by Isovar
1.62.2. Validate its sample/source, read/fragment set identities, counts and
occurrence scope before using it. Check each measured set against the separate
candidate union and, when present, its original observation ledger. Retain
measured zero and the entire producer evidence without merging alternatives.

Historical v4 exports retain their current conservative behavior: no witnesses
means zero, a single occurrence can use the union, and multiple unresolved
occurrences stay held out. A malformed explicit field fails instead of falling
back to a broader union. Ambiguity, partial-sequence, minimum-support and
tumor-specificity gates remain unchanged.

Test unequal, shared and empty sets, scope/count corruption, raw-vs-export
parity, older files and unchanged admission gates. Update the dependency floor
and documentation, run lint/full tests/CI, then merge and deploy.

Evidence units follow the SAM specification (templates vs segments vs alignment
records): https://samtools.github.io/hts-specs/SAMv1.pdf, and the released producer
contract: https://github.com/openvax/isovar/blob/v1.62.2/docs/sv-rna.md.

Release order: follow the existing #711 → #712 → #715 series as version
3.47.8. This dependency coordinates release versions; the SV reader change is
independent of report/export behavior. Retarget to main after #715 ships.
