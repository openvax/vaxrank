# Ordered consistency fixes

Requested order: #871, #836, #837 and #840, then #793.

1. Restore legacy protein-hypothesis reporting when external source evidence is
   missing or ambiguous. Extract only available origin facts; retain strict
   producer-input validation and unchanged saved antigen identities. Cover
   complete and partial evidence through rendered reports. Release this fix first.
2. Replace CTA and tissue-risk representative-transcript ranking with Topiary's
   public rule, preserving each caller's candidate filters and all-isoform mode.
   Verify the biological semantics against primary references and raise the
   dependency floor to the first compatible public API when needed.
3. Share expressed-allele validation through mhcgnomes/Topiary, then prevent
   explicit genome pins from being silently overridden by --ensembl-release.
   Exercise CLI and library inputs, including nonhuman alleles and conflicts.
4. Finish #793 by replacing the remaining hand-written FASTA parser with the
   shared reader while preserving the caller's documented behavior.

For each PR: include a version bump, run lint, the full suite and applicable
CLI/report smokes, require green CI, merge, deploy from clean main, and verify
PyPI and release artifacts. Preserve user files and the shared Python environment.
