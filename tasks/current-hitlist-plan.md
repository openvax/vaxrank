# Current Hitlist compatibility (#676)

The `hitlist` extra excludes all Hitlist 2.x releases, preventing a current
OpenVax stack from installing together. Validate the public portable evidence
contract with current Hitlist 2.0.35 and retain the existing archived 1.65.1
producer fixture. Broaden the extra to the validated major series only after
the current writer/verifier, CLI capture, metadata mismatch rejection and
offline replay tests pass. No scientific rules, cache paths, saved artifact
bytes or original producer provenance should change.

Validate a tiny synthetic evidence bundle produced with current Hitlist;
exercise the real writer and verifier rather than duplicating their
implementation. Keep the existing lower-bound compatibility coverage. Run the
full test suite, lint, package checks and the real combined installer with
PirlyGenes' current-OncoRef candidate. Include current Hitlist in CI so the
extra cannot silently resolve to an older major. Bump the patch version, merge
only after green CI, deploy from clean main and verify public artifacts.

This is coordinated with pirl-unc/pirlygenes#652. Both releases must retain the
same validated OncoRef pin. Existing user environments and caches are tested
through isolated environments and copy-on-write copies.

## Validation before publication

- The 14 archived-bundle, capture and replay tests pass with Hitlist 1.65.1
  and 2.0.34; the full candidate suite uses current Hitlist 2.0.35.
- A disposable integration test against Hitlist 2.0.35 uses its real writer
  and verifier, then Vaxrank capture. It verifies exact peptide capture,
  annotation mismatch rejection and artifact tamper rejection.
- The combined current OpenVax environment resolves and passes `pip check`.
  PirlyGenes keeps its existing 6.0.13 data bundle and both libraries use
  OncoRef 1.19.2.
- Integration tests use copied caches and a filesystem guard denying writes
  to existing user caches. OncoRef shard-cache configuration is tracked
  separately in pirl-unc/oncoref#626.
