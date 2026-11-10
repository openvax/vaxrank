# Ranked selection with explicit antigen-type limits

Remove category-count balancing from new peptide and mRNA designs. With no
explicit type limits, preserve the existing score/allele-coverage order after
sequence deduplication. Minimum goals and hard maximum caps remain opt-in; a
candidate does not gain priority merely because its type is underrepresented.

Remove the `--antigen-allocation` CLI selector. Default configuration and the
Python policy resolve to ranked selection; YAML no longer needs a mode field.
Keep historical mode values readable for checksum-verified native replay, but
reject attempts to assemble a new design with `mode: balanced`. Do not retain
the old balancing algorithm or silently reinterpret a saved balanced design.
An explicit ranked policy can be used to redesign an archived candidate pool.

Verify that a weak viral candidate cannot displace stronger candidates by
default in either modality; explicit type goals/caps still work. Capture a
small synthetic archive with the released implementation before changing it,
then prove offline replay preserves its exact products and audit. Exercise
CLI defaults, obsolete-option rejection, native redesign and the mixed-design
report smoke. Run lint, full tests and GitHub CI before merge and release.
