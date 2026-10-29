# Shared Isovar antigen read controls

Coordinate openvax/isovar#640 with the released Isovar 1.76 API. Permit the
existing `options.reads` contract for viral and antigen-isoform reconstruction,
using the producer's shared RNA-path defaults. Record the resolved settings;
saved-result replay must still reject requests to reconstruct with new controls.
Keep typed data/index identities out of scalar reconstruction options. Reject
an input ReadCollector instead of silently replacing it: Vaxrank's adapter uses
`options.reads`; custom callbacks/profiles can be applied in Isovar before replay.

Coordinate the exact OncoRef owner pin with PirlyGenes at published 1.22.0 so
the forthcoming Isovar catalogue API can be installed in the same stack. This
package upgrade preserves the legacy CTA admission/self-exclusion APIs and
does not adopt the distinct `cta_selection` policy (vaxrank#687).

Verify both live families, default and overridden policy, inappropriate input
rejection, saved replay, CTA/peptide exclusions, joint installation, the full
suite and an equivalent CLI/report smoke before release.
