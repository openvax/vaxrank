# Shared Isovar antigen read controls

Coordinate openvax/isovar#640 with the released Isovar 1.76 API. Permit the
existing `options.reads` contract for viral and antigen-isoform reconstruction,
using the producer's shared RNA-path defaults. Record the resolved settings;
saved-result replay must still reject requests to reconstruct with new controls.
Keep typed data/index identities out of scalar reconstruction options. Reject
an input ReadCollector instead of silently replacing it: Vaxrank's adapter uses
`options.reads`; custom callbacks/profiles can be applied in Isovar before replay.

Accept the reviewed OncoRef compatibility interval 1.19.2 through 1.22.0.
Exact simultaneous consumer pins create a circular release constraint:
PirlyGenes tests the latest published Vaxrank, while Vaxrank reaches PirlyGenes
through Topiary. This bounded interval allows Vaxrank to ship against the old
consumer and the jointly tested 1.22.0 stack before PirlyGenes moves its pin.
It preserves legacy CTA admission/self-exclusion APIs and records actual
owner/data/snapshot versions; it does not adopt the distinct `cta_selection`
policy (vaxrank#687). Verify both endpoint environments explicitly.

Verify both live families, default and overridden policy, inappropriate input
rejection, saved replay, CTA/peptide exclusions, joint installation, the full
suite and an equivalent CLI/report smoke before release.
