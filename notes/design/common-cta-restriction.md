# Common experimental CTA restriction policy

The shared expression input must evaluate the same evidence in the same way for
human, dog and mouse. Taxon scopes the reference; it does not decide admission.
Keep the existing human default and saved decisions intact. Select the new policy
explicitly in the existing `cta_expression` API/manifest.

Declare the common policy's OncoRef floor in the `cta` install extra. Resolve all
CI features together and run `pip check`; a base-only install may select the
historical owner version through PirlyGenes' exact pin. Keep that lower endpoint
usable for existing workflows, with a capability error before source access when
the common policy is requested. The owner-pin coordination is already tracked in
pirl-unc/pirlygenes#656; do not weaken its materialized-data contract here.

Use a named, versioned policy with explicit native RNA units, required somatic and
reproductive tissues, quantitative bounds and independent sample/donor minima.
Do not supply a purportedly validated biological threshold as a library default.
The input bundle binds its exact OncoRef snapshot and full ORFdata host inventory.
It retains source checksums, study lineage, global donor/sample identities,
pooling, exclusions and independently allocated transcript RNA bounds. Gene RNA
in the patient table remains gene RNA and selects only a reference protein.

Enumerate every host occurrence of the selected full protein, including extra
ORFs. Require coverage for every occurrence; unknown or ambiguous allocation
cannot become zero. Pooled samples and technical replicates cannot satisfy an
independent-donor minimum. Retain all supplied observations and all OncoRef facts,
including excluded observations. Positive somatic evidence is not removable by a
source exclusion. Qualitative, unmatched-unit or unresolved-reference facts remain
unknown. A target passes only when the declared evidence policy and patient RNA
threshold both pass. Keep all other loci and extra ORFs in downstream self checks.

Reference adapters depend on declared identifier namespaces and source roles.
Versioned Ensembl and RefSeq transcripts are matched exactly. A candidate-only
reference subset cannot stand in for the complete host protein inventory.
Use the annotation provider's metadata labels for Ensembl host scopes and verify
the exact owner assembly accession and release before acquisition. Assembly
display names and accession strings must not produce a false mismatch against a
matching variant annotation; source names are compared without case differences.

Freeze the policy, evidence, all equivalent occurrences and assessment reasons in
native results and antigen evidence. Validate the same assessment on replay. Reuse
existing target masks, native intervals, MHCtools/Topiary capability checks and
final-product safety audits. Research export remains draft; restriction evidence
does not establish binding, presentation, processing or immunogenicity.

Validation uses released snapshots for held/unknown cases and explicitly
synthetic complete evidence for a positive case in each species. Include the
other-locus heart counterexample, missing RNA/tissues, duplicate and pooled donors,
wrong sequence/reference identities, frozen replay, native constructs, unsupported
DLA predictions and boundary windows. Run the full suite and real CLI smoke.

Primary evidence: [Hindriks 2025](https://pubmed.ncbi.nlm.nih.gov/40359694/)
reports canine CTCFL/DDX53 kidney expression; [Sarma 1999](https://pubmed.ncbi.nlm.nih.gov/10049945/)
reports P1A expression in normal mouse tissues. These observations must survive
selection and remain distinct from tumor RNA and antigen presentation.
