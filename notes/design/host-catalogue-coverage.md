# Host ORF catalogue coverage

Report the actual reference searched by the experimental shared-input workflow.
An annotated-only dog inventory must not imply that dog ncORFs were searched or
found absent. Use ORFdata's offline catalogue status API for species/assembly
availability; never infer ncORF coverage from quantitative RNA/protein datasets.

Keep sequence identity and self-match behavior unchanged. Attach a frozen coverage
record to HostReference: reference identity, retained annotated/extra-ORF counts,
exact source layers and their existing selections, and the ORFdata version and
compatible catalogue releases observed when the reference was prepared. Separate
registered availability from installation and from inclusion in the inventory.
Preserve opaque merged/custom source labels; unresolved catalogue attribution is
unknown, not an assertion that another catalogue was unsearched.

Allow an optional `catalogues` declaration on the existing host-reference mapping
and Python wrapper. It pins the exact source labels recorded in the supplied
snapshot; an empty list explicitly selects annotated-only coverage. A mismatch,
including requested extras absent from a dog reference, fails before prediction.
Composition/acquisition remains in ORFdata; this declaration never downloads.

Native reference transport retains the coverage record and its digest. Reload
validates it against the actual snapshot without querying the current registry.
Historical files without this field report unrecorded historical availability.
Combining references must not silently replace a frozen coverage record.

Show the same coverage in the Python API, native outputs, run summary and
HTML/PDF/text reports, including zero-candidate runs and relocated offline replay.
Use existing report rendering and resource transport. Validate real registry dog,
mouse and human cases; included, filtered-to-zero, custom/merged and unrecorded
sources; declaration errors before acquisition; tamper detection; and actual CLI
export/replay after disabling registry access. Run lint, full tests and CI before
merge, then the mandatory clean-main deployment and published-wheel acceptance.

Scope: remaining coverage declaration/reporting in #625. A full mouse CTA host
reference, new biological datasets and dependency-pin coordination are separate.
