# Unified antigen inputs and target reports

## Outcome

Run ordinary somatic variants, fusion/SV reconstructions, ORFdata-backed ncORF
mutations, and CTA expression through the same prediction, ranking, design and
report pipeline. Source switches control eligibility, including a separate switch
for ncORF mutations. Existing source adapters retain their biological evidence;
no source needs a bespoke design script. Native exports preserve the information
needed to reproduce both selection and the report without live reference queries.

## Implementation

1. Remove the remaining prohibition on combining direct VCF/BAM and native
   sequence sources in a patient-scoped manifest. Validate the complete scope
   before predictions. Keep source-specific reconstruction and evidence controls.
2. Retain compact biological context on each candidate: source event, coding
   model, altered residues/junctions, RNA support, catalogued translation evidence,
   and available reference/patient-background comparisons. Ordinary local RNA
   windows do not require observation of an entire protein. Unresolved frames and
   absent somatic provenance remain explicit rather than silently asserted.
3. Consume OncoRef's verified current sequence bundle for CTA peptide targeting.
   Exact containment in an Unsafe/exclude segment blocks that peptide as a CTA
   target, retaining every matching occurrence and the policy/reference identity.
   Review and absent coverage remain distinguishable from exclusion. Never join
   across ambiguous residues, use a blanket CTA-gene exemption, or describe an
   unmatched peptide as clinically safe. Saved decisions replay offline.
4. Add a compact biological overview to normal HTML/PDF/text reports and retain
   corresponding fields in native/tabular output. Show fusion/SV boundaries,
   ncORF type and translation evidence, tumor/reference differences, and CTA
   eligible/excluded peptide intervals with matching Unsafe proteins. Detailed
   source evidence remains available below the overview.
5. Exercise mixed inputs through CLI selection, peptide and mRNA products,
   native replay and report rendering. Check each source switch independently,
   repeated sequence occurrences, incomplete evidence and unsafe shared peptides.

## Verification and delivery

Use an isolated repository-local environment; preserve the shared environment and
untracked work. Run focused regression tests, lint, the complete test suite and
an equivalent CLI/report smoke with visual inspection. Bump the version, open a
PR, require green CI, then merge and deploy under the repository release gates.

## Scientific contracts

- OncoRef `docs/sequence-bundle.md` and `docs/cta-selection.md`: current verified
  sequence roles, all source occurrences, interrupted segments, coverage gaps.
- Mudge et al. (2022), https://pmc.ncbi.nlm.nih.gov/articles/PMC9757701/:
  ncORF classes and the distinction between Ribo-seq translation evidence and
  evidence of a stable functional protein. Reference evidence and patient RNA
  support are reported separately.
