# #355: validate allele evidence at input boundaries

- Use the existing Topiary-backed `is_allele_scoped_kind` definition; do not
  mistake allele-free processing or half-life observations for malformed MHC
  evidence. Parse stated allele identifiers with mhcgnomes, including non-human
  alleles and class-II pairs, without hand-written HLA prefix rules.
- Reject malformed allele-scoped leaves before filtering/scoring in native
  flat and serialized inputs, LENS and pVACseq. Include the input filename,
  row and offending kind/comparator where those are available.
- Keep rendering defensive for directly constructed legacy objects: show an
  unavailable allele/score, never invent an allele, a zero or patient evidence.
  Missing scores remain numeric missing values in tabular output.
- Test nested non-WT comparators, blank/null spellings, valid class-II and mouse
  alleles, allele-free kinds, report rendering and actual load paths. No ranking
  policy changes. Review, lint, full tests, B16 smoke, green CI and deployment
  are required before closing the issue.
