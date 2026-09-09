# Processing and final-construct audit PR series

The immediate correctness fix is #424: a processing result must identify the
peptide's resolved occurrence, not just its sequence and source. Repeated
peptides can have different local cleavage scores. The canonical record and
report join will use `(peptide_sequence, source_sequence, peptide_offset,
predictor_name)`, where the offset is zero-based in the scored source. Existing
offset relocation must be shared by annotation and lookup. There will be no
coordinate-free fallback. Ranking and model inference remain unchanged.

Acceptance: input-order-independent repeated-occurrence results and report
values, relocation parity, same-occurrence sharing across alleles, terminal
positions, one inference per source, full tests, report smoke and release.

Preferred delivery order after that fix:

1. #422: validated site-level processing evidence and complete-sequence audits,
   with explicit context, coverage, predictor identity and target-ligand overlap.
2. #421: native-versus-manufactured sequence and modification provenance, then
   reassessment of the actual modified sequence and new boundaries.
3. #311: integrate complete-window non-CTA self-risk assessment into product
   audit/reporting before its separately reviewed enforcement policy.
4. Extracellular/serum processing evidence, after primary validation and model
   availability review in openvax/mhctools#278/#279 and integration under
   openvax/topiary#288. A whole-peptide half-life model
   is not a cleavage-site predictor or a pMHC dissociation-stability model.
5. #414 and #423: independent final-selection and manufactured-construct
   comparisons against the documented Sid vaccine sequences.

RNA-source expansion stays in openvax/isovar#218; cell/UMI support units and
cell-stratified reconstruction are tracked separately in openvax/isovar#226.
Consume reviewed releases without modifying the active upstream worktree.

Across the series, use target antigen/target ligand terminology. Mutation, CTA
and viral provenance and targetable masks remain distinct. Scores are evidence,
not proof of antigen destruction, presentation, TCR recognition or safety.

Primary processing scope:
https://doi.org/10.1093/bioinformatics/btab628
https://doi.org/10.1371/journal.pone.0178943
