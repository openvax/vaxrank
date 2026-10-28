# Interpretable selection, replay and reports

The October 2026 example audit found successful executions whose outputs were
misleading or difficult to reuse. This series addresses the underlying contracts,
with separate PRs and releases for independently reviewable changes.

1. **Empty results are results (#705).** A valid native result with no candidate
   rows must replay its scope, policies and empty products without acquiring
   predictions or inventing evidence. Empty arbitrary external files must still
   fail validation. Exercise both the public Python loader and CLI replay.
2. **Products and sources have different identities (#703).** Apply design caps
   to distinct selected products, retain every supporting source and evidence
   mapping, and make overlapping alternatives from one event explicit. Cover
   peptide and mRNA outputs, minimal epitopes, user-requested alternatives, and
   deterministic replay. Do not sum support across alternative reconstructions.
3. **Reports describe what was emitted (#704, #700, #702).** Lead with final
   products and exact sequences, then distinguish ranked candidates and held-out
   evidence. Use source-neutral labels, local styling and shared report data.
   Remove redundant copies of the same candidate report. Check rendered PDFs
   against the actual manifests and FASTAs.
4. **Export and replay should be predictable.** Review source-dependent report
   flags, duplicate serialization, optional decision outputs, reference resource
   bundles and documentation. Give the same option the same meaning across input
   routes. Keep required resources discoverable and verifiable; remove redundant
   artifacts rather than hiding required evidence or adding another wrapper API.

Initial architectural critique:

- `prepare_reports` recognizes some kinds of zero-row evidence but rejects an
  equally valid empty native run. Validity should come from the input contract.
- Distinct candidate IDs currently consume separate product slots even when the
  emitted sequences are identical. Grouping must preserve source multiplicity.
- `emit_outputs` ignores `--output-xlsx-report` on native input, while
  `--output-csv` changes table meaning by input route. These are surprising public
  contracts, not merely missing help text.
- Default report naming and repeated per-modality candidate reports obscure the
  difference between a source window and a selected manufactured sequence.
- Native reference resources serve real portability and checksum requirements.
  Their necessity must be evaluated separately from optional policy reports and
  repeated object payloads; deleting them indiscriminately would lose evidence.

Acceptance uses the preserved Sid, COLO829 and synthetic source controls, plus
small focused regression fixtures. Each PR bumps the version, passes lint/full
tests and CI, and is deployed after merging. Scientific admission and prediction
policies remain explicit; these usability changes do not assert clinical safety.
