# Sid policy comparison, 2026-10-02

The frozen `openvax-v1` epitope policy reproduced the existing scores, target
sets and selected sequences in all seven pinned Sid cases. Six passed RNA
filters; H1-2 did not. Four of the six eligible baseline sequences matched the
documented historical sequence. These are behavioral checks, not evidence of
clinical superiority or an attempt to fit undocumented provider settings.

[comparison-summary.json](comparison-summary.json) contains measured scores,
windows, self-content changes, serum-weight sensitivity, rejected-window counts,
HLA discrepancies, package versions and the full comparison's SHA-256.
[baseline.tsv](baseline.tsv) is the compact baseline comparison. The experiment
uses real pinned BAM reads and cached NetMHCpan-4.2c outputs, plus fresh MHCflurry
2.3.11 presentation predictions and Pepsickle 0.1.3 inference. The checked-in
[20-row evidence subset](../../../tests/data/selection_policy/sid-evidence.tsv)
provides an offline regression for all four epitope policies.

## Observed window changes

The first two experimental policies (presentation and equal-weight
presentation/BA) agreed on windows, changing three of six eligible cases:

| Case | Frozen baseline | Presentation / BA blend |
| --- | --- | --- |
| DYNC1H1 minimal | `KRFHATISF` | Same |
| DYNC1H1 long | `LKHGKRFHATISFDTDTGLKQALETVNDYN` | Same |
| DYNC1H1 JLF | `LKHGKRFHATISFDTDTGLKQA` | Same |
| DYNC1H1 CeGaT | `KHGKRFHATISFDTDTGL` | `HGKRFHATISFDTDTGLK` |
| EXOC4 JLF | `SVIRTLSTIDDVEDRENEKGR` | `ISVIRTLSTIDDVEDRENEKG` |
| EXOC4 CeGaT | `SVIRTLSTIDDVEDREN` | `ISVIRTLSTIDDVEDRE` |

Adding the Pepsickle C-terminal multiplier changed only the minimal DYNC1H1
case, to `HGKRFHATI`. Just **5.9%** of raw prediction rows in that short source
context had an unpadded, internal C-terminal score. Other cases had 66.3–77.4%
coverage. The extra short-window change is substantially a missing-context
effect; it is not evidence for promoting the cleavage policy to the default.
Coverage here is over raw rows (different kinds/models repeat occurrences),
not an independent biological sample count.

With self-weight 0.25, the 30-aa and 22-aa DYNC1H1 windows both trimmed to
`KHGKRFHATISFDTDTGLKQA` (21 aa). In the BA/presentation blend, both retained
100% of the original target score, while self scores fell from 0.214315 and
0.088136, respectively, to 0.084404. No window became demonstrably free of
normal-tissue risk: the reference is only a small transcript subset.

Serum weights 0.25 and 0.5 produced the same emitted windows as the self-only
peptide refinement where the full requested enzyme panel was assessable.
They penalized 0.54033 target-score units at risk in the DYNC1H1 CeGaT blend
window. At serum weight zero, cleavage inference is disabled, so the recorded
risk of zero means **not evaluated**, not resistance. The minimal `KRFHATISF`
window was excluded at nonzero serum weights because its DPP4 qPISA evidence
was unassessed. These are first-cut susceptibility heuristics; motif rules do
not establish serum degradation rates.

## Reproduce and interpret

From the repository with its dependencies, MHCflurry presentation assets and
Pepsickle installed:

```sh
python -m examples.osteosarc_test_data.compare_policies \
  --processing --output NEW_DIRECTORY
```

The experiment overrides peptide length per historical record (9, 30, 22, 18,
21 or 17 aa); it does not claim those lengths are the bundled 25-aa default.
Self-trimming explores 15 aa through each case's original length, with the
9-aa case kept at 9 aa. Serum refinement explores within each already selected
window; shared self trimming can slide within the full available RNA context.
No sequences, affinity values, support counts or model predictions are invented.

A complete run retains raw model outputs, native candidate datasets, Topiary
policy evidence/decisions/criteria, detailed window audits, the small annotation
reference and SHA-256 manifests. The saved typed policy evidence is independently
replayable without predictors. Native mutation reconstruction currently needs
the saved reference at its original location; relocation is tracked in
[#563](https://github.com/openvax/vaxrank/issues/563).

The clinical typing includes HLA-A*01:11N; no expressed allele is substituted
for that null call. Computational typing differs, and class-II evidence is not
included. The run uses Isovar 1.39.15, matching current Vaxrank's dependency
range; consuming the already released flanking-indel fix is tracked in
[#562](https://github.com/openvax/vaxrank/issues/562). The broader audit of actual
final Sid vaccine constructs remains [#423](https://github.com/openvax/vaxrank/issues/423).

See [policy definitions and primary scientific sources](../../../docs/selection-policies.md)
for score formulas, enzyme evidence and interpretation limits.
