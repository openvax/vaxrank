# MS-supported prediction scores

Enable the opt-in `ms-evidence-v1` policy with `--ms-evidence`. Hitlist is a
regular dependency; installed, current data are optional. Missing or stale data
remain explicitly unavailable and add no support. `--ms-evidence-package DIR`
also enables the policy and uses a verified Hitlist MHC evidence package
(using the producer's package verifier). A multi-modality EvidenceStore requires its own typed relation adapter
and is currently rejected instead of being treated as the older package format.

```yaml
evidence:
  ms:
    lookup: required             # auto (default), off, required
    package: /path/to/mhc-package
    parameters:
      gain: 1
      max_relaxation: null       # uncapped; an explicit cap must be >= 1
      reference_discount: 0.25
      prediction_ic50_cutoff: 50
      direct_admission: true
      weights:
        sequence: 0.25
        predicted_sequence: 0.5
        typed: 1
        recorded: 1
        predicted: 3
        experimental: 6
        monoallelic: 8
```

`EpitopeConfig(ms_evidence={...})` accepts the same object as `evidence.ms`.
CLI config overrides work, for example
`--config-value evidence.ms.parameters.gain=2`.

For each exact peptide and patient allele, use the strongest applicable facet
of each eligible observation. Exact peptide sequence, chemical form, MHC species
and class constrain transfer. Unresolved class can contribute sequence support
only. A typed multi-allelic source is not an experimentally restricted ligand.
A qualifying class-I prediction is the selected affinity model's IC50 <= 50 nM;
it contributes the prediction-supported facet only for the compatible allele.
The current native package adapter does not assert source-prediction qualification
without its model provenance. Unknown chemical forms and contaminants are held out.

Patient prediction alleles now use Topiary's exact-allele boundary consistently:
allele groups, null alleles, and ambiguous class-I pairs are rejected at input.
Source evidence can retain those unresolved declarations without promoting them
to an exact restriction or assigning them to every patient allele.

Resolved native donors or engineered origins contribute at most once each.
`R` sums those unit maxima, `P` is the largest presence contribution, and `J`
is discounted support across connected references. Known shared units and reuse
identities connect references before quality exclusion. Descriptive labels and
row counts never establish independent donors. Packages lacking reviewed lineage
can contribute P and J; they do not acquire invented recurrence counts.

Topiary evaluates `M=max(R,P,J)`, `S=log2(1+M)`, `G=1+gain*S`, and `C=B^(1/G)`.
An optional cap applies to G. B is the bounded prediction score before its soft
cutoff. The default is normalized logistic affinity. A custom `base_score` needs
an explicit `soft_threshold`; `log_base_score` may supply a matching log value.
These weights and scores are heuristic ranking choices, not probabilities.

C can pass the soft gate; the user's `epitopes.filter_expr` remains a hard gate.
Exact experimental/verified monoallelic support can admit a candidate directly
only with positive support and G>1. Missing predictions remain missing unless an
explicit `ms_only_score` in [0,1] is configured. A true zero stays zero.
Native protein scans enumerate the models' declared peptide lengths and alleles
independently of returned predictions; absent measurements remain blank rows.

Each admitted allele retains its prediction-only baseline contribution b0.
An epitope adds the single increment `max(C-b0)` to its baseline allele sum.
A window adds this increment once per unique peptide sequence. It is never
multiplied by the number of matching alleles or repeated sequence positions.

Raw Hitlist captures, normalized observations, availability, per-allele baseline
and adjusted values, R/P/J, M/S/G, and Topiary policy decisions travel with native
datasets. Replay uses those captured observations without querying installed data.
Report calculation details show the baseline allele contributions and MS increment.
Keep native datasets to retain the complete acquisition and replay contract.

See [the full design](https://github.com/openvax/vaxrank/blob/main/notes/design/antigen-ms-evidence.md) for richer
EvidenceStore relationships, source quality calibration, and evaluation work.
