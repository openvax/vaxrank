# Named selection policies

Use the existing repeated `--config` option. `builtin:openvax-v1` is the frozen
Vaxrank 3.35 baseline; overlay YAML files derive new bundles. There is no second
`--ruleset` precedence system.

```bash
vaxrank --config builtin:openvax-v1 ...
vaxrank --config builtin:openvax-v1 --config study.yaml ...
```

Mappings merge recursively, scalars and lists replace earlier values, and an
explicit null clears an inherited value. CLI flags take precedence. Put a new
top-level `name` on a derived bundle. Names are labels; content hashes identify
definitions. Never edit the released `openvax-v1.yaml` to change a policy.

The bundle includes RNA-context, epitope, window, manufacturability and construct
settings. `epitopes.selection_policy` is Topiary's portable `SelectionPolicy`,
including named criteria, minimum scores, unknown handling and representative
tie-breaks. This broader bundle is intentionally not called `RankingProfile`.
For example, an overlay can contain:

```yaml
name: study-binding-v1
epitopes:
  selection_policy:
    name: study-binding-v1
    criteria:
      - name: binding
        role: score
        expression: affinity.value.logistic_normalized(350, 150)
    score_by: criterion("binding")
    filter_by: null
    score_fill: 0.0
    min_score: 0.00001
```

Authoring YAML omits Topiary's derived `expanded` definitions. Vaxrank resolves
them after composition and saves the complete canonical definition. Native
reload uses that saved definition. A policy replaces the legacy scalar scoring
formula: do not also set `filter_expr` or `score_expr`. Use policy fields to edit
its formula; `--min-epitope-score` explicitly overrides its minimum. Vaxrank
requires higher-is-better scores because it aggregates epitope contributions.

`selection_policy.json` records resolved settings, ordered configuration sources
and hashes. `policy_evidence/` contains Topiary typed evidence, occurrence
decisions and named-criterion audit tables, including rejected observations.
Topiary's `replay_selection_policy(read_tsv(path))` reproduces those decisions
without predictors. Native `--output-epitopes` datasets embed the policy evidence
and retain original predictions. Direct VCF/BAM runs also write complete policy
evidence alongside their legacy selected-epitope export. If only
`--output-epitopes FILE` is requested, the decision sidecars go in `FILE.policy/`.

The frozen definition fixes selection rules, not the measurements supplied to
them. Predictor versions, resolved defaults and source-specific evaluation
contexts are recorded at runtime. Reproducing a result requires the same input
evidence and configuration. Adding another affinity predictor can change which
method an unqualified `affinity` reference resolves to; use `default_methods`
or a qualified reference when the experiment requires a particular model.

## Experimental alternatives

Compose these after `builtin:openvax-v1`; they are opt-in and class-I only:

| Overlay | Epitope score |
| --- | --- |
| `presentation-v1` | MHCflurry presentation, with NetMHCpan-4.2 BA < 5,000 nM eligibility |
| `presentation-ba-v1` | Equal-weight MHCflurry presentation and normalized NetMHCpan-4.2 BA |
| `presentation-ba-cterm-v1` | Previous blend × (0.75 + 0.25 × Pepsickle C-terminal score) |
| `self-trim-v1` | Shared 15–25 aa window enumeration and non-CTA self-content penalty |
| `peptide-serum-v1` | Peptide SLP window refinement, self penalty and serum susceptibility penalty |

Weights are exploratory design choices, not fitted probabilities. Presentation
already incorporates processing; multiplying it again by MHCflurry's processing
head would double-count correlated information. The presentation-only and
binding/presentation alternatives provide ablations for the additional Pepsickle
weight. We have not added standalone TAP/ERAP filters: their extra contribution
and context assumptions need evaluation before introducing more correlated
scores. RNA support remains in the existing combined-score expression.

MHCflurry's integrated predictor was evaluated against held-out ligand datasets
([O'Donnell et al., 2020](https://doi.org/10.1016/j.cels.2020.06.010)). Pepsickle
predicts proteasomal processing, with different epitope-trained and digestion
models ([Weeder et al., 2021](https://doi.org/10.1093/bioinformatics/btab628)).
These results do not validate this particular ensemble or predict T-cell
recognition. NetMHCpan BA and EL are distinct outputs; these alternatives use
BA as binding evidence, rather than treating EL scores as nM affinities
([Reynisson et al., 2020](https://doi.org/10.1093/nar/gkaa379)).

Selection does not silently run new predictors. Supply real MHCflurry and
NetMHCpan observations, with `pepsickle_cterm_score` for the C-terminal overlay.
The included Sid experiment explicitly runs MHCflurry through Topiary and
Pepsickle through mhctools, using complete available source sequences. The
feature is the score at the epitope's C-terminal bond; it is missing at a source
endpoint or where Pepsickle had to pad missing context. It is not mhctools'
default C-terminal × anti-internal-cleavage composite. Missing evidence remains
missing and the policy's minimum-score gate excludes it.

```bash
python -m examples.osteosarc_test_data.compare_policies \
  --processing --output sid-policy-comparison
```

This needs installed MHCflurry presentation assets and Pepsickle. It uses the
pinned real NetMHCpan-4.2 cache; no licensed NetMHCpan executable is needed for
replay. Without `--processing`, it runs the frozen-baseline comparison only.
Original artifacts are never overwritten. Full evidence, native datasets,
per-policy decisions and SHA-256 manifests are written to the new directory.
The experiment retains the custom reference beside its native mutation records.
Keep that directory in place when reloading native constructs: reference-path
relocation is tracked in [#563](https://github.com/openvax/vaxrank/issues/563).
The typed policy evidence itself can be moved and replayed independently.

## Serum susceptibility and self-content windows

The initial human panel is deliberately narrow:

| Enzyme/model | Evidence used | Interpretation |
| --- | --- | --- |
| DPP4 / `dpp4-qpisa` | Published substrate-depletion model for the exposed N-terminal triplet | Flag bond 2 at predicted log2 depletion ≥ 1; the threshold is a heuristic |
| FAP / `fap-endo-gp`, `fap-dipeptidyl` | Reviewed recognition rules for internal Gly-Pro and exposed N-terminal X-Pro | Motif evidence, without a rate estimate |
| ACE / `ace-dipeptidyl` | Ordinary free-C-terminal dipeptide-removal rule | Does not cover all ACE routes or substrates |
| CPN1 / `cpn-basic` | Free-C-terminal Lys/Arg removal | Sequence context affects rates |

DPP4-mediated degradation is observed in human blood specimens, with different
behavior across peptides and handling conditions
([Yi et al., 2015](https://doi.org/10.1371/journal.pone.0134427)). The qPISA
coefficients describe an in-vitro substrate-depletion assay, not circulating
peptide half-life ([Gudipati et al., 2024](https://doi.org/10.1038/s44320-024-00071-4)).
Soluble FAP activity was isolated from human plasma
([Lee et al., 2006](https://pubmed.ncbi.nlm.nih.gov/16223769/)); its cleavage
preferences depend on more than a two-residue motif
([Lee et al., 2009](https://pubmed.ncbi.nlm.nih.gov/19402713/)). ACE and CPN
contributions were measured for bradykinin in human plasma
([Kuoppala et al., 2000](https://pubmed.ncbi.nlm.nih.gov/10749699/)). These are
reasons to screen susceptibility, not evidence that every matching vaccine
peptide will be rapidly degraded.

CPB2/TAFI requires an activation assumption and is excluded from the default
panel. MME, ANPEP, ENPEP and XPNPEP2 are candidates for separate extracellular
or tissue-exposure scenarios. Cytosolic peptidases and constitutively active
coagulation/fibrinolysis assumptions do not belong in this serum panel. PREP's
plasma activity is also less straightforward than a generic post-proline motif
would suggest ([Lee et al., 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC4711262/)).

For each candidate window let T be its unique target sequence/allele score,
R the score whose every copy has a supported internal cut, and S its unique
non-CTA exact-self sequence/allele score. The experimental utility is:

```
utility = (T - serum_weight * R) / (1 + self_weight * S)
```

First retain windows with at least `min_target_fraction` (default 0.95) of the
best available T; then maximize utility. Exact ties prefer the configured
length and existing tie-breaks. This permits an explicit target/self tradeoff;
it does not guarantee preservation of every target or every HLA allele.
The shared `self-trim-v1` formula substitutes `window_epitope_score` for the
legacy target score in the RNA-weighted combined score.

A cut must be strictly inside a target epitope to damage that copy. Cuts outside
it can release an intact epitope and are not penalized. A protected repeated
copy preserves its sequence/allele content. These are first-cut assessments;
the algorithm does not assume a chain of subsequent cleavages. Enzymes are not
treated as independent events, and the score is not a survival probability.

Serum refinement runs independently during peptide assembly, before allele
coverage selection, only for one-SLP-per-construct products. It examines shifts
and trims within the available ranked antigen context, using the emitted free,
acetylated or amidated termini. An unassessed requested enzyme excludes that
window from this experimental selection; it is not treated as resistance.
`window_selection.json` retains rejected alternatives even if nothing is
emitted. Selected constructs include target-cut details and enzyme provenance
in their manifest. The mRNA branch retains its original inputs.

The self term uses the existing `occurs_in_non_CTA_reference` annotation. CTA
classification remains with oncoref, and a peptide shared with a non-CTA gene
still counts. Missing self provenance is reported; an absent match in a partial
reference is not proof of tumor specificity. This heuristic does not replace
normal-tissue expression or TCR cross-reactivity assessment, nor close the
broader [self-risk orchestration issue #311](https://github.com/openvax/vaxrank/issues/311).

The Sid experiment is a small behavioral check, not a clinical validation.
Its pinned reference is a transcript subset; its prediction HLA panel differs
from parts of the reported typing. The report retains those discrepancies,
RNA rejection, missing cleavage context and disagreements with historical
selections. See the [checked-in comparison](../examples/osteosarc_test_data/policy_comparison/README.md)
for measured results.
