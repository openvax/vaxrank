# Core delivery across candidate SLP windows

`vaxrank-window-exposure` connects saved epitope scores and antigen contexts to
MHCtools' core-preserving exposure model. It enumerates **all admitted,
target-containing windows** in the requested length range before ordinary
best-score and construct-count pruning. This opt-in analysis does not change the
existing vaccine ranking or products.

Start with an enriched file from `vaxrank --output-epitopes native.tsv`. It must
contain explicit antigen contexts and saved per-allele scores. Missing context,
missing scores and antigen admission exclusions appear as coverage gaps; the
command does not reconstruct unrecorded source context or run new predictors.
Combined observations use the existing Topiary representative policy.

## Enumerate and assess

A minimal profile deliberately leaves kinetics unassessed:

```json
{
  "scenario": "IM SLP window comparison; kinetics not yet supplied",
  "times_hours": [0, 1, 2, 4, 8, 24, 48],
  "administration": "IM"
}
```

```sh
vaxrank-window-exposure --input-epitopes native.tsv --profile profile.json \
  --min-length 20 --max-length 30 --n-term free --c-term free \
  --output-dir window-exposure
```

The new directory contains `request.json`, `exposure.json`, `cores.csv` and
`report.md`. The request is directly usable by `mhctools vaccine-windows`.
Choose the actual terminal chemistry; amidated, acetylated or otherwise modified
inputs remain explicit and unassessed by the current sequence-only trajectory
engine. They are never silently replaced with free peptide.

Free linear SLPs default to 1 µg **of each actual window**. Each contributing
scored occurrence/allele keeps its original score, coordinates and source
provenance. Gene, transcript, protein, species and input-scope metadata accompany
every window, together with the native antigen kind, targetability and admission
record. Repeated sequences at distinct positions remain distinct targets.

## Supply kinetics or run sensitivity scenarios

`kinetics` in the profile maps MHC class (`I` or `II`) to the mhctools trajectory
rates, length bound and exact fragment cleavage assessments. To use different
kinetics or a different dose for a particular window, put those fields under
`candidates[window_id]`; IDs are listed in the initial request. Unknown window
IDs are rejected. Scores, ranks and protease motifs are not converted to rates.

For a complete explicitly synthetic sensitivity workflow, run the example from
a mhctools source checkout:

```sh
python examples/vaccine_sensitivity_scenarios.py \
  --input /path/to/window-exposure/request.json --output-dir sensitivity
```

This replaces the profiles with declared synthetic transport, processing and
production assumptions, then varies all unproductive loss/clearance rates
together. It is an illustrative scenario range, not fitted human physiology.
For all parameters and their units, see the
[mhctools window guide](https://openvax.github.io/mhctools/vaccine-windows/).

## Scored ligands and binding cores

By default, the tracked interval is the **saved scored ligand**, not the full
SLP. When a shorter binding core is known, supply it explicitly in the profile:

```json
{
  "core_intervals": [
    {
      "prediction_id": "the-original-scored-occurrence-id",
      "allele": "HLA-DRB1*04:01",
      "start": 3,
      "end": 12,
      "source": "Identifier of the binding-core annotation or assay"
    }
  ]
}
```

These coordinates are zero-based, half-open **within the scored ligand**.
They move correctly as the containing SLP window moves. The original ligand and
its score remain recorded alongside the supplied core and source. When different scored ligands document the same core occurrence and allele,
the adapter combines their existing score contributions and retains every ligand,
score and annotation source as a component. It runs that physical core once.
Repeated evidence for the same original ligand must be resolved by the existing
representative policy before assessment. An unknown
shorter MHC-II core is not inferred from a 15-mer's affinity alone. Delivery of
a supplied core does not establish that its isolated peptide is a loadable
ligand; class-II loading still uses the explicitly supplied precursor-length
bound and kinetics.

The report distinguishes nodal APC arrival, blood persistence and pMHC exposure.
An IM blood curve is not an intrinsic circulation elimination half-life. Separate
core marginals are not summed as a unique-molecule dose fraction. A weighted mean
is available only when all contributing cores have numerical assessments.

## mRNA constructs

Use the same mhctools batch schema for **actual full encoded constructs**, with
translated-core coordinates, delivery type, explicit signal boundary when
applicable, and a transcript dose. The SLP enumerator does not invent an RNA
sequence or its molecular weight from a protein window. mhctools supports direct
APC expression/processing, secretion and uptake, and MHC-I cross-dressing; see
[mechanisms and evidence](https://openvax.github.io/mhctools/vaccine-core-exposure/).
