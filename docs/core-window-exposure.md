# Experimental window exposure

**Experimental:** explore delivery assumptions with the MHCtools exposure model.
Results do not establish clinical safety or efficacy and do not change vaccine
ranking or products.

## Run

Start with an enriched `--output-epitopes` file containing native antigen contexts
and saved per-allele scores. The command evaluates all admitted, target-containing
SLP windows in the length range, before ranking or count limits.

Save this minimal `profile.json`; kinetics remain unassessed until supplied:

```json
{
  "scenario": "IM SLP comparison; kinetics not supplied",
  "times_hours": [0, 1, 2, 4, 8, 24, 48],
  "administration": "IM"
}
```

```sh
vaxrank-window-exposure --input-epitopes native.tsv --profile profile.json \
  --min-length 20 --max-length 30 --n-term free --c-term free \
  --output-dir window-exposure
```

Outputs: `request.json`, `exposure.json`, `cores.csv` and `report.md`.
The request also works with `mhctools vaccine-windows`.

## Interpret

- Free linear peptides default to **1 µg per window**. Supply the actual terminal
  chemistry; modified peptides remain unassessed by the current engine.
- Missing context, scores or admission appear as coverage gaps. No new prediction
  or source reconstruction occurs; provenance and admission records are retained.
- Supply explicit kinetic rates and cleavage assessments under `kinetics`.
  Per-window overrides use `candidates[window_id]`. Binding scores do not become
  kinetic rates; synthetic sensitivity scenarios are illustrative.
- The tracked interval is the saved scored ligand. A shorter documented core
  uses `core_intervals` with prediction ID, allele, source and zero-based,
  half-open coordinates **within that ligand**. Affinity alone does not identify
  a shorter core. Original ligand scores and provenance remain recorded.
- APC arrival, blood persistence and pMHC exposure are different outputs. Core
  marginals cannot be summed as unique dose fractions; weighted means require
  numerical assessments for every contributing core.

For rates, units and sensitivity examples, see the
[MHCtools window guide](https://openvax.github.io/mhctools/vaccine-windows/).
For mRNA, use the [full-construct model](https://openvax.github.io/mhctools/vaccine-core-exposure/)
with the actual encoded construct, transcript dose and delivery assumptions.
The SLP command does not infer an RNA sequence from a protein window.
