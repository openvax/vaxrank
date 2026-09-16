# Osteosarcoma assembly figures

These figures use the repository's small, checksum-pinned RNA fixtures to show
why variant annotation is not a substitute for RNA assembly around a mutation.
They are software regression examples, not clinical recommendations.

The six curated cases cover distinct assembly outcomes:

- **DYNC1H1 p.Val314Ile:** short-read RNA assembly confirms the annotation-only
  protein sequence.
- **EXOC4 p.Ser34Ile:** short-read RNA assembly confirms a second substitution
  in a different transcript context.
- **H1-2 p.Ala197_Lys201del:** long-read RNA assembly refines the junction in a
  lysine/alanine/proline-rich repeat where sequence placement is ambiguous.
- **GTF3C5 p.Glu503_Glu506del:** long-read RNA resolves another low-complexity
  deletion context and retains its assembled transcript provenance.
- **MAP2 p.Leu867fs:** annotation predicts a novel frameshift tail, but the
  selected RNA fixture has no alternate fragments, so Vaxrank withholds an
  RNA-assembled protein sequence rather than presenting the prediction as
  RNA-supported.
- **PIP5K1A p.Gly474fs:** a second annotation-predicted frameshift is withheld
  when the selected RNA fixture has no alternate fragments.

Rebuild the source CSV and a new UTC-stamped SVG/PDF/PNG run from the
repository root:

```bash
python examples/osteosarc_mutation_figures/generate.py
```

For a byte-for-byte organized run name, pass an explicit timestamp:

```bash
python examples/osteosarc_mutation_figures/generate.py \
  --timestamp 2026-09-15T230000Z
```

Each run contains one directory per variant, with a white-background SVG, PDF,
3600×2100 PNG, and the compact input record used for the figure. `manifest.json`
records the Vaxrank version, source CSV checksum, timestamp, raster dimensions,
and every generated file.
The source CSV is retained under `source/`; the timestamped outputs live under
`runs/`.

The website uses VCF-style anchored deletion coordinates. Isovar removes the
shared anchor before displaying its HGVS-like genomic deletion, so the H1-2 and
MAP2 positions in these figures begin one base later while representing the
same alleles linked in each footer.
