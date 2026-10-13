# Analysis outputs and figures

Analysis outputs are requested with their own flags, with or without construct
files.

## Analysis reports

Tables of ranked candidates, predicted epitopes and applicable
manufacturability metrics.

| Flag | Output |
|---|---|
| `--output-ascii-report` | Plain-text summary |
| `--output-html-report` | HTML report |
| `--output-pdf-report` | PDF report (wkhtmltopdf or WeasyPrint backend) |
| `--output-xlsx-report` | Excel workbook with one sheet per variant |
| `--output-csv` | Flat CSV |
| `--output-json-file` | Full ranked-vaccine-peptides graph as JSON |

## Neoepitope report

A per-(peptide, allele) table in XLSX or CSV, chosen by file extension. It is
the default output of the LENS / pVACseq input path and is also available on
the full pipeline.

| Flag | Output |
|---|---|
| `--output-neoepitope-report` | XLSX (default) or CSV (by extension) |

## Mutation and assembly figures

`vaxrank-mutation-figure` turns the Isovar table produced by
`--output-isovar-csv` into publication-ready, white-background transcript and
protein context figures. When the input includes the optional `*_cdna_sequence`
and `*_cdna_variant_{start,end}` columns, each variant also gets a
5-prime-to-3-prime nucleotide comparison. Both views separate the reference
transcript, the annotation-only mutant prediction and the RNA-assembled
sequence, so unsupported or assembly-refined claims stay explicit.

The `random` predictor below is for testing the figure workflow only.

```bash
vaxrank --vcf variants.vcf --bam tumor-rna.bam \
  --mhc-predictor random --mhc-alleles 'HLA-A*02:01' \
  --output-isovar-csv patient-isovar.csv

vaxrank-mutation-figure patient-isovar.csv \
  --output-root figures \
  --variant DYNC1H1 \
  --format svg --format pdf --format png \
  --png-scale 3
```

Outputs are grouped under `figures/YYYY-MM-DDTHHMMSSZ/`. Every variant gets an
SVG, a PDF, a high-resolution PNG (3600×2100 at the default 3× scale) and a
compact JSON record. The figure and JSON distinguish annotation
gene/transcript provenance from RNA-assembly provenance. The run manifest
records the input checksum, Vaxrank version, formats, raster dimensions and
file inventory. Omit `--variant` to render every row. Use `--timestamp` for
reproducible builds or tests; an existing run is never overwritten.

## Evidence figures for structural variants

For structural variants and fusion junctions, use a cross-platform evidence
JSON rather than an Isovar table:

```bash
vaxrank-evidence-figure evidence.json \
  --output-root evidence-figures \
  --format svg --format pdf --format png
```

This view keeps DNA, short-read RNA, short-read assembly and long-read RNA as
separate evidence stages. It displays the junction nucleotide sequence and
withholds a protein unless the input supplies an established coding sequence.
See `examples/osteosarc_evidence_figures/` for DNA-only, cross-platform,
matched-sample long-read rescue, and unresolved long-read-only examples.

## Complex-result figures

To show the downstream Vaxrank decision rather than the evidence alone, use a
complex-result record. This view separates DNA interpretation, RNA support,
translated protein, target epitopes and the final vaccine-construct selection:

```bash
vaxrank-complex-result-figure results.json \
  --output-root result-figures \
  --combined-output vaxrank-complex-variant-results.pdf
```

The timestamped run contains a compact combined PDF and one high-resolution
PNG/SVG page per result. See `examples/osteosarc_complex_results/` for a real
assembly-dependent frameshift selection, an assembled indel with no passing
binder, an ambiguous fusion held out of ranking, and an unresolved long-read
fusion whose protein is deliberately withheld.
