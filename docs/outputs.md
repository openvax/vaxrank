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
| `--output-xlsx-report` | Ranked candidate windows: All sheet and one sheet per source |
| `--output-csv` | The same ranked candidate rows and columns as the workbook |
| `--output-json-file` | Explicit legacy ranked-object JSON export |
| `--output-policy-dir` | Optional resolved configuration and policy decision/evidence tables |

CSV/XLSX candidate rows use one schema for VCF/BAM, external and native inputs.
`source_rank` orders the ranked source groups; `window_rank` orders candidate
windows within each group. Source kind, identifiers, genes, transcripts, species
and zero-based, end-exclusive targetable intervals accompany each sequence.
Mutation coordinates and expression scores are blank when no mutation fragment
is available. Scores retain their numeric precision. These are candidate windows;
manifests and FASTAs describe the emitted products. Empty results still write
headers and an empty `All` worksheet.

With `--output-dir`, the default bundle is `ranked_candidates.csv`,
`candidate_predictions.tsv`, `vaccine_report.txt`, `vaccine_report.pdf` and the
modality construct files. Native evidence retains source provenance, configuration
and policy decisions. `--output-policy-dir DIR` additionally writes
`DIR/selection_policy.json` and `DIR/policy_evidence/` for inspecting individual
criteria. Removing those optional exports does not prevent native replay.
Checksum-verified `candidate_predictions.tsv.references/` resources, when present,
are required for portable annotation/reference access; keep them with the dataset.

Migration from earlier output bundles:

- `--output-csv` now always exports ranked windows. To export prediction evidence
  from external/native input, use `--output-neoepitope-report predictions.csv`.
- `source_rank` / `window_rank` replace the mutation-specific `variant_rank` /
  `peptide_rank` column names. The workbook no longer repeats rank in an unnamed column.
- Default runs no longer write `input_provenance.json`, a sibling `.policy/`
  directory or legacy ranked JSON. Explicit `--output-json-file` still works.
- Legacy `--input-json-file` supports reports/products but cannot reconstruct a
  complete native prediction dataset. It does not auto-export native evidence;
  an explicit `--output-epitopes` request with that input raises an explanatory error.

## Joint protein collections

The Python API exports retained source proteins, evidence and Vaxrank's
selection state to ORFdata without running prediction or ranking again:

```python
from orfdata import ReferenceContext
from vaxrank import to_proteome

reference = ReferenceContext(9606, "human", "GRCh38", "Ensembl", "115")
patient = to_proteome(
    results,
    sample_id="tumor-1",
    source="vaxrank-run-1",
    reference=reference,
)
joint = reference_transcriptome.proteome(nc_orf_sources="all").merge(patient)
joint.write_fasta("joint.fasta")
joint.save("joint.proteome.json.gz")
```

To convert RNA-backed results again without the original reads, save Isovar's
portable hypotheses with the same sample and run IDs:

```python
from pathlib import Path
from vaxrank.vaxrank_results import VaxrankResults

Path("results.json").write_text(
    results.to_json(sample_id="tumor-1", source="vaxrank-run-1")
)
restored = VaxrankResults.from_json(Path("results.json").read_text())
replayed = to_proteome(
    restored, sample_id="tumor-1", source="vaxrank-run-1", reference=reference
)
assert replayed.identity == patient.identity
```

This JSON retains the source hypotheses and evidence for ORFdata replay. It
does not reconstruct live Isovar read objects.

Every retained Isovar alternative comes through Isovar's own adapter. Vaxrank
adds an export path and its recorded fragment/selection state without creating
a second source protein. DNA-predicted and supplied fragments keep their exact
window and mutation intervals; they are not presented as reconstructed
full-length proteins. Unknown source references remain unknown and require
ORFdata's explicit compatibility policy when composed with a known host.

## Reviewing a vaccine design

With `--output-dir`, one `vaccine_report.txt` and `vaccine_report.pdf` summarize
all emitted modalities. The product section shows exact sequences, lengths,
antigen placements and source mappings; the candidate table links ranked source
windows to the products that use them. Historical files without source mappings
say that the information is unavailable. Successful empty designs explicitly
report zero products.

`--report-candidate-details` puts annotated source windows and paired epitope/WT
comparisons first, after a short patient and predictor overview. Exact emitted
products and raw run settings follow as linked appendices. Without that flag,
the report retains its compact, annotation-free ranking view. Candidate evidence
and manufactured product coverage remain separately labeled.

Target highlighting follows recorded intervals (including disjoint regions and
deletion junctions). WT differences are marked for equal-length recorded pairs;
unequal-length pairs are not given an invented alignment. WT is a comparator,
not evidence of an observed patient-normal protein. Product highlights require
unambiguous placement and agreement among the linked source masks.

The predictor overview describes evidence actually available in the reported
candidates, including imported results. Missing MHCflurry results are explicit;
disabled Pepsickle is distinguished from missing or partial cleavage evidence.
NetMHCpan presentation is labeled **Score_EL**, while MHCflurry presentation and
processing retain their model names. Presentation scores run higher-is-stronger;
presentation percentile ranks run lower-is-stronger. These model outputs are
separate from Vaxrank's configured ranking score and from Pepsickle cleavage
annotations. The report's score guide defines each quantity.

The per-row WT affinity must match the target's allele, predictor and version;
a missing matching value is reported as unavailable. Gene/antigen names lead
product headings, with stable product IDs retained for manifest/order-form lookup.
Report-only commands without `--output-dir` show candidate results and state
that no construct files were emitted. HTML/PDF styling is bundled locally.

## Neoepitope report

A per-(peptide, allele) table in XLSX or CSV, chosen by file extension. Request
it explicitly for VCF/BAM, external or native input; it supplements the ranked
candidate table and native prediction dataset.

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
