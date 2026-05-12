[![Tests](https://github.com/openvax/vaxrank/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/vaxrank/actions/workflows/tests.yml)
[![Coverage Status](https://coveralls.io/repos/github/openvax/vaxrank/badge.svg?branch=master)](https://coveralls.io/github/openvax/vaxrank?branch=master)
[![Docs](https://github.com/openvax/vaxrank/actions/workflows/docs.yml/badge.svg)](https://github.com/openvax/vaxrank/actions/workflows/docs.yml)
[![GitHub Pages](https://github.com/openvax/vaxrank/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/openvax/vaxrank/actions/workflows/pages/pages-build-deployment)
<a href="https://pypi.python.org/pypi/vaxrank/">
    <img src="https://img.shields.io/pypi/v/vaxrank.svg?maxAge=1000" alt="PyPI" />
</a>

# vaxrank

Vaxrank is the neoantigen ranking component of the
[OpenVax](https://www.openvax.org/) pipeline for designing personalized
cancer vaccines.  Given either (a) a patient's somatic mutations + tumor
RNA-seq + HLA type, or (b) a pre-computed neoepitope report from
[LENS](https://github.com/openvax/lens) or pVACseq, Vaxrank selects and
ranks the mutant antigens most likely to elicit a T-cell response and
emits them as the vaccine type(s) the user requests — peptide pools,
mRNA constructs, or analysis reports for review.

## Contents

- [Features](#features)
- [Quick Start](#quick-start)
- [Overview](#overview)
- [Vaccine designs](#vaccine-designs) — 8 design points across two modalities
- [Vaccine types and output modes](#vaccine-types-and-output-modes) —
  reports, peptide constructs, mRNA constructs, external inputs
- [Clinical Use](#clinical-use)
- [Installation](#installation)
- [Configuration](#configuration) — YAML config + topiary DSL
- [MHC Binding Predictors](#mhc-binding-predictors)
- [How It Works](#how-it-works) — pipeline internals + data model
- [Papers & Citations](#papers--citations)
- [Dependencies](#dependencies)
- [Legacy flags](#legacy-flags)
- [Development](#development)

## Features

- **Two input paths.** Full pipeline from a tumor VCF + RNA-seq BAM,
  or drive ranking from a pre-computed neoepitope report produced by
  [LENS](https://github.com/openvax/lens) or pVACseq via
  `--input-lens` / `--input-pvacseq`.
- **Multiple vaccine types from one run.** `--vaccine-type` is
  multi-valued: pass `peptide`, `mrna`, or both. The same ranked
  candidates render to peptide-pool FASTAs, mRNA constructs, or both
  in parallel.
- **Eight vaccine designs across two modalities.** Two orthogonal
  axes — `--antigen-content` (mutation-spanning SLP vs minimal MHC
  ligand) × `--antigens-per-construct` (1 vs N) — give four designs
  per modality, including canonical SLPs (PGV-001) and BioNTech-style
  multi-antigen mRNA (FixVac / iNeST).
- **Analysis reports** in CSV, XLSX, ASCII, HTML, PDF, and JSON, plus
  a per-(peptide, allele) neoepitope report. Reports run independently
  of vaccine-type dispatch.
- **MHC predictor integration** via [mhctools](https://github.com/openvax/mhctools):
  MHCflurry, NetMHCpan family, BigMHC, NetMHCIIpan, Pepsickle, and IEDB
  web predictors.
- **Shared linker library** with primary-source citations and a
  compositional grammar — `(G4S)2`, `AAY`, 2A self-cleaving peptides,
  furin motifs — usable interchangeably across peptide and mRNA
  designs.

## Quick Start

Full pipeline from a tumor VCF + RNA-seq BAM, emitting ranked vaccine
peptides as text, HTML, and PDF reports:

```sh
vaxrank \
    --vcf tests/data/b16.f10/b16.vcf \
    --bam tests/data/b16.f10/b16.combined.bam \
    --mhc-predictor netmhc \
    --mhc-alleles H2-Kb,H2-Db \
    --output-ascii-report vaccine-peptides.txt \
    --output-pdf-report vaccine-peptides.pdf
```

Required inputs:
- `--vcf` — somatic variants (VCF from any variant caller)
- `--bam` — tumor RNA-seq alignments (used by Isovar to assemble mutant transcripts)
- `--mhc-alleles` — patient HLA alleles (e.g. `HLA-A*02:01,HLA-B*07:02`)
- `--mhc-predictor` — which MHC binding predictor to use (see
  [MHC Binding Predictors](#mhc-binding-predictors))

Drive design from a pre-computed neoepitope report (LENS or pVACseq)
when upstream MHC prediction has already been done:

```sh
vaxrank --input-lens patient.lens.tsv \
        --vaccine-type mrna --output-dir mrna_out/ \
        --ensembl-release 102
```

Emit both peptide and mRNA constructs in one run — outputs land in
per-modality subdirs:

```sh
vaxrank --vcf v.vcf --bam r.bam \
        --vaccine-type peptide mrna --output-dir vaccines/
# → vaccines/peptide/{vaccine.fasta, manifest.json, order_form.csv}
# → vaccines/mrna/{cds.fasta, no_polyA.fasta, full.fasta, manifest.json, layers.csv}
```

## Overview

Personalized cancer vaccines (also called neoantigen vaccines) work by
training the immune system to recognise peptides that arise from somatic
mutations unique to a patient's tumor.  Designing such a vaccine requires
a computational pipeline that bridges raw sequencing data and the
peptide synthesiser:

1. **Variant calling** — Whole-exome or whole-genome sequencing of the
   tumor and matched normal identifies somatic mutations.  This is
   typically done with tools such as MuTect or Strelka, upstream of
   Vaxrank.
2. **Mutant transcript assembly** — Tumor RNA-seq reads overlapping each
   mutation are assembled by [Isovar](https://github.com/openvax/isovar)
   to determine the true mutant protein sequence.  This step phases
   nearby germline variants and captures any mutation-associated splicing
   differences, producing a more accurate reading frame than DNA-only
   prediction.
3. **MHC binding prediction** — Candidate epitopes (short peptide
   subsequences spanning the mutation) are scored for predicted binding
   to the patient's HLA class I molecules using
   [mhctools](https://github.com/openvax/mhctools), which wraps
   predictors such as MHCflurry, NetMHCpan, and BigMHC.
4. **Vaccine peptide selection** — Vaxrank assembles longer synthetic long
   peptides (SLPs, typically 25-mers) around the mutation, scores them by
   the number and strength of their predicted MHC-binding epitopes,
   filters out peptides that appear in the reference proteome, annotates
   known cancer hotspot mutations, and ranks candidates by a combined
   immunogenicity and manufacturability score.
5. **Vaccine-type dispatch** — the ranked candidates are written out as
   one or more of the vaccine types selected via `--vaccine-type`: a
   peptide pool ready for synthesis, an mRNA construct ready for IVT, or
   both. Analysis reports are emitted independently. Steps 1-3 are
   skipped when an external neoepitope report is supplied via
   `--input-lens` or `--input-pvacseq`; the ranking and dispatch steps
   are identical.

## Vaccine designs

Vaxrank's vaccine design space is two **orthogonal axes** (shared
across vaccine types) plus the type itself:

| Axis | Values | What it controls |
|---|---|---|
| `--vaccine-type` | `peptide` / `mrna` (multi-valued) | The platform(s); pass multiple for parallel design |
| `--antigen-content` | `mutation_spanning` / `minimal_epitope` | What each antigen *is* |
| `--antigens-per-construct` | `1` / `N` | How many antigens to concatenate per construct |

Combined, the matrix yields 8 distinct designs — 4 per vaccine type:

| Type | Content | Per-construct | Design name | Reference |
|---|---|---|---|---|
| peptide | mutation_spanning | 1 | **SLP** (default) | [PGV-001 (Saxena 2025)](https://pubmed.ncbi.nlm.nih.gov/40094414/) |
| peptide | mutation_spanning | N | Multi-SLP / multi-epitope long peptide | |
| peptide | minimal_epitope | 1 | Minimal-ligand peptide | |
| peptide | minimal_epitope | N | Concatenated minimal-ligand peptide | |
| mrna | mutation_spanning | N | **BioNTech FixVac / iNeST** (default for mRNA) | [Sahin 2017](https://doi.org/10.1038/nature23003) / [Rojas 2023](https://doi.org/10.1038/s41586-023-06063-y) |
| mrna | mutation_spanning | 1 | Single-antigen mRNA | |
| mrna | minimal_epitope | N | "String of beads" mRNA | [Whitton 1993](https://pubmed.ncbi.nlm.nih.gov/7677954/) |
| mrna | minimal_epitope | 1 | Single-ligand mRNA | |

A third knob, `--epitopes-per-antigen`, controls how many top MHC
ligands to take *per ranked vaccine peptide* when content is
`minimal_epitope`. The default `1` is the "single top ligand"
semantics; `>1` packs multiple top ligands from the same variant as
separate antigens.

### Peptide designs

**SLP (default).** Mutation-spanning long peptide, one antigen per
construct — the PGV-001 canonical design.

```sh
vaxrank --vcf v.vcf --bam r.bam --output-dir vaccine_out/
# → vaccine_out/{vaccine.fasta, manifest.json, order_form.csv}
```

**Multi-epitope concatenated peptide.** Several mutation-spanning
antigens linked into one longer peptide. Use `--peptide-linker` to
pick the spacer; `AAY` is the proteasome-friendly default.

```sh
vaxrank --vcf v.vcf --bam r.bam \
        --output-dir vaccine_out/ \
        --peptide-antigens-per-construct 5 --peptide-linker AAY
```

**Minimal-epitope peptide.** A single short MHC ligand per construct
— useful when minimum-length manufacturability matters more than
flanking context.

```sh
vaxrank --vcf v.vcf --bam r.bam \
        --output-dir vaccine_out/ \
        --antigen-content minimal_epitope
```

### mRNA designs

**BioNTech FixVac / iNeST canonical.** Multi-antigen mutation-spanning
mRNA — the default for `--vaccine-type mrna`. Antigens are linked with
`(G4S)2` and emitted as CDS, no-polyA, and full (with polyA) FASTAs
plus a structured manifest.

```sh
vaxrank --vcf v.vcf --bam r.bam --vaccine-type mrna --output-dir mrna_out/
# → mrna_out/{cds.fasta, no_polyA.fasta, full.fasta, manifest.json, layers.csv}
```

**String-of-beads mRNA.** Concatenated minimal-epitope antigens —
short MHC ligands linked together rather than mutation-spanning
windows.

```sh
vaxrank --vcf v.vcf --bam r.bam --vaccine-type mrna --output-dir out/ \
        --mrna-antigen-content minimal_epitope --mrna-antigens-per-construct 8 \
        --mrna-linker AAY
```

**Top-N ligands per variant in a string-of-beads.** Pack multiple top
MHC ligands from each ranked vaccine peptide as separate antigens.

```sh
vaxrank --vcf v.vcf --bam r.bam --vaccine-type mrna --output-dir out/ \
        --mrna-antigen-content minimal_epitope \
        --mrna-epitopes-per-antigen 2 --mrna-antigens-per-construct 16
```

### Both modalities in one run

Multi-valued `--vaccine-type` writes per-modality subdirs in
`--output-dir`.

```sh
vaxrank --vcf v.vcf --bam r.bam --vaccine-type peptide mrna --output-dir vaccines/
# → vaccines/peptide/{vaccine.fasta, manifest.json, order_form.csv}
# → vaccines/mrna/{cds.fasta, no_polyA.fasta, full.fasta, manifest.json, layers.csv}
```

## Vaccine types and output modes

Vaxrank always ranks. The vaccine-type writer fires only when both
`--vaccine-type` and `--output-dir` are set. `--vaccine-type` is
multi-valued (default `peptide`): pass one or more of `peptide` / `mrna`.
Single-mode runs write canonical files directly in `--output-dir`;
multi-mode runs scope into per-modality subdirs (`DIR/peptide/`,
`DIR/mrna/`, …). Analysis reports use their own `--output-*` flags
and are independent of the vaccine-type dispatch.

```sh
# Peptide pool (default vaccine type)
vaxrank --vcf v.vcf --bam r.bam --output-dir vaccine_out/

# mRNA construct
vaxrank --vcf v.vcf --bam r.bam --vaccine-type mrna --output-dir mrna_out/

# Both at once (per-modality subdirs in mixed_out/)
vaxrank --vcf v.vcf --bam r.bam --vaccine-type peptide mrna --output-dir mixed_out/

# Reports only (no vaccine constructs)
vaxrank --vcf v.vcf --bam r.bam --output-pdf-report report.pdf

# Drive vaccine design from a pre-computed LENS report
vaxrank --input-lens patient.lens.tsv --vaccine-type mrna \
        --output-dir mrna_out/ \
        --ensembl-release 102

# Full ASCII summary report from a LENS file (transcripts resolved)
vaxrank --input-lens patient.lens.tsv --output-ascii-report report.txt \
        --ensembl-release 102
```

### Analysis reports

Per-variant tables of ranked vaccine peptide candidates, predicted
epitopes, and manufacturability scores. Independent of vaccine-type
dispatch — runs whenever any report flag is set.

| Flag | Output |
|---|---|
| `--output-ascii-report` | Plain-text summary |
| `--output-html-report` | HTML report |
| `--output-pdf-report` | PDF report (wkhtmltopdf or WeasyPrint backend) |
| `--output-xlsx-report` | Excel workbook with one sheet per variant |
| `--output-csv` | Flat CSV |
| `--output-json-file` | Full ranked-vaccine-peptides graph as JSON |

### Neoepitope report

Per-(peptide, allele) report (XLSX/CSV). Default output of the
LENS / pVACseq input path; also available on the full pipeline.

| Flag | Output |
|---|---|
| `--output-neoepitope-report` | XLSX (default) or CSV (by extension) |

### Peptide constructs

`vaccine.fasta` + `manifest.json` + `order_form.csv` written into
`--output-dir` (or `--output-dir/peptide/` in multi-mode). The
peptide design comes from `--antigen-content` and
`--peptide-antigens-per-construct` (see [Vaccine designs](#vaccine-designs)
above).

| Flag | Purpose |
|---|---|
| `--output-dir` | Where to write the construct files |
| `--peptide-linker` | Inter-antigen spacer (e.g. `AAY`, `(G4S)2`); default `G4S3` |
| `--peptide-max-antigen-length-aa` | Truncate antigens longer than this |
| `--peptide-n-terminal-acetyl` | Add N-terminal acetylation note to the manifest |
| `--peptide-c-terminal-amide` | Add C-terminal amide note to the manifest |

### mRNA constructs

A directory containing three FASTAs (`cds.fasta`, `no_polyA.fasta`,
`full.fasta`) plus `manifest.json` (per-element view) and `layers.csv`
(long-format per-element table with AA + nt). Codon optimization uses
[DnaChisel](https://github.com/Edinburgh-Genome-Foundry/DnaChisel); 2A
self-cleaving peptides preserve their published codon usage
automatically.

Flags are grouped by what they configure:

**Construct anatomy**

| Flag | Purpose |
|---|---|
| `--output-dir` | Where to write construct files (or `--output-dir/mrna/` in multi-mode) |
| `--mrna-signal-peptide` | Leader peptide: `HLA-A`, `HLA-B`, `tPA`, `IgK`, `CD8A`, `CD28` |
| `--mrna-5p-utr` | 5' UTR (e.g. `HBB`, `HBB_FI` tandem) |
| `--mrna-3p-utr` | 3' UTR |
| `--mrna-include-mitd` / `--mrna-no-mitd` | Include the BioNTech MITD trafficking domain |
| `--mrna-mitd` | Which MITD variant (`HLA-A` / `HLA-B`) |
| `--mrna-max-length-nt` | Hard cap on construct length (nt) |
| `--mrna-antigens-per-construct` | Antigens per CDS |
| `--mrna-max-constructs` | Stop emitting after this many constructs |

**PolyA tail**

| Flag | Purpose |
|---|---|
| `--mrna-poly-a-length` | Length of polyA tail (default `120`) |
| `--mrna-poly-a-segmented` | Use BNT162b2-style segmented pattern (A30 + linker + A70) |
| `--mrna-poly-a-first-segment` | Length of the first segment when segmented |
| `--mrna-poly-a-segment-linker` | Inter-segment linker sequence |

**Linker optimization**

Per-junction MHC-aware linker swap minimizes predicted presentation
of chimeric k-mers spanning antigen junctions.

| Flag | Purpose |
|---|---|
| `--mrna-linker` | Default inter-antigen spacer (e.g. `(G4S)2`) |
| `--mrna-optimize-linkers` / `--mrna-no-optimize-linkers` | Per-junction MHC-aware swap (on by default) |
| `--mrna-junction-candidates` | Candidate linkers considered at each junction |
| `--mrna-junction-rank-strong` | Strong-binder %-rank threshold |
| `--mrna-junction-rank-mild` | Mild-binder %-rank threshold |

**Codon optimization**

| Flag | Purpose |
|---|---|
| `--mrna-codon-species` | Target organism for codon usage (default `h_sapiens`) |
| `--mrna-codon-method` | DnaChisel optimization strategy |
| `--mrna-csv-no-full-rows` | Skip the `full` polyA'd rows in `layers.csv` (saves disk space) |

### External-input mode

Drive vaccine design from a pre-computed neoepitope report instead of
VCF + BAM. Same downstream dispatch — peptide and mRNA construct
outputs work identically.

| Flag | Input format |
|---|---|
| `--input-lens` | LENS report TSV |
| `--input-pvacseq` | pVACseq aggregated TSV (`*all_epitopes.aggregated.tsv`) |

### Manifest schema

The peptide and mRNA construct JSON manifests share a back-compat
schema (`modality`, `name`, `length`, `length_unit`, `antigen_names`,
`components`, `manufacturability`). The mRNA manifest additionally
exposes `cds`, `no_polya_nt`, `full_nt`, per-antigen `antigens` (each
with AA + nt), and a structured `elements` dict with one entry per
layer (5' UTR, signal peptide, antigens, linkers per junction, MITD,
stop codon, 3' UTR, polyA) — every layer carrying both AA (where
applicable) and nt forms for direct inspection.

### Shared linker library and grammar

Both vaccine types consume the same set of linker names so a single
construct design can be ported between peptide and mRNA backbones.

**Static entries:**

| Name | Type | Use |
|---|---|---|
| `G2S`, `G3S`, `G4S`, `G5S` | flexible (Gly_n_Ser) | The (Gly4Ser)n family ([Huston *PNAS* 1988](https://doi.org/10.1073/pnas.85.16.5879)); used clinically in BioNTech FixVac / iNeST as `(G4S)2` |
| `EAAAK` | rigid α-helical | When fused antigens need separation rather than flex ([Arai *Protein Eng* 2001](https://doi.org/10.1093/protein/14.8.529)) |
| `RKRR`, `RVKR`, `RKRKR` | furin cleavage | R-X-(K/R)-R motif ([Hosaka *J Biol Chem* 1991](https://pubmed.ncbi.nlm.nih.gov/1905715/)); preclinical in DNA vaccines, no clinical vaccine use as of 2025 |
| `AAY` | proteasome-friendly | Empirical foundation: [Livingston *Vaccine* 2001](https://doi.org/10.1016/S0264-410X(01)00233-X); see citation in `vaxrank/vaccine_library.py` for the AAY-vs-GGGS empirical landscape ([Yang 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4514284/) vs [Aguilar-Gurrieri 2023](https://doi.org/10.1007/s00262-023-03409-3)) |
| `AAA` | alanine spacer | [Aguilar-Gurrieri *Cancer Immunol Immunother* 2023](https://doi.org/10.1007/s00262-023-03409-3) — strongest empirical alanine spacer for MHC-I presentation |
| `GPGPG` | helper-T spacer | Between MHC-II epitopes ([Livingston *J Immunol* 2002](https://doi.org/10.4049/jimmunol.168.11.5499)) |
| `P2A`, `T2A`, `F2A`, `E2A` | self-cleaving 2A | Co-translational ribosomal skipping for mRNA constructs ([Donnelly *J Gen Virol* 2001](https://doi.org/10.1099/0022-1317-82-5-1027); [Kim *PLoS ONE* 2011](https://doi.org/10.1371/journal.pone.0018556)). In peptide mode these are functionally inert and the manifest annotates them as such. |

**Compositional grammar** (parsed at lookup time):

| Form | Meaning | Example |
|---|---|---|
| `(BASE)N` / `(BASE)xN` / `BASExN` | Repeat N times | `(G4S)2` → `GGGGSGGGGS`, `G4Sx2` → same |
| `GnSm` | Literal n glycines + m serines (single unit, **not** a repeat) | `G6S` → `GGGGGGS`, `G4S2` → `GGGGSS` |
| `AnY` | n alanines + tyrosine | `A3Y` → `AAAY` |
| `An` | n alanines (no Y) | `A4` → `AAAA` |
| `Gn` | n glycines (no S) | `G4` → `GGGG` |

Repeat counts are capped at 100. 2A entries (codon-frozen, positional)
are rejected in repeat forms — use the base linker once.

Every name resolves through `vaccine_library.get_linker(name)` and
returns a `Linker` with primary-source citations attached. The
default mRNA inter-antigen linker is `(G4S)2` (BioNTech FixVac
canonical, [Sahin *Nature* 2017](https://www.nature.com/articles/nature23003));
the default peptide linker is `G4S3`. Per-junction MHC-aware linker
swap (`--mrna-optimize-linkers`, on by default) considers `G3S`,
`G4S`, `(G3S)2`, `(G4S)2`, `AAA` per junction and substitutes
whichever minimizes predicted presentation of chimeric k-mers
spanning the junction.

All sequences carry primary-source citations in `vaxrank/vaccine_library.py`.

## Clinical Use

Vaxrank is the ranking engine behind the OpenVax neoantigen vaccine
pipeline, which has been used in several clinical trials of personalized
cancer vaccines at Mount Sinai:

- **PGV001** ([NCT02721043](https://clinicaltrials.gov/study/NCT02721043)) —
  A phase I study of personalised neoantigen vaccines in patients with
  solid and haematologic malignancies.  All 11 treated patients developed
  neoantigen-specific T-cell responses
  ([Saxena et al., Cancer Discovery 2025](https://pubmed.ncbi.nlm.nih.gov/40094414/)).
- **PGV001 + atezolizumab in urothelial cancer**
  ([NCT03359239](https://clinicaltrials.gov/study/NCT03359239)) —
  A phase I trial combining PGV001 with checkpoint inhibition.
  The combination was safe and induced neoantigen-specific CD4+ and CD8+
  T-cell responses in all evaluated patients
  ([Saxena et al., Nature Cancer 2025](https://pubmed.ncbi.nlm.nih.gov/40346292/)).
- **PGV001 + TTFields in newly diagnosed glioblastoma**
  ([NCT03223103](https://clinicaltrials.gov/study/NCT03223103)) —
  A phase I trial combining PGV001 with tumor treating fields and
  standard-of-care temozolomide (paper in preparation).

The computational pipeline used in these trials is described in
[Kodysh & Rubinsteyn, Methods Mol. Biol. 2020](https://link.springer.com/protocol/10.1007/978-1-0716-0327-7_10).

## Installation

```
pip install vaxrank
```

**Requirements:** Python 3.9+

Vaxrank uses [PyEnsembl](https://github.com/openvax/pyensembl) for
reference genome annotation.  Install an Ensembl release matching your
reference genome:

```sh
# GRCh38
pyensembl install --release 113 --species human
# GRCh37 (legacy)
pyensembl install --release 75 --species human
```

PDF report generation uses [wkhtmltopdf](http://wkhtmltopdf.org/) by default:

```
brew install --cask wkhtmltopdf
```

Alternatively, pass `--pdf-backend=weasyprint` to use
[WeasyPrint](https://weasyprint.org/) (experimental), which has no external
binary dependency:

```
pip install weasyprint
# macOS also needs: brew install pango
```

On Apple Silicon, WeasyPrint loads Pango via dyld, which doesn't search
Homebrew's `/opt/homebrew/lib` by default. Add this to your shell profile:

```sh
export DYLD_FALLBACK_LIBRARY_PATH="/opt/homebrew/lib:$DYLD_FALLBACK_LIBRARY_PATH"
```

(Intel macOS doesn't need this — Homebrew's `/usr/local/lib` is in dyld's
default fallback path.)

## Configuration

### YAML config file

Common parameters can be stored in a YAML file to avoid repeating them
on every run:

```sh
vaxrank --config my_config.yaml --vcf variants.vcf --bam tumor.bam
```

Example `my_config.yaml`:

```yaml
epitopes:
  min_score: 0.00001                        # drop epitopes below this score
  scoring_mode: affinity                    # "affinity" or "percentile_rank"
  logistic_midpoint: 350.0                  # IC50 (nM) at which score = 0.5
  logistic_width: 150.0                     # steepness of logistic curve
  affinity_cutoff: 5000.0                   # IC50 >= this → score 0
  percentile_rank_cutoff: 10.0              # rank >= this → score 0 (percentile mode)
  top_epitopes_per_candidate: 1000          # 0 = keep all

vaccine_peptides:
  preferred_length: 25                      # target amino acids per vaccine peptide
  min_length: 25                            # minimum vaccine peptide length
  max_length: 25                            # maximum vaccine peptide length
  padding_around_mutation: 5                # off-centre windows to consider
  per_mutation: 1                           # peptides to keep per variant
  max_epitopes_per_candidate: 1000          # 0 = keep all
  score_fraction_of_best: 0.99              # drop candidates scoring < 99% of best
  manufacturability:                        # GRAVY = mean hydropathy
    max_c_terminal_hydropathy: 1.5          # max GRAVY of C-terminal 7-mer
    min_kmer_hydropathy: 0.0                # min max-7mer GRAVY (floor)
    max_kmer_hydropathy_low_priority: 1.5   # low-priority max-7mer GRAVY cap
    max_kmer_hydropathy_high_priority: 2.5  # high-priority max-7mer GRAVY cap
```

### Custom filtering and scoring with the topiary DSL

For anything beyond the scalar logistic / percentile-rank defaults, set
`epitopes.filter_expr` and/or `epitopes.score_expr` to a topiary DSL
string. Both accept the full topiary 5.0 expression grammar (kind
accessors like `affinity` / `presentation`, arithmetic, `&` / `|`,
`.logistic(...)` / `.clip(...)` transforms, `column(col_name)` for raw
DataFrame columns, etc.).

```yaml
epitopes:
  # Drop rows wholesale before scoring
  filter_expr: "affinity <= 500 & affinity.rank <= 2.0"
  # Compute a per-(peptide, allele) score in [0, 1] (binder-quality score)
  score_expr:  "affinity.logistic_normalized(350, 150)"
```

When `filter_expr` is omitted, no rows are dropped up-front; the default
`score_expr` is synthesized from the scalar fields above
(`binding_affinity_cutoff`, `logistic_midpoint`, `logistic_width`, etc.)
and masked so `ic50 >= affinity_cutoff → 0`, reproducing the pre-5.0
behavior byte-for-byte.

Use `affinity.logistic_normalized(m, w)` for a `[0, 1]` binder-quality
score (the topiary 5.1+ primitive); the plain `affinity.logistic(m, w)`
is the raw sigmoid and caps below 1 (≈0.912 at default `m=350, w=150`).

Invalid DSL strings are rejected at config load (not mid-pipeline), so
typos in the YAML surface before any predictions run.

### CLI overrides

CLI arguments override YAML values.  You can also use `--config-value` to
override individual keys without editing the file:

```sh
vaxrank --config my_config.yaml \
  --config-value vaccine_peptides.score_fraction_of_best=0.95 \
  --config-value epitopes.percentile_rank_cutoff=5.0
```

Use `--config-text` when the right-hand side should be kept as a raw
string instead of being YAML-parsed.

### Resolution order

Config values are resolved in order (later wins):

1. Compiled-in defaults (see `vaxrank/config/defaults.py`)
2. YAML config file (`--config`)
3. `--config-value` / `--config-text` overrides
4. Dedicated CLI flags (e.g. `--vaccine-peptide-length`)

### Config reference

#### `EpitopeConfig` — epitope scoring and filtering

| Field | Default | Description |
|-------|---------|-------------|
| `logistic_epitope_score_midpoint` | 350.0 | IC50 (nM) at which epitope score = 0.5 |
| `logistic_epitope_score_width` | 150.0 | Steepness of logistic scoring curve |
| `min_epitope_score` | 0.00001 | Epitopes scoring below this are dropped |
| `binding_affinity_cutoff` | 5000.0 | IC50 >= this → score 0 |
| `scoring_mode` | `"affinity"` | `"affinity"` (IC50-based) or `"percentile_rank"` |
| `percentile_rank_cutoff` | 10.0 | Rank >= this → score 0 (percentile mode) |
| `filter_expr` | `None` | Topiary DSL string; drops rows where the expression is false. Parsed eagerly at config load. |
| `score_expr` | `None` | Topiary DSL string; overrides the default per-`(peptide, allele)` score. |

#### `VaccineConfig` — peptide assembly and manufacturability

| Field | Default | Description |
|-------|---------|-------------|
| `preferred_peptide_length` | 25 | Preferred amino acids per vaccine peptide |
| `min_peptide_length` | 25 | Minimum vaccine peptide length |
| `max_peptide_length` | 25 | Maximum vaccine peptide length |
| `padding_around_mutation` | 5 | Off-centre window positions to consider |
| `max_vaccine_peptides_per_variant` | 1 | Peptides to keep per variant |
| `num_target_epitopes_to_keep` | 1000 | Max epitope predictions per peptide (0 = all) |
| `score_fraction_of_best` | 0.99 | Drop candidates scoring below this fraction of the best |
| `max_c_terminal_hydropathy` | 1.5 | Max GRAVY score of the C-terminal 7-mer |
| `min_kmer_hydropathy` | 0.0 | Minimum max-7mer GRAVY (floor) |
| `max_kmer_hydropathy_low_priority` | 1.5 | Low-priority max-7mer GRAVY cap |
| `max_kmer_hydropathy_high_priority` | 2.5 | High-priority max-7mer GRAVY cap |

The four `*_hydropathy*` fields control the manufacturability tie-breaking
in vaccine peptide ranking.  See `VaccinePeptide.peptide_synthesis_difficulty_score_tuple`
for details on how each threshold is applied.

## MHC Binding Predictors

Vaxrank integrates with MHC binding predictors via
[mhctools](https://github.com/openvax/mhctools).
Use `--mhc-predictor <name>` to select one:

| `--mhc-predictor` | Tool | MHC Class | Notes |
|--------------------|------|-----------|-------|
| `mhcflurry` | [MHCflurry](https://github.com/openvax/mhcflurry) | I | Open-source neural network; installed with mhctools |
| `bigmhc` | [BigMHC](https://github.com/KarchinLab/bigmhc) | I | Auto-detects EL or IM model |
| `bigmhc-el` | [BigMHC](https://github.com/KarchinLab/bigmhc) EL | I | Presentation (eluted ligand) model |
| `bigmhc-im` | [BigMHC](https://github.com/KarchinLab/bigmhc) IM | I | Immunogenicity model |
| `pepsickle` | [Pepsickle](https://github.com/pdxgx/pepsickle) | I | Proteasomal cleavage predictor |
| `netmhc` | [NetMHC](https://services.healthtech.dtu.dk/services/NetMHC-4.0/) | I | Auto-detects NetMHC3 or NetMHC4 |
| `netmhc3` | NetMHC 3.x | I | Requires local install |
| `netmhc4` | [NetMHC 4.0](https://services.healthtech.dtu.dk/services/NetMHC-4.0/) | I | Requires local install |
| `netmhcpan` | [NetMHCpan](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) | I | Auto-detects installed version |
| `netmhcpan28` | NetMHCpan 2.8 | I | Requires local install |
| `netmhcpan3` | NetMHCpan 3.x | I | Requires local install |
| `netmhcpan4` | [NetMHCpan 4.0](https://services.healthtech.dtu.dk/services/NetMHCpan-4.0/) | I | Default mode (EL + BA) |
| `netmhcpan4-ba` | NetMHCpan 4.0 | I | Binding affinity mode only |
| `netmhcpan4-el` | NetMHCpan 4.0 | I | Eluted ligand mode only |
| `netmhcpan41` | [NetMHCpan 4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) | I | Default mode (EL + BA) |
| `netmhcpan41-ba` | NetMHCpan 4.1 | I | Binding affinity mode only |
| `netmhcpan41-el` | NetMHCpan 4.1 | I | Eluted ligand mode only |
| `netmhcpan42` | NetMHCpan 4.2 | I | Default mode (EL + BA) |
| `netmhcpan42-ba` | NetMHCpan 4.2 | I | Binding affinity mode only |
| `netmhcpan42-el` | NetMHCpan 4.2 | I | Eluted ligand mode only |
| `netmhccons` | [NetMHCcons](https://services.healthtech.dtu.dk/services/NetMHCcons-1.1/) | I | Requires local install |
| `netmhcstabpan` | [NetMHCstabpan](https://services.healthtech.dtu.dk/services/NetMHCstabpan-1.0/) | I | Stability predictor; requires local install |
| `netchop` | [NetChop](https://services.healthtech.dtu.dk/services/NetChop-3.1/) | -- | Proteasomal cleavage predictor |
| `netmhciipan` | [NetMHCIIpan](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/) | II | Auto-detects installed version |
| `netmhciipan3` | NetMHCIIpan 3.x | II | Requires local install |
| `netmhciipan4` | [NetMHCIIpan 4.0](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.0/) | II | Default mode (EL + BA) |
| `netmhciipan4-ba` | NetMHCIIpan 4.0 | II | Binding affinity mode only |
| `netmhciipan4-el` | NetMHCIIpan 4.0 | II | Eluted ligand mode only |
| `netmhciipan43` | [NetMHCIIpan 4.3](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/) | II | Default mode (EL + BA) |
| `netmhciipan43-ba` | NetMHCIIpan 4.3 | II | Binding affinity mode only |
| `netmhciipan43-el` | NetMHCIIpan 4.3 | II | Eluted ligand mode only |
| `mixmhcpred` | [MixMHCpred](https://github.com/GfellerLab/MixMHCpred) | I | Requires local install |
| `netmhcpan-iedb` | NetMHCpan via IEDB | I | Uses IEDB web API |
| `netmhccons-iedb` | NetMHCcons via IEDB | I | Uses IEDB web API |
| `netmhciipan-iedb` | NetMHCIIpan via IEDB | II | Uses IEDB web API |
| `smm-iedb` | SMM via IEDB | I | Uses IEDB web API |
| `smm-pmbec-iedb` | SMM-PMBEC via IEDB | I | Uses IEDB web API |
| `random` | Random | -- | Returns random scores; for testing only |

## How It Works

### Upstream inputs

Vaxrank accepts two distinct input shapes, both producing the same
ranked-vaccine-peptides intermediate:

**Full pipeline** (VCF + BAM): Vaxrank does not perform variant
calling or read alignment itself.  Those steps happen upstream,
typically as part of a larger bioinformatics pipeline (e.g.
[neoantigen-vaccine-pipeline](https://github.com/openvax/neoantigen-vaccine-pipeline)):

1. Tumor and matched-normal DNA are sequenced and aligned; a variant
   caller (MuTect, Strelka, etc.) produces a VCF of somatic mutations.
2. Tumor RNA is sequenced and aligned to produce a BAM file.
3. The patient's HLA class I alleles are typed (from sequencing data or
   clinical records).

Vaxrank takes these three inputs — the VCF, the tumor RNA BAM, and the
HLA alleles — runs Isovar transcript assembly + MHC binding prediction
+ ranking, and produces vaccine peptide candidates.

**External-input mode** (`--input-lens` or `--input-pvacseq`): when an
upstream tool (e.g. [LENS](https://github.com/openvax/lens) or pVACseq)
has already produced a per-(peptide, allele) neoepitope report, Vaxrank
skips Isovar + MHC prediction and consumes the report directly. The
per-row `pep_context` (LENS) or `Best Peptide` (pVACseq aggregate) is
used as the SLP-style antigen window. Downstream dispatch — reports +
peptide constructs + mRNA constructs — is identical to the full
pipeline.

### Mutant transcript assembly (Isovar)

For each somatic variant, [Isovar](https://github.com/openvax/isovar)
extracts RNA-seq reads overlapping the mutant locus and assembles them
into a mutant protein fragment.  This is more accurate than simply
applying the DNA variant to the reference transcript because it:

- **Phases** adjacent germline and somatic variants that fall on the same
  read, producing the true amino acid sequence
- **Captures splicing differences** such as intron retention events that
  may alter the reading frame near the mutation
- **Confirms expression** — variants with no supporting RNA reads are
  filtered out

### CandidateEpitope scoring

Each mutant protein fragment is sliced into overlapping subsequences of
epitope length (typically 8–15 amino acids).  These candidate epitopes
are scored for predicted MHC binding affinity using the selected
predictor.  Binding predictions are converted to a score between 0 and 1
via a logistic function parameterised by the `EpitopeConfig` settings.

### Vaccine peptide ranking

Candidate vaccine peptides (longer SLPs, typically 25-mers) are
constructed around each mutation.  Each candidate is scored by the
combined immunogenicity of the epitopes it contains.  Candidates are
then filtered and ranked by:

1. **CandidateEpitope content** — total predicted immunogenicity score
2. **Reference proteome filtering** — peptides matching the human
   reference proteome are removed to ensure only truly novel sequences
   are selected
3. **Cancer hotspot annotation** — variants at known recurrently mutated
   positions (bundled data from
   [cancerhotspots.org](https://www.cancerhotspots.org/), ~2,700
   mutations across cancer types) are flagged
4. **Manufacturability** — tie-breaking by hydropathy-based synthesis
   difficulty (C-terminal and 7-mer window GRAVY scores)

### Data model

Vaxrank's central data unit is the **VaccinePeptide** (VP) — one ranked
candidate of "this is a vaccine peptide we should consider for this
variant." A VP bundles:

- a **`MutantProteinFragment`** — the SLP-style amino-acid sequence
  with mutation positions, gene name, source variant, and the
  ranking-driving expression metrics (`n_alt_reads`, etc.);
- a list of **`EpitopePrediction`** records — per-(k-mer, HLA-allele)
  MHC binding predictions, sorted into a mutant set (overlapping the
  mutation, drives ranking) and a wildtype set (cross-reactivity
  candidates).

The pipeline output is a list of `(varcode.Variant, [VaccinePeptide, ...])`
tuples — each variant has 1 or more VPs depending on
`max_vaccine_peptides_per_variant`:

```
ranked_variants_with_vaccine_peptides = [
    (Variant_A, [VP_A1, VP_A2, ...]),    # multiple windows around variant A's mutation
    (Variant_B, [VP_B1]),                 # single SLP for variant B
    ...
]
```

For each variant, vaxrank can emit multiple alternate constructs:

- `--vaccine-peptide-length` + `--padding-around-mutation` — control
  how the SLP window slides over the mutation site.
- `max_vaccine_peptides_per_variant` (config) — controls how many
  alternate windows per variant make it into the ranked output.
- `--peptide-candidates-per-slot` / `--mrna-candidates-per-slot`
  (CLI) — controls how many VP alternates per variant slot the
  *construct assembler* renders into FASTAs.

Reports render one section per variant; within a section, each VP gets
its own per-epitope sub-table — column counts can differ per VP
(e.g. when pepsickle credibility tagging succeeded for one VP and
failed for another, only the successful VP's table shows the
processing columns).

### Key modules

**Shared upstream:**

- `core_logic.py`: Main vaccine peptide selection algorithm
- `epitope_logic.py`: CandidateEpitope scoring and filtering
- `epitope_io.py`: LENS / pVACseq / vaxrank-native I/O for epitope predictions
- `external_input.py`: Synthesize the canonical ranked-vaccine-peptides shape from a LENS / pVACseq report so external-input runs reach the same dispatch as VCF + BAM
- `reference_proteome.py`: Set-based kmer index for reference proteome filtering (O(1) lookup, built once and cached)
- `cancer_hotspots.py`: Cancer mutation hotspot annotation
- `vaccine_peptide.py`: Vaccine peptide scoring and manufacturability
- `vaccine_library.py`: Shared linker vocabulary + compositional grammar (`(BASE)N`, `GnSm`, `AnY`, `An`, `Gn`) with primary-source citations

**Vaccine-type-specific (downstream):**

- `peptide.py`: Peptide construct assembly + FASTA / JSON manifest / vendor order-form CSV writers; sub-modes `slp` / `minimal_epitope` / `multi_epitope`
- `mrna.py`: mRNA construct assembly + three-FASTA / structured manifest / long-format CSV writers. DnaChisel codon optimization, 2A frozen-codon handling, configurable polyA tail (default A120, optional segmented BNT162b2 pattern), per-junction MHC-aware linker swap (issue #247)
- `mrna_library.py`: mRNA-specific elements (5'/3' UTRs incl. tandem 2× HBB FI; signal peptides HLA-A / HLA-B / tPA / IgK / CD8A / CD28; MITD HLA-A / HLA-B)
- `junction_swap.py`: Per-junction linker optimizer that minimizes predicted MHC presentation of chimeric k-mers spanning antigen junctions

**Reports:**

- `report.py`: Analysis-report generation (ASCII, HTML, PDF, XLSX, CSV, JSON)

## Papers & Citations

**Vaxrank algorithm:**

> Rubinsteyn, A., Hodes, I., Kodysh, J. & Hammerbacher, J.
> [Vaxrank: A Computational Tool For Designing Personalized Cancer Vaccines.](https://doi.org/10.1101/142919)
> *bioRxiv* (2017).

**OpenVax pipeline (methods):**

> Kodysh, J. & Rubinsteyn, A.
> [OpenVax: An Open-Source Computational Pipeline for Cancer Neoantigen Prediction.](https://link.springer.com/protocol/10.1007/978-1-0716-0327-7_10)
> *Methods Mol. Biol.* 2120, 147–160 (2020).

**PGV001 clinical results:**

> Saxena, Marron, Kodysh, et al.
> [PGV001, a Multi-Peptide Personalized Neoantigen Vaccine Platform: Phase I Study in Patients with Solid and Hematologic Malignancies in the Adjuvant Setting.](https://pubmed.ncbi.nlm.nih.gov/40094414/)
> *Cancer Discovery* 15(5), 930–947 (2025).

> Saxena, Anker, Kodysh, et al.
> [Atezolizumab plus personalized neoantigen vaccination in urothelial cancer: a phase 1 trial.](https://pubmed.ncbi.nlm.nih.gov/40346292/)
> *Nature Cancer* 6(6), 988–999 (2025).

BibTeX for the Vaxrank paper:

    @article {Rubinsteyn142919,
        author = {Rubinsteyn, Alex and Hodes, Isaac and Kodysh, Julia and Hammerbacher, Jeffrey},
        title = {Vaxrank: A Computational Tool For Designing Personalized Cancer Vaccines},
        year = {2017},
        doi = {10.1101/142919},
        publisher = {Cold Spring Harbor Laboratory},
        URL = {https://www.biorxiv.org/content/early/2017/05/27/142919},
        journal = {bioRxiv}
    }

## Dependencies

Vaxrank is built on the [OpenVax](https://github.com/openvax) ecosystem:

- [pyensembl](https://github.com/openvax/pyensembl): Reference genome annotation
- [varcode](https://github.com/openvax/varcode): Variant effect prediction from DNA
- [isovar](https://github.com/openvax/isovar): RNA-based mutant transcript assembly and variant phasing
- [mhctools](https://github.com/openvax/mhctools): Unified interface to MHC binding predictors

Other key dependencies:
- `msgspec`: Configuration serialization (YAML/JSON)
- `pandas`, `numpy`: Data processing
- `jinja2`, `pdfkit`/`weasyprint`: Report generation

## Legacy flags

For back-compat with older scripts:

- `--peptide-mode {slp, minimal_epitope, multi_epitope}` is a
  shorthand for the orthogonal axes:
  - `slp` ≡ `--antigen-content mutation_spanning --peptide-antigens-per-construct 1`
  - `minimal_epitope` ≡ `--antigen-content minimal_epitope --peptide-antigens-per-construct 1`
  - `multi_epitope` ≡ `--antigen-content mutation_spanning --peptide-antigens-per-construct N`

The orthogonal axes are preferred for new designs.

## Development

To install Vaxrank for local development:

```bash
git clone git@github.com:openvax/vaxrank.git
cd vaxrank
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
pip install -e .
# Examples; adjust release to match your reference
pyensembl install --release 113 --species human
pyensembl install --release 113 --species mouse
```

Run linting and tests:

```bash
./lint.sh && ./test.sh
```

The first run of the tests may take a while to build the reference proteome kmer index, but subsequent runs will use the cached index.

### Scripts

- `develop.sh`: installs the package in editable mode and sets `PYTHONPATH` to the repo root.
- `lint.sh`: runs ruff on `vaxrank` and `tests`.
- `test.sh`: runs pytest with coverage.
- `deploy.sh`: runs lint/tests, builds a distribution with `build`, uploads via `twine`, and tags the release (`vX.Y.Z`). Deploy is restricted to the `main`/`master` branch.
