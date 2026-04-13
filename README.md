[![Tests](https://github.com/openvax/vaxrank/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/vaxrank/actions/workflows/tests.yml)
[![Coverage Status](https://coveralls.io/repos/github/openvax/vaxrank/badge.svg?branch=master)](https://coveralls.io/github/openvax/vaxrank?branch=master)
[![Docs](https://github.com/openvax/vaxrank/actions/workflows/docs.yml/badge.svg)](https://github.com/openvax/vaxrank/actions/workflows/docs.yml)
[![GitHub Pages](https://github.com/openvax/vaxrank/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/openvax/vaxrank/actions/workflows/pages/pages-build-deployment)
<a href="https://pypi.python.org/pypi/vaxrank/">
    <img src="https://img.shields.io/pypi/v/vaxrank.svg?maxAge=1000" alt="PyPI" />
</a>

# vaxrank

Vaxrank is the vaccine peptide ranking component of the
[OpenVax](https://www.openvax.org/) pipeline for designing personalized
cancer vaccines.  Given a patient's somatic mutations, tumor RNA sequencing
data, and HLA type, Vaxrank selects and ranks the mutant peptides most
likely to elicit a T-cell response, producing a report suitable for
guiding vaccine manufacture.

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

Vaxrank outputs ranked reports in ASCII, HTML, PDF, and XLSX formats.
Each report lists the top vaccine peptide candidates per variant, their
predicted epitopes, and supporting evidence from the RNA data.

## Clinical Use

Vaxrank is the ranking engine behind the OpenVax neoantigen vaccine
pipeline, which has been used in several clinical trials of personalized
cancer vaccines at Mount Sinai:

- **PGV001** ([NCT02721043](https://clinicaltrials.gov/study/NCT02721043)) —
  A phase I study of personalised neoantigen vaccines in patients with
  solid and haematologic malignancies.  All 11 treated patients developed
  neoantigen-specific T-cell responses
  ([Bortman et al., Cancer Discovery 2025](https://pubmed.ncbi.nlm.nih.gov/40094414/)).
- **PGV001 + atezolizumab in urothelial cancer**
  ([NCT03359239](https://clinicaltrials.gov/study/NCT03359239)) —
  A phase I trial combining PGV001 with checkpoint inhibition.
  The combination was safe and induced neoantigen-specific CD4+ and CD8+
  T-cell responses in all evaluated patients
  ([Galsky et al., Nature Cancer 2025](https://www.nature.com/articles/s43018-025-00966-7)).
- **PGV001 + TTFields in newly diagnosed glioblastoma**
  ([NCT03223103](https://clinicaltrials.gov/study/NCT03223103)) —
  A phase I trial combining PGV001 with tumor treating fields and
  standard-of-care temozolomide (paper in preparation).

The computational pipeline used in these trials is described in
[Kodysh & Rubinsteyn, Methods Mol. Biol. 2020](https://link.springer.com/protocol/10.1007/978-1-0716-0327-7_10).

## Quick Start

```sh
vaxrank \
    --vcf tests/data/b16.f10/b16.vcf \
    --bam tests/data/b16.f10/b16.combined.bam \
    --vaccine-peptide-length 25 \
    --mhc-predictor netmhc \
    --mhc-alleles H2-Kb,H2-Db \
    --padding-around-mutation 5 \
    --output-ascii-report vaccine-peptides.txt \
    --output-pdf-report vaccine-peptides.pdf \
    --output-html-report vaccine-peptides.html
```

Inputs:
- `--vcf` — Somatic variants (VCF from any variant caller)
- `--bam` — Tumor RNA-seq alignments (used by Isovar to assemble mutant transcripts)
- `--mhc-alleles` — Patient HLA alleles (e.g. `HLA-A*02:01,HLA-B*07:02`)
- `--mhc-predictor` — Which MHC binding predictor to use (see table below)

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
  # Compute a per-(peptide, allele) score (numeric DSL node)
  score_expr:  "affinity.logistic(350, 150)"
```

When `filter_expr` is omitted, no rows are dropped up-front; the default
`score_expr` is synthesized from the scalar fields above
(`binding_affinity_cutoff`, `logistic_midpoint`, `logistic_width`, etc.)
and masked so `ic50 >= affinity_cutoff → 0`, reproducing the pre-5.0
behavior byte-for-byte.

> **Note:** a user-supplied `score_expr` of `affinity.logistic(m, w)` is
> the *raw* topiary sigmoid (maximum ≈ 0.912 at default params), not the
> normalized `[0, 1]` score the defaults produce. Divide by
> `1/(1+exp(-m/w))` if you need the normalized form. Follow-up node
> [openvax/topiary#116](https://github.com/openvax/topiary/issues/116)
> adds a `logistic_normalized` sibling so the constant doesn't have to
> live in config.

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
| `num_mutant_epitopes_to_keep` | 1000 | Max epitope predictions per peptide (0 = all) |
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

Vaxrank does not perform variant calling or read alignment itself.
Those steps happen upstream, typically as part of a larger
bioinformatics pipeline (e.g.
[neoantigen-vaccine-pipeline](https://github.com/openvax/neoantigen-vaccine-pipeline)):

1. Tumor and matched-normal DNA are sequenced and aligned; a variant
   caller (MuTect, Strelka, etc.) produces a VCF of somatic mutations.
2. Tumor RNA is sequenced and aligned to produce a BAM file.
3. The patient's HLA class I alleles are typed (from sequencing data or
   clinical records).

Vaxrank takes these three inputs — the VCF, the tumor RNA BAM, and the
HLA alleles — and produces a ranked list of vaccine peptide candidates.

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

### Epitope scoring

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

1. **Epitope content** — total predicted immunogenicity score
2. **Reference proteome filtering** — peptides matching the human
   reference proteome are removed to ensure only truly novel sequences
   are selected
3. **Cancer hotspot annotation** — variants at known recurrently mutated
   positions (bundled data from
   [cancerhotspots.org](https://www.cancerhotspots.org/), ~2,700
   mutations across cancer types) are flagged
4. **Manufacturability** — tie-breaking by hydropathy-based synthesis
   difficulty (C-terminal and 7-mer window GRAVY scores)

### Key modules

- `core_logic.py`: Main vaccine peptide selection algorithm
- `epitope_logic.py`: Epitope scoring and filtering
- `reference_proteome.py`: Set-based kmer index for reference proteome filtering (O(1) lookup, built once and cached)
- `cancer_hotspots.py`: Cancer mutation hotspot annotation
- `vaccine_peptide.py`: Vaccine peptide scoring and manufacturability
- `report.py`: Report generation (ASCII, HTML, PDF, XLSX)

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

> Bortman et al.
> [PGV001, a Multi-Peptide Personalized Neoantigen Vaccine Platform: Phase I Study in Patients with Solid and Hematologic Malignancies in the Adjuvant Setting.](https://pubmed.ncbi.nlm.nih.gov/40094414/)
> *Cancer Discovery* 15(5), 930–945 (2025).

> Galsky et al.
> [Atezolizumab plus personalized neoantigen vaccination in urothelial cancer: a phase 1 trial.](https://www.nature.com/articles/s43018-025-00966-7)
> *Nature Cancer* (2025).

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
