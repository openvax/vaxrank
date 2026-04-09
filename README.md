[![Tests](https://github.com/openvax/vaxrank/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/vaxrank/actions/workflows/tests.yml)
[![Coverage Status](https://coveralls.io/repos/github/openvax/vaxrank/badge.svg?branch=master)](https://coveralls.io/github/openvax/vaxrank?branch=master)
[![Docs](https://github.com/openvax/vaxrank/actions/workflows/docs.yml/badge.svg)](https://github.com/openvax/vaxrank/actions/workflows/docs.yml)
[![GitHub Pages](https://github.com/openvax/vaxrank/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/openvax/vaxrank/actions/workflows/pages/pages-build-deployment)
<a href="https://pypi.python.org/pypi/vaxrank/">
    <img src="https://img.shields.io/pypi/v/vaxrank.svg?maxAge=1000" alt="PyPI" />
</a>

# vaxrank

Selection of mutated protein fragments for therapeutic personalized cancer vaccines.

## Usage

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

### Using a YAML Configuration File

You can specify common parameters in a YAML configuration file to avoid
repeating them on every run:

```sh
vaxrank --config my_config.yaml --vcf variants.vcf --bam tumor.bam
```

Example `my_config.yaml`:

```yaml
epitope_config:
  logistic_epitope_score_midpoint: 350.0  # IC50 (nM) at which score = 0.5
  logistic_epitope_score_width: 150.0     # steepness of logistic curve
  min_epitope_score: 0.00001              # drop epitopes below this score
  binding_affinity_cutoff: 5000.0         # IC50 >= this → score 0
  scoring_mode: affinity                  # "affinity" or "percentile_rank"
  percentile_rank_cutoff: 10.0            # rank >= this → score 0 (percentile mode)

vaccine_config:
  vaccine_peptide_length: 25              # amino acids per vaccine peptide
  padding_around_mutation: 5              # off-centre windows to consider
  max_vaccine_peptides_per_variant: 1     # peptides to keep per variant
  num_mutant_epitopes_to_keep: 1000       # 0 = keep all
  score_fraction_of_best: 0.99           # drop candidates scoring < 99% of best

  # Manufacturability thresholds (GRAVY = mean hydropathy)
  max_c_terminal_hydropathy: 1.5          # max GRAVY of C-terminal 7-mer
  min_kmer_hydropathy: 0.0                # min max-7mer GRAVY (floor)
  max_kmer_hydropathy_low_priority: 1.5   # low-priority max-7mer GRAVY cap
  max_kmer_hydropathy_high_priority: 2.5  # high-priority max-7mer GRAVY cap
```

CLI arguments override values from the config file.  You can also use
`--set` to override any config value without editing the file:

```sh
vaxrank --config my_config.yaml \
  --set vaccine_config.score_fraction_of_best=0.95 \
  --set epitope_config.percentile_rank_cutoff=5.0
```

## Installation

Vaxrank can be installed using [pip](https://packaging.python.org/installing/#use-pip-for-installing):

```
pip install vaxrank
```

**Requirements:** Python 3.9+

Note: to generate PDF reports, you first need to install [wkhtmltopdf](http://wkhtmltopdf.org/), which you can do (on macOS) like so:

```
brew install --cask wkhtmltopdf
```

Vaxrank uses [PyEnsembl](https://github.com/openvax/pyensembl) for accessing information about the reference genome. You must install an Ensembl release corresponding to the reference genome associated with the mutations provided to Vaxrank.

Example for GRCh38 (adjust release to match your reference):
```
pyensembl install --release 113 --species human
```

Example for GRCh37 (legacy):
```
pyensembl install --release 75 --species human
```

If your variants were called from alignments against hg19 then you can still use GRCh37 but should ignore mitochondrial variants.

## Features

### Reference Proteome Filtering

Vaxrank filters out peptides that exist in the reference proteome to focus on truly novel mutant sequences. This uses a set-based kmer index for O(1) membership testing. The index is built once and cached locally for subsequent runs.

### Cancer Hotspot Annotation

Vaxrank annotates variants that occur at known cancer mutation hotspots using bundled data from [cancerhotspots.org](https://www.cancerhotspots.org/) (Chang et al. 2016, 2017). This helps identify clinically relevant mutations. The hotspot data includes ~2,700 recurrently mutated positions across cancer types.

### MHC Binding Prediction

Vaxrank integrates with multiple MHC binding predictors via [mhctools](https://github.com/openvax/mhctools), including:
- NetMHC / NetMHCpan
- MHCflurry (open source, installed by default)

## Paper & Citation

There is a Vaxrank paper on biorxiv called [Vaxrank: A Computational Tool For Designing Personalized Cancer Vaccines](https://www.biorxiv.org/content/early/2017/05/27/142919) which can be cited as:

    @article {Rubinsteyn142919,
        author = {Rubinsteyn, Alex and Hodes, Isaac and Kodysh, Julia and Hammerbacher, Jeffrey},
        title = {Vaxrank: A Computational Tool For Designing Personalized Cancer Vaccines},
        year = {2017},
        doi = {10.1101/142919},
        publisher = {Cold Spring Harbor Laboratory},
        abstract = {Therapeutic vaccines targeting mutant tumor antigens ({\textquotedblleft}neoantigens{\textquotedblright}) are an increasingly popular form of personalized cancer immunotherapy. Vaxrank is a computational tool for selecting neoantigen vaccine peptides from tumor mutations, tumor RNA data, and patient HLA type. Vaxrank is freely available at www.github.com/hammerlab/vaxrank under the Apache 2.0 open source license and can also be installed from the Python Package Index.},
        URL = {https://www.biorxiv.org/content/early/2017/05/27/142919},
        eprint = {https://www.biorxiv.org/content/early/2017/05/27/142919.full.pdf},
        journal = {bioRxiv}
    }


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

## Architecture

### Configuration

Vaxrank uses [msgspec](https://jcristharif.com/msgspec/) frozen Struct objects
for configuration, with all defaults centralised in
`vaxrank/config/defaults.py`.  Config values are resolved in order:

1. Compiled-in defaults
2. YAML config file (`--config`)
3. `--set` / `--expr` CLI overrides
4. Dedicated CLI flags (e.g. `--vaccine-peptide-length`)

#### `EpitopeConfig` — epitope scoring and filtering

| Field | Default | Description |
|-------|---------|-------------|
| `logistic_epitope_score_midpoint` | 350.0 | IC50 (nM) at which epitope score = 0.5 |
| `logistic_epitope_score_width` | 150.0 | Steepness of logistic scoring curve |
| `min_epitope_score` | 0.00001 | Epitopes scoring below this are dropped |
| `binding_affinity_cutoff` | 5000.0 | IC50 >= this → score 0 |
| `scoring_mode` | `"affinity"` | `"affinity"` (IC50-based) or `"percentile_rank"` |
| `percentile_rank_cutoff` | 10.0 | Rank >= this → score 0 (percentile mode) |

#### `VaccineConfig` — peptide assembly and manufacturability

| Field | Default | Description |
|-------|---------|-------------|
| `vaccine_peptide_length` | 25 | Amino acids per vaccine peptide |
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

### Key Modules

- **`reference_proteome.py`**: Set-based kmer index for checking if peptides exist in the reference proteome
- **`cancer_hotspots.py`**: Lookup for known cancer mutation hotspots
- **`epitope_logic.py`**: Epitope scoring and filtering logic
- **`core_logic.py`**: Main vaccine peptide selection algorithm
- **`report.py`**: Report generation (ASCII, HTML, PDF, XLSX)

## Dependencies

Key dependencies:
- `pyensembl`: Reference genome annotation
- `varcode`: Variant effect prediction
- `isovar`: RNA-based variant calling
- `mhctools`: MHC binding prediction
- `msgspec`: Configuration serialization (YAML/JSON)
- `pandas`, `numpy`: Data processing
- `jinja2`, `pdfkit`: Report generation

## Scripts

Helper scripts included in the repo:
- `develop.sh`: installs the package in editable mode and sets `PYTHONPATH` to the repo root.
- `lint.sh`: runs ruff on `vaxrank` and `tests`.
- `test.sh`: runs pytest with coverage.
- `deploy.sh`: runs lint/tests, builds a distribution with `build`, uploads via `twine`, and tags the release (`vX.Y.Z`). Deploy is restricted to the `main`/`master` branch.
