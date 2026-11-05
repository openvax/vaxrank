[![Tests](https://github.com/openvax/vaxrank/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/vaxrank/actions/workflows/tests.yml)
[![Coverage Status](https://coveralls.io/repos/github/openvax/vaxrank/badge.svg?branch=main)](https://coveralls.io/github/openvax/vaxrank?branch=main)
[![Docs](https://github.com/openvax/vaxrank/actions/workflows/docs.yml/badge.svg)](https://github.com/openvax/vaxrank/actions/workflows/docs.yml)
[![GitHub Pages](https://github.com/openvax/vaxrank/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/openvax/vaxrank/actions/workflows/pages/pages-build-deployment)
<a href="https://pypi.python.org/pypi/vaxrank/">
    <img src="https://img.shields.io/pypi/v/vaxrank.svg?maxAge=1000" alt="PyPI" />
</a>

# Vaxrank

Vaxrank selects mutated protein fragments for personalized therapeutic cancer
vaccines. It is the antigen ranking and vaccine design component of the
[OpenVax](https://www.openvax.org/) pipeline.

It takes somatic variants, tumor RNA-seq and MHC alleles, or previously scored
candidate epitopes, and returns ranked antigen windows,
peptide or mRNA construct files, and reports for reviewing the evidence behind
each design. Unified mutation/CTA/viral inputs and the expanded peptide/mRNA
design workflows are experimental. Vaxrank builds on the openvax stack:
[varcode](https://github.com/openvax/varcode) annotates variants,
[isovar](https://github.com/openvax/isovar) assembles mutant protein sequence
from tumor RNA, and [topiary](https://github.com/openvax/topiary) with
[mhctools](https://github.com/openvax/mhctools) predicts and ranks MHC ligands.

## Install

```sh
pip install vaxrank
pyensembl install --release 113 --species human   # match your reference genome
```

Vaxrank needs Python 3.11 or newer. PDF reports additionally need wkhtmltopdf or
WeasyPrint's system libraries; see [installation](https://openvax.github.io/vaxrank/installation/).

## Quick start

Choose the vaccine type first: `peptide`, `mrna`, or both. This smoke example
uses the mouse B16 melanoma files in a repository checkout and a **random test
predictor**. It exercises the workflow; its scores are not biological predictions.

```sh
vaxrank \
    --vaccine-type peptide \
    --download-reference-genome-data \
    --vcf tests/data/b16.f10/b16.vcf \
    --bam tests/data/b16.f10/b16.combined.bam \
    --mhc-predictor random \
    --mhc-alleles H2-Kb,H2-Db \
    --output-ascii-report vaccine-peptides.txt
```

The report starts with what was read and what survived each stage:

```
Vaccine Design Report
Vaccine type: peptide
MHC alleles: H-2-Kb H-2-Db
...
Parsed variants loaded: 5
Variants with annotation results: 5
Coding-effect variants: 4
Coding-effect variants with an RNA-derived protein: 2
```

It continues with the selection rules, each ranked candidate window with its
scored epitopes, per-source outcomes and provenance. Add `--output-pdf-report
vaccine-peptides.pdf` for a PDF (see [PDF reports](https://openvax.github.io/vaxrank/installation/#pdf-reports)).

For a real VCF + RNA run, supply:

- `--vcf`: somatic variants, from any variant caller
- `--bam`: tumor RNA-seq alignments, used by Isovar to assemble mutant transcripts
- `--mhc-alleles`: the patient's MHC alleles, as separate quoted arguments or a
  quoted comma- or space-separated list
- `--mhc-predictor`: which binding predictor to run
  ([available predictors](https://openvax.github.io/vaxrank/mhc-predictors/))

These allele-list forms work:

```sh
--mhc-alleles 'HLA-A*02:01' 'HLA-B*07:02'
--mhc-alleles 'HLA-A*02:01,HLA-B*07:02'
--mhc-alleles 'HLA-A*02:01 HLA-B*07:02'
```

Keep the quotes: an unquoted `*` can expand to filenames, or cause zsh to
reject the command before Vaxrank starts.

To emit vaccine constructs, add `--output-dir` (both modalities by default).
Use `--vaccine-type peptide` or `--vaccine-type mrna` to select one modality.
[Vaccine designs](https://openvax.github.io/vaxrank/vaccine-designs/) describes the peptide and mRNA
designs and the files each one writes. For an offline canine DNA/RNA replay, see
the [WLOS_0001 TP53 example](https://github.com/openvax/vaxrank/tree/main/examples/wlos0001_tp53).

## Command line

`vaxrank` is the main command; run `vaxrank --help` for every option and
`vaxrank --version` for the installed release. Besides VCF + BAM input as
above, it designs directly from a LENS or pVACseq report when MHC prediction
has happened upstream:

```sh
vaxrank --input-lens patient.lens.tsv \
        --vaccine-type mrna --output-dir mrna_out/ \
        --ensembl-release 102
```

Use the Ensembl release the report was annotated with. Original predictions are
reused by default. [Input workflows](https://openvax.github.io/vaxrank/input-workflows/) covers pVACseq
reports, saved Vaxrank evidence, combining several reports in one run, and when
to compute new predictions instead.

Companion commands share the installed package:

| Command | Purpose | Documented in |
|---|---|---|
| `vaxrank-mutation-figure` | Transcript and protein context figures from `--output-isovar-csv` | [Analysis outputs](https://openvax.github.io/vaxrank/outputs/#mutation-and-assembly-figures) |
| `vaxrank-evidence-figure` | Evidence summary figure from an evidence JSON | [Analysis outputs](https://openvax.github.io/vaxrank/outputs/) |
| `vaxrank-complex-result-figure` | Figures for complex-variant results JSON | [Analysis outputs](https://openvax.github.io/vaxrank/outputs/) |
| `vaxrank-window-exposure` | Delivery exposure of admitted SLP windows | [Window exposure](https://openvax.github.io/vaxrank/core-window-exposure/) |
| `vaxrank-test-data` | Download the shared Sid test data | [Test data](https://github.com/openvax/vaxrank/blob/main/TEST_DATA.md) |

## Documentation

Full documentation: <https://openvax.github.io/vaxrank/>. In reading order:

1. [Installation](https://openvax.github.io/vaxrank/installation/): reference annotation and PDF backends.
2. [Input workflows](https://openvax.github.io/vaxrank/input-workflows/): VCF/BAM, LENS, pVACseq, Topiary tables, CTA expression, antigen sequences and saved evidence.
3. [Vaccine designs](https://openvax.github.io/vaxrank/vaccine-designs/): the design axes, peptide and mRNA construct files, linker library, manifest schema.
4. [Analysis outputs](https://openvax.github.io/vaxrank/outputs/): report formats and the figure commands.
5. [Configuration](https://openvax.github.io/vaxrank/configuration/): YAML config, override order, named policy bundles, RNA reconstruction context, Topiary scoring DSL.
6. [How Vaxrank works](https://openvax.github.io/vaxrank/how-it-works/): pipeline stages, ranking criteria, data model, module map.

The site also covers [MHC binding predictors](https://openvax.github.io/vaxrank/mhc-predictors/), the
[canine workflow](https://openvax.github.io/vaxrank/canine-vcf-workflow/),
[which library owns which decision](https://openvax.github.io/vaxrank/library-responsibilities/), and the
experimental [antigen](https://openvax.github.io/vaxrank/antigen-inputs/),
[Isovar](https://openvax.github.io/vaxrank/isovar-inputs/),
[CTA expression](https://openvax.github.io/vaxrank/design/expression-cta-inputs/),
[mixed-antigen](https://openvax.github.io/vaxrank/mixed-antigen-design/) and
[window exposure](https://openvax.github.io/vaxrank/core-window-exposure/) workflows.

## Development

```sh
git clone https://github.com/openvax/vaxrank.git && cd vaxrank
python -m venv venv && source venv/bin/activate
pip install -e "." pytest pytest-cov ruff
./lint.sh
./test.sh
```

Always run tests through `./test.sh`, not bare `pytest`. Every pull request
bumps `vaxrank/version.py`. See [Development](https://openvax.github.io/vaxrank/development/) for test
data, reference downloads and releases.

<a id="clinical-use"></a>

## Clinical studies

The OpenVax pipeline has used Vaxrank in these Mount Sinai trials:

- [PGV001 / NCT02721043](https://clinicaltrials.gov/study/NCT02721043): personalized vaccines for solid and hematologic malignancies.
- [NCT03359239](https://clinicaltrials.gov/study/NCT03359239): PGV001 with atezolizumab in urothelial cancer.
- [NCT03223103](https://clinicaltrials.gov/study/NCT03223103): PGV001 with tumor treating fields and temozolomide in glioblastoma.

Methods and published results are cited below.

<a id="papers--citations"></a>

## Citation

Vaxrank algorithm:

> Rubinsteyn, A., Hodes, I., Kodysh, J. & Hammerbacher, J.
> [Vaxrank: A Computational Tool For Designing Personalized Cancer Vaccines.](https://doi.org/10.1101/142919)
> *bioRxiv* (2017).

OpenVax pipeline (methods):

> Kodysh, J. & Rubinsteyn, A.
> [OpenVax: An Open-Source Computational Pipeline for Cancer Neoantigen Prediction.](https://link.springer.com/protocol/10.1007/978-1-0716-0327-7_10)
> *Methods Mol. Biol.* 2120, 147–160 (2020).

PGV001 clinical results:

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

## License

Apache 2.0. Contributions are welcome; see
[CONTRIBUTING.md](https://github.com/openvax/vaxrank/blob/main/CONTRIBUTING.md)
and the
[code of conduct](https://github.com/openvax/vaxrank/blob/main/code-of-conduct.md).
