[![Tests](https://github.com/openvax/vaxrank/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/vaxrank/actions/workflows/tests.yml)
[![Coverage Status](https://coveralls.io/repos/github/openvax/vaxrank/badge.svg?branch=master)](https://coveralls.io/github/openvax/vaxrank?branch=master)
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
design workflows are experimental.

## Install

```sh
pip install vaxrank
pyensembl install --release 113 --species human   # match your reference genome
```

Vaxrank needs Python 3.11 or newer. PDF reports additionally need wkhtmltopdf or
WeasyPrint; see [installation](docs/installation.md).

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
    --output-ascii-report vaccine-peptides.txt \
    --output-pdf-report vaccine-peptides.pdf \
    --pdf-backend weasyprint
```

For an offline canine DNA/RNA replay, see the
[WLOS_0001 TP53 example](examples/wlos0001_tp53/README.md). It uses original
tumor RNA and a published matched-normal somatic call with a random test
predictor; the dog's DLA genotype is unknown.

For a real VCF + RNA run, supply:

- `--vcf`: somatic variants, from any variant caller
- `--bam`: tumor RNA-seq alignments, used by Isovar to assemble mutant transcripts
- `--mhc-alleles`: the patient's MHC alleles, e.g. `HLA-A*02:01,HLA-B*07:02`
- `--mhc-predictor`: which binding predictor to run
  ([available predictors](docs/mhc-predictors.md))

To get vaccine constructs rather than reports, add `--output-dir` and, for mRNA,
`--vaccine-type mrna`. [Vaccine designs](docs/vaccine-designs.md) describes the
peptide and mRNA designs and the files each one writes.

### Start from an already-scored report

When MHC prediction has happened upstream, design directly from a LENS or
pVACseq report instead of a VCF and BAM:

```sh
vaxrank --input-lens patient.lens.tsv \
        --vaccine-type mrna --output-dir mrna_out/ \
        --ensembl-release 102
```

Use the Ensembl release the report was annotated with. Original predictions are
reused by default. [Input workflows](docs/input-workflows.md) covers pVACseq
reports, saved Vaxrank evidence, combining several reports in one run, and when
to compute new predictions instead.

## Documentation

Full documentation: <https://openvax.github.io/vaxrank/>

- [Installation](docs/installation.md): reference annotation and PDF backends.
- [Input workflows](docs/input-workflows.md): VCF/BAM, LENS, pVACseq, Topiary
  tables, CTA expression, antigen sequences and saved evidence.
- [Vaccine designs](docs/vaccine-designs.md): the design axes, peptide and mRNA
  construct files, linker library, manifest schema.
- [Analysis outputs](docs/outputs.md): report formats and the figure commands.
- [Configuration](docs/configuration.md): YAML config, override order, named
  policy bundles, RNA reconstruction context, Topiary scoring DSL.
- [MHC binding predictors](docs/mhc-predictors.md): supported predictor names.
- [How Vaxrank works](docs/how-it-works.md): pipeline stages, ranking criteria,
  data model, module map.
- [Varcode, Isovar, and Vaxrank](docs/library-responsibilities.md): which
  library owns which decision.
- [Development](docs/development.md): local setup, test data, release scripts.

Experimental workflows: [unified antigen inputs](docs/antigen-inputs.md),
[CTA expression inputs](docs/design/expression-cta-inputs.md),
[mixed mutation/CTA/viral selection](docs/mixed-antigen-design.md) and
[window exposure](docs/core-window-exposure.md).

<a id="clinical-use"></a>

## Clinical studies

The OpenVax pipeline has used Vaxrank in these Mount Sinai trials:

- [PGV001 / NCT02721043](https://clinicaltrials.gov/study/NCT02721043): personalized vaccines for solid and hematologic malignancies.
- [NCT03359239](https://clinicaltrials.gov/study/NCT03359239): PGV001 with atezolizumab in urothelial cancer.
- [NCT03223103](https://clinicaltrials.gov/study/NCT03223103): PGV001 with tumor treating fields and temozolomide in glioblastoma.

Methods and published results are cited below.

<a id="papers--citations"></a>

## Papers and citations

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
