[![Tests](https://github.com/openvax/vaxrank/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/vaxrank/actions/workflows/tests.yml)
[![Coverage Status](https://coveralls.io/repos/github/openvax/vaxrank/badge.svg?branch=master)](https://coveralls.io/github/openvax/vaxrank?branch=master)
[![Docs](https://github.com/openvax/vaxrank/actions/workflows/docs.yml/badge.svg)](https://github.com/openvax/vaxrank/actions/workflows/docs.yml)
[![GitHub Pages](https://github.com/openvax/vaxrank/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/openvax/vaxrank/actions/workflows/pages/pages-build-deployment)
<a href="https://pypi.python.org/pypi/vaxrank/">
    <img src="https://img.shields.io/pypi/v/vaxrank.svg?maxAge=1000" alt="PyPI" />
</a>

# vaxrank

Vaxrank selects mutated protein fragments for personalized therapeutic cancer
vaccines. It is the antigen ranking and vaccine design component of the
[OpenVax](https://www.openvax.org/) pipeline.

It takes somatic variants, tumor RNA-seq and HLA types — or candidate epitopes
that another tool has already scored — and returns ranked vaccine peptides,
peptide or mRNA construct files, and reports for reviewing the evidence behind
each design.

## Install

```sh
pip install vaxrank
pyensembl install --release 113 --species human   # match your reference genome
```

Vaxrank needs Python 3.10 or newer. PDF reports additionally need wkhtmltopdf or
WeasyPrint; see [installation](docs/installation.md).

## Quick start

Ranked vaccine peptides from a tumor VCF plus an RNA-seq BAM, as text and PDF
reports (this command uses the mouse B16 melanoma data in this repository):

```sh
vaxrank \
    --vcf tests/data/b16.f10/b16.vcf \
    --bam tests/data/b16.f10/b16.combined.bam \
    --mhc-predictor netmhc \
    --mhc-alleles H2-Kb,H2-Db \
    --output-ascii-report vaccine-peptides.txt \
    --output-pdf-report vaccine-peptides.pdf
```

Four inputs are required:

- `--vcf`: somatic variants, from any variant caller
- `--bam`: tumor RNA-seq alignments, used by Isovar to assemble mutant transcripts
- `--mhc-alleles`: the patient's HLA alleles, e.g. `HLA-A*02:01,HLA-B*07:02`
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
  tables, saved evidence, and how to combine them.
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

Experimental research workflows: [unified antigen inputs](docs/antigen-inputs.md),
[CTA expression inputs](docs/design/expression-cta-inputs.md),
[mixed mutation/CTA/viral selection](docs/mixed-antigen-design.md) and
[window exposure](docs/core-window-exposure.md). These do not establish clinical
safety or efficacy.

## Clinical use

Vaxrank is the ranking engine behind the OpenVax neoantigen vaccine pipeline,
which has been used in several clinical trials of personalized cancer vaccines
at Mount Sinai:

- [NCT02721043](https://clinicaltrials.gov/study/NCT02721043) (PGV001): phase I
  study in solid and hematologic malignancies. All 11 treated patients developed
  neoantigen-specific T-cell responses
  ([Saxena et al., Cancer Discovery 2025](https://pubmed.ncbi.nlm.nih.gov/40094414/)).
- [NCT03359239](https://clinicaltrials.gov/study/NCT03359239): phase I of PGV001
  plus atezolizumab in urothelial cancer. The combination was safe and induced
  neoantigen-specific CD4+ and CD8+ T-cell responses in all evaluated patients
  ([Saxena et al., Nature Cancer 2025](https://pubmed.ncbi.nlm.nih.gov/40346292/)).
- [NCT03223103](https://clinicaltrials.gov/study/NCT03223103): phase I of PGV001
  plus tumor treating fields and standard-of-care temozolomide in newly
  diagnosed glioblastoma (paper in preparation).

The computational pipeline used in these trials is described in
[Kodysh & Rubinsteyn, Methods Mol. Biol. 2020](https://link.springer.com/protocol/10.1007/978-1-0716-0327-7_10).

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
