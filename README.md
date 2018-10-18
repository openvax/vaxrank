[![Build Status](https://travis-ci.org/openvax/vaxrank.svg?branch=master)](https://travis-ci.org/openvax/vaxrank) [![Coverage Status](https://coveralls.io/repos/github/openvax/vaxrank/badge.svg?branch=master)](https://coveralls.io/github/openvax/vaxrank?branch=master)
<a href="https://pypi.python.org/pypi/vaxrank/">
    <img src="https://img.shields.io/pypi/v/vaxrank.svg?maxAge=1000" alt="PyPI" />
</a>

# vaxrank

Selection of mutated protein fragments for therapeutic personalized cancer vaccines.

## Usage

```sh

vaxrank \
    --vcf test/data/b16.f10/b16.vcf \
    --bam test/data/b16.f10/b16.combined.bam \
    --vaccine-peptide-length 25 \
    --mhc-predictor netmhc \
    --mhc-alleles H2-Kb,H2-Db \
    --padding-around-mutation 5 \
    --output-ascii-report vaccine-peptides.txt \
    --output-pdf-report vaccine-peptides.pdf \
    --output-html-report vaccine-peptides.html
```

## Installation

Vaxrank can be installed using [pip](https://packaging.python.org/installing/#use-pip-for-installing):

```
pip install vaxrank
```

Note: to generate PDF reports, you first need to install [wkhtmltopdf](http://wkhtmltopdf.org/), which you can do (on OS X) like so:

```
brew install Caskroom/cask/wkhtmltopdf
```

Vaxrank uses [PyEnsembl](https://github.com/openvax/pyensembl) for accessing information about the reference genome. You must install an Ensembl release corresponding to the reference genome associated with the mutations provided to Vaxrank.

The latest release for GRCh38 is Ensembl 93:
```
pyensembl install --release 93 --species human
```

The last release for GRCh37 is Ensembl 75:
```
pyensembl install --release 75 --species human
```

If your variants were called from alignments against hg19 then you can still use GRCh37 but should ignore mitochondrial variants.

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


# Development

To install Vaxrank for local development, you may do the below:

```
git clone git@github.com:openvax/vaxrank.git
conda create -q -n vaxrank-dev-env python=3.5.2 numpy scipy nose pandas pylint
source activate vaxrank-dev-env
pip install -r requirements.txt
pip install .
pyensembl install --release 87 --species human
pyensembl install --release 87 --species mouse
```

You should run the linter and the test suite as you work on Vaxrank (and these will be run automatically by our continuous integration server up on a PR being made).

```
./lint.sh
nosetests test
```

The first run of the tests may take a while (8 minutes on a 2016 Macbook Pro) to create the FM index of the proteome, but subsequent tests should take only a few seconds.

