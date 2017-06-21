[![Build Status](https://travis-ci.org/hammerlab/vaxrank.svg?branch=master)](https://travis-ci.org/hammerlab/vaxrank) [![Coverage Status](https://coveralls.io/repos/github/hammerlab/vaxrank/badge.svg?branch=master)](https://coveralls.io/github/hammerlab/vaxrank?branch=master)

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

Vaxrank uses [PyEnsembl](https://github.com/hammerlab/pyensembl) for accessing information about the reference genome. You must install an Ensembl release corresponding to the reference genome associated with the mutations provided to Vaxrank. 

The latest release for GRCh38 is Ensembl 87:
```
pyensembl install --release 87 --species human
```

The last release for GRCh37 is Ensembl 75:
```
pyensembl install --release 75 --species human
```

If your variants were called from alignments against hg19 then you can still use GRCh37 but should ignore mitochondrial variants.  

# Development

To install Vaxrank for local development, you may do the below:

```
git clone git@github.com:hammerlab/vaxrank.git
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
