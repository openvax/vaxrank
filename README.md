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
````
