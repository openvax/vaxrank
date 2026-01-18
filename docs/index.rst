.. vaxrank documentation master file

Vaxrank Documentation
=====================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started
   configuration
   api

Getting Started With Vaxrank
============================

Overview
--------
Vaxrank is a tool for selecting mutated peptides for use in personalized therapeutic cancer vaccination. Vaxrank determines which peptides should be used in a vaccine from tumor-specific somatic mutations, tumor RNA sequencing data, and a patient's HLA type. Additionally, Vaxrank considers surrounding non-mutated residues in a peptide to prioritize vaccine peptide candidates and improve the odds of successful synthesis.

Questions, Bug Reporting, and Issue Tracking
--------------------------------------------
Questions, bug reporting and issue tracking are provided by GitHub. Please report all bugs by creating a new issue. You can ask questions by creating a new issue with the question tag.

Installation
============

Vaxrank can be installed using `pip <https://packaging.python.org/installing/#use-pip-for-installing>`_:

.. code-block:: bash

    pip install vaxrank

**Requirements:** Python 3.9+

Note: to generate PDF reports, you first need to install `wkhtmltopdf <http://wkhtmltopdf.org/>`_, which you can do (on OS X) like so:

.. code-block:: bash

    brew install Caskroom/cask/wkhtmltopdf

Vaxrank uses `PyEnsembl <https://github.com/openvax/pyensembl>`_ for accessing information about the reference genome. You must install an Ensembl release corresponding to the reference genome associated with the mutations provided to Vaxrank.

The latest supported release for GRCh38 is Ensembl 93:

.. code-block:: bash

    pyensembl install --release 93 --species human

The latest release for GRCh37 is Ensembl 75:

.. code-block:: bash

    pyensembl install --release 75 --species human


Running Vaxrank
===============

Basic Usage
-----------

.. code-block:: bash

    vaxrank \
        --vcf somatic-variants.vcf \
        --bam tumor-rna.bam \
        --mhc-predictor netmhc \
        --mhc-alleles A*02:01,A*02:03 \
        --mhc-epitope-lengths 8 \
        --padding-around-mutation 5 \
        --vaccine-peptide-length 25 \
        --output-ascii-report vaccine-peptides-report.txt

This tells Vaxrank to:

- consider each variant from the input VCF file against the RNA evidence in the input BAM file;
- predict MHC binding of each resulting mutant protein sequence using the NetMHC prediction algorithm with the A*02:01 and A*02:03 MHC alleles, evaluating sequences of length 8 for purposes of MHC binding prediction;
- choose protein vaccine candidates, each composed of 25 amino acids; and
- generate a report written to vaccine-peptides-report.txt, containing the top ranked variants with their associated vaccine proteins.

Using a YAML Configuration File
-------------------------------

You can specify common parameters in a YAML configuration file:

.. code-block:: bash

    vaxrank --config my_config.yaml --vcf variants.vcf --bam tumor.bam

Example configuration file:

.. code-block:: yaml

    epitope_config:
      min_epitope_score: 0.001
      logistic_epitope_score_midpoint: 350.0
      logistic_epitope_score_width: 150.0

    vaccine_config:
      vaccine_peptide_length: 25
      padding_around_mutation: 5
      max_vaccine_peptides_per_variant: 1

CLI arguments override values from the config file.


Variant Parameters
------------------
Vaxrank starts with a set of candidate genomic variants and considers each for inclusion in the vaccine. There are several ways to specify a set of variants for Vaxrank to consider:

--vcf VCF_FILE
    Genomic variants in `VCF <https://samtools.github.io/hts-specs/VCFv4.2.pdf>`_ format.
--maf MAF_FILE
    Genomic variants in `MAF <https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification>`_ format.
--json-variants JSON_VARIANTS
    Path to Varcode.VariantCollection object serialized as a JSON
    file. To learn more about Varcode, see `docs <https://github.com/openvax/varcode>`_.

MHC Prediction Parameters
-------------------------

Vaxrank uses a patient's HLA type information to predict which of the candidate vaccine peptides are most likely to be seen and targeted by the patient's immune system.

--mhc-alleles-file MHC_ALLELES_FILE
  File with one HLA allele per line
--mhc-alleles MHC_ALLELES
  Comma-separate or space-separated list of MHC alleles, e.g. "HLA-A*02:01,HLA-A*02:03".
--mhc-peptide-lengths MHC_PEPTIDE_LENGTHS
  Comma-separated list of epitope lengths to consider for MHC binding prediction, e.g. "8,9,10,11".
--mhc-predictor MHC_PREDICTOR
  MHC predictor to use. MHCFlurry is an open-source predictor installed by default.

RNA Parameters
--------------

Vaxrank uses input tumor RNA data to see whether the input somatic variants are sufficiently expressed.

--bam BAM
  BAM file containing tumor RNA reads.
--min-alt-rna-reads MIN_ALT_RNA_READS
  Minimum number of RNA reads supporting the variant allele. Default: 2.
--min-variant-sequence-coverage MIN_VARIANT_SEQUENCE_COVERAGE
  Minimum number of reads supporting a variant sequence. Default: 2.

Vaccine Peptide Parameters
--------------------------

--vaccine-peptide-length VACCINE_PEPTIDE_LENGTH
  Number of amino acids in the resulting vaccine peptides. Default: 25.
--padding-around-mutation PADDING_AROUND_MUTATION
  Number of off-center windows around the mutation to consider. Default: 0.
--min-epitope-score MIN_EPITOPE_SCORE
  Ignore epitopes whose normalized score falls below this threshold. Default: 0.001.

Output Formats
--------------

Vaxrank can generate many types of outputs:

--output-ascii-report OUTPUT_ASCII_REPORT
    Path to ASCII vaccine peptide report
--output-html-report OUTPUT_HTML_REPORT
    Path to HTML vaccine peptide report
--output-pdf-report OUTPUT_PDF_REPORT
    Path to PDF vaccine peptide report
--output-xlsx-report OUTPUT_XLSX_REPORT
    Path to XLSX vaccine peptide report worksheet
--output-json-file OUTPUT_JSON_FILE
    Path to JSON vaccine peptide data


Features
========

Reference Proteome Filtering
----------------------------

Vaxrank filters out peptides that exist in the reference proteome to focus on truly novel mutant sequences. This uses a set-based kmer index for O(1) membership testing. The index is built once and cached locally for subsequent runs.

Cancer Hotspot Annotation
-------------------------

Vaxrank annotates variants that occur at known cancer mutation hotspots using bundled data from `cancerhotspots.org <https://www.cancerhotspots.org/>`_ (Chang et al. 2016, 2017). This helps identify clinically relevant mutations. The hotspot data includes ~2,700 recurrently mutated positions across cancer types.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
