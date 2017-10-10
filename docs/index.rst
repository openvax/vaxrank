.. vaxrank documentation master file, created by
   sphinx-quickstart on Tue Oct 10 16:59:03 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Getting Started With Vaxrank
============================

Overview
--------
Vaxrank is a tool for selecting mutated peptides for use in personalized therapeutic cancer vaccination. Vaxrank determines which peptides should be used in a vaccine from tumor-specific somatic mutations, tumor RNA sequencing data, and a patient's HLA type. Additionally, Vaxrank considers surrounding non-mutated residues in a peptide to prioritize vaccine peptide candidates and improve the odds of successful synthesis.

Vaxrank is being actively developed at the Icahn School of Medicine at Mount Sinai.

Questions, Bug Reporting, and Issue Tracking
--------------------------------------------
Questions, bug reporting and issue tracking are provided by GitHub. Please report all bugs by creating a new issue. You can ask questions by creating a new issue with the question tag.

Installation
============

Vaxrank can be installed using `pip <https://packaging.python.org/installing/#use-pip-for-installing>`_:

.. code-block:: bash

    pip install vaxrank

Note: to generate PDF reports, you first need to install `wkhtmltopdf <http://wkhtmltopdf.org/>`_, which you can do (on OS X) like so:

.. code-block:: bash

    brew install Caskroom/cask/wkhtmltopdf

Vaxrank uses `PyEnsembl <https://github.com/hammerlab/pyensembl>`_ for accessing information about the reference genome. You must install an Ensembl release corresponding to the reference genome associated with the mutations provided to Vaxrank.

The latest release for GRCh38 is Ensembl 87:

.. code-block:: bash

    pyensembl install --release 87 --species human

The latest release for GRCh37 is Ensembl 75:

.. code-block:: bash

    pyensembl install --release 75 --species human


Running Vaxrank
===============

Basic Vaxrank usage involves these parameters:

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

For a complete description of parameters supported by Vaxrank, keep on reading.


Variant Parameters
------------------
Vaxrank starts with a set of candidate genomic variants and considers each for inclusion in the vaccine. There are several ways to specify a set of variants for Vaxrank to consider:

--vcf VCF_FILE
    Genomic variants in `VCF <https://samtools.github.io/hts-specs/VCFv4.2.pdf>`_ format.
--maf MAF_FILE
    Genomic variants in `MAF <https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification>`_ format.
--json-variants JSON_VARIANTS
    Path to Varcode.VariantCollection object serialized as a JSON
    file. To learn more about Varcode, see `docs <https://github.com/hammerlab/varcode>`_. 

MHC Prediction Parameters
-------------------------

Vaxrank uses a patient's HLA type information to predict which of the candidate vaccine peptides are most likely to be seen and targeted by the patient's immune system. The MHC alleles can be passed in either in a file or as a comma-separated list of inputs.

--mhc-alleles-file MHC_ALLELES_FILE
  File with one HLA allele per line
--mhc-alleles MHC_ALLELES
  Comma-separate or space-separated list of MHC alleles, e.g. "HLA-A*02:01,HLA-A*02:03".
--mhc-peptide-lengths MHC_PEPTIDE_LENGTHS
  Comma-separated list of epitope lengths to consider for MHC binding prediction, e.g. "8,9,10,11". This can also take a range of values, e.g. "8-11".

In addition, the user can specify different MHC binding predictors for Vaxrank to use:

--mhc-predictor MHC_PREDICTOR
  MHC predictor to use. MHCFlurry is an open-source predictor installed by default. Note that to use NetMHC predictors, you need to have locally installed the NetMHC suite software, with binaries like NetMHCpan as executable files on your path. See a list of all supported predictors `here <https://github.com/hammerlab/mhctools>`_.

RNA Parameters
--------------

Vaxrank uses input tumor RNA data to see whether the input somatic variants are sufficiently expressed. 

--bam BAM
  BAM file containing tumor RNA reads.

Each variant's effect on a resulting protein is predicted and matched against what we see in the input RNA. There are many options available to the power user, but the only actual required argument is the location of the tumor RNA BAM; all values listed below come with reasonable defaults.

--min-alt-rna-reads MIN_ALT_RNA_READS
  Minimum number of RNA reads supporting the variant allele. Default: 2.
--min-variant-sequence-coverage MIN_VARIANT_SEQUENCE_COVERAGE
  Minimum number of reads supporting a variant sequence. Variant sequences will be trimmed to positions supported by at least this number of RNA reads. Default: 2.
--disable-variant-sequence-assembly
  By default, variant cDNA sequences are assembled from overlapping reads. Include this argument to disable the assembly behavior.
--protein-sequence-length
  Vaxrank will try to translate protein sequences of this length, though sometimes the resulting sequence may be shorter (depending on the RNA data, presence of stop codons, etc.). Default: 20.
--max-reference-transcript-mismatches MAX_REFERENCE_TRANSCRIPT_MISMATCHES
  Maximum number of mismatches between the variant sequence being constructed and the reference sequence before the variant sequence gets dropped from consideration. Default: 2.
--include-mismatches-after-variant
  By default, only mismatches that occur before the actual variant locus count against --max-reference-transcript-mismatches. Set this value to True if you also want to count mismatches after the variant locus towards the total. Default: false.
--min-transcript-prefix-length MIN_TRANSCRIPT_PREFIX_LENGTH
  Number of nucleotides before the variant we try to match against a reference transcript. Default: 10.
--min-mapping-quality MIN_MAPPING_QUALITY
  Minimum MAPQ value to allow for a read. Default: 1.
--use-duplicate-reads
  Use a read even if it's been marked as a duplicate. Default: false.
--drop-secondary-alignments
  If true, Vaxrank will use a read even at a location that isn't its primary alignment. Default: false.

Vaccine Peptide Parameters
--------------------------
There are some more options to specify the desired characteristics of the output vaccine peptides, which will contain shorter sequences that contain the mutation and are predicted to be strong MHC binders.

--vaccine-peptide-length VACCINE_PEPTIDE_LENGTH
  Number of amino acids in the resulting vaccine peptides. Default: 25.
--padding-around-mutation PADDING_AROUND_MUTATION
  Number of off-center windows around the mutation to consider as vaccine peptides. Default: 0.
--min-epitope-score MIN_EPITOPE_SCORE
  Ignore epitopes whose normalized score falls below this threshold. Default: 0.001. 

Output Parameters
-----------------

By default, the report will contain all high-confidence vaccine peptides, but the report can be made more restrictive using the following parameters:

--max-vaccine-peptides-per-mutation MAX_VACCINE_PEPTIDES_PER_MUTATION
                        Number of vaccine peptides to generate for each
                        mutation
--max-mutations-in-report MAX_MUTATIONS_IN_REPORT
                        Number of mutations to report

Output Formats
^^^^^^^^^^^^^^

Vaxrank can generate many types of outputs. The most basic output is an ASCII-formatted report, listing each high-scoring variant and its associated vaccine peptides. However, the user can also generate a PDF report and two types of Excel reports.

Options related to report generation:
  --output-ascii-report OUTPUT_ASCII_REPORT
                        Path to ASCII vaccine peptide report
  --output-html-report OUTPUT_HTML_REPORT
                        Path to HTML vaccine peptide report
  --output-pdf-report OUTPUT_PDF_REPORT
                        Path to PDF vaccine peptide report
  --output-xlsx-report OUTPUT_XLSX_REPORT
                        Path to XLSX vaccine peptide report worksheet, one
                        sheet per variant. This is meant for use by the
                        vaccine manufacturer.
  --output-neoepitope-report OUTPUT_NEOEPITOPE_REPORT
                        Path to XLSX neoepitope report, containing information
                        focusing on short peptide sequences.

Vaxrank can also output all variants and vaccine sequences in a JSON file, which can be used for further programmatic processing if necessary. The file output location should be specified by:

--output-json-file OUTPUT_JSON_FILE
                    Path to JSON vaccine peptide data
