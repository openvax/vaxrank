Getting Started
===============

This guide will help you get started with Vaxrank for personalized cancer vaccine design.

Prerequisites
-------------

Before using Vaxrank, you'll need:

1. **Somatic variants**: A VCF or MAF file containing tumor-specific mutations
2. **Tumor RNA-seq data**: A BAM file with aligned RNA-seq reads from the tumor
3. **HLA typing**: The patient's HLA alleles (e.g., HLA-A*02:01)
4. **Reference genome data**: Ensembl annotation data for the reference genome

Installation
------------

Install Vaxrank using pip:

.. code-block:: bash

    pip install vaxrank

Install reference genome data:

.. code-block:: bash

    # For GRCh38 (human)
    pyensembl install --release 93 --species human

    # For GRCh37 (human)
    pyensembl install --release 75 --species human

For PDF report generation, install wkhtmltopdf:

.. code-block:: bash

    # macOS
    brew install Caskroom/cask/wkhtmltopdf

    # Ubuntu/Debian
    apt-get install wkhtmltopdf


Quick Start
-----------

Run Vaxrank with minimal options:

.. code-block:: bash

    vaxrank \
        --vcf somatic_variants.vcf \
        --bam tumor_rna.bam \
        --mhc-alleles HLA-A*02:01,HLA-B*07:02 \
        --mhc-predictor mhcflurry \
        --output-ascii-report report.txt

This will:

1. Load variants from the VCF file
2. Check RNA evidence for each variant in the BAM file
3. Predict MHC binding for mutant peptides
4. Rank vaccine peptide candidates
5. Generate a text report

Understanding the Output
------------------------

The output report contains:

- **Ranked variants**: Mutations ordered by their vaccine potential
- **Vaccine peptides**: Suggested peptide sequences for each variant
- **Epitope information**: Predicted MHC-binding epitopes within each peptide
- **Expression data**: RNA read support for each variant
- **Database annotations**: Links to cancer hotspot databases if applicable

Example Output
^^^^^^^^^^^^^^

::

    Variant: chr17 g.7577538C>T (TP53 p.R248Q)
    Gene: TP53
    RNA reads: 45 alt / 120 total
    Cancer Hotspot: Yes (https://www.cancerhotspots.org/#/gene/TP53)

    Vaccine Peptide #1:
      Sequence: SSCMGGMNRRPILTIITLEDS
      Length: 21 aa
      Mutation position: 11

      Top Epitopes:
        MGGMNRRPI (HLA-A*02:01): IC50=45nM, score=0.89
        RRPILTII (HLA-B*07:02): IC50=120nM, score=0.72


Next Steps
----------

- See :doc:`configuration` for detailed configuration options
- See :doc:`api` for programmatic usage
- Check the `GitHub repository <https://github.com/openvax/vaxrank>`_ for examples
