Configuration
=============

Vaxrank can be configured through command-line arguments or a YAML configuration file.
Configuration objects are defined using `msgspec <https://jcristharif.com/msgspec/>`_ Struct classes for type safety and serialization.

YAML Configuration File
-----------------------

Use the ``--config`` argument to specify a YAML configuration file:

.. code-block:: bash

    vaxrank --config my_config.yaml --vcf variants.vcf --bam tumor.bam

CLI arguments always override values from the configuration file.

Example Configuration
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

    # Epitope scoring and filtering parameters
    epitope_config:
      # IC50 value (nM) at which epitope score equals 0.5
      logistic_epitope_score_midpoint: 350.0

      # Width parameter for logistic scoring function
      logistic_epitope_score_width: 150.0

      # Minimum normalized score threshold
      min_epitope_score: 0.00001

      # Maximum IC50 to consider (nM)
      binding_affinity_cutoff: 5000.0

    # Vaccine peptide assembly parameters
    vaccine_config:
      # Length of vaccine peptides (amino acids)
      vaccine_peptide_length: 25

      # Off-center window positions to consider
      padding_around_mutation: 5

      # Maximum peptides per variant
      max_vaccine_peptides_per_variant: 1

      # Epitopes to retain per variant
      num_mutant_epitopes_to_keep: 1000


EpitopeConfig
-------------

The ``EpitopeConfig`` class controls epitope scoring and filtering.

.. autoclass:: vaxrank.epitope_config.EpitopeConfig
   :members:
   :undoc-members:
   :show-inheritance:

Epitope Scoring Algorithm
^^^^^^^^^^^^^^^^^^^^^^^^^

Epitopes are scored using a logistic function that transforms IC50 binding
affinity values into a normalized score between 0 and 1:

.. math::

    rescaled = \\frac{ic50 - midpoint}{width}

    score = \\frac{1}{1 + e^{rescaled}} \\cdot \\frac{1}{normalizer}

Where normalizer ensures score ≈ 1.0 when IC50 ≈ 0:

.. math::

    normalizer = \\frac{1}{1 + e^{-midpoint / width}}

Parameters:

- **midpoint**: IC50 value (nM) around which scores transition (default: 350 nM)
- **width**: Controls steepness of the scoring curve (default: 150)

Lower IC50 values indicate stronger MHC binding and result in higher scores.


VaccineConfig
-------------

The ``VaccineConfig`` class controls vaccine peptide assembly.

.. autoclass:: vaxrank.vaccine_config.VaccineConfig
   :members:
   :undoc-members:
   :show-inheritance:


Programmatic Configuration
--------------------------

Configuration objects can be created programmatically:

.. code-block:: python

    from vaxrank.epitope_config import EpitopeConfig
    from vaxrank.vaccine_config import VaccineConfig

    # Create with defaults
    epitope_config = EpitopeConfig()

    # Create with custom values
    strict_config = EpitopeConfig(
        min_epitope_score=0.01,
        binding_affinity_cutoff=1000.0
    )

    vaccine_config = VaccineConfig(
        vaccine_peptide_length=30,
        max_vaccine_peptides_per_variant=3
    )

Configuration objects are immutable (frozen) and can be safely shared.
