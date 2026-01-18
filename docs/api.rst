API Reference
=============

This section documents the public API of Vaxrank.

Configuration Classes
---------------------

EpitopeConfig
^^^^^^^^^^^^^

.. automodule:: vaxrank.epitope_config
   :members:
   :undoc-members:
   :show-inheritance:

VaccineConfig
^^^^^^^^^^^^^

.. automodule:: vaxrank.vaccine_config
   :members:
   :undoc-members:
   :show-inheritance:


Reference Proteome
------------------

The reference proteome module provides functionality for checking if peptide
sequences exist in a reference proteome. This is used to filter out peptides
that are not truly novel (i.e., they exist in the normal/wildtype proteome).

.. automodule:: vaxrank.reference_proteome
   :members:
   :undoc-members:
   :show-inheritance:


Cancer Hotspots
---------------

The cancer hotspots module provides lookup functionality for known cancer
mutation hotspots based on data from cancerhotspots.org (Chang et al. 2016, 2017).

.. automodule:: vaxrank.cancer_hotspots
   :members:
   :undoc-members:
   :show-inheritance:


Core Logic
----------

The core logic module contains the main vaccine peptide selection algorithm.

.. automodule:: vaxrank.core_logic
   :members:
   :undoc-members:
   :show-inheritance:


Epitope Logic
-------------

The epitope logic module handles epitope scoring and filtering.

.. automodule:: vaxrank.epitope_logic
   :members:
   :undoc-members:
   :show-inheritance:


Data Classes
------------

MutantProteinFragment
^^^^^^^^^^^^^^^^^^^^^

.. automodule:: vaxrank.mutant_protein_fragment
   :members:
   :undoc-members:
   :show-inheritance:

VaccinePeptide
^^^^^^^^^^^^^^

.. automodule:: vaxrank.vaccine_peptide
   :members:
   :undoc-members:
   :show-inheritance:


CLI Configuration
-----------------

Epitope Config Args
^^^^^^^^^^^^^^^^^^^

.. automodule:: vaxrank.cli.epitope_config_args
   :members:
   :undoc-members:
   :show-inheritance:

Vaccine Config Args
^^^^^^^^^^^^^^^^^^^

.. automodule:: vaxrank.cli.vaccine_config_args
   :members:
   :undoc-members:
   :show-inheritance:
