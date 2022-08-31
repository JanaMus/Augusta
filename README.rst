Augusta
==========

Python package: From RNA-Seq to the Boolean Network through the Gene Regulatory Network

Documentation and tutorials are available at https://augusta.readthedocs.io.

Quick Guide
----------------

**Installation:**

.. code-block::

   $ pip install Augusta

Dependencies:
docker

**Usage:**

.. code-block:: python

   >>> import Augusta
   >>> RNASeq_to_SBML(count_table_input = 'Ecoli_DREAM4.csv', promoter_length = 1000, genbank_file_input = 'Ecoli.gb', normalization_type = 'TPM')


Example data files are available in the "data" folder.

*Note: run time for C. Beijerinckii example data approximates 2 days.*
