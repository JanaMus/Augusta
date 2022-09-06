Augusta
==========

Python package: From RNA-Seq to the Boolean Network through the Gene Regulatory Network

Documentation and tutorials are available at `augusta.readthedocs.io <https://augusta.readthedocs.io>`_

Quick Guide
----------------

**Installation:**

.. code-block::

   $ pip3 install Augusta

Dependencies:
- Python 3, up to 3.8
- Docker

**Usage:**

.. code-block:: 

   $ python3
   >>> import Augusta
   
GRN and BN inference using RNA-Seq:

.. code-block:: 

   >>> Augusta.RNASeq_to_SBML(count_table_input, promoter_length, genbank_file_input, normalization_type)

GRN inference using RNA-Seq:

.. code-block:: 

   >>> GRN = Augusta.RNASeq_to_GRN(count_table_input, promoter_length, genbank_file_input, normalization_type)


BN inference using GRN:

.. code-block:: 

   >>> Augusta.GRNtoBN(GRN_input, promoter_length, genbank_file_input, add_dbs_info)


