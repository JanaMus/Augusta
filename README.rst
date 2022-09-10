Augusta
==========

Python package: From RNA-Seq to the Boolean Network through the Gene Regulatory Network

Documentation and tutorials are available at `augusta.readthedocs.io <https://augusta.readthedocs.io>`_.

Quick Guide
----------------

Dependencies:

- Python 3, versions 3.7 and 3.8
- Docker

**Installation:**

We highly recomment installing and using Augusta in a virtual environment.

.. code-block::

   $ conda create -n Augusta_venv python=3.7 anaconda
   $ conda activate Augusta_venv
   

.. code-block::

   $ pip install Augusta


**Usage:**

See `Inputs <https://augusta.readthedocs.io/en/latest/User%20guide.html>`_ for details about input files and variables.

.. code-block:: 

   $ python
   >>> import Augusta
   
GRN and BN inference using RNA-Seq:

.. code-block:: 

   >>> Augusta.RNASeq_to_SBML(count_table_input = 'MyCT_file.csv', promoter_length = My_number, genbank_file_input = 'MyGB_file.gb', normalization_type = 'My_string')

GRN inference using RNA-Seq:

.. code-block:: 

   >>> GRN = Augusta.RNASeq_to_GRN(count_table_input = 'MyCT_file.csv', promoter_length = My_number, genbank_file_input = 'MyGB_file.gb', normalization_type = 'My_string')


BN inference using GRN:

.. code-block:: 

   >>> Augusta.GRNtoBN(GRN_input = 'MyGRN_file.csv', promoter_length = My_number, genbank_file_input = 'MyGB_file.gb', add_dbs_info = 'My_string')
   



