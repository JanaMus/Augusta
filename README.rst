Augusta
==========

Python package: From RNA-Seq to the Boolean Network through the Gene Regulatory Network

Documentation and tutorials are available at `augusta.readthedocs.io <https://augusta.readthedocs.io>`_.

Credits
----------------
The Augusta project is based on research detailed in the following paper. Please cite this paper when using or referencing our work:

Augusta: From RNA‐Seq to gene regulatory networks and Boolean models. Jana Musilova, Zdenek Vafek, Bhanwar Lal Puniya, Ralf Zimmer, Tomas Helikar, and Karel Sedlar. *Computational and Structural Biotechnology Journal*, 2024. DOI: `10.1016/j.csbj.2024.01.013 <https://doi.org/10.1016/j.csbj.2024.01.013>`_.


Contributors
----------------
- Jana Musilova, musilovaj22@gmail.com
- Zdenek Vafek
- Karel Sedlar, sedlar@vut.cz


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

   >>> Augusta.RNASeq_to_BN(count_table_input = 'MyCT_file.csv', promoter_length = My_number, genbank_file_input = 'MyGB_file.gb', normalization_type = 'My_string', motifs_max_time = My_seconds)

GRN inference using RNA-Seq:

.. code-block:: 

   >>> Augusta.RNASeq_to_GRN(count_table_input = 'MyCT_file.csv', promoter_length = My_number, genbank_file_input = 'MyGB_file.gb', normalization_type = 'My_string', motifs_max_time = My_seconds)


BN inference using GRN:

.. code-block:: 

   >>> Augusta.GRN_to_BN(GRN_input = 'MyGRN_file.csv', promoter_length = My_number, genbank_file_input = 'MyGB_file.gb', add_dbs_info = 'My_string')


GRN refinement:

.. code-block:: 

   >>> Augusta.refineGRN(GRN_input = 'MyGRN_file.csv', genbank_file_input = 'MyGB_file.gb', count_table_input = 'MyCT_file.csv', promoter_length = My_number, motifs_max_time = My_seconds)

   



