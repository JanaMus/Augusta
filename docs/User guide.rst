User guide
----------

Inputs
^^^^^^
Below is the list of all input files and parameters. See section :ref:`Usage` for particular functions inputs.

* **count table** file

 * *count_table_input*
 * matrix NxM; N = genes, M = time points
 * CSV file format
 * examples: Table 1 or `"data" directory <https://github.com/JanaMus/Augusta/tree/master/data>`_ on GitHub

.. list-table:: Table 1: Count table example
   :widths: 20 20 20 20
   :header-rows: 1
   :stub-columns: 1
   :align: center

   * -
     - Time 1
     - Time 2
     - Time 3
   * - Gene 1
     - 255
     - 1,596
     - 80
   * - Gene 2
     - 112
     - 63
     - 0
   * - Gene 3
     - 56
     - 3,582
     - 27
   * - Gene 4
     - 559
     - 865
     - 91


* **normalization type** parameter

 * *normalization_type*
 * optional, default: None
 * options: RPKM, CPM, TPM
 * Please note that the count table must be in a normalized form in order to infer GRN correctly!

* **GenBank** file

 * *genbank_file_input*
 * optional but several parts are skipped if not provided
 * Gene names in the count table must match locus tags (i.e. "locus_tag") or names (i.e. "name") in the GenBank file!
 * required format: GenBank full (i.e. containing nucleotide sequence in the ORIGIN section)
 * examples: `"data" directory <https://github.com/JanaMus/Augusta/tree/master/data>`_ on GitHub or :ref:`Examples`)


* **promoter length** parameter

 * *promoter_length*
 * optional, default: 1000 bp

* **Gene Regulatory Network** file

 * *GRN_input*
 * Adjacency matrix NxN; N = genes
 * CSV file format
 * examples: Table 2 or `"data" directory <https://github.com/JanaMus/Augusta/tree/master/data>`_ on GitHub

.. list-table:: Table 2: GRN example
   :widths: 20 20 20 20
   :header-rows: 1
   :stub-columns: 1
   :align: center

   * -
     - Gene 1
     - Gene 2
     - Gene 3
   * - Gene 1
     - 0
     - 1
     - -1
   * - Gene 2
     - 1
     - 0
     - 0
   * - Gene 3
     - 1
     - -1
     - 0


* **Add database information** parameter

 * *add_dbs_info*
 * optional, default: None
 * options to add DBs info: 'yes', 'Yes', 1; every other parameter´s value results in skipping Cell Collective DB search


Usage
^^^^^^
Below is the list of Augusta's functions along with the inputs. See :ref:`Examples` for further description and tutorials.

Import Augusta:

.. code-block:: python

   > python3
   >>> import Augusta
   
   
GRN and BN inference using RNA-Seq
""""""""""""""""""""""""""""""""""""""""""""""""""""""""
`RNASeq_to_SBML` is the main function for inferring both networks using RNA-Seq dataset as an input.

Usage:

.. code-block:: python

   >>> Augusta.RNASeq_to_SBML(count_table_input, promoter_length, genbank_file_input, normalization_type)


*Note: count_table_input is the only indispensable input, the remaining ones are optional.*
*Not providing GenBank file results in only inferring GRN by computing mutual information. Further steps such as verification and BN inference would be skipped.*


GRN inference using RNA-Seq
""""""""""""""""""""""""""""
`RNASeq_to_GRN` is the function for inferring only the Gene Regulatory Network using RNA-Seq dataset as an input.

Usage:

.. code-block:: python

   >>> GRN = Augusta.RNASeq_to_GRN(count_table_input, promoter_length, genbank_file_input, normalization_type)

*Note: count_table_input is the only indispensable input, the remaining ones are optional.*
*Not providing GenBank file results in only inferring GRN by computing mutual information. Further steps such as verification and BN inference would be skipped.*


BN inference using GRN
"""""""""""""""""""""""
`GRN_toBN` is the function for inferring the Boolean Network (BN) using the Gene Regulatory Network (GRN) file as an input.

Usage:

.. code-block:: python

   >>> Augusta.GRNtoBN(GRN_input, promoter_length, genbank_file_input, add_dbs_info)


*Note: GRN_input is the only indispensable input, the remaining ones are optional. Not providing GenBank file and/or not setting add_dbs_info only results in a GRN to BN conversion. Cell Collective database would not be searched.*



Outputs
^^^^^^^^
All output files are stored in generated "output" directory.
During motif search is moreover generated temporary file "temporary_coreg_seq.fasta" which is deleted at the end of the verification process.

* Gene Regulatory Network

 * adjancency matrix in CSV file format
 * "GRN.csv"

* Boolean Network

 * SBML-qual file format
 * "BN.sbml"

* transcription motifs

 * all motifs discovered in the genome assigned to their transcription factor
 * Stockholm file format
 * "discovered_motifs.sto"

* genes interactions

 * all interactions searched across databases stored as "DBs_interactions_list.csv"
 * uncertain interactions stored as "DBs_interactions_uncertain.csv" (i.e. the same gene pair has different interaction type in different DBs)
