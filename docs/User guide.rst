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
 * The length of a sequence in which a TFBM (transcription factor binding motif) is searched. E.g. promoter_length=1000 means that the sequence 1000 bp upstream of a gene is taken.


* **maximum time for a TFMB search** parameter

 * *motifs_max_time*
 * optional, default: 180 s
 * Maximum time in seconds to search TFBM for individual TF using MEME Suite. The recommended time by MEME Suite is 180 seconds, but it may take longer for large genomes. If the search is terminated due to timeout, a message will be displayed and the parameter should be extended.

* **Gene Regulatory Network** file

 * *GRN_input*
 * Adjacency matrix NxN; N = number of genes; TF = transcription factor (regulator); TG = target (regulated) gene
 * CSV file format
 * examples: Table 2 or `"data" directory <https://github.com/JanaMus/Augusta/tree/master/data>`_ on GitHub

.. list-table:: Table 2: GRN example
   :widths: 20 20 20 20
   :header-rows: 1
   :stub-columns: 1
   :align: center

   * -
     - TF 1
     - TF 2
     - TF 3
   * - TG 1
     - 0
     - 1
     - -1
   * - TG 2
     - 1
     - 0
     - 0
   * - TG 3
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
`RNASeq_to_BN` is the main function for inferring both networks (GRN and BN) using RNA-Seq dataset as an input.

Usage:

.. code-block:: python

   >>> Augusta.RNASeq_to_BN(count_table_input, promoter_length, genbank_file_input, normalization_type, motifs_max_time)


*Note: count_table_input is the only indispensable input, the remaining ones are optional.*
*Not providing GenBank file results in only inferring GRN by computing mutual information. Further steps such as count table normalization, GRN validation (TFBM and DBs search), and Cell Collective DB search would be skipped.*


GRN inference using RNA-Seq
""""""""""""""""""""""""""""
`RNASeq_to_GRN` is the function for inferring only the Gene Regulatory Network using RNA-Seq dataset as an input.

Usage:

.. code-block:: python

   >>> Augusta.RNASeq_to_GRN(count_table_input, promoter_length, genbank_file_input, normalization_type, motifs_max_time)

*Note: count_table_input is the only indispensable input, the remaining ones are optional.*
*Not providing GenBank file results in only inferring GRN by computing mutual information. Further steps such as count table normalization, GRN validation (TFBM and DBs search) would be skipped.*


BN inference using GRN
"""""""""""""""""""""""
`GRN_to_BN` is the function for inferring the Boolean Network (BN) using the Gene Regulatory Network (GRN) file as an input.

Usage:

.. code-block:: python

   >>> Augusta.GRN_to_BN(GRN_input, promoter_length, genbank_file_input, add_dbs_info)


*Note: GRN_input is the only indispensable input, the remaining ones are optional. Not providing GenBank file and/or not setting add_dbs_info only results in a GRN to BN conversion. CC DB would not be searched.*



Outputs
^^^^^^^^
All output files are stored in generated "output" directory.
During motif search, the temporary file "temporary_coreg_seq.fasta" is generated and deleted at the end of the verification process.

* Gene Regulatory Network

 * adjancency matrix in CSV file format
 * rows: TFs (trascription factors / regulators), cols: TGs (target / regulated genes)
 * "GRN.csv"

* Boolean Network

 * SBML-qual file format
 * "BN.sbml"
 * *Note: GRN is primarily converted to the temporary file "BN.txt". If memory is sufficient, the "BN.txt" is converted to "BN.sbml". Otherwise, "BN.txt" is the final output.*

* motifs

 * all TFBM discovered in the genome assigned to their transcription factor
 * Stockholm file format
 * "discovered_motifs.sto"

* genes interactions

 * all interactions searched across databases stored as "DBs_interactions_list.csv"
 * uncertain interactions stored as "DBs_interactions_uncertain.csv" (i.e. the same gene pair has different interaction type in different DBs)
