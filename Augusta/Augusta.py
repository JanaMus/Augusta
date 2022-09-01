from .core import *

def RNASeq_to_GRN(count_table_input, promoter_length = 1000, genbank_file_input = None, normalization_type = None):
   count_table = import_CountTable(count_table_input)  # data import
   gene_names = count_table.index
   output_directory = os.path.join(os.getcwd(), r'output')  # data export
   if not os.path.exists(output_directory):
      os.makedirs(output_directory)
   if genbank_file_input:
      gene_lengths, gene_promoters, organism, taxon, genes_IDs = genbank_process(genbank_file_input, gene_names, promoter_length)
      if normalization_type:  # count table normalization
         normalized_count_table = normalize_CountTable(count_table, gene_lengths, normalization_type)
      else:
         normalized_count_table = count_table
      GRNmi, count_table_differences = RNASeq_to_GRNmi(normalized_count_table)  # GRN inference using MI
      GRNmotifs = GRNmi_to_GRNmotifs(GRNmi, gene_promoters, count_table_differences)  # GRN verification using motif search
      GRNdb = get_DBs_data(organism, genes_IDs, taxon, GRNmotifs)  # GRN verification using databases
      if GRNdb is None:
         print('No databases interaction data available; GRN verification skipped.')
         GRNdb = GRNmotifs.copy()
      export_GRN(GRNdb)
      return GRNdb
   elif genbank_file_input is None:
      if normalization_type:
         print('Count table normalization not available - GenBank missing; skipped.')
      GRNmi, count_table_differences = RNASeq_to_GRNmi(count_table)
      export_GRN(GRNmi)
      return GRNmi


def GRN_to_BN(GRN_input, promoter_length = 1000, genbank_file_input = None, add_dbs_info = None):
   try:
      GRN_imported = import_GRN(GRN_input) # GRN input from file
   except TypeError:
      GRN_imported = GRN_input.copy() # GRN input from variable

   output_directory = os.path.join(os.getcwd(), r'output') # data export
   if not os.path.exists(output_directory):
      os.makedirs(output_directory)

   if (add_dbs_info == 'yes' or add_dbs_info == 'Yes' or add_dbs_info == 1):
      if genbank_file_input:
         gene_names = GRN_imported.index
         gene_lengths, gene_promoters, organism, taxon, genes_IDs = genbank_process(genbank_file_input, gene_names, promoter_length)
         GRNdb = get_DBs_data(organism, genes_IDs, taxon, GRN_imported)
         if GRNdb is None:
            print('No databases interaction data available; GRN verification skipped.')
            GRNdb = GRN_imported.copy()
         export_GRN(GRNdb)
         GRNcc, CC_conditions = get_CC_data(genes_IDs, GRNdb)  # GRN verification using Cell Colective db
         GRN_to_SBML(GRNcc, CC_conditions) # GRN to BN
      else:
         print('Databases info add not available - GenBank missing; skipped.')
         GRN_to_SBML(GRN_imported) # GRN to BN
   else:
      GRN_to_SBML(GRN_imported) # GRN to BN


def RNASeq_to_SBML(count_table_input, promoter_length=1000, genbank_file_input=None, normalization_type=None):
   GRN = RNASeq_to_GRN(count_table_input, promoter_length, genbank_file_input, normalization_type)
   GRN_to_BN(GRN, promoter_length, genbank_file_input, add_dbs_info=None)
