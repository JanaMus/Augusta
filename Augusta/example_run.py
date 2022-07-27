from Augusta import *

## import data
parent_directory = os.path.abspath('..')
# RNA-Seq count table
count_table_input = parent_directory + '\data\Cbeijerinckii.csv'
# promoter length - optional (default is 1000)
promoter_length = 1000
# GenBank file
genbank_file_input = parent_directory + '\data\Cbeijerinckii.gb'
# normalization - optional (default is 'None')
normalization_type = 'RPKM'

## run:
GRN = RNASeq_to_GRN(count_table_input, promoter_length, genbank_file_input, normalization_type) # all data
#GRN = RNASeq_to_GRN(count_table_input) # only indispesable data

#GRN_input = 'output\GRN.csv'
#GRN_to_BN(GRN_input, promoter_length, genbank_file_input, add_dbs_info = 1) # all data
#GRN_to_BN(GRN, genbank_file_input, add_dbs_info = 0) # test input GRN variable ###
#GRN_to_BN(GRN_input, add_dbs_info = 1) # test no genbank ###
#GRN_to_BN(GRN_input) # only indispesable data

#RNASeq_to_SBML(count_table_input, promoter_length, genbank_file_input, normalization_type)  # all data
#RNASeq_to_SBML(count_table_input, genbank_file_input) # test ###
#RNASeq_to_SBML(count_table_input) # only indispesable data
