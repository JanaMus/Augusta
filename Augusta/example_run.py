from Augusta import *

## import data
parent_directory = os.path.abspath('..')

# RNA-Seq count table
count_table_input = parent_directory + '\data\Cbeijerinckii.csv' # Windows
#count_table_input = parent_directory + '/data/Cbeijerinckii.csv' # Linux

# promoter length - optional (default is 1000)
promoter_length = 1000

# GenBank file
genbank_file_input = parent_directory + '\data\Cbeijerinckii.gb'# Windows
#genbank_file_input = parent_directory + '/data/Cbeijerinckii.gb'# Linux

# normalization - optional (default is 'None')
normalization_type = 'RPKM'
print(count_table_input)

# GRN
GRN_input = 'output\GRN_example.csv' # Windows
#GRN_input = 'output/GRN_example.csv' # Linux

## RUN (uncomment selected line):
#GRN = RNASeq_to_GRN(count_table_input, promoter_length, genbank_file_input, normalization_type) # all data
#GRN = RNASeq_to_GRN(count_table_input) # only indispesable data

#GRN_to_BN(GRN_input, promoter_length, genbank_file_input, add_dbs_info = 1) # all data
#GRN_to_BN(GRN_input) # only indispesable data

#RNASeq_to_SBML(count_table_input, promoter_length, genbank_file_input, normalization_type)  # all data
#RNASeq_to_SBML(count_table_input) # only indispesable data
