# Augusta
Python package: From RNA-Seq to the Boolean Network through the Gene Regulatory Network

## Introduction

## Installation
From PyPi repository:
```
pip install Augusta
```
From GitHub:
```
git clone https://github.com/JanaMus/Augusta.git
cd Augusta
python setup.py install
```
Dependencies:
docker

## Usage
### Generate Boolean network from RNA-Seq data
Inputs:
count table file (matrix NxM; N = genes, M =  time points; see examples in data folder for correct formating)
promoter length (optional, default: 1000)
GenBank file (optional but several parts are skipped if not provided)
normalization type (optional, default: None, options: RPKM, CPM, TPM)

Example:
```
import Augusta
RNASeq_to_SBML(count_table_input = 'Ecoli_DREAM4.csv', promoter_length = 500, genbank_file_input = 'Ecoli.gb', normalization_type = 'RPKM')
```
Example data files are available in the data folder.
Note: run time for C. Beijerinckii example data approximates 2 days.

### Generate Gene Regulatory network from RNA-Seq data
Inputs:
count table file (matrix NxM; N = genes, M =  time points; see examples in data folder for correct formating)
promoter length (optional, default: 1000)
GenBank file (optional but several parts are skipped if not provided)
normalization type (optional, default: None, options: RPKM, CPM, TPM)

Example:
```
import Augusta
GRN = RNASeq_to_GRN(count_table_input = 'Ecoli_DREAM4.csv', promoter_length = 500, genbank_file_input = 'Ecoli.gb', normalization_type = 'RPKM')
```
Example data files are available in the data folder.

### Generate Boolean network from Gene Regulatory network
Inputs:
Gene regulatory network file (Adjacency matrix NxN; N = genes; see examples in data folder for correct formating)
promoter length (optional, default: 1000; needed only if databases info is added to verify GRN)
GenBank file (optional but several parts are skipped if not provided)
Databases info add (optional, default: None; to search databases and use found information to verify GRN, use inputs 'yes' or 1)

Example:
```
import Augusta
GRN_toBN(GRN_input = 'GRN_example.csv', promoter_length = 500, genbank_file_input = 'Ecoli.gb', add_dbs_info = 1)
```
Example data files are available in the data folder.
