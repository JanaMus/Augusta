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
Dependencies:<br />
docker

## Usage
### Generate Boolean network from RNA-Seq data
**Inputs:**<br />
- count table file (matrix NxM; N = genes, M =  time points; see examples in data folder for correct formating)<br />
- promoter length (optional, default: 1000)<br />
- GenBank file (optional but several parts are skipped if not provided)<br />
- normalization type (optional, default: None, options: RPKM, CPM, TPM)<br />
Please note that count table must be in a normalized form in order to infer GRN (using MI computation) correctly.

**Example:**
```
import Augusta
RNASeq_to_SBML(count_table_input = 'Ecoli_DREAM4.csv', promoter_length = 500, genbank_file_input = 'Ecoli.gb', normalization_type = 'RPKM')
```
Example data files are available in the data folder.<br />
*Note: run time for C. Beijerinckii example data approximates 2 days.*

### Generate Gene Regulatory network from RNA-Seq data
**Inputs:**<br />
- count table file (matrix NxM; N = genes, M =  time points; see examples in data folder for correct formating)<br />
- promoter length (optional, default: 1000)<br />
- GenBank file (optional but several parts are skipped if not provided)<br />
- normalization type (optional, default: None, options: RPKM, CPM, TPM)<br />
Please note that count table must be in a normalized form in order to infer GRN (using MI computation) correctly.

**Example:**
```
import Augusta
GRN = RNASeq_to_GRN(count_table_input = 'Ecoli_DREAM4.csv', promoter_length = 500, genbank_file_input = 'Ecoli.gb', normalization_type = 'RPKM')
```
Example data files are available in the data folder.

### Generate Boolean network from Gene Regulatory network
**Inputs:**<br />
- Gene regulatory network file (Adjacency matrix NxN; N = genes; see examples in data folder for correct formating)<br />
- promoter length (optional, default: 1000; needed only if databases info is added to verify GRN)<br />
- GenBank file (optional but several parts are skipped if not provided)<br />
- Databases info add (optional, default: None; to search databases and use found information to verify GRN, use inputs 'yes' or 1)<br />

**Example:**
```
import Augusta
GRN_toBN(GRN_input = 'GRN_example.csv', promoter_length = 500, genbank_file_input = 'Ecoli.gb', add_dbs_info = 1)
```
Example data files are available in the data folder.
