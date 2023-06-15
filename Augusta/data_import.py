from pandas import read_table
from numpy import shape
from Bio import SeqIO

### import count table (RNA-Seq) data: matrix MxN (M = gene locus tag / name; N = time)
def import_CountTable(count_table):
    input_matrix = read_table(count_table, index_col=[0], sep = None, engine = 'python')
    if (shape(input_matrix)[1]) < 3:
        print('Expression data error: Not enough time points included. Provide at least 3 time points.')
        return
    elif (shape(input_matrix)[0]) < 2:
        print ('Expression data error: Not enough genes included. Provide at least 2 genes.')
        return
    elif (len(input_matrix.columns) != shape(input_matrix)[1]) or (len(input_matrix.index) != shape(input_matrix)[0]):
        print('Expression data error: Provide names of rows and columns.')
        return
    else:
        print('Count table uploaded.')
        return input_matrix


### import GenBank file
def genbank_process(gb_file, gene_names, promoter_length):
    gene_lengths = [0] * shape(gene_names)[0]
    gene_promoters = []
    qualifier = find_qualifier(gb_file, gene_names)

    for record in SeqIO.parse(gb_file, 'genbank'):
        sequence = record.seq
        indexes = index_genbank_features(record, 'gene', qualifier)
        genes_IDs = dict.fromkeys(gene_names)
        for i in range(0, shape(gene_names)[0]):
            # get gene promoter
            try:
                feature = record.features[indexes[gene_names[i]]]
                gene_start = feature.location._start.position + 1 # (exact gene start - python marks gene start at position - 1)
                gene_end = feature.location._end.position
            except:
                raise Exception('Gene name ', gene_names[i], ' does not correspond to Genbank IDs. Needed gene names in count table: Locus tag or Gene name.')
            gene_lengths[i] = gene_end - gene_start + 1
            promoter_start = gene_start - promoter_length
            promoter_end = gene_start - 1
            # verification if promoter is not outside of the sequence and if it isn´t part of a previous gene (if yes, promoter is shorten)
            if promoter_start <= 0:
                promoter_start = 1
            else:
                feature_prev = record.features[indexes[gene_names[i]]-1]
                try:
                    gene_end_prev = feature_prev.location._end.position
                except AttributeError: # feature has more locations (join ...)
                    gene_end_prev = feature_prev.location.parts[-1]._end.position
                if gene_end_prev >= promoter_start:
                    promoter_start = promoter_start + (gene_end_prev - promoter_start) + 1
            # append gene_promoters with promoter sequence based on genes strand
            promoter = sequence[promoter_start:promoter_end]
            if feature.strand == -1:
                promoter = promoter.reverse_complement()
            gene_promoters.append(promoter)

            # get gene IDs
            for f in feature.qualifiers:
                if f != str(qualifier) and (f == "gene" or f == "locus_tag" or f == "gene_synonym"):
                    if genes_IDs[gene_names[i]]:  # if value exists (is not None)
                        genes_IDs[gene_names[i]].extend(feature.qualifiers[f])
                    else:
                        genes_IDs[gene_names[i]] = feature.qualifiers[f]

        # get organism name
        try: # catch organism name is not defined
            organism = record.annotations['organism']
        except KeyError: # organism line in not defined in GB file
            organism = input('Enter organism as: genus species (e.g. Homo sapiens): ')
        while len(organism.split()) == 0: # organism is not defined in GB file or user added blank value
            organism = input('Organism not defined. Enter organism as: genus species (e.g. Homo sapiens): ')

        # get taxon
        try:
            taxon = record.features[0].qualifiers['db_xref'][0]
            taxon = taxon.replace('taxon:', '')
        except KeyError:
            taxon = input('Enter taxon (e.g. 9606): ') # taxon line is not defined in GB file
        while True: # check if taxon was defined correctly (using numbers)
            try:
                taxon_int = int(taxon)
            except ValueError:
                taxon = input('Taxon not defined. Enter taxon (e.g. 9606): ')
            else:
                break
    print('GenBank uploaded.')
    return gene_lengths, gene_promoters, organism, taxon_int, genes_IDs


# search type of qualifier in input GRN (gene / locus tag)
def find_qualifier(gb_file, gene_names):
    qualifier = None
    for record in SeqIO.parse(gb_file, 'genbank'):
        for name in gene_names:
            for f in record.features:
                if f.type == 'gene':
                    if 'locus_tag' in f.qualifiers and f.qualifiers['locus_tag'][0] == str(name):
                        qualifier = 'locus_tag'
                    elif 'gene' in f.qualifiers and f.qualifiers['gene'][0] == str(name):
                        qualifier = 'gene'
                if qualifier:
                    break
            if qualifier:
                break
    if qualifier is None:
        raise Exception('Input gene IDs does not correspond to Genbank IDs. Needed identifier: Locus tag or Gene name.')
    return qualifier

# get positions of each gene in genbank file
def index_genbank_features(gb_record, feature_type, qualifier):
    answer = dict()
    for (index, feature) in enumerate(gb_record.features):
        if feature.type == feature_type:
            if qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier]:
                    if value in answer:
                        print('WARNING - Duplicate key %s for %s features %i and %i' \
                           % (value, feature_type, answer[value], index))
                    else:
                        answer[value] = index
    return answer

### import gene regulatory network
def import_GRN(input_GRN):
    GRN = read_table(input_GRN, index_col=[0], sep=None, engine='python') # import table with unknown separator
    print('GRN uploaded.')
    return GRN
