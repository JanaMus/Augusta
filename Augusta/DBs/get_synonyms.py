from EcoNameTranslator import to_common, synonyms
from numpy import unique
from mygene import MyGeneInfo

### translate scientific organism name into common name
def organism_synonyms(organism):
    if len(organism.split()) > 2: # remove strain from organism name and save as list
        organism = [' '.join(organism.split()[0:2])]
    else:
        organism = [organism]
    try: common_names = to_common(organism)
    except: common_names = []
    else: common_names = list(common_names.values())[0][1]

    try: synonym_names = synonyms(organism)
    except: synonym_names = []
    else: synonym_names = list(synonym_names.values())[0][1]

    organism_all_names = organism + common_names + synonym_names
    return organism_all_names

### get gene name synonyms
def gene_synonyms(genes_IDs_input, taxon):
    genes_IDs = genes_IDs_input.copy()
    mg = MyGeneInfo()
    for gene_key in genes_IDs_input:
        if genes_IDs_input[gene_key]:
            genes = [gene_key] + genes_IDs_input[gene_key]
        else:
            genes = [gene_key]
        all_IDs = genes
        try:
            new_IDs = mg.querymany(genes, scopes='symbol, locus_tag, homologene, alias, accession', fields='symbol, alias', as_dataframe=True, species=taxon)
        except:
            new_IDs = None
        if new_IDs is not None:
            if 'notfound' in new_IDs:
                new_IDs = new_IDs[new_IDs['notfound'] != True]  # filter searched names
                new_IDs = new_IDs.drop('notfound', axis=1)
            if not new_IDs.empty:
                if '_score' in new_IDs:
                    new_IDs = new_IDs.drop('_score', axis=1)
                new_IDs_list = sum([new_IDs[i].tolist() for i in new_IDs.columns], [])
                new_IDs_extracted = []  # unlist listed values
                for value in new_IDs_list:
                    if isinstance(value, list):
                        new_IDs_extracted = value + new_IDs_extracted
                    elif isinstance(value, str):
                        new_IDs_extracted = [value] + new_IDs_extracted
                all_IDs += new_IDs_extracted
        all_IDs_unique = unique(all_IDs).tolist()
        genes_IDs[gene_key] = all_IDs_unique
    return genes_IDs
