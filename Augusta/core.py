from pandas import DataFrame, concat
from numpy import where, sum
import sys
from os import path, remove, devnull

from .GRN_MI import GRNmi_infer, preprocessing
from .motif_discovery import find_motifs
from .DBs import DBs_search, DBs_filter, CC_search_filter, get_synonyms
from .data_export import export_interactions, export_GRN


# infer GRN from count table using mutual information
def RNASeq_to_GRNmi(count_table):
    print('Mutual information computation...')
    gene_names = count_table.index
    count_table_differences, highest_difference = GRNmi_infer.preprocess(count_table)
    MI_matrix = GRNmi_infer.calc_MI(count_table, highest_difference)
    CoeN_np = where(MI_matrix > 0, 1, MI_matrix).astype(dtype='int8')
    CoeN = DataFrame(CoeN_np, columns=gene_names, index=gene_names, dtype='int8')
    GRNmi = preprocessing.find_interaction_type(CoeN, count_table_differences)
    return GRNmi, count_table_differences

# verify GRN using motif search
def GRNmi_to_GRNmotifs(GRNmi, gene_promoters, count_table_differences, motifs_max_time):
    n_of_reg_genes = sum(GRNmi > 0)
    if sum(n_of_reg_genes >= 5) == 0:  # if no TF regulates at least 5 genes, motif search is not possible
        print('Motif search not available; skipped.')
        return GRNmi
    else:
        CoeNmotifs = find_motifs(GRNmi, gene_promoters, motifs_max_time)
        if path.exists('temporary_coreg_seq.fasta'):
            remove('temporary_coreg_seq.fasta')
        if path.getsize('output/discovered_motifs.sto') == 0:
            remove('output/discovered_motifs.sto')
            print('No motif found; skipped.')
            return GRNmi
        elif sum(CoeNmotifs.sum(axis=1)) == 0:
            print('No motifs reversely searched in genome. Discovered motifs have been saved as "discovered_motifs.sto".')
            return GRNmi
        else:
            GRNmotifs = preprocessing.find_interaction_type(CoeNmotifs, count_table_differences)
            return GRNmotifs

# search genes interactions in databases
def get_DBs_data(organism, genes_IDs_gb, taxon, GRN = None):
    global organism_all_names
    print('Synonym organism names search...')
    organism_all_names = get_synonyms.organism_synonyms(organism)
    print('Synonym genes names search...')
    sys.stdout = open(devnull, 'w')  # block print
    genes_IDs_all = get_synonyms.gene_synonyms(genes_IDs_gb, taxon)
    sys.stdout = sys.__stdout__  # enable print
    print('Synonym names search done.')

    # search interaction databases
    op_interactions = DBs_search.search_OmniPath(taxon)
    signor_interactions = DBs_search.search_Signor(taxon)
    sl_interactions = DBs_search.search_SignaLink(taxon)
    trrust_interactions = DBs_search.search_TRRUST(organism_all_names)
    try:
        allDBs_interactions = concat([signor_interactions, op_interactions, sl_interactions, trrust_interactions])
    except ValueError:
        print('No data searched in interaction databases.')
        return GRN, genes_IDs_all
    else:
        interactions_all, interactions_uncertain, interactions_filtered = DBs_filter.filter_interactions(allDBs_interactions, genes_IDs_all)
        export_interactions(interactions_all, interactions_uncertain)
        if len(interactions_filtered) > 0:
            print('%s interactions found in databases (%s genes in GRN).' % (len(interactions_filtered), len(genes_IDs_all)))
            if GRN is None: # GRN only from databases interactions data - GRN was not inferred / imported
                GRNdb = DBs_filter.DBsInteractions_to_GRN(interactions_filtered)
            else:
                gene_names = GRN.index
                GRNdb = DBs_filter.fill_GRN(GRN, gene_names, interactions_filtered) # add databases interactions to inferred / imported GRN
        else: GRNdb = GRN # no interactions found, input GRN is returned
    export_GRN(GRNdb)
    return GRNdb, genes_IDs_all

# search Cell Collective database for interactions and logical rules
def get_CC_data(genes_IDs, GRN):
    CC_models = CC_search_filter.CC_search(organism_all_names)
    if len(CC_models) > 0:
        GRNcc, CC_conditions = CC_search_filter.CC_filter(genes_IDs, GRN, CC_models)
    else:
        GRNcc = GRN.copy(); CC_conditions = {}
    return GRNcc, CC_conditions
