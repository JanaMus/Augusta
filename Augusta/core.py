import pandas as pd
import numpy as np
import sys, os

from data_import import import_CountTable, genbank_process, import_GRN
import GRN_MI
from RNASeq_normalization import normalize_CountTable
from motif_discovery import find_motifs
import DBs
from SBML_generation import GRN_to_SBML
from data_export import export_interactions, export_GRN


# infer GRN from count table using mutual information
def RNASeq_to_GRNmi(count_table):
    print('Mutual information computation...')
    gene_names = count_table.index
    count_table_differences, highest_difference = GRN_MI.preprocess(count_table)
    MI_matrix = GRN_MI.calc_MI(count_table, highest_difference)
    CoeN_np = np.where(MI_matrix > 0, 1, MI_matrix).astype(int)
    CoeN = pd.DataFrame(CoeN_np, columns=gene_names, index=gene_names)
    GRNmi = GRN_MI.find_interaction_type(CoeN, count_table_differences)
    return GRNmi, count_table_differences

# verify GRN using motif search
def GRNmi_to_GRNmotifs(GRNmi, gene_promoters, count_table_differences):
    n_of_reg_genes = np.sum(GRNmi > 0)
    if np.sum(n_of_reg_genes >= 5) == 0:  # if no TF regulates at least 5 genes, motif search is not possible
        print('Motif search not available; skipped.')
        return GRNmi
    else:
        CoeNmotifs = find_motifs(GRNmi, gene_promoters)
        if os.path.getsize('output/discovered_motifs.sto') == 0:
            os.remove('output/discovered_motifs.sto')
            print('No motif found; skipped.')
            return GRNmi
        elif sum(CoeNmotifs.sum(axis=1)) == 0:
            print('No motifs reversely searched in genome. Discovered motifs have been saved as "discovered_motifs.sto".')
            return GRNmi
        else:
            GRNmotifs = GRN_MI.find_interaction_type(CoeNmotifs, count_table_differences)
            return GRNmotifs

# search genes interactions in databases
def get_DBs_data(organism, genes_IDs_gb, taxon, GRN = None):
    global organism_all_names
    print('Synonym organism names search...')
    organism_all_names = DBs.organism_synonyms(organism)
    print('Synonym genes names search...')
    sys.stdout = open(os.devnull, 'w')  # block print
    genes_IDs_all = DBs.gene_synonyms(genes_IDs_gb, taxon)
    sys.stdout = sys.__stdout__  # enable print
    print('Synonym names search done.')

    # search interaction databases
    op_interactions = DBs.search_OmniPath(taxon)
    signor_interactions = DBs.search_Signor(taxon)
    sl_interactions = DBs.search_SignaLink(taxon)
    trrust_interactions = DBs.search_TRRUST(organism_all_names)
    try:
        allDBs_interactions = pd.concat([signor_interactions, op_interactions, sl_interactions, trrust_interactions])
    except ValueError:
        print('No data searched in interaction databases.')
        return GRN
    else:
        interactions_all, interactions_uncertain, interactions_filtered = DBs.filter_interactions(allDBs_interactions, genes_IDs_all)
        export_interactions(interactions_all, interactions_uncertain)
        if len(interactions_filtered) > 0:
            print('%s interactions found in databases (%s genes in GRN).' % (len(interactions_filtered), len(genes_IDs_all)))
            if GRN is None: # GRN only from databases interactions data - GRN was not inferred / imported
                GRNdb = DBs.DBsInteractions_to_GRN(interactions_filtered)
            else:
                gene_names = GRN.index
                GRNdb = DBs.fill_GRN(GRN, gene_names, interactions_filtered) # add databases interactions to inferred / imported GRN
        else: GRNdb = GRN # no interactions found, input GRN is returned
    export_GRN(GRNdb)
    return GRNdb

# search Cell Collective database for interactions and logical rules
def get_CC_data(genes_IDs, GRN):
    CC_models = DBs.CC_search(organism_all_names)
    if len(CC_models) > 0:
        GRNcc, CC_conditions = DBs.CC_filter(genes_IDs, GRN, CC_models)
    else:
        GRNcc = GRN.copy(); CC_conditions = {}
    return GRNcc, CC_conditions