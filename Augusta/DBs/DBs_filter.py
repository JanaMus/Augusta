import numpy as np
import pandas as pd

### filter databases duplications and unclear interactions
def filter_interactions(allDBs_interactions, genes_IDs_all):
    allDBs_interactions_merged = allDBs_interactions['references'].groupby([allDBs_interactions['source_GeneSymbol'], allDBs_interactions['target_GeneSymbol'], allDBs_interactions['interaction_type']]).apply(set).reset_index()
    for i in range(0, len(allDBs_interactions_merged)):
        allDBs_interactions_merged.loc[i, 'NofRefs'] = len(allDBs_interactions_merged['references'][i])

    uncertain_interactions = allDBs_interactions_merged[allDBs_interactions_merged.duplicated(['source_GeneSymbol', 'target_GeneSymbol'], keep=False) == True]
    uncertain_interactions_positions = uncertain_interactions.index
    DBs_interactions = allDBs_interactions_merged.drop(uncertain_interactions_positions)  # delete rows with uncertain interactions
    uncertain_interactions.index = np.array(range(len(uncertain_interactions)))
    DBs_interactions.index = np.array(range(len(DBs_interactions)))
    position = 0
    while position < len(uncertain_interactions):
        if uncertain_interactions.loc[position, 'NofRefs'] > uncertain_interactions.loc[position + 1, 'NofRefs']:
            DBs_interactions.loc[len(DBs_interactions.index)] = uncertain_interactions.loc[position, :]
        elif uncertain_interactions.loc[position, 'NofRefs'] < uncertain_interactions.loc[position + 1, 'NofRefs']:
            DBs_interactions.loc[len(DBs_interactions.index)] = uncertain_interactions.loc[position + 1, :]
        position = position + 2
    DBs_interactions_drop = DBs_interactions.drop(['references', 'NofRefs'], axis=1)
    DBs_interactions_names = filter_GRNGeneMatches(DBs_interactions_drop, genes_IDs_all)
    DBs_interactions_drop['source_GeneID'] = DBs_interactions_names['source_GeneID']
    DBs_interactions_drop['target_GeneID'] = DBs_interactions_names['target_GeneID']
    DBs_interactions_out = DBs_interactions_drop.where(pd.notnull(DBs_interactions_drop), None)

    DBs_interactions_GRNGeneMatches = DBs_interactions_names[['source_GeneID', 'target_GeneID', 'interaction_type']]
    DBs_interactions_GRNGeneMatches.index = np.array(range(len(DBs_interactions_GRNGeneMatches)))

    return DBs_interactions_out, uncertain_interactions, DBs_interactions_GRNGeneMatches

def filter_GRNGeneMatches(DBs_interactions, genes_IDs_all):
    genes_IDs_all_values = sum(genes_IDs_all.values(),[]) # select only DBs interactions (genes) available in the input GRN
    DBs_interactions_genesNames = DBs_interactions.loc[(DBs_interactions['source_GeneSymbol'].isin(genes_IDs_all_values) & DBs_interactions['target_GeneSymbol'].isin(genes_IDs_all_values))]
    DBs_interactions_genesNames['source_GeneID'] = None
    DBs_interactions_genesNames['target_GeneID'] = None
    position = DBs_interactions_genesNames.index
    for key in genes_IDs_all.keys():
        for i in position: #range(0, len(DBs_interactions_genesNames['interaction_type'])):
            if DBs_interactions_genesNames.loc[i, 'source_GeneSymbol'] in genes_IDs_all[key]:
                DBs_interactions_genesNames.loc[i, 'source_GeneID'] = key
            if DBs_interactions_genesNames.loc[i, 'target_GeneSymbol'] in genes_IDs_all[key]:
                DBs_interactions_genesNames.loc[i, 'target_GeneID'] = key
    return DBs_interactions_genesNames

### add databases interactions to existing GRN
def fill_GRN(GRN, names, DBs_interactions_GRNGeneMatches):
    GRN_filled = GRN.copy()
    for i in range(0, len(DBs_interactions_GRNGeneMatches)):
        r = int(np.where(names == DBs_interactions_GRNGeneMatches.loc[i, 'source_GeneID'])[0])
        c = int(np.where(names == DBs_interactions_GRNGeneMatches.loc[i, 'target_GeneID'])[0])
        GRN_filled.iloc[r, c] = DBs_interactions_GRNGeneMatches.loc[i, 'interaction_type']
    return GRN_filled

### infer GRN from databases intractions
def DBsInteractions_to_GRN(DBs_interactions_GRNGeneMatches):
    DB_GRN_names = np.unique(DBs_interactions_GRNGeneMatches[['source_GeneID', 'target_GeneID']].values)
    DBsGRN_np = np.zeros((len(DB_GRN_names), len(DB_GRN_names)), dtype='int8')
    DBsGRN = pd.DataFrame(DBsGRN_np, columns=DB_GRN_names, index=DB_GRN_names, dtype='int8')
    DBsGRN_filled = fill_GRN(DBsGRN, DB_GRN_names, DBs_interactions_GRNGeneMatches)
    return DBsGRN_filled
