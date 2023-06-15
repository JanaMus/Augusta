from pandas import DataFrame

### export interactions found in databases
def export_interactions(interactions_all, interactions_uncertain):
    if len(interactions_all) > 0:
        DataFrame(interactions_all).to_csv('output/DBs_interactions_list.csv', index=False)
        print('All interactions searched across databases stored as "DBs_interactions_list.csv".')
    if len(interactions_uncertain) > 0:
        DataFrame(interactions_uncertain).to_csv('output/DBs_interactions_uncertain.csv', index=False)
        print('Uncertain interactions stored as "DBs_interactions_uncertain.csv". The more prevaled interaction type was used for the GRN inference.')

### export GRN
def export_GRN(GRN):
    DataFrame(GRN).to_csv('output/GRN.csv')
    print('GRN stored as "GRN.csv".')
