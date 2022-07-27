import omnipath as op
import requests
from bs4 import BeautifulSoup
import lxml
import html5lib
import pandas as pd
from io import BytesIO

###search OmniPath database
def search_OmniPath(taxon):
    op_all_interactions = None
    op_organisms = dict((o.code, o.value) for o in op.constants.Organism)
    i = 0
    while True:
        if taxon == list(op_organisms.keys())[i]:
            try:
                op_all_interactions = op.interactions.AllInteractions.get(organism=list(op_organisms.values())[i], genesymbols = True)
                op_all_interactions = op_all_interactions.dropna(subset=['is_stimulation', 'is_inhibition', 'source_genesymbol', 'target_genesymbol'])
            except:
                print('OmniPath database not available; skipped.')
                return op_all_interactions
            break
        i += 1
        if i == len(op_organisms): break

    # filter OmniPath data
    if (op_all_interactions is not None) and len(op_all_interactions) > 0:
        op_all_interactions['interaction'] = 0
        op_all_interactions.loc[(op_all_interactions['is_stimulation'] == True) & (op_all_interactions['is_inhibition'] == False), 'interaction'] = 1
        op_all_interactions.loc[(op_all_interactions['is_stimulation'] == False) & (op_all_interactions['is_inhibition'] == True), 'interaction'] = -1
        op_interactions = op_all_interactions[['source_genesymbol', 'target_genesymbol', 'interaction', 'references']]
        op_interactions = op_interactions[op_interactions['interaction'] != 0]
        op_interactions.columns = ['source_GeneSymbol', 'target_GeneSymbol', 'interaction_type', 'references']
        op_interactions['references'] = op_interactions['references'].fillna('NOREFID') # replace nan values
        op_interactions['references'] = (op_interactions['references'].str.upper()).str.split(';')
        op_interactions = op_interactions.explode('references')
        return op_interactions

    else: return op_all_interactions

### search Signor database
def search_Signor(taxon):
    signor_interactions = None
    signor_organisms_url = 'https://signor.uniroma2.it'
    try:
        signor_source_code = requests.get(signor_organisms_url).text
        signor_soup = BeautifulSoup(signor_source_code, features='lxml')
        # extract organisms names in Signor
        select = signor_soup.find('select', {'name': 'organism'})
        signor_organisms = [option['value'] for option in select.find_all('option')]
    except:
        print('Signor database not available; skipped.')
        return signor_interactions

    taxon_str = '(%a)'%(taxon)
    for value in signor_organisms:
        if (taxon_str in value):
            signor_url = 'https://signor.uniroma2.it/getData.php'
            res = requests.post(signor_url, data={'submit': 'Download'})
            if res.status_code == 200:
                rc = res.content
                signor_all_interactions = pd.read_csv(BytesIO(rc), header=None, delimiter="\t")
                if len(signor_all_interactions.columns) == 29:
                    signor_all_interactions.columns = ['ENTITYA','TYPEA','IDA','DATABASEA','ENTITYB','TYPEB','IDB','DATABASEB','EFFECT',
                                                                  'MECHANISM','RESIDUE','SEQUENCE','TAX_ID','CELL_DATA','TISSUE_DATA','MODULATOR_COMPLEX',
                                                                  'TARGET_COMPLEX','MODIFICATIONA','MODASEQ','MODIFICATIONB','MODBSEQ','PMID','DIRECT','NOTES',
                                                                  'ANNOTATOR','SENTENCE','SIGNOR_ID','SCORE','interaction']
                    # filter Signor data
                    signor_all_interactions['interaction'] = 0
                    signor_all_interactions = signor_all_interactions.dropna(subset=['TAX_ID', 'DIRECT', 'EFFECT', 'ENTITYA', 'ENTITYB'])
                    signor_all_interactions = signor_all_interactions[(signor_all_interactions['TAX_ID'] == taxon) & (signor_all_interactions['DIRECT'] == 't')]
                    if len(signor_all_interactions) > 0:
                        signor_all_interactions.loc[signor_all_interactions['EFFECT'].str.contains('up-regulates'), 'interaction'] = 1
                        signor_all_interactions.loc[signor_all_interactions['EFFECT'].str.contains('down-regulates'), 'interaction'] = -1
                        signor_interactions = signor_all_interactions[['ENTITYA', 'ENTITYB', 'interaction', 'PMID']]
                        signor_interactions = signor_interactions[signor_interactions['interaction'] != 0]
                        signor_interactions.columns = ['source_GeneSymbol', 'target_GeneSymbol', 'interaction_type', 'references']
                        signor_interactions['references'] = ['PMID:'] + signor_interactions['references']
                        signor_interactions.loc[signor_interactions['references'] == 'PMID:Other', 'references'] = 'NOREFID'

            else:
                print('Signor database not available; skipped.')
    return signor_interactions

### search SignaLink database
def search_SignaLink(taxon):
    sl_interactions = None
    try:
        sl_organisms = pd.read_table('http://korcsmaroslab.org/slk3/slk3_current_speices.tsv')
    except:
        print('SignaLink database not available; skipped.')
        return sl_interactions
    for value in sl_organisms["Taxon.ID"]:
        if (int(taxon) == value):
            sl_organism_split = sl_organisms.loc[sl_organisms['Taxon.ID'] == taxon, 'Taxon.name'].item().lower().split()
            sl_organism = sl_organism_split[0] + "_" + sl_organism_split[1]
            try:
                sl_all_interactions = pd.read_table('http://korcsmaroslab.org/slk3/' + sl_organism + '_layer_0_1_3_filtered.tsv', low_memory=False)
                sl_all_interactions = sl_all_interactions.dropna(subset=['source_name', 'source_speciesID', 'target_name', 'target_speciesID', 'directness', 'interaction_type'])
            except:
                print('SignaLink database not available; skipped.')
                return sl_interactions

            # filter SignaLink data
            sl_all_interactions = sl_all_interactions[(sl_all_interactions['directness'] == 'directed')]
            if len(sl_all_interactions) > 0:
                sl_all_interactions['interaction'] = 0
                sl_all_interactions.loc[sl_all_interactions['interaction_type'].str.contains('MI:0624'), 'interaction'] = 1
                sl_all_interactions.loc[sl_all_interactions['interaction_type'].str.contains('MI:0623'), 'interaction'] = -1
                sl_interactions = sl_all_interactions[['source_name', 'target_name', 'interaction', 'references']]
                sl_interactions = sl_interactions[sl_interactions['interaction'] != 0]
                sl_interactions['references'] = sl_interactions['references'].fillna('NOREFID')
                sl_interactions.columns = ['source_GeneSymbol', 'target_GeneSymbol', 'interaction_type', 'references']
                sl_interactions['references'] = (sl_interactions['references'].str.upper()).str.split('|')
                sl_interactions = sl_interactions.explode('references')
    return sl_interactions

### search TRRUST database
def search_TRRUST(input_organism_all_names):
    trrust_interactions = None
    trrust_url = 'https://www.grnpedia.org/trrust/downloadnetwork.php'
    try:
        trrust_source_code = requests.get(trrust_url).text
        trrust_soup = BeautifulSoup(trrust_source_code, features='html5lib')
        # extract organisms names in TRRUST
        trrust_table = trrust_soup.find(lambda tag: tag.name == 'table' and tag.find(lambda ttag: ttag.name == 'th' and ttag.text == 'Species'))
        trrust_organisms_list = [row.td.text for row in trrust_table.find_all('tr')[1:]]
    except:
        print('TRRUST database not available; skipped.')
        return trrust_interactions

    # download TRRUST if it contains input organism´s data
    trrust_downloaded = False
    for trrust_organism in trrust_organisms_list:
        for input_organism in input_organism_all_names:
            if trrust_organism.lower() == input_organism.lower():
                try:
                    trrust_all_interactions = pd.read_table('https://www.grnpedia.org/trrust/data/trrust_rawdata.' + trrust_organism.lower() + '.tsv',header=None)
                    trrust_all_interactions.columns = ['source_name', 'target_name', 'interaction_type', 'references (PMID)']
                    trrust_downloaded = True
                    break
                except:
                    print('TRRUST database not available; skipped.')
                    return trrust_interactions
        if trrust_downloaded:
            trrust_all_interactions = trrust_all_interactions.dropna(subset=['source_name', 'target_name', 'interaction_type'])
            if len(trrust_all_interactions) > 0:
                trrust_all_interactions['interaction'] = 0
                trrust_all_interactions.loc[trrust_all_interactions['interaction_type'].str.contains('Activation'), 'interaction'] = 1
                trrust_all_interactions.loc[trrust_all_interactions['interaction_type'].str.contains('Repression'), 'interaction'] = -1
                trrust_interactions = trrust_all_interactions[['source_name', 'target_name', 'interaction', 'references (PMID)']]
                trrust_interactions.columns = ['source_GeneSymbol', 'target_GeneSymbol', 'interaction_type', 'references']
                trrust_interactions = trrust_interactions[trrust_interactions['interaction_type'] != 0]
                trrust_interactions['references'] = trrust_interactions['references'].str.split(';')  # str to list
                trrust_interactions = trrust_interactions.explode('references')
                trrust_interactions['references'] = ['PMID:'] + trrust_interactions['references']
            break
    return trrust_interactions