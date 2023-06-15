from ccapi import Client
import logging

### get interactions and logical rules from Cell Collective (CC) database
def CC_search(organism_all_names):
    models = []
    try:
        logging.disable(logging.INFO) # disable print log info messages
        logging.disable(logging.WARNING)
        client = Client()
        fetch_models = client.get('model', size = 1000, filters = {'modelTypes': ['boolean']})
        logging.disable(logging.NOTSET)
    except:
        print('Cell Collective database not available; skipped.')
        return models

    # filter CC models based on input organism
    pos = 0
    while pos < len(fetch_models):
        one_model = fetch_models[pos]
        for name in organism_all_names:
            if (one_model.name and (name in one_model.name)) or (one_model.description and (name in one_model.description)):
                models.append(one_model.default_version)
                pos += 1
                break
            elif name == organism_all_names[-1]:
                pos += 1
    print('Cell Collective database search done.')
    return models

### filter CC data based on input organism´s genes
def CC_filter(genes_IDs_all, GRN, models):
    GRN_CC = GRN.copy()
    global genes_IDs_all_values
    genes_IDs_all_values = sum(genes_IDs_all.values(), [])
    CC_conditions = {}
    for model in models:
        for component in model.internal_components: # check if component and regulator genes are present in the input data
            if component.name in genes_IDs_all_values:
                helper = 0 # helper for adding unique key to conditions dictionary
                for regulator in component.regulators:
                    if regulator.species.name in genes_IDs_all_values: # add edges to the GRN
                        reg_gene_name = get_keys_from_value(genes_IDs_all, regulator.species.name)
                        comp_gene_name = get_keys_from_value(genes_IDs_all, component.name)
                        if regulator.type == 'positive':
                            GRN_CC.loc[reg_gene_name, comp_gene_name] = 1
                        elif regulator.type == 'negative':
                            GRN_CC.loc[reg_gene_name, comp_gene_name] = -1

                        if len(regulator.conditions) > 0: # add conditions
                            for condition in regulator.conditions:
                                condition_values = get_conditions(condition, genes_IDs_all)
                                subcondition_values = []
                                if len(condition.sub_conditions) > 0: # add subconditions
                                    for subcondition in condition.sub_conditions:
                                        subcondition_values.append(get_conditions(subcondition, genes_IDs_all))
                                values = [regulator.species.name, condition_values, subcondition_values]
                                if len(values) > 0:
                                    CC_conditions[component.name+'_'+str(helper)] = values
                                    helper += 1
    return GRN_CC, CC_conditions


def get_keys_from_value(d, val):
    name = [key for key, value in d.items() if val in value]
    return name[0]

# get individual conditions, check if CC nodes are in input data
def get_conditions(condition, genes_IDs_all):
    should_continue = True
    names = []
    values = []
    for cComp in condition.components:
        if cComp.name in genes_IDs_all_values:  # check if condition genes in present in the input data
            name = get_keys_from_value(genes_IDs_all, cComp.name)
            names.append(name)
        else:
            should_continue = False
            break
    if should_continue:
        names.append(cComp.name)
        values = [names, condition.type, condition.relation, condition.state]
    return values
