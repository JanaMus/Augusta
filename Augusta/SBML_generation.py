from .cno.io import cnograph
from numpy import all
from os import remove
from gc import collect

### convert GRN to BN in SBML-qual file format
def GRN_to_SBML(GRN, CC_conditions = {}):
    lonely_nodes = GRN_to_BNequations(GRN, CC_conditions)
    try:
        c_sbml = BNequations_to_SBML('output/BN.txt')
    except MemoryError:
        print('Low memory for exporting BN to SBML-qual. BN is saved in a form of logical equations as "BN.txt".')
        return
    else:
        SBML_to_file(c_sbml, lonely_nodes)
        remove('output/BN.txt')

### convert GRN to Boolean netwok in logical equations
def GRN_to_BNequations(GRN, CC_conditions):
    lonely_nodes = []
    with open('output/BN.txt', 'w') as BN_file:  # motif sequences in Stockholm format
        for source in GRN.columns: # get column name
            neg_values = GRN[source].index[GRN[source] < 0].tolist() # extract negative regulators
            pos_values = GRN[source].index[GRN[source] > 0].tolist() # extract positive regulators
            react = str()
            if (len(CC_conditions) > 0) and ((source + '_0') in CC_conditions): # get conditions if exist for the source and target nodes
                source_conditions = get_conditions(source, CC_conditions)
            else: source_conditions = []

            if (len(neg_values) > 0) and (len(pos_values) > 0):
                for pos_regulator in pos_values:
                    pos_regulator_number = 0
                    more_pos_regulator_conditions = True
                    while more_pos_regulator_conditions == True:
                        add_pos_regulator, more_pos_regulator_conditions, pos_regulator_number = search_condition(pos_regulator, pos_regulator_number, source_conditions)
                        if len(react) > 0:
                            react = react + '+'
                        react = react + add_pos_regulator
                        for neg_regulator in neg_values:
                            neg_regulator_number = 0
                            more_neg_regulator_conditions = True
                            while more_neg_regulator_conditions == True:
                                add_neg_regulator, more_neg_regulator_conditions, neg_regulator_number = search_condition(neg_regulator, neg_regulator_number, source_conditions)
                                react = react + '^!' + add_neg_regulator

            elif (len(neg_values) > 0) and (len(pos_values) == 0):
                for neg_regulator in neg_values:
                    neg_regulator_number = 0
                    more_neg_regulator_conditions = True
                    while more_neg_regulator_conditions == True:
                        add_neg_regulator, more_neg_regulator_conditions, neg_regulator_number = search_condition(neg_regulator, neg_regulator_number, source_conditions)
                        if len(react) > 0:
                            react = react + '^'
                        react = react + '!' + add_neg_regulator

            elif (len(neg_values) == 0) and (len(pos_values) > 0):
                for pos_regulator in pos_values:
                    pos_regulator_number = 0
                    more_pos_regulator_conditions = True
                    while more_pos_regulator_conditions == True:
                        add_pos_regulator, more_pos_regulator_conditions, pos_regulator_number = search_condition(pos_regulator, pos_regulator_number, source_conditions)
                        if len(react) > 0:
                            react = react + '+'
                        react = react + add_pos_regulator

            if len(react) > 0:
                react = react + '=' + source
                BN_file.write(react + '\n')
            elif all((GRN.loc[source]) == 0): # check for nodes without edges
                lonely_nodes.append(source)
            collect()
    return lonely_nodes

### convert Boolean netwok in logical equations into BN is SBML-qual file format
def BNequations_to_SBML(filename):
    c = cnograph.CNOGraph() # generate SBML-qual object
    with open(filename) as file:
        for line in file:
            c.add_reaction(line)
    c_sbml = c.to_sbmlqual()
    return c_sbml

### write SBML-qual file; add nodes without edges (lonely nodes)
def SBML_to_file(c_sbml, lonely_nodes):
    with open('output/BN.sbml','w') as f_out:
        for line_no, line in enumerate(c_sbml.split("\n")):
            f_out.write(line + "\n")
            if "<qual:listOfQualitativeSpecies" in line:
                for lonely_node in lonely_nodes:
                    f_out.write("""<qual:qualitativeSpecies qual:constant="false" qual:compartment="main" qual:id="{0}"/>\n""".format(lonely_node))
    print('Boolean network stored as "BN.sbml".')

# add reaction´s equation in case conditions are available
def get_conditions(source, CC_conditions):
    source_conditions = {}
    condition_number = 0  # helper variable for getting condition in case one component has more conditions
    regulator_number = 0 # helper variable for getting condition in case one regulator has more conditions
    more_conditions = True
    while more_conditions == True:
        if (source + '_' + str(condition_number)) in CC_conditions:  # get conditions if exist for the source and target nodes
            condition_values = CC_conditions[source + '_' + str(condition_number)]
            condition_number += 1
            condition_react = condition_values[0]
            for pos, condition in enumerate(condition_values[1][0]): # get condition nodes
                    if (pos > 0) and (condition_values[1][2] == 'independent'): # add next condition nodes
                        condition_react = condition_react + '+' + condition_values[0]
                    condition_react = get_conditions_parameters(condition_values[1], condition_react, condition)
                    if ((pos+1) == len(condition_values[1][0])) or (condition_values[1][2] == 'independent'):
                        if len(condition_values[2]) > 0:  # get subconditions
                            subcond_react = condition_react
                            for sub_pos, subcondition in enumerate(condition_values[2][0][0]):
                                if (condition_values[2][0][2] == 'cooperative'):
                                    condition_react = get_conditions_parameters(condition_values[2][0], condition_react, subcondition)
                                elif (condition_values[2][0][2] == 'independent'):
                                    if sub_pos >0:
                                        condition_react = condition_react + '+' + get_conditions_parameters(condition_values[2][0], subcond_react, subcondition)
                                    else:
                                        condition_react = get_conditions_parameters(condition_values[2][0], subcond_react, subcondition)
            while condition_values[0]+'_'+str(regulator_number) in source_conditions: # get unique key in conditions dictionary
                regulator_number += 1
            source_conditions[condition_values[0]+'_'+str(regulator_number)] = condition_react
            regulator_number = 0
        else:
            more_conditions = False
    return source_conditions

# add condition´s parameter to reaction equation
def get_conditions_parameters(condition_value, condition_react, condition):
    if (condition_value[1] == 'if' and condition_value[3] == 'on') or (condition_value[1] == 'unless' and condition_value[3] == 'off'):
        condition_react = condition_react + '^' + condition
    elif (condition_value[1] == 'if' and condition_value[3] == 'off') or (condition_value[1] == 'unless' and condition_value[3] == 'on'):
        condition_react = condition_react + '^!' + condition
    return condition_react

def search_condition(regulator, regulator_number, source_conditions):
    if (str(regulator) + '_' + str(regulator_number)) in source_conditions:  # add conditions
        add_regulator = source_conditions[regulator + '_' + str(regulator_number)]
        regulator_number += 1
        more_regulator_conditions = True
    else:
        add_regulator = regulator
        more_regulator_conditions = False
    return add_regulator, more_regulator_conditions, regulator_number
