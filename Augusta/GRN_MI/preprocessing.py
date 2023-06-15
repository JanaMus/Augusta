from numpy import zeros, shape
from pandas import DataFrame

### compute times differences from count table data
def expression_differences(input_matrix, gene_names):
    dif_matrix_np = zeros((shape(input_matrix)[0], shape(input_matrix)[1] - 1))
    column_names = list(range(1, shape(dif_matrix_np)[1] + 1))
    difference_matrix = DataFrame(dif_matrix_np, columns=column_names, index=gene_names)

    for i in range(0, shape(difference_matrix)[0]):
        for j in range(0, shape(difference_matrix)[1]):
            difference_matrix.iloc[i, j] = (input_matrix.iloc[i, j+1] - input_matrix.iloc[i, j])
    return difference_matrix


### select highest difference in times differences
def select_expression(difference_matrix, gene_names):
    highest_expr_np = zeros((shape(difference_matrix)[0], 2))
    column_names = ['highest expr value', 'n of diff']
    highest_expression = DataFrame(highest_expr_np, columns = column_names, index=gene_names)

    position = difference_matrix.abs().values.argmax(axis=1)
    highest_expression['n of diff'] = position  # column name in dif_matrix where the highest difference occures (i.e. 1 = first difference, column 0)
    highest_expression['highest expr value'] = ([difference_matrix.values[i][position[i]] for i in range(len(difference_matrix.values))])
    return highest_expression

### get type of interactions (positive or negative)
def find_interaction_type(coexpression_matrix, dif_matrix):
    dif_matrix.columns = list(range(0,shape(dif_matrix)[1]))
    gene_names = coexpression_matrix.index
    GRN_matrix = coexpression_matrix.copy()
    max_expr_pos = dif_matrix.abs().values.argmax(axis=1)

    for col in range(0, shape(GRN_matrix)[1]):
        pos_target = max_expr_pos[col]
        if pos_target == 0: # max expression change in zero time - target couldn´t be regulated
            GRN_matrix[gene_names[col]] = GRN_matrix[gene_names[col]].replace(1, 0)
            continue

        expr_value_target = dif_matrix.loc[gene_names[col], pos_target]
        for row in range(0, shape(GRN_matrix)[0]):
            if GRN_matrix.iloc[row, col] == 1:
                expr_value_source = dif_matrix.loc[gene_names[row], pos_target-1]
                if ((expr_value_source > 0 and expr_value_target < 0) or (expr_value_source < 0 and expr_value_target > 0)):
                    GRN_matrix.iloc[row, col] = -1
    return GRN_matrix
