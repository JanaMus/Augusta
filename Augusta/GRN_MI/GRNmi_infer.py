from math import floor, sqrt
from pandas import DataFrame
from numpy import histogram2d, zeros, shape, arange, delete
from .preprocessing import expression_differences, select_expression
from sklearn.metrics import mutual_info_score

# preprocess count table for MI computation
def preprocess(count_table):
    gene_names = count_table.index
    difference_matrix = expression_differences(count_table, gene_names)
    highest_expression = select_expression(difference_matrix, gene_names)
    return difference_matrix, highest_expression

# compute MI
def MI_value(x, y, bins):
    bin_xy = histogram2d(x, y, bins)[0]
    MI = mutual_info_score(None, None, contingency=bin_xy)
    return MI

# create matrix of MIs
def calc_MI(input_matrix, highest_expression):
    gene_names = input_matrix.index
    MI_matrix_np = zeros((shape(input_matrix)[0], shape(input_matrix)[0]))
    MI_matrix = DataFrame(MI_matrix_np, columns=gene_names, index=gene_names)

    # set number of bins to discretize normalized RNA-seq
    bins = floor(sqrt(shape(input_matrix)[0]/5))
    if bins > 10:
        bins = 10

    # compute MI
    vector_cols = arange(1, shape(input_matrix)[0])
    for i in range(0, shape(input_matrix)[0]):
        #print(f'MI for gene: {i+1} / {shape(input_matrix)[0]}')
        for j in vector_cols:
            if highest_expression.iloc[i, 1] != highest_expression.iloc[j, 1]: # compute MI only if highest expression is in different time points
                if highest_expression.iloc[i, 1] < highest_expression.iloc[j, 1]: # set direction of edges
                    MI_matrix.iloc[i, j] = MI_value(input_matrix.iloc[i, :], input_matrix.iloc[j, :], bins)
                else:
                    MI_matrix.iloc[j, i] = MI_value(input_matrix.iloc[i, :], input_matrix.iloc[j, :], bins)
        if len(vector_cols) > 0:
            vector_cols = delete(vector_cols, 0) # delete first item in vector in order to fill only triangle in MI_matrix
    print('Mutual information computation done.')
    return MI_matrix
