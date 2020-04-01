import pandas as pd

#Minimizes distance between observed values and regression line from given data
def calc_sum_of_squares(target_file, new_file, num_genes):
    """
    Calculates the sum of squares of the target input file and the simulated output file to determine if the sum of squares decreased with a new mutation or not.
    SSE = sum of (y - y-hat)^2
    """

    sos = 0

    df = new_file
    df1 = target_file
    df_diff = df - df1
    df_squared = df_diff ** 2
    sum_df_squared = df_squared.sum()
    for index in range(num_genes):
        sos+=sum_df_squared[index]
    return sos
