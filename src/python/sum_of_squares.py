import pandas as pd

#Minimizes distance between observed values and regression line from given data
def calc_sum_of_squares(target_file, new_file):
    """
    Calculates the sum of squares of the target input file and the simulated output file to determine if the sum of squares decreased with a new mutation or not.
    SSE = sum of (y - y-hat)^2
    """

    df = new_file
    df1 = target_file
    df_diff = df - df1
    df_squared = df_diff ** 2
    sum_df_squared = df_squared.sum()
    sos = sum_df_squared[0] + sum_df_squared[1] + sum_df_squared[2]
    return sos
