import pandas as pd

#Minimizes distance between observed values and regression line from given data
def calc_sse(target_file, new_file, num_genes):
    """
    Calculates the sum of squares of the target input file and the simulated output file to determine if the sum of squares decreased with a new mutation or not.
    SSE = sum of (y - y-hat)^2
    """

    sse = 0

    df = new_file
    df1 = target_file
    df_diff = df - df1
    df_squared = df_diff ** 2
    sum_df_squared = df_squared.sum()
    for index in range(num_genes):
        sse+=sum_df_squared[index]
    return sse

def calc_accepted_sse_range(target_file, genome_tracker_new):
    """
    Calculates the highest SSE value allowed before the program has found the a suitable architecture to reproduce the target data.
    """

    altered_file = target_file[['protein1', 'protein2', 'protein3']] * 0.9
    max_sse = calc_sse(target_file, altered_file, genome_tracker_new['num_genes'])
    return max_sse
