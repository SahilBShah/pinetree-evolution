import pandas as pd

#Minimizes distance between observed values and regression line from given data
def calc_sse(target_file, new_file, num_genes):
    """
    Calculates the sum of squares of the target input file and the simulated output file to determine if the sum of squares decreased with a new mutation or not.
    SSE = sum of (y - y-hat)^2
    Input(s):
    target_file is the user-inputted tsv file containing transcript abundances for each gene.
    new_file is the simulator-generated tsv file containng transcript abundances for each gene.
    num_genes refers to the number of genes in the genome.
    """

    sse = 0

    #The two dataframes are compared to one another
    df = new_file
    df1 = target_file
    #Each row is subtracted from its corresponding row in the other dataframe
    df_diff = df - df1
    #The subtracted values are squared
    df_squared = df_diff ** 2
    #The squared values are added together
    sum_df_squared = df_squared.sum()
    #The added values for each gene are added together for an overall sum of squared value
    for index in range(num_genes):
        sse+=sum_df_squared[index]
    return sse

def calc_accepted_sse_range(target_file, genome_tracker_new):
    """
    Calculates the highest SSE value allowed before the program has found the a suitable architecture to reproduce the target data.
    Input(s):
    target_file is the user-inputted tsv file containing transcript abundances for each gene.
    genome_tracker_new is the dataframe containing the most recent edited genomic data.
    """

    #Decreases the transcript abundances of the target file by 10%
    altered_file = target_file[['protein1', 'protein2', 'protein3']] * 0.9
    #Calculates the sum of squared error between the reduced values and the target values and returns it
    max_sse = calc_sse(target_file, altered_file, genome_tracker_new['num_genes'])
    return max_sse
