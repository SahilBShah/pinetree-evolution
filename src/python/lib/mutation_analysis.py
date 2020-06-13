#Common imports
import numpy as np
import pandas as pd

#lib imports
import genome_simulator
import file_setup
import sum_of_squares

def analyze_mutation(genome_tracker_new, output_dir, df, mutation_number, deg_rate, sse_range=False):
    """
    Calls the pinetree script to simulate the genome x number of times and calculate the average of transcript abundances.
    Returns the sum of squares value to compare the fitness of the new and previously accepted genome.
    """

    dfs = []
    sse_list = []

    #Creates test files from pinetree to find average number of transcripts at each time
    for i in range(1, mutation_number+1):
        if deg_rate:
            #If each individual RNase site's binding strength IS being altered
            genome_simulator.pt_call_alt(output_dir, genome_tracker_new)
        else:
            #If each individual RNase site's binding strength is NOT being altered
            genome_simulator.pt_call(output_dir, genome_tracker_new)
        if sse_range:
            nf = pd.read_csv(output_dir+"expression_pattern.tsv", header=0, sep='\t')
            nf = file_setup.rearrange_file(nf, genome_tracker_new)
            sse_list.append(sum_of_squares.calc_sse(df, nf, genome_tracker_new['num_genes']))
        save_df = pd.read_csv(output_dir+"expression_pattern.tsv", header=0, sep='\t')
        save_df['time'] = save_df['time'].round().astype(int)
        dfs.append(save_df)

    if sse_range:
        return sse_list

    #Averages all the values in each file and creates a new file with those averages
    df_concat = pd.concat(dfs)
    df_gb = df_concat.groupby(['time', 'species'], as_index=False)
    df_mean = df_gb.sum()
    df_mean[['protein', 'transcript', 'ribo_density']] = df_mean[['protein', 'transcript', 'ribo_density']] / mutation_number
    df_mean.to_csv(output_dir+'expression_pattern.tsv', sep='\t', index=False)
    #New file is read in as a dataframe
    nf = pd.read_csv(output_dir+"expression_pattern.tsv", header=0, sep='\t')
    nf = file_setup.rearrange_file(nf, genome_tracker_new)

    return sum_of_squares.calc_sse(df, nf, genome_tracker_new['num_genes'])
