#Common imports
import pandas as pd

#lib imports
from lib import genome_simulator
from lib import file_setup
from lib import root_mean_square_error

def analyze_mutation(genome_tracker_new, output_dir, df, mutation_number):
    """
    Calls the pinetree script to simulate the genome x number of times and calculate the average of transcript abundances.
    Returns the sum of squares value to compare the fitness of the new and previously accepted genome.
    Input(s):
    genome_tracker_new is the dataframe containing the most recently edited genomic data.
    output_dir is a string containing information of the path to the directory in which all the saved files are stored by the program.
    df is the dataframe containing information regarding the target transcript abundances.
    mutation_number refers to the number of times to run a simulation on the same architecture to reduce noise.
    Output(s):
    Returns a floating point number that refers to the sum of squared error for the recently simulated genome.
    """

    dfs = []

    #Creates test files from pinetree to find average number of transcripts at each time
    for i in range(1, mutation_number+1):
        #If each individual RNase site's binding strength IS being altered
        genome_simulator.pt_call(output_dir, genome_tracker_new, df.index[-1])
        save_df = pd.read_csv(output_dir+"expression_pattern.tsv", header=0, sep='\t')
        save_df['time'] = save_df['time'].round().astype(int)
        dfs.append(save_df)

    #Averages all the values in each file and creates a new file with those averages
    df_concat = pd.concat(dfs)
    df_gb = df_concat.groupby(['time', 'species'], as_index=False)
    df_mean = df_gb.sum()
    df_mean[['protein', 'transcript', 'ribo_density']] = df_mean[['protein', 'transcript', 'ribo_density']] / mutation_number
    df_mean.to_csv(output_dir+'expression_pattern.tsv', sep='\t', index=False)
    #New file is read in as a dataframe
    nf = pd.read_csv(output_dir+"expression_pattern.tsv", header=0, sep='\t')
    nf = file_setup.rearrange_file(nf, df.index[-1], genome_tracker_new['num_genes'])

    return root_mean_square_error.calc_nrmse(df, nf)