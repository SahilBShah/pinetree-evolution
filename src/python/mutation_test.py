import glob
import file_setup
from os_genome import call_pt
import numpy as np
import pandas as pd
import sum_of_squares


# def test_mutation(df, output_dir):
#     """
#     Mutation is tested to either be accepted or rejected based on calculated fitness values or based on a probability to allow for some neutral mutations to be accepted.
#     """
#
#     call_pt.pt_call(output_dir)
#     #genome_simulator.call_pt()
#     nf = pd.read_csv(output_dir+"three_genes_replicated.tsv", header=0, sep='\t')
#     nf = file_setup.rearrange_file(nf)
#     return(sum_of_squares.calc_sum_of_squares(df, nf))


def test_mutation(df, output_dir, mutation_number):

    ss_list = []

    #Creates test files from pinetree to find average number of transcripts at each time
    for i in range(1, mutation_number+1):
        call_pt.pt_call(output_dir)
        save_df = pd.read_csv(output_dir+"three_genes_replicated.tsv", header=0, sep='\t')
        save_df.to_csv(output_dir+"test_{}.tsv".format(i), sep='\t', index=False)
    #Calculates the sum of squares for each test file
    for replicate_file in glob.glob(output_dir + 'test*'):
        nf = pd.read_csv(replicate_file, header=0, sep='\t')
        nf = file_setup.rearrange_file(nf)
        ss_ind = sum_of_squares.calc_sum_of_squares(df, nf)
        ss_list.append(ss_ind)
    #Returns the average sum of squares for all test files
    return np.mean(ss_list)
