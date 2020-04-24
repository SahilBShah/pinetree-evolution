import numpy as np
import os
import pandas as pd
import pinetree as pt
import random


def calc_average_run(mutation_number, element, beg, end, increment):
    """
    Averages across pinetree simulations and saves file.
    """
    for strength in np.arange(beg, end, increment):

        dfs = []

        #Creates test files from pinetree to find average number of transcripts at each time
        for i in range(1, mutation_number+1):
            pt_call(strength, element)
            save_df = pd.read_csv("three_genes_replicated.tsv", header=0, sep='\t')
            save_df['time'] = save_df['time'].round().astype(int)
            if 'promoter' in element:
                save_df.to_csv('../../data/parameter_testing/promoter/{}_test_{}_{}.tsv'.format(element, strength, i), sep='\t', index=False)
            elif 'terminator' in element:
                save_df.to_csv('../../data/parameter_testing/terminator/{}_test_{}_{}.tsv'.format(element, strength, i), sep='\t', index=False)
            elif 'rnase' in element:
                save_df.to_csv('../../data/parameter_testing/rnase/{}_test_{}_{}.tsv'.format(element, strength, i), sep='\t', index=False)
            dfs.append(save_df)

        #Averages all the values in each file and creates a new file with those averages
        df_concat = pd.concat(dfs)
        df_gb = df_concat.groupby(['time', 'species'], as_index=False)
        df_mean = df_gb.sum()
        df_mean[['protein', 'transcript', 'ribo_density']] = df_mean[['protein', 'transcript', 'ribo_density']] / mutation_number
        if 'promoter' in element:
            df_mean.to_csv('../../data/parameter_testing/promoter/{}_average_test_{}.tsv'.format(element, strength), sep='\t', index=False)
        elif 'terminator' in element:
            df_mean.to_csv('../../data/parameter_testing/terminator/{}_average_test_{}.tsv'.format(element, strength), sep='\t', index=False)
        elif 'rnase' in element:
            df_mean.to_csv('../../data/parameter_testing/rnase/{}_average_test_{}.tsv'.format(element, strength), sep='\t', index=False)
        os.remove("three_genes_replicated.tsv")

def pt_call(strength, element):

    #Creating starting three genes sequence
    sim = pt.Model(cell_volume=8e-16)
    sim.seed(random.randint(0, 10e6))
    sim.add_polymerase(name="rnapol", copy_number=4, speed=40, footprint=10)
    sim.add_ribosome(copy_number=100, speed=30, footprint=10)

    plasmid = pt.Genome(name="plasmid", length=450,
                        transcript_degradation_rate_ext=1e-2,
                        rnase_speed=20,
                        rnase_footprint=10)

    plasmid.add_promoter(name="p1", start=1, stop=10,
                         interactions={"rnapol": 10e7})
    if 'promoter' in element:
        plasmid.add_promoter(name="p2", start=1, stop=10,
                             interactions={"rnapol": float('10e{}'.format(int(strength)))})

    if 'terminator' in element:
        plasmid.add_terminator(name="t1", start=134, stop=135,
                               efficiency={"rnapol": strength})
    if 'rnase' in element:
        plasmid.add_rnase_site(name='r1', start=134, stop=144, rate=strength)

    plasmid.add_gene(name="protein1", start=26, stop=121,
                     rbs_start=11, rbs_stop=26, rbs_strength=1e7)

    plasmid.add_gene(name="protein2", start=159, stop=280,
                     rbs_start=144, rbs_stop=159, rbs_strength=1e7)
    plasmid.add_gene(name="protein3", start=319, stop=449,
                     rbs_start=303, rbs_stop=318, rbs_strength=1e7)
    sim.register_genome(plasmid)
    sim.simulate(time_limit=240, time_step=1,
                 output = "three_genes_replicated.tsv")

def main():
    mutation_number = int(input("Please enter the number to average the simulations over: "))
    element = input("Please enter the element you wish to modify: \'promoter\', \'terminator\', or \'rnase\': ")
    while element != 'promoter' and element != 'terminator' and element != 'rnase':
        element = input("Please enter the element you wish to modify: \'promoter\', \'terminator\', or \'rnase\': ")
    start = float(input("Please enter the desired strength to start with: "))
    stop = float(input("Please enter the desired strength to end with: "))
    increment = float(input("Please enter the increment in which to increase the strength by: "))
    calc_average_run(mutation_number, element, start, stop, increment)

main()
