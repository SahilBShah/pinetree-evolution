#Common imports
import sys

sys.path.insert(1, './lib/')

import argparse
import copy
import datetime
import math
import numpy as np
import os
import pandas as pd
import random

#lib imports
import file_setup
import fitness_score
import genome_simulator
import initialize_yaml
import mutation_choices
import mutation_analysis
import sum_of_squares
import yaml


class Command_line_args(object):

    def __init__(self):

        #Command line arguments
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument('target_transcript_data', help='Input target tsv file name with transcript information.')
        self.parser.add_argument('genome_config', help='Input genome parameters file name.')
        self.parser.add_argument('run_number', type=int, default=1, nargs='?', help='Input the file number so that folder names are unique.')
        self.parser.add_argument('generation_number', type=int, default=1, nargs='?', help='Input number of generations to run evolutionary program.')
        self.parser.add_argument('replicate_mutation_number', type=int, default=1, nargs='?', help='Input number of times to simulate proposed mutation.')
        self.parser.add_argument('N', type=float, default=1.0, nargs='?', help='Input desired effective population size value.')
        self.parser.add_argument('beta', type=float, default=0.001, nargs='?', help='Input desired beta value based on starting sum of squares value for fitness calculation.')
        self.parser.add_argument('dynamic_deg_rate', type=bool, default='', nargs='?', help='Input \'True\' if rnase site degredation rate should be an option to be modified or leave blank if not.')
        self.args = self.parser.parse_args()


def organize_output_dir(arguments):

    #Output file directory structure
    year = datetime.date.today().year
    month = datetime.date.today().month
    day = datetime.date.today().day
    output_dir = '../../results/{}_{}_{}/'.format(year, month, day)

    output_dir = output_dir + '{}_nf{}_rep{}_nmut{}/'.format(arguments.args.target_transcript_data.strip('.tsv'), arguments.args.N, arguments.args.run_number, arguments.args.replicate_mutation_number)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    final_data_dir = output_dir+'/final'
    if not os.path.exists(final_data_dir):
        os.makedirs(final_data_dir)

    return output_dir

def setup_configuration_files(output_dir, arguments):

    #Opens yaml files containing genome coordinates
    starting_genome = output_dir + 'config.yml'
    with open('../../data/gene_parameters/'+arguments.args.genome_config, 'r') as gene_parameters:
        genome_config = yaml.safe_load(gene_parameters)
    initialize_yaml.create_yaml(starting_genome, genome_config)
    with open(starting_genome, 'r') as gene_elements:
        genome_tracker_new = yaml.safe_load(gene_elements)
    genome_tracker_old = copy.deepcopy(genome_tracker_new)
    with open(output_dir+'gene_0.yml', 'w') as save_yaml:
        yaml.dump(genome_tracker_old, save_yaml)

    return (genome_tracker_new, genome_tracker_old)

def sim_initial_genome(output_dir, genome_tracker_new, target_file, arguments):

    dfs = []

    for x in range(1, arguments.args.replicate_mutation_number+1):
        if arguments.args.dynamic_deg_rate:
            genome_simulator.pt_call_alt(output_dir, genome_tracker_new)
        else:
            genome_simulator.pt_call(output_dir, genome_tracker_new)
        save_df = pd.read_csv(output_dir+"expression_pattern.tsv", header=0, sep='\t')
        save_df['time'] = save_df['time'].round().astype(int)
        dfs.append(save_df)
    #Averages all the values in each file and creates a new file with those averages
    df_concat = pd.concat(dfs)
    df_gb = df_concat.groupby(['time', 'species'], as_index=False)
    df_mean = df_gb.sum()
    df_mean[['protein', 'transcript', 'ribo_density']] = df_mean[['protein', 'transcript', 'ribo_density']] / arguments.args.replicate_mutation_number
    df_mean.to_csv(output_dir+'expression_pattern_0.tsv', sep='\t', index=False)
    orig_file = pd.read_csv(output_dir+'expression_pattern_0.tsv', header=0, sep='\t')
    orig_file = file_setup.rearrange_file(orig_file, genome_tracker_new)
    ss_old = sum_of_squares.calc_sse(target_file, orig_file, genome_tracker_new['num_genes'])

    return ss_old

def enumerate_mutation_options(genome_tracker_new, dynamic_deg_rate):
    """
    All the mutation possibilities are added to a dictionary.
    """

    possibilities = {}
    modify_possibilities = {}

    #Possible addition and removal mutation choices are enumerated
    for gene in range(genome_tracker_new['num_genes']+1):
        promoter = 'promoter_{}'.format(gene)
        terminator = 'terminator_{}'.format(gene)
        rnase = 'rnase_{}'.format(gene)
        if gene != 0 and gene != genome_tracker_new['num_genes']:
            if genome_tracker_new[promoter]['start'] == 0:
                possibilities[promoter+'.add'] = 'add'
            elif genome_tracker_new[promoter]['start'] > 0:
                possibilities[promoter+'.remove'] = 'remove'
        if gene != 0:
            if genome_tracker_new[terminator]['start'] == 0:
                possibilities[terminator+'.add'] = 'add'
            elif genome_tracker_new[terminator]['start'] > 0:
                possibilities[terminator+'.remove'] = 'remove'
        if gene != genome_tracker_new['num_genes']:
            if genome_tracker_new[rnase]['start'] == 0:
                possibilities[rnase+'.add'] = 'add'
            elif genome_tracker_new[rnase]['start'] > 0:
                possibilities[rnase+'.remove'] = 'remove'
    #Enumerate modification possibilities
    for gene in range(genome_tracker_new['num_genes']+1):
        promoter = 'promoter_{}'.format(gene)
        terminator = 'terminator_{}'.format(gene)
        rnase = 'rnase_{}'.format(gene)
        if gene != 0:
            if genome_tracker_new[terminator]['start'] > 0:
                for term in range(len(possibilities)*13):
                    modify_possibilities[terminator+'.modify{}'.format(term)] = 'modify'
        if gene != genome_tracker_new['num_genes']:
            if genome_tracker_new[promoter]['start'] > 0:
                for prom in range(len(possibilities)*13):
                    modify_possibilities[promoter+'.modify{}'.format(prom)] = 'modify'
            if genome_tracker_new[rnase]['start'] > 0:
                for rna in range(len(possibilities)*13):
                    if dynamic_deg_rate:
                        modify_possibilities[rnase+'.modify{}'.format(rna)] = 'modify'
    possibilities.update(modify_possibilities)

    return possibilities


def run_evolution(output_dir, genome_tracker_new, genome_tracker_old, target_file, arguments, all_sse_list, max_sse):

    #General setup
    sse_iter_list = []
    is_accepted = []
    i = 1
    count = 0
    generation_number = arguments.args.generation_number
    ss_old = all_sse_list[0]

    #Start of evolution program
    while i <= generation_number:

        possibilities = enumerate_mutation_options(genome_tracker_new, arguments.args.dynamic_deg_rate)

        #Mutation is chosen and performed on best genome
        item = random.choice(list(possibilities.keys()))
        if possibilities[item] == 'add':
            genome_tracker_new = mutation_choices.add_element(genome_tracker_new, output_dir, genome_tracker_new['num_genes'], arguments.args.dynamic_deg_rate, item.split('.')[0])
        elif possibilities[item] == 'remove':
            genome_tracker_new = mutation_choices.remove_element(genome_tracker_new, output_dir, genome_tracker_new['num_genes'], arguments.args.dynamic_deg_rate, item.split('.')[0])
        else:
            genome_tracker_new = mutation_choices.modify_element(genome_tracker_new, output_dir, genome_tracker_new['num_genes'], arguments.args.dynamic_deg_rate, item.split('.')[0])

        #Sum of squared error is calculated and the mutation is accepted or rejected based off of its calculated fitness value

        #pinetree called in test_mutation
        ss_new = mutation_analysis.analyze_mutation(genome_tracker_new, output_dir, target_file, arguments.args.replicate_mutation_number, arguments.args.dynamic_deg_rate)

        accept_prob = fitness_score.calc_fitness(ss_new, ss_old, arguments.args.N, arguments.args.beta)
        all_sse_list.append(ss_new)
        sse_iter_list.append(i)
        if accept_prob > random.random():
            #If accepted the old genome is replace by the new genome and files for the new genome are saved
            genome_tracker_old = copy.deepcopy(genome_tracker_new)
            with open(output_dir+'gene_{}.yml'.format(i), 'w') as save_yaml:
                yaml.dump(genome_tracker_old, save_yaml)
            save_df = pd.read_csv(output_dir+"expression_pattern.tsv", header=0, sep='\t')
            save_df.to_csv(output_dir+"expression_pattern_{}.tsv".format(i), sep='\t', index=False)
            ss_old = ss_new
            is_accepted.append("yes")
        else:
            #If mutation is rejected then the new genome becomes the last accepted genome
            genome_tracker_new = copy.deepcopy(genome_tracker_old)
            is_accepted.append("no")

        if ss_old <= max_sse and count == 0:
            count = 1
            generation_number+=500
        elif count != 0 and count != 500:
            count+=1
        elif count == 500:
            break

        i+=1
        if i <= generation_number:
            print('generation =', i)

    return (genome_tracker_new, all_sse_list, sse_iter_list, is_accepted)

def save_final_data(output_dir, genome_tracker_new, arguments, target_file, all_sse_list, sse_iter_list, is_accepted):

    all_sse_list = all_sse_list[1:]
    sse_dataframe = pd.DataFrame(data=zip(sse_iter_list, all_sse_list, is_accepted), columns=["Iteration", "SSE", "Accepted"])
    export_csv = sse_dataframe.to_csv(output_dir + 'final/sse_data.tsv', index=False, sep='\t')
    print('Cleaning up genome architecture...')
    mutation_choices.cleanup_genome(output_dir, target_file, sse_dataframe, genome_tracker_new['num_genes'], arguments.args.dynamic_deg_rate)
    os.remove(output_dir+'expression_pattern.tsv')


def main():

    arguments = Command_line_args()
    output_dir = organize_output_dir(arguments)
    genome_trackers = setup_configuration_files(output_dir, arguments)
    genome_tracker_new = genome_trackers[0]
    genome_tracker_old = genome_trackers[1]

    #Target file inputted as dataframe
    target_file = pd.read_csv('../../data/'+arguments.args.target_transcript_data, header=0, sep='\t')
    target_file = file_setup.rearrange_file(target_file, genome_tracker_new)
    #Determines the highest SSE value allowed before finding and accepting a suitable architecture
    max_sse = sum_of_squares.calc_accepted_sse_range(target_file, genome_tracker_new)

    #Creates test files from pinetree to find average number of transcripts at generation 0
    print('generation =', 0)
    ss_old = sim_initial_genome(output_dir, genome_tracker_new, target_file, arguments)
    all_sse_list = [ss_old]

    print('generation =', 1)
    evo_vars = run_evolution(output_dir, genome_tracker_new, genome_tracker_old, target_file, arguments, all_sse_list, max_sse)
    genome_tracker_new = evo_vars[0]
    all_sse_list = evo_vars[1]
    sse_iter_list = evo_vars[2]
    is_accepted = evo_vars[3]

    save_final_data(output_dir, genome_tracker_new, arguments, target_file, all_sse_list, sse_iter_list, is_accepted)


if __name__ == '__main__':
    main()
