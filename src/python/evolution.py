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
import yaml

#lib imports
import file_setup
import fitness_score
import genome_simulator
import initialize_yaml
import mutation_choices
import mutation_analysis
import root_mean_square_error


class Command_line_args(object):
    """
    This class contains all the arguments the user inputs for the class to run.
    Input(s):
    No other inputs needed.
    Output(s):
    None.
    """

    def __init__(self):

        #Command line arguments
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument('target_transcript_data', help='Input target tsv file name with transcript information.')
        self.parser.add_argument('genome_config', help='Input genome parameters file name.')
        self.parser.add_argument('run_number', type=int, default=1, nargs='?', help='Input the file number so that folder names are unique.')
        self.parser.add_argument('generation_number', type=int, default=1, nargs='?', help='Input number of generations to run evolutionary program.')
        self.parser.add_argument('replicate_mutation_number', type=int, default=1, nargs='?', help='Input number of times to simulate proposed mutation.')
        self.parser.add_argument('dynamic_deg_rate', type=bool, default='', nargs='?', help='Input \'True\' if rnase site degredation rate should be an option to be modified or leave blank if not.')
        self.parser.add_argument('progress_bar_out', type=bool, default='', nargs='?', help='Input \'True\' if the progress bar should be outputted or leave blank if not.')
        self.args = self.parser.parse_args()


def organize_output_dir(arguments):
    """
    Set up output directories for where the genome files can be stored.
    Input(s):
    The only argument passes is the class containing the command line arguments the user previously inputted when running the program.
    Output(s):
    output_dir is a string containing information of the path to the directory in which all the saved files are stored by the program.
    """

    #Output file directory structure
    year = datetime.date.today().year
    month = datetime.date.today().month
    day = datetime.date.today().day
    output_dir = '../../results/{}_{}_{}/'.format(year, month, day)

    output_dir = output_dir + '{}_rep{}_nmut{}/'.format(arguments.args.target_transcript_data.strip('.tsv'), arguments.args.run_number, arguments.args.replicate_mutation_number)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    final_data_dir = output_dir+'/final'
    if not os.path.exists(final_data_dir):
        os.makedirs(final_data_dir)

    return output_dir

def setup_configuration_files(output_dir, arguments):
    """
    Sets up all the files the program needs to run, such as the yaml file containing all the genome's information.
    Input(s):
    output_dir is a string containing information of the path to the directory in which all the saved files are stored by the program.
    arguments is the class containing all the command line arguments the user previously inputted when running the program.
    Output(s):
    genome_tracker_new is the dataframe containing the most recently edited genomic data.
    genome_tracker_old is the dataframe containing the previously edited genomic data.
    """

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

def sim_initial_genome(output_dir, genome_tracker_new, target_df, arguments):
    """
    Simulates the template genome containng only one promoter upstream of the first gene and all the user specified genes.
    Input(s):
    output_dir is a string containing information of the path to the directory in which all the saved files are stored by the program.
    genome_tracker_new is the dataframe containing the most recently edited genomic data.
    target_df is the user-inputted tsv dataframe containing transcript abundances for each gene.
    arguments is the class containing all the command line arguments the user previously inputted when running the program.
    Output(s):
    rmse_old is the float that refers to the root mean square error value related to genome_tracker_old.
    """

    dfs = []

    #Template genome is simulated and its gene expression pattern is compared to the target
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
    orig_df = file_setup.rearrange_file(orig_file, target_df.index[-1], genome_tracker_new['num_genes'])
    rmse_old = root_mean_square_error.calc_nrmse(target_df, orig_df)

    return rmse_old

def enumerate_mutation_options(genome_tracker_new, dynamic_deg_rate):
    """
    All the mutation possibilities are added to a dictionary.
    Input(s):
    genome_tracker_new is the dataframe containing the most recently edited genomic data.
    dynamic_deg_rate is a command line argument that specifies if rnase degredation rates should be individually specified or not.
    Output(s):
    possibilities is a dictionary containing all the information regarding which mutations are possible to occur.
    """

    possibilities = {}
    modify_possibilities = {}

    #Possible addition and removal mutation choices are enumerated
    for gene in range(genome_tracker_new['num_genes']+1):
        promoter = 'promoter_{}'.format(gene)
        terminator = 'terminator_{}'.format(gene)
        rnase = 'rnase_{}'.format(gene)
        #If the region upstream of the first gene and region downstream of the last gene are not selected
        if gene != 0 and gene != genome_tracker_new['num_genes']:
            if genome_tracker_new[promoter]['start'] == 0:
                possibilities[promoter+'.add'] = 'add'
            elif genome_tracker_new[promoter]['start'] > 0:
                possibilities[promoter+'.remove'] = 'remove'
        #If the region upstream of the first gene is not selected
        if gene != 0:
            if genome_tracker_new[terminator]['start'] == 0:
                possibilities[terminator+'.add'] = 'add'
            elif genome_tracker_new[terminator]['start'] > 0:
                possibilities[terminator+'.remove'] = 'remove'
        #If the region downstream of the last gene is not selected
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
        #If the region upstream of the first gene is not selected
        if gene != 0:
            if genome_tracker_new[terminator]['start'] > 0:
                for term in range(len(possibilities)*13):
                    modify_possibilities[terminator+'.modify{}'.format(term)] = 'modify'
        #If the region downstream of the last gene is not selected
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

def progress_bar(count, total, replicates, status=''):
    """
    Displays progress bar within terminal.
    Adapted from Vladimir Ignatev and given permission to use through a free software license.
    Input(s):
    count is the integer that refers to the generation that the simulation is on.
    total is the integer that refers to the total number of generation dictated by the user.
    replicates is the integer that refers to the number of times each architecture is simulated.
    status is the string displayed next to the progress bar.
    Output(s):
    Prints out a progress bar in the terminal.
    """

    #Setting progress bar length
    bar_len = int(os.get_terminal_size()[0] * 0.6)
    filled_len = int(round(bar_len * count / float(total)))

    #percents = round(100.0 * count / float(total), 1)
    #Bar displays equals sign when progressing
    bar = '=' * filled_len + '-' * (bar_len - filled_len)
    #Calculates how much time is left until the simulation is complete
    time_left = round(((total*replicates*0.47) - (count*replicates*0.43)) / 60, 2)

    #Writes out the progress bar information
    sys.stdout.write('%s: [%s] %.2fmin\r' % (status, bar, time_left))
    sys.stdout.flush()

    return

def run_evolution(output_dir, genome_tracker_new, genome_tracker_old, target_df, arguments, all_rmse_list):
    """
    The main part of the evolution program where the genome is modified with each generation and determine if the mutation should be acccepted.
    Input(s):
    output_dir is a string containing information of the path to the directory in which all the saved files are stored by the program.
    genome_tracker_new is the dataframe containing the most recently edited genomic data.
    genome_tracker_old is the dataframe containing the previously edited genomic data.
    target_df is the user-inputted tsv dataframe containing thranscript abundances for each gene.
    arguments is the class containing all the command line arguments the user previously inputted when running the program.
    all_rmse_list is a list containing each root mean square error value calculated.
    Output(s):
    genome_tracker_new is the dataframe containing the most recently edited genomic data.
    all_rmse_list is a list containing each root mean square error value calculated.
    rmse_iter_list is a list containing the number of each generation that corresponds to each root mean square error value.
    is_accepted is a list containing information regarding if a mutation was accepted or not.
    """

    #General setup
    rmse_iter_list = [0]
    is_accepted = ['yes']
    i = 1
    generation_number = arguments.args.generation_number
    rmse_old = all_rmse_list[0]

    #Start of evolution program
    while i <= generation_number:

        #Progress bar is printed out as each generation progresses
        if arguments.args.progress_bar_out:
            progress_bar(i, generation_number, arguments.args.replicate_mutation_number, 'Simulation running')
        elif i == 1:
            print('Running simulation...\nEstimated time:', round(((generation_number*arguments.args.replicate_mutation_number*0.47) - (i*arguments.args.replicate_mutation_number*0.47)) / 60, 2))

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

        #pinetree called in analyze_mutation
        rmse_new = mutation_analysis.analyze_mutation(genome_tracker_new, output_dir, target_df, arguments.args.replicate_mutation_number, 1, arguments.args.dynamic_deg_rate)

        accept_prob = fitness_score.calc_fitness(rmse_new, rmse_old, generation_number, i)
        all_rmse_list.append(rmse_new)
        rmse_iter_list.append(i)
        if accept_prob > random.random():
            #If accepted, the old genome is replace by the new genome and files for the new genome are saved
            genome_tracker_old = copy.deepcopy(genome_tracker_new)
            with open(output_dir+'gene_{}.yml'.format(i), 'w') as save_yaml:
                yaml.dump(genome_tracker_old, save_yaml)
            save_df = pd.read_csv(output_dir+"expression_pattern.tsv", header=0, sep='\t')
            save_df.to_csv(output_dir+"expression_pattern_{}.tsv".format(i), sep='\t', index=False)
            rmse_old = copy.copy(rmse_new)
            is_accepted.append("yes")
        else:
            #If mutation is rejected then the new genome becomes the last accepted genome
            genome_tracker_new = copy.deepcopy(genome_tracker_old)
            is_accepted.append("no")

        i+=1

    return (genome_tracker_old, all_rmse_list, rmse_iter_list, is_accepted)

def save_final_data(output_dir, genome_tracker_old, arguments, target_df, all_rmse_list, rmse_iter_list, is_accepted, max_rmse):
    """
    Saves the sum of squared error data, final/best genomic architecture it found, and removes any unnecessary files unrelated to the output.
    Input(s):
    output_dir is a string containing information of the path to the directory in which all the saved files are stored by the program.
    genome_tracker_old is the dataframe containing the most recently accepted genomic data.
    arguments is the class containing all the command line arguments the user previously inputted when running the program.
    target_df is the user-inputted tsv dataframe containing thranscript abundances for each gene.
    all_rmse_list is a list containing each root mean square error value calculated.
    rmse_iter_list is a list containing the number of each generation that corresponds to each root mean square error value.
    is_accepted is a list containing information regarding if a mutation was accepted or not.
    max_rmse is a float that refers to the highest root means square error value the program deems as a successfully found genomic architecture.
    Output(s):
    found_arch is a boolean that that flags when a suitable genome architecture has been found.
    Saves all the relevant files to the output directory.
    """

    found_arch = False

    #Genome is tested to determine if any elements on the genome significantly alter the expression pattern produced when deleted
    print('\nCleaning up genome architecture...')
    rmse_best = mutation_choices.cleanup_genome(output_dir, genome_tracker_old, target_df, arguments.args.replicate_mutation_number, arguments.args.dynamic_deg_rate)
    
    #Sum of squared error data is saved to output directory
    rmse_iter_list.append(arguments.args.generation_number+1)
    all_rmse_list.append(rmse_best)
    is_accepted.append('yes')
    rmse_dataframe = pd.DataFrame(data=zip(rmse_iter_list, all_rmse_list, is_accepted), columns=["Iteration", "NRMSE", "Accepted"])
    export_csv = rmse_dataframe.to_csv(output_dir + 'final/rmse_data.tsv', index=False, sep='\t')
    os.remove(output_dir+'expression_pattern.tsv')

    #If the genome architecture pattern produces an expression at least 90% similar to the target then only run for 500 more generations
    if rmse_best <= max_rmse:
        found_arch = True


    return found_arch

def main():
    """
    Main function that contains all other functions in order so that the program can run.
    Input(s):
    No inputs necessary, however command line arguments are needed for the program to run.
    Output(s):
    Prints out if the simulation was successful or not.
    """

    print ("\r\x1b[8;25;90t\r")

    #Class containing command line arguments are called
    arguments = Command_line_args()
    print("Configuring simulation...")
    #Output directory structure is organized
    output_dir = organize_output_dir(arguments)
    #Yaml files containing genome information are setup
    genome_trackers = setup_configuration_files(output_dir, arguments)
    genome_tracker_new = genome_trackers[0]
    genome_tracker_old = genome_trackers[1]

    #Target file inputted as dataframe
    target_file = pd.read_csv('../../data/targets/'+arguments.args.target_transcript_data, header=0, sep='\t')
    target_df = file_setup.rearrange_file(target_file, target_file.iloc[-1]['time'], genome_tracker_new['num_genes'])
    #Determines the highest RMSE value allowed before finding and accepting a suitable architecture
    max_rmse = root_mean_square_error.calc_accepted_rmse_range(target_df)

    #Creates test files from pinetree to find average number of transcripts at generation 0
    #Template genome is simulated and compared against target dataframe
    rmse_old = sim_initial_genome(output_dir, genome_tracker_new, target_df, arguments)
    all_rmse_list = [rmse_old]

    #Evolution program run and mutations are tested
    evo_vars = run_evolution(output_dir, genome_tracker_new, genome_tracker_old, target_df, arguments, all_rmse_list)
    genome_tracker_old = evo_vars[0]
    all_rmse_list = evo_vars[1]
    rmse_iter_list = evo_vars[2]
    is_accepted = evo_vars[3]

    #Files containing sum of squared data and genome information are saved in output directory
    found = save_final_data(output_dir, genome_tracker_old, arguments, target_df, all_rmse_list, rmse_iter_list, is_accepted, max_rmse)
    if found:
        print('Simulation successfully found an architecture!')
    else:
        print('Simulation did not find a sound architecture.')


if __name__ == '__main__':
    main()
