import argparse
import copy
import datetime
import file_setup
import fitness_score
import initialize_yaml
import mutation_choices
import mutation_test
import numpy as np
import os
import sum_of_squares
import pandas as pd
import random
import yaml

#Command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('filename', help='Input target tsv file name with transcript information.')
parser.add_argument('gene_file', help='Input target gene parameters file name.')
parser.add_argument('run_number', type=int, default=1, nargs='?', help='Input the file number so that folder names are unique.')
parser.add_argument('generation_number', type=int, default=1, nargs='?', help='Input number of generations to run evolutionary program.')
parser.add_argument('replicate_mutation_number', type=int, default=1, nargs='?', help='Input number of times to simulate proposed mutation.')
parser.add_argument('initial_sum_of_squares', type=int, default=1000000, nargs='?', help='Input the desired initial sum of squares value')
parser.add_argument('N', type=int, default=10, nargs='?', help='Input desired effective population size value.')
parser.add_argument('beta', type=float, default=0.001, nargs='?', help='Input desired beta value based on starting sum of squares value for fitness calculation.')
parser.add_argument('dynamic_deg_rate', type=bool, default='', nargs='?', help='Input \'True\' if rnase site degredation rate should be an option to be modified or leave blank if not.')
args = parser.parse_args()

#General setup
ss_old = args.initial_sum_of_squares
all_sos_list = [ss_old]
sos_iter_list = []
is_accepted = []
i = 1

#Output file directory structure
year = datetime.date.today().year
month = datetime.date.today().month
day = datetime.date.today().day
output_dir = '../../results/{}_{}_{}/'.format(year, month, day)

output_dir = output_dir + '{}_nf{}_rep{}_nmut{}/'.format(args.filename.strip('.tsv'), args.N, args.run_number, args.replicate_mutation_number)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

#Opens yaml files containing genome coordinates
starting_file = output_dir + 'config.yml'
with open('../../data/gene_parameters/'+args.gene_file, 'r') as gene_parameters:
    gene_file = yaml.safe_load(gene_parameters)
initialize_yaml.create_yaml(starting_file, gene_file)
with open(starting_file, 'r') as gene_elements:
    genome_tracker_new = yaml.safe_load(gene_elements)
genome_tracker_old = copy.deepcopy(genome_tracker_new)

#Target file inputted as dataframe
df = pd.read_csv('../../data/'+args.filename, header=0, sep='\t')
df = file_setup.rearrange_file(df, genome_tracker_new)

#Start of evolution program
while i <= args.generation_number:

    #Mutation is chosen and performed on best genome
    possibilities = [mutation_choices.add_element, mutation_choices.remove_element, mutation_choices.modify_element]
    genome_tracker_new = random.choice(possibilities)(genome_tracker_new, starting_file, genome_tracker_new['num_genes'], args.dynamic_deg_rate)

    #Sum of squares is calculated and the mutation is accepted or rejected based off of its calculated fitness value

    #PineTree called in test_mutation
    ss_new = mutation_test.test_mutation(df, genome_tracker_new, output_dir, args.replicate_mutation_number, args.dynamic_deg_rate)

    accept_prob = fitness_score.calc_fitness(ss_new, ss_old, args.N, args.beta)
    all_sos_list.append(ss_new)
    sos_iter_list.append(i)
    if accept_prob > random.random():
        #If accepted the old genome is replace by the new genome and files for the new genome are saved
        genome_tracker_old = copy.deepcopy(genome_tracker_new)
        with open(output_dir+'gene_{}.yml'.format(i), 'w') as save_yaml:
            yaml.dump(genome_tracker_old, save_yaml)
        save_df = pd.read_csv(output_dir+"three_genes_replicated.tsv", header=0, sep='\t')
        save_df.to_csv(output_dir+"three_genes_replicated_{}.tsv".format(i), sep='\t', index=False)
        ss_old = ss_new
        is_accepted.append("yes")
    else:
        #If mutation is rejected then the new genome becomes the last accepted genome
        genome_tracker_new = copy.deepcopy(genome_tracker_old)
        is_accepted.append("no")

    i+=1


all_sos_list = all_sos_list[1:]
sos_dataframe = pd.DataFrame(data=zip(sos_iter_list, all_sos_list, is_accepted), columns=["Iteration", "Sum_of_Squares", "Accepted"])
export_csv = sos_dataframe.to_csv(output_dir + 'sos_data.tsv', index=False, sep='\t')
