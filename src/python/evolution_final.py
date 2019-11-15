import argparse
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


#General setup
ss_old = 1000000
all_sos_list = [ss_old]
sos_iter_list = []
is_accepted = []
accepted = 0
i = 1

#Output file directory structure
year = datetime.date.today().year
month = datetime.date.today().month
day = datetime.date.today().day
output_dir = '../../results/{}_{}_{}/'.format(year, month, day)


#Command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('filename', help='Input target tsv file name with transcript information.')
parser.add_argument('replicate_number', type=int, default=1, nargs='?', help='Input the file number so that folder names are unique.')
parser.add_argument('generation_number', type=int, default=1, nargs='?', help='Input number of generations to run evolutionary program.')
parser.add_argument('replicate_mutation_number', type=int, default=1, nargs='?', help='Input number of times to simulate proposed mutation.')
parser.add_argument('N', type=int, default=10, nargs='?', help='Input desired effective population size value.')
parser.add_argument('beta', type=float, default=0.001, nargs='?', help='Input desired beta value based on starting sum of squares value for fitness calculation.')
args = parser.parse_args()

output_dir = output_dir + '{}_nf{}_rep{}_nmut{}/'.format(args.filename.strip('.tsv'), args.N, args.replicate_number, args.replicate_mutation_number)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

#Target file inputted as dataframe
df = pd.read_csv('../../data/'+args.filename, header=0, sep='\t')
df = file_setup.rearrange_file(df)

#Opens yaml files containing genome coordinates
starting_file = output_dir + 'new_gene.yml'
initialize_yaml.create_yaml(starting_file)
with open(starting_file, 'r') as gene_parameters:
    genome_tracker_new = yaml.safe_load(gene_parameters)
    genome_tracker_old = yaml.safe_load(gene_parameters)

#Start of evolution program
while i <= args.generation_number:

    #Mutation is chosen and performed on best genome
    possibilities = [mutation_choices.modify_promoter, mutation_choices.modify_rnase, mutation_choices.modify_terminator]
    random.choice(possibilities)(genome_tracker_new, starting_file)
    with open(starting_file, 'r') as current:
        genome_tracker_new = yaml.safe_load(current)

    #Sum of squares is calculated and the mutation is accepted or rejected based off of its calculated fitness value

    #PineTree called in test_mutation
    ss_new = mutation_test.test_mutation(df, output_dir, args.replicate_mutation_number)

    accept_prob = fitness_score.calc_fitness(ss_new, ss_old, args.N, args.beta)
    all_sos_list.append(ss_new)
    sos_iter_list.append(i)
    if accept_prob > random.random():
        #If accepted the old genome is replace by the new genome
        genome_tracker_old = genome_tracker_new
        with open(output_dir+'gene_{}.yml'.format(i), 'w') as save_yaml:
            yaml.dump(genome_tracker_old, save_yaml)
        save_df = pd.read_csv(output_dir+"three_genes_replicated.tsv", header=0, sep='\t')
        save_df.to_csv(output_dir+"three_genes_replicated_{}.tsv".format(i), sep='\t', index=False)
        ss_old = ss_new
        is_accepted.append("yes")
    else:
        #If mutation is rejected then the new genome becomes the last accepted genome
        genome_tracker_new = genome_tracker_old
        is_accepted.append("no")

    i+=1


all_sos_list = all_sos_list[1:]
sos_dataframe = pd.DataFrame(data=zip(sos_iter_list, all_sos_list, is_accepted), columns=["Iteration", "Sum_of_Squares", "Accepted"])
export_csv = sos_dataframe.to_csv(output_dir + 'sos_data.tsv', index=False, sep='\t')
