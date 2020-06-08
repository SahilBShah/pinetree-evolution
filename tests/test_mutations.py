import sys

sys.path.insert(1, '../src/python/lib/')

import copy
import filecmp
import os
import pandas as pd
import pytest
import mutation_choices
import yaml

def test_expand_genome():

    with open('./inputs/testing.yml', 'r') as gene_parameters:
        genome_tracker = yaml.safe_load(gene_parameters)
    genome_tracker_old = copy.deepcopy(genome_tracker)
    region_choice = 2
    element_choice = 'promoter_2'
    genome_shift = 35
    num_genes = 3
    genome_tracker_new = mutation_choices.expand_genome(genome_tracker, num_genes, region_choice, genome_shift, element_choice)
    assert genome_tracker_new['length_of_genome'] > genome_tracker_old['length_of_genome']
    assert genome_tracker_new['length_of_genome'] - 35 == genome_tracker_old['length_of_genome']
    assert not (genome_tracker_new['length_of_genome'] < genome_tracker_old['length_of_genome'])
    assert genome_tracker_new['rnase_2']['start'] > genome_tracker_old['rnase_2']['start']

def test_shrink_genome():

    with open('./inputs/testing.yml', 'r') as gene_parameters:
        genome_tracker = yaml.safe_load(gene_parameters)
    genome_tracker_old = copy.deepcopy(genome_tracker)
    region_choice = 2
    element_choice = 'promoter_2'
    genome_shift = 35
    num_genes = 3
    genome_tracker_new = mutation_choices.shrink_genome(genome_tracker, num_genes, region_choice, genome_shift, element_choice)
    assert genome_tracker_new['length_of_genome'] < genome_tracker_old['length_of_genome']
    assert genome_tracker_new['length_of_genome'] + 35 == genome_tracker_old['length_of_genome']
    assert not (genome_tracker_new['length_of_genome'] > genome_tracker_old['length_of_genome'])
    assert genome_tracker_new['rnase_2']['start'] < genome_tracker_old['rnase_2']['start']

def test_cleanup_genome():

    target_file = pd.read_csv('./inputs/test_compare.tsv', header=0, sep='\t')
    sse_df = pd.read_csv('./inputs/sse_data.tsv', header=0, sep='\t')
    output_dir = './inputs/'
    num_genes = 3
    deg_rate = True
    mutation_choices.cleanup_genome(target_file, sse_df, output_dir, num_genes, deg_rate)
    assert not filecmp.cmp('./inputs/gene_2.yml', './inputs/final/gene_best.yml')
    os.remove('./inputs/expression_pattern.tsv')
    os.remove('./inputs/final/expression_pattern_best.tsv')

def test_add_element():

    with open('./inputs/testing.yml', 'r') as gene_parameters:
        genome_tracker = yaml.safe_load(gene_parameters)
    genome_tracker_old = copy.deepcopy(genome_tracker)
    output_dir = './inputs/'
    num_genes = 3
    deg_rate = True
    element_choice = 'promoter_1'
    genome_tracker_new = mutation_choices.add_element(genome_tracker, output_dir, num_genes, deg_rate, element_choice)
    with open('./inputs/testing_tmp_comp.yml', 'w') as save_yaml:
        yaml.dump(genome_tracker_new, save_yaml)
    assert not filecmp.cmp('./inputs/testing_tmp_comp.yml', './inputs/testing.yml')
    assert genome_tracker_new['promoter_1']['start'] > 0
    assert genome_tracker_new['promoter_1']['stop'] > 0
    assert genome_tracker_new['promoter_1']['current_strength'] > 0
    assert genome_tracker_new['length_of_genome'] > genome_tracker_old['length_of_genome']
    assert genome_tracker_new['rnase_2']['start'] != genome_tracker_old['rnase_2']['start']
    assert genome_tracker_new['rnase_2']['start'] <= genome_tracker_old['rnase_2']['start'] + 35
    os.remove('./inputs/testing_tmp_comp.yml')

def test_remove_element():

    with open('./inputs/testing.yml', 'r') as gene_parameters:
        genome_tracker = yaml.safe_load(gene_parameters)
    genome_tracker_old = copy.deepcopy(genome_tracker)
    output_dir = './inputs/'
    num_genes = 3
    deg_rate = True
    element_choice = 'promoter_2'
    genome_tracker_new = mutation_choices.remove_element(genome_tracker, output_dir, num_genes, deg_rate, element_choice)
    with open('./inputs/testing_tmp_comp.yml', 'w') as save_yaml:
        yaml.dump(genome_tracker_new, save_yaml)
    assert not filecmp.cmp('./inputs/testing_tmp_comp.yml', './inputs/testing.yml')
    assert genome_tracker_new['promoter_2']['start'] == 0
    assert genome_tracker_new['promoter_2']['stop'] == 0
    assert genome_tracker_new['promoter_2']['current_strength'] == 0
    assert genome_tracker_new['promoter_2']['previous_strength'] > 0
    assert genome_tracker_new['length_of_genome'] < genome_tracker_old['length_of_genome']
    assert genome_tracker_new['rnase_2']['start'] != genome_tracker_old['rnase_2']['start']
    assert genome_tracker_new['rnase_2']['start'] == genome_tracker_old['rnase_2']['start'] - 35
    os.remove('./inputs/testing_tmp_comp.yml')

def test_modify_element():

    with open('./inputs/testing.yml', 'r') as gene_parameters:
        genome_tracker = yaml.safe_load(gene_parameters)
    genome_tracker_old = copy.deepcopy(genome_tracker)
    output_dir = './inputs/'
    num_genes = 3
    deg_rate = True
    element_choice = 'promoter_0'
    genome_tracker_new = mutation_choices.modify_element(genome_tracker, output_dir, num_genes, deg_rate, element_choice)
    with open('./inputs/testing_tmp_comp.yml', 'w') as save_yaml:
        yaml.dump(genome_tracker_new, save_yaml)
    assert not filecmp.cmp('./inputs/testing_tmp_comp.yml', './inputs/testing.yml')
    assert genome_tracker_new['promoter_0']['current_strength'] != genome_tracker_old['promoter_0']['current_strength']
    assert genome_tracker_new['promoter_0']['current_strength'] != genome_tracker_new['promoter_0']['previous_strength']
    assert genome_tracker_new['promoter_2']['current_strength'] == genome_tracker_old['promoter_2']['current_strength']
    assert genome_tracker_new['rnase_2']['start'] == genome_tracker_old['rnase_2']['start']
    os.remove('./inputs/testing_tmp_comp.yml')
