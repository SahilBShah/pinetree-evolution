#Common imports
import copy
import filecmp
import os
import pandas as pd
import pytest
import yaml

#lib imports
from ..mutation_choices import *

def test_expand_genome():

    #Verify that increasing the size of the genome is working as intended
    with open(os.path.dirname(os.path.abspath(__file__))+'/inputs/testing.yml', 'r') as gene_parameters:
        genome_tracker = yaml.safe_load(gene_parameters)
    genome_tracker_old = copy.deepcopy(genome_tracker)
    region_choice = 2
    element_choice = 'promoter_2'
    genome_shift = 35
    num_genes = 3
    genome_tracker_new = expand_genome(genome_tracker, num_genes, region_choice, genome_shift, element_choice)
    assert genome_tracker_new['length_of_genome'] > genome_tracker_old['length_of_genome']
    assert genome_tracker_new['length_of_genome'] - 35 == genome_tracker_old['length_of_genome']
    assert not (genome_tracker_new['length_of_genome'] < genome_tracker_old['length_of_genome'])
    assert genome_tracker_new['rnase_2']['start'] > genome_tracker_old['rnase_2']['start']

def test_shrink_genome():

    #Verify decreasing the size of the genome is working as intended
    with open(os.path.dirname(os.path.abspath(__file__))+'/inputs/testing.yml', 'r') as gene_parameters:
        genome_tracker = yaml.safe_load(gene_parameters)
    genome_tracker_old = copy.deepcopy(genome_tracker)
    region_choice = 2
    element_choice = 'promoter_2'
    genome_shift = 35
    num_genes = 3
    genome_tracker_new = shrink_genome(genome_tracker, num_genes, region_choice, genome_shift, element_choice)
    assert genome_tracker_new['length_of_genome'] < genome_tracker_old['length_of_genome']
    assert genome_tracker_new['length_of_genome'] + 35 == genome_tracker_old['length_of_genome']
    assert not (genome_tracker_new['length_of_genome'] > genome_tracker_old['length_of_genome'])
    assert genome_tracker_new['rnase_2']['start'] < genome_tracker_old['rnase_2']['start']

def test_add_element():

    #Verify the insertion mutation is accurate
    with open(os.path.dirname(os.path.abspath(__file__))+'/inputs/testing.yml', 'r') as gene_parameters:
        genome_tracker = yaml.safe_load(gene_parameters)
    genome_tracker_old = copy.deepcopy(genome_tracker)
    output_dir = os.path.dirname(os.path.abspath(__file__))+'/inputs/'
    num_genes = 3
    element_choice = 'promoter_1'
    genome_tracker_new = add_element(genome_tracker, output_dir, num_genes, element_choice)
    with open(os.path.dirname(os.path.abspath(__file__))+'/inputs/testing_tmp_comp.yml', 'w') as save_yaml:
        yaml.dump(genome_tracker_new, save_yaml)
    assert not filecmp.cmp(os.path.dirname(os.path.abspath(__file__))+'/inputs/testing_tmp_comp.yml', os.path.dirname(os.path.abspath(__file__))+'/inputs/testing.yml')
    assert genome_tracker_new['promoter_1']['start'] > 0
    assert genome_tracker_new['promoter_1']['stop'] > 0
    assert genome_tracker_new['promoter_1']['current_strength'] > 0
    assert genome_tracker_new['length_of_genome'] > genome_tracker_old['length_of_genome']
    assert genome_tracker_new['rnase_2']['start'] != genome_tracker_old['rnase_2']['start']
    assert genome_tracker_new['rnase_2']['start'] <= genome_tracker_old['rnase_2']['start'] + 35
    os.remove(os.path.dirname(os.path.abspath(__file__))+'/inputs/testing_tmp_comp.yml')

def test_remove_element():

    #Verify the deletion mutation is accurate
    with open(os.path.dirname(os.path.abspath(__file__))+'/inputs/testing.yml', 'r') as gene_parameters:
        genome_tracker = yaml.safe_load(gene_parameters)
    genome_tracker_old = copy.deepcopy(genome_tracker)
    output_dir = os.path.dirname(os.path.abspath(__file__))+'/inputs/'
    num_genes = 3
    element_choice = 'promoter_2'
    genome_tracker_new = remove_element(genome_tracker, output_dir, num_genes, element_choice)
    with open(os.path.dirname(os.path.abspath(__file__))+'/inputs/testing_tmp_comp.yml', 'w') as save_yaml:
        yaml.dump(genome_tracker_new, save_yaml)
    assert not filecmp.cmp(os.path.dirname(os.path.abspath(__file__))+'/inputs/testing_tmp_comp.yml', os.path.dirname(os.path.abspath(__file__))+'/inputs/testing.yml')
    assert genome_tracker_new['promoter_2']['start'] == 0
    assert genome_tracker_new['promoter_2']['stop'] == 0
    assert genome_tracker_new['promoter_2']['current_strength'] == 0
    assert genome_tracker_new['promoter_2']['previous_strength'] > 0
    assert genome_tracker_new['length_of_genome'] < genome_tracker_old['length_of_genome']
    assert genome_tracker_new['rnase_2']['start'] != genome_tracker_old['rnase_2']['start']
    assert genome_tracker_new['rnase_2']['start'] == genome_tracker_old['rnase_2']['start'] - 35
    os.remove(os.path.dirname(os.path.abspath(__file__))+'/inputs/testing_tmp_comp.yml')

def test_modify_element():

    #Verify the modificaiton of an element's strength is accurate
    with open(os.path.dirname(os.path.abspath(__file__))+'/inputs/testing.yml', 'r') as gene_parameters:
        genome_tracker = yaml.safe_load(gene_parameters)
    genome_tracker_old = copy.deepcopy(genome_tracker)
    output_dir = os.path.dirname(os.path.abspath(__file__))+'/inputs/'
    num_genes = 3
    element_choice = 'promoter_0'
    genome_tracker_new = modify_element(genome_tracker, output_dir, num_genes, element_choice)
    with open(os.path.dirname(os.path.abspath(__file__))+'/inputs/testing_tmp_comp.yml', 'w') as save_yaml:
        yaml.dump(genome_tracker_new, save_yaml)
    assert not filecmp.cmp(os.path.dirname(os.path.abspath(__file__))+'/inputs/testing_tmp_comp.yml', os.path.dirname(os.path.abspath(__file__))+'/inputs/testing.yml')
    assert genome_tracker_new['promoter_0']['current_strength'] != genome_tracker_old['promoter_0']['current_strength']
    assert genome_tracker_new['promoter_0']['current_strength'] != genome_tracker_new['promoter_0']['previous_strength']
    assert genome_tracker_new['promoter_2']['current_strength'] == genome_tracker_old['promoter_2']['current_strength']
    assert genome_tracker_new['rnase_2']['start'] == genome_tracker_old['rnase_2']['start']
    os.remove(os.path.dirname(os.path.abspath(__file__))+'/inputs/testing_tmp_comp.yml')
