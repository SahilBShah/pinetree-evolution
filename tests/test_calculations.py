#Common imports
import sys

sys.path.append('../src/python/')

import math
import os
import pandas as pd
import pytest
import yaml

#lib imports
from lib.fitness_score import calc_fitness
from lib.root_mean_square_error import calc_nrmse
from lib.mutation_analysis import analyze_mutation


def test_calc_nrme():

    #Verify nRMSE is calculated correctly
    assert round(calc_nrmse(pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]]), pd.DataFrame([[9, 8, 7], [6, 5, 4], [3, 2, 1]])), 2) == 1.06
    assert round(calc_nrmse(pd.DataFrame([[10, 20, 30], [40, 50, 60], [70, 80, 90]]), pd.DataFrame([[90, 80, 70], [60, 50, 40], [30, 20, 10]])), 2) == 1.06
    assert round(calc_nrmse(pd.DataFrame([[10, 20, 30], [40, 50, 60], [70, 80, 90]]), pd.DataFrame([[90, 80, 70], [60, 50, 40], [30, 20, 10]])), 2) != 0.55
    assert round(calc_nrmse(pd.DataFrame([[35, 99, 106], [32, 76, 83], [22, 5, 68]]), pd.DataFrame([[49, 53, 110], [12, 88, 91], [70, 80, 55]])), 2) == 0.67

def test_calc_fitness():

    #Verify fitness score is calculated correctly
    assert calc_fitness(1.05, 2.00, 5000, 1000) == 1.0
    assert calc_fitness(3.15, 1.25, 5000, 1000) != 1.0

def test_mutation_analysis():

    #Testing if nRMSE value is returned from mutation_analysis
    with open('./inputs/testing.yml', 'r') as gene_parameters:
        genome_tracker = yaml.safe_load(gene_parameters)
    df = pd.read_csv('./inputs/pt_expression_data_compare.tsv', header=0, sep='\t')
    df = df.set_index([pd.Index([i for i in range(1, 251)])])
    mutation_number = 1
    nrmse = analyze_mutation(genome_tracker, './', df, mutation_number)
    assert nrmse != 3.10
    os.remove('expression_pattern.tsv')
