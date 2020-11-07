#Common imports
import sys

sys.path.append('../src/python/')

import filecmp
import os
import pytest
import yaml

#lib imports
from lib.initialize_yaml import create_yaml

def test_initialize_yaml():

	#Verify that the yaml files are properly configured
    with open('./inputs/initial_test.yml', 'r') as gene_parameters:
        gene_file = yaml.safe_load(gene_parameters)
    starting_file = './inputs/initial_test_new.yml'
    create_yaml(starting_file, gene_file)
    assert filecmp.cmp('./inputs/initial_test_new.yml', './inputs/initial_test_compare.yml')
    os.remove('./inputs/initial_test_new.yml')
