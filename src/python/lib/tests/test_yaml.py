#Common imports
import filecmp
import os
import pytest
import yaml

#lib imports
from ..initialize_yaml import create_yaml

def test_initialize_yaml():

	#Verify that the yaml files are properly configured
    with open(os.path.dirname(os.path.abspath(__file__))+'/inputs/initial_test.yml', 'r') as gene_parameters:
        gene_file = yaml.safe_load(gene_parameters)
    starting_file = os.path.dirname(os.path.abspath(__file__))+'/inputs/initial_test_new.yml'
    create_yaml(starting_file, gene_file)
    assert filecmp.cmp(os.path.dirname(os.path.abspath(__file__))+'/inputs/initial_test_new.yml', os.path.dirname(os.path.abspath(__file__))+'/inputs/initial_test_compare.yml')
    os.remove(os.path.dirname(os.path.abspath(__file__))+'/inputs/initial_test_new.yml')
