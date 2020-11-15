#Common imports
import filecmp
import os
import pandas as pd
import pytest
import yaml

#lib imports
from ..file_setup import rearrange_file

def test_rearrange_file():

    #Verify dataframe formatting is correct
    file_comp = pd.read_csv(os.path.dirname(os.path.abspath(__file__))+"/inputs/pt_expression_data_compare.tsv", header=0, sep='\t')
    file_new = pd.read_csv(os.path.dirname(os.path.abspath(__file__))+"/inputs/pt_expression_data.tsv", header=0, sep='\t')
    num_genes = 3
    max_time = 250
    with open(os.path.dirname(os.path.abspath(__file__))+'/inputs/testing.yml', 'r') as gene_parameters:
        genome_tracker = yaml.safe_load(gene_parameters)
    file_new = rearrange_file(file_new, max_time, num_genes)
    file_new.to_csv(os.path.dirname(os.path.abspath(__file__))+"/inputs/pt_expression_data_new.tsv", sep='\t', index=False)
    assert filecmp.cmp(os.path.dirname(os.path.abspath(__file__))+'/inputs/pt_expression_data_new.tsv', os.path.dirname(os.path.abspath(__file__))+'/inputs/pt_expression_data_compare.tsv')
    os.remove(os.path.dirname(os.path.abspath(__file__))+'/inputs/pt_expression_data_new.tsv')
