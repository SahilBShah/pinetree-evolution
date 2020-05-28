import sys

sys.path.insert(1, '../src/python/lib/')

import filecmp
import os
import pandas as pd
import pytest
from file_setup import rearrange_file
import yaml

def test_rearrange_file():

    file_comp = pd.read_csv("./inputs/test_compare.tsv", header=0, sep='\t')
    file_new = pd.read_csv("./inputs/test.tsv", header=0, sep='\t')
    with open('./inputs/testing.yml', 'r') as gene_parameters:
        genome_tracker = yaml.safe_load(gene_parameters)
    file_new = rearrange_file(file_new, genome_tracker)
    file_new.to_csv("test_new.tsv", sep='\t', index=False)
    assert filecmp.cmp('test_new.tsv', './inputs/test_compare.tsv')
    os.remove('test_new.tsv')
