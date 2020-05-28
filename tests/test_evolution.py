import sys

sys.path.insert(1, '../src/python/lib/')

import copy
import datetime
import filecmp
import pandas as pd
import os
import pytest
import yaml
import evolution_ignore


def test_enumerate_mutation_options():

    with open('./inputs/testing.yml', 'r') as gene_parameters:
        genome_tracker = yaml.safe_load(gene_parameters)
    dynamic_deg_rate = True
    mutation_possibilities = evolution_ignore.enumerate_mutation_options(genome_tracker, dynamic_deg_rate)
    assert mutation_possibilities != {}
    assert 'promoter_2.remove' in mutation_possibilities
    assert mutation_possibilities['promoter_2.remove'] == 'remove'
    assert mutation_possibilities['promoter_2.remove'] != 'add'
    assert 'terminator_1.add' in mutation_possibilities
    assert mutation_possibilities['terminator_1.add'] == 'add'
    assert 'rnase_2.modify0' in mutation_possibilities
    assert mutation_possibilities['rnase_2.modify0'] == 'modify'
