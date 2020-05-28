import sys

sys.path.insert(1, '../src/python/lib/')

from os import path
from os import remove
import pytest
import pandas as pd
from genome_simulator import pt_call_alt
import yaml

def test_pt_call():

    with open('testing.yml', 'r') as gene_parameters:
        genome_tracker = yaml.safe_load(gene_parameters)
    pt_call_alt('./', genome_tracker)
    assert path.exists("expression_pattern.tsv")
    remove("expression_pattern.tsv")
