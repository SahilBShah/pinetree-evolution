#Common imports
from os import path
from os import remove
import pytest
import pandas as pd
import yaml

#lib imports
from ..genome_simulator import pt_call

def test_pt_call():

	#Verify the pinetree call works
	max_time = 300
	with open(path.dirname(path.abspath(__file__))+'/inputs/testing.yml', 'r') as gene_parameters:
		genome_tracker = yaml.safe_load(gene_parameters)
	pt_call(path.dirname(path.abspath(__file__))+'/', genome_tracker, max_time)
	assert path.exists(path.dirname(path.abspath(__file__))+"/expression_pattern.tsv")
	remove(path.dirname(path.abspath(__file__))+"/expression_pattern.tsv")
