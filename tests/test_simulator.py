#Common imports
import sys

sys.path.append('../src/python/')

from os import path
from os import remove
import pytest
import pandas as pd
import yaml

#lib imports
from lib.genome_simulator import pt_call

def test_pt_call():

	#Verify the pinetree call works
	max_time = 300
	with open('./inputs/testing.yml', 'r') as gene_parameters:
		genome_tracker = yaml.safe_load(gene_parameters)
	pt_call('./', genome_tracker, max_time)
	assert path.exists("expression_pattern.tsv")
	remove("expression_pattern.tsv")
