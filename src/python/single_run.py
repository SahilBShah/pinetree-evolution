#Common imports

import sys

sys.path.insert(1, './lib/')

import os
import pandas as pd
import pinetree as pt
import yaml

#lib imports
from genome_simulator import pt_call_alt
import initialize_yaml

def main():
    """
    Runs a single simulation using pinetree.
    """

    #Sets output directory
    output_dir = "../../tests/single_run_outputs/"
    #Opens yaml files containing genome coordinates
    starting_genome = output_dir + 'config.yml'
    #Opens yaml file
    with open('../../results/positive_control3/rep5/final/gene_best_alt.yml', 'r') as gene_parameters:
        genome_config = yaml.safe_load(gene_parameters)
    #Create a configuration file from yaml file
    initialize_yaml.create_yaml(starting_genome, genome_config)
    #Open yaml file containing the information for the full genome architecture
    with open(starting_genome, 'r') as gene_elements:
        genome_tracker_new = yaml.safe_load(gene_elements)
    #Simulates genome from yaml file
    pt_call_alt(output_dir, genome_tracker_new, 300)
    print('Simulation successful.')


if __name__ == '__main__':
    main()
