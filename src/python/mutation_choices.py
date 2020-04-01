import numpy as np
import random
import pandas as pd
import yaml

starting_promoter_strength = 1e10
max_promoter_strength = 3e15
min_promoter_strength = 1e9
promoter_offset = 9
promoter_min_space = 18
max_rnase_strength = 10e17
min_rnase_strength = 1e-2
rnase_offset = 9
rnase_min_space = 12
max_terminator_strength = 1.0
min_terminator_strength = 0.0
terminator_starting_strength = 0.85
terminator_offset = 1
terminator_min_space = 3


def add_element(genome_tracker_new, output_dir, num_genes, deg_rate):

    elements_list = []
    genome_elements = []
    spaces_dict = {}
    possible_spaces = {}
    is_added = False
    #Possible elements that can be added to the genome
    for num_prom in range(1, num_genes):
        elements_list.append('promoter{}'.format(num_prom))
    for num_term in range(1, num_genes+1):
        elements_list.append('terminator{}'.format(num_term))
    for num_rnase in range(num_genes):
        elements_list.append('rnase{}'.format(num_rnase))

    element_choice = random.choice(elements_list)
    region_choice = 'region{}'.format(element_choice[-1:])
    genome_tracker_saved = genome_tracker_new
    #If the element choice is already present on genome then remove it from the genome and shrink the genome
    if genome_tracker_new[element_choice]['start'] > 0:
        genome_tracker_new = remove_element(genome_tracker_new, output_dir, num_genes, deg_rate, element_choice)
    #If promoter is chosen then the required amount of space to place a promoter is set
    if 'promoter' in element_choice:
        required_spaces = promoter_min_space
    #If terminator is chosen then the required amount of space to place a promoter is set
    elif 'terminator' in element_choice:
        required_spaces = terminator_min_space
    #If rnase is chosen then the required amount of space to place a promoter is set
    else:
        required_spaces = rnase_min_space
    #Appends elements that are already on the genome on to a list
    promoter = 'promoter{}'.format(int(region_choice[-1:]))
    terminator = 'terminator{}'.format(int(region_choice[-1:]))
    rnase = 'rnase{}'.format(int(region_choice[-1:]))
    gene = 'gene{}'.format(int(region_choice[-1:]))
    gene_offsetted = 'gene{}'.format(int(region_choice[-1:])+1)
    if int(region_choice[-1:]) != num_genes:
        if genome_tracker_new[promoter]['start'] > 0:
            if int(region_choice[-1:]) != 0:
                genome_elements.append(genome_tracker_new[promoter]['start'])
            genome_elements.append(genome_tracker_new[promoter]['stop'])
        if genome_tracker_new[rnase]['start'] > 0:
            genome_elements.append(genome_tracker_new[rnase]['start'])
            genome_elements.append(genome_tracker_new[rnase]['stop'])
        genome_elements.append(genome_tracker_new[gene_offsetted]['start'])
    elif int(region_choice[-1:]) == num_genes and genome_tracker_new[terminator]['stop'] != genome_tracker_new['length_of_genome']:
        genome_elements.append(genome_tracker_new['length_of_genome'])
    if region_choice != 'region0':
        genome_elements.append(genome_tracker_new[gene]['stop'])
        if genome_tracker_new[terminator]['start'] > 0 and genome_tracker_new[terminator]['start'] != genome_tracker_new[gene]['stop']:
            genome_elements.append(genome_tracker_new[terminator]['start'])
            if int(region_choice[-1:]) != num_genes or genome_tracker_new[terminator]['stop'] != genome_tracker_new['length_of_genome']:
                genome_elements.append(genome_tracker_new[terminator]['stop'])
    #Determines areas within the intergenic regions that are available
    space_index = 0
    highest_position = max(genome_elements)
    while True:
        starting_position = min(genome_elements)
        genome_elements.remove(starting_position)
        ending_position = min(genome_elements)
        genome_elements.remove(ending_position)
        spaces_dict.update({'space{}'.format(space_index): dict(start=starting_position, stop=ending_position)})
        if ending_position == highest_position or genome_elements == []:
            break
        space_index+=1
    space_index = 0
    #Determines if there is a sufficient amount of space needed for the chosen element
    for item in spaces_dict:
        if spaces_dict[item]['stop'] - spaces_dict[item]['start'] >= required_spaces:
            possible_spaces.update({'space{}'.format(space_index): dict(start=spaces_dict[item]['start'], stop=spaces_dict[item]['stop'])})
            space_index+=1
    #Determines the position of the chosen element
    if possible_spaces != {}:
        key = random.choice(list(possible_spaces.keys()))
        starting_position = possible_spaces[key]['start']
        ending_position = possible_spaces[key]['stop']
        if 'promoter' in element_choice:
            #If promoter is chosen then the binding strength, and starting and endging positions are determined and the genome length is changed appropriately
            genome_shift = 10
            if ending_position - starting_position >= promoter_min_space:
                prom_start = genome_tracker_new[element_choice]['start'] = random.randint(starting_position+1, ending_position-17)
                prom_stop = genome_tracker_new[element_choice]['stop'] = prom_start + promoter_offset
                genome_tracker_new[element_choice]['previous_strength'] = genome_tracker_new[element_choice]['current_strength']
                genome_tracker_new[element_choice]['current_strength'] = starting_promoter_strength
                genome_tracker_new = expand_genome(genome_tracker_new, num_genes, int(region_choice[-1:]), genome_shift, element_choice)
            else:
                #If the promoter is not added then revert back to previous accepted genome before adjustments were made
                return genome_tracker_saved
        elif 'terminator' in element_choice:
            #If terminator is chosen then the binding strength, and starting and endging positions are determined and the genome length is changed appropriately
            genome_shift = 2
            if ending_position - starting_position >= terminator_min_space and int(region_choice[-1:]) != num_genes:
                term_start = genome_tracker_new[element_choice]['start'] = random.randint(starting_position+1, ending_position-1)
            elif ending_position - starting_position >= terminator_min_space-1 and int(region_choice[-1:]) == num_genes:
                term_start = genome_tracker_new[element_choice]['start'] = random.randint(starting_position+1, ending_position)
            else:
                #If the terminator is not added then revert back to previous accepted genome before adjustments were made
                return genome_tracker_saved
            term_stop = genome_tracker_new[element_choice]['stop'] = term_start + terminator_offset
            genome_tracker_new[element_choice]['previous_strength'] = genome_tracker_new[element_choice]['current_strength']
            genome_tracker_new[element_choice]['current_strength'] = terminator_starting_strength
            genome_tracker_new = expand_genome(genome_tracker_new, num_genes, int(region_choice[-1:]), genome_shift, element_choice)
        else:
            #If RNase is chosen then the binding strength, and starting and endging positions are determined and the genome length is changed appropriately
            genome_shift = 10
            if ending_position - starting_position >= rnase_min_space:
                rnase_start = genome_tracker_new[element_choice]['start'] = random.randint(starting_position+1, ending_position-8)
                rnase_stop = genome_tracker_new[element_choice]['stop'] = rnase_start + rnase_offset
                if deg_rate == 'yes':
                    genome_tracker_new[element_choice]['previous_strength'] = genome_tracker_new[element_choice]['current_strength']
                    genome_tracker_new[element_choice]['current_strength'] = min_rnase_strength
                genome_tracker_new = expand_genome(genome_tracker_new, num_genes, int(region_choice[-1:]), genome_shift, element_choice)
            else:
                #If the rnase is not added then revert back to previous accepted genome before adjustments were made
                return genome_tracker_saved
    else:
        #If there are no possible spaces to add an element then revert back to previous accepted genome before adjustments were made
        return genome_tracker_saved

    with open(output_dir+'config.yml', 'w') as outfile:
        yaml.dump(genome_tracker_new, outfile, default_flow_style=False)
    return genome_tracker_new

def remove_element(genome_tracker_new, output_dir, num_genes, deg_rate, specific_element=None):

    if specific_element == None:
        possibilities = []
        #Determining possible elements on the genome that can be removed
        for num_prom in range(1, num_genes):
            if genome_tracker_new['promoter{}'.format(num_prom)]['start'] > 0:
                possibilities.append('promoter{}'.format(num_prom))
        for num_term in range(1, num_genes+1):
            if genome_tracker_new['terminator{}'.format(num_term)]['start'] > 0:
                possibilities.append('terminator{}'.format(num_term))
        for num_rnase in range(num_genes):
            if genome_tracker_new['rnase{}'.format(num_rnase)]['start'] > 0:
                possibilities.append('rnase{}'.format(num_rnase))
        #Removes the selected element from the genome
        if possibilities != []:
            element_choice = random.choice(possibilities)
            region_choice = 'region{}'.format(element_choice[-1:])
            if 'promoter' in element_choice or 'rnase' in element_choice:
                genome_shift = 10
            elif 'terminator' in element_choice:
                genome_shift = 2
            genome_tracker_new = shrink_genome(genome_tracker_new, num_genes, int(region_choice[-1:]), genome_shift, element_choice)
            genome_tracker_new[element_choice]['start'] = 0
            genome_tracker_new[element_choice]['stop'] = 0
            #If modifying each individual RNase site, then the strength is set to 0
            if (deg_rate == 'yes') or (deg_rate == 'no' and 'rnase' not in element_choice):
                genome_tracker_new[element_choice]['previous_strength'] = genome_tracker_new[element_choice]['current_strength']
                genome_tracker_new[element_choice]['current_strength'] = 0
        else:
            #If there are no elements to remove, then that means there are no elements to add
            genome_tracker_new = add_element(genome_tracker_new, output_dir, num_genes, deg_rate)
    else:
        region_choice = 'region{}'.format(specific_element[-1:])
        if 'promoter' in specific_element or 'rnase' in specific_element:
            genome_shift = 10
        elif 'terminator' in specific_element:
            genome_shift = 2
        genome_tracker_new = shrink_genome(genome_tracker_new, num_genes, int(region_choice[-1:]), genome_shift, specific_element)
        genome_tracker_new[specific_element]['start'] = 0
        genome_tracker_new[specific_element]['stop'] = 0
        #If modifying each individual RNase site, then the strength is set to 0
        if (deg_rate == 'yes') or (deg_rate == 'no' and 'rnase' not in specific_element):
            genome_tracker_new[specific_element]['previous_strength'] = genome_tracker_new[specific_element]['current_strength']
            genome_tracker_new[specific_element]['current_strength'] = 0

    return genome_tracker_new

def modify_element(genome_tracker_new, output_dir, num_genes, deg_rate):

    possibilities = []
    #Determining which elements are present on the genome so that they can potentially be modified
    for num_prom in range(num_genes):
        if genome_tracker_new['promoter{}'.format(num_prom)]['start'] > 0:
            possibilities.append('promoter{}'.format(num_prom))
    for num_term in range(1, num_genes+1):
        if genome_tracker_new['terminator{}'.format(num_term)]['start'] > 0:
            possibilities.append('terminator{}'.format(num_term))
    if deg_rate == 'yes':
        for num_rnase in range(num_genes):
            if genome_tracker_new['rnase{}'.format(num_rnase)]['start'] > 0:
                possibilities.append('rnase{}'.format(num_rnase))
    #Strengths of current elements on genome are modified
    element_choice = random.choice(possibilities)
    if 'promoter' in element_choice:
        genome_tracker_new[element_choice]['previous_strength'] = genome_tracker_new[element_choice]['current_strength']
        genome_tracker_new[element_choice]['current_strength'] = genome_tracker_new[element_choice]['current_strength'] * np.random.normal(1, 0.1)
        if genome_tracker_new[element_choice]['current_strength'] < min_promoter_strength or genome_tracker_new[element_choice]['current_strength'] > max_promoter_strength:
            genome_tracker_new[element_choice]['current_strength'] = starting_promoter_strength
    elif 'terminator' in element_choice:
        genome_tracker_new[element_choice]['previous_strength'] = genome_tracker_new[element_choice]['current_strength']
        genome_tracker_new[element_choice]['current_strength'] = genome_tracker_new[element_choice]['current_strength'] * np.random.normal(1, 0.1)
        if genome_tracker_new[element_choice]['current_strength'] <= min_terminator_strength or genome_tracker_new[element_choice]['current_strength'] >= max_terminator_strength:
            genome_tracker_new[element_choice]['current_strength'] = terminator_starting_strength
    else:
        genome_tracker_new[element_choice]['previous_strength'] = genome_tracker_new[element_choice]['current_strength']
        genome_tracker_new[element_choice]['current_strength'] = genome_tracker_new[element_choice]['current_strength'] * np.random.normal(1, 0.1)
        if genome_tracker_new[element_choice]['current_strength'] <= min_rnase_strength or genome_tracker_new[element_choice]['current_strength'] >= max_rnase_strength:
            genome_tracker_new[element_choice]['current_strength'] = min_rnase_strength

    return genome_tracker_new

def expand_genome(genome_tracker_new, num_genes, region_choice, genome_shift, element_choice):

    #Increases the genome size if an element is added
    for beg_point in range(region_choice, num_genes+1):

        region = 'region{}'.format(beg_point)
        terminator = 'terminator{}'.format(beg_point)
        promoter = 'promoter{}'.format(beg_point)
        rnase = 'rnase{}'.format(beg_point)
        gene = 'gene{}'.format(beg_point)

        if beg_point != region_choice:
            genome_tracker_new[region]['start']+=genome_shift
        genome_tracker_new[region]['stop']+=genome_shift
        if beg_point != 0 and int(region[-1]) != num_genes:
             if element_choice != promoter:
                if genome_tracker_new[promoter]['start'] >= genome_tracker_new[element_choice]['start'] and genome_tracker_new[promoter]['start'] > 0:
                    genome_tracker_new[promoter]['start']+=genome_shift
                    genome_tracker_new[promoter]['stop']+=genome_shift
        if beg_point != 0:
            if element_choice != terminator:
                if genome_tracker_new[terminator]['start'] >= genome_tracker_new[element_choice]['start'] and genome_tracker_new[terminator]['start'] > 0:
                   genome_tracker_new[terminator]['start']+=genome_shift
                   genome_tracker_new[terminator]['stop']+=genome_shift
        if int(region[-1]) != num_genes:
            if element_choice != rnase:
                if genome_tracker_new[rnase]['start'] >= genome_tracker_new[element_choice]['start'] and genome_tracker_new[rnase]['start'] > 0:
                    genome_tracker_new[rnase]['start']+=genome_shift
                    genome_tracker_new[rnase]['stop']+=genome_shift
        if region != 'region0':
            genome_tracker_new[gene]['start']+=genome_shift
            genome_tracker_new[gene]['stop']+= genome_shift
    genome_tracker_new['length_of_genome']+=genome_shift

    return genome_tracker_new

def shrink_genome(genome_tracker_new, num_genes, region_choice, genome_shift, element_choice):

    #Decreases the genome size if an element is added
    for beg_point in range(region_choice, num_genes+1):

        region = 'region{}'.format(beg_point)
        terminator = 'terminator{}'.format(beg_point)
        promoter = 'promoter{}'.format(beg_point)
        rnase = 'rnase{}'.format(beg_point)
        gene = 'gene{}'.format(beg_point)

        if beg_point != region_choice:
            genome_tracker_new[region]['start']-=genome_shift
        genome_tracker_new[region]['stop']-=genome_shift
        if beg_point != 0 and int(region[-1]) != num_genes:
            if element_choice != promoter:
                if genome_tracker_new[promoter]['start'] >= genome_tracker_new[element_choice]['start'] and genome_tracker_new[promoter]['start'] > 0:
                    genome_tracker_new[promoter]['start']-=genome_shift
                    genome_tracker_new[promoter]['stop']-=genome_shift
        if beg_point != 0:
            if element_choice != terminator:
                if genome_tracker_new[terminator]['start'] >= genome_tracker_new[element_choice]['start'] and genome_tracker_new[terminator]['start'] > 0:
                    genome_tracker_new[terminator]['start']-=genome_shift
                    genome_tracker_new[terminator]['stop']-=genome_shift
        if int(region[-1]) != num_genes:
            if element_choice != rnase:
                if genome_tracker_new[rnase]['start'] >= genome_tracker_new[element_choice]['start'] and genome_tracker_new[rnase]['start'] > 0:
                    genome_tracker_new[rnase]['start']-=genome_shift
                    genome_tracker_new[rnase]['stop']-=genome_shift
        if region != 'region0':
            genome_tracker_new[gene]['start']-=genome_shift
            genome_tracker_new[gene]['stop']-= genome_shift
    genome_tracker_new['length_of_genome']-=genome_shift

    return genome_tracker_new
