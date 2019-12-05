import numpy as np
import random
import pandas as pd
import yaml

starting_promoter_strength = 10e8
max_promoter_strength = 3e10
min_promoter_strength = 10e5
promoter_offset = 9
max_rnase_strength = 10e17
min_rnase_strength = 1e-2
rnase_offset = 10
max_terminator_strength = 1.0
min_terminator_strength = 0.0
terminator_starting_strength = 0.85
terminator_offset = 1


def modify_promoter(genome_tracker_new, output_dir):
    """
    Promoters are either added with a starting polymerase strength, removed from the genome all together, or has its polymerase strength modified.
    """

    promoter_modification = ['add', 'remove', 'modify']
    chosen_prom_modification = random.choice(promoter_modification)
    promoter_slots = ['A', 'B']

    if chosen_prom_modification == 'add':
        promoter_possibilities = ["promoter2", "promoter3"]
        chosen_promoter = random.choice(promoter_possibilities)
        if chosen_promoter == 'promoter2':
            regionA = 'region2a'
            regionB = 'region2b'
            rnase = 'rnase1'
            terminator = 'terminator1'
            region_endA = genome_tracker_new['region2a']['start']
            region_endB = genome_tracker_new['region2b']['start']
        elif chosen_promoter == 'promoter3':
            regionA = 'region3a'
            regionB = 'region3b'
            rnase = 'rnase2'
            terminator = 'terminator2'
            region_endA = genome_tracker_new['region3a']['start']
            region_endB = genome_tracker_new['region3b']['start']

        #Adding in a promoter between genes 1 and 2 or genes 2 and 3
        items = [rnase, terminator]
        for item in items:
            if (genome_tracker_new[item]['start'] >= genome_tracker_new[regionA]['start']) and (genome_tracker_new[item]['start'] <= genome_tracker_new[regionA]['stop']):
                promoter_slots.remove('A')
            if (genome_tracker_new[item]['start'] >= genome_tracker_new[regionB]['start']) and (genome_tracker_new[item]['start'] <= genome_tracker_new[regionB]['stop']):
                promoter_slots.remove('B')

        if promoter_slots == []:
            genome_tracker_new[chosen_promoter]['start'] = 0
            genome_tracker_new[chosen_promoter]['stop'] = 0
        else:
            available_slot = random.choice(promoter_slots)
            if available_slot == 'A':
                region = regionA
                region_end = region_endA
            elif available_slot == 'B':
                region = regionB
                region_end = region_endB
            prom_start = genome_tracker_new[chosen_promoter]['start'] = random.randint(genome_tracker_new[region]['start'], region_end)
            prom_stop = genome_tracker_new[chosen_promoter]['stop'] = prom_start + promoter_offset

        genome_tracker_new[chosen_promoter]['previous_strength'] = genome_tracker_new[chosen_promoter]['current_strength']
        genome_tracker_new[chosen_promoter]['current_strength'] = starting_promoter_strength

    #Altering the polymerase binding strength
    if chosen_prom_modification == 'modify':
        promoter_possibilities = ['promoter1']
        for promoter in ['promoter2', 'promoter3']:
            if genome_tracker_new[promoter]['start'] > 0:
                promoter_possibilities.append(promoter)
        chosen_promoter = random.choice(promoter_possibilities)
        genome_tracker_new[chosen_promoter]['previous_strength'] = genome_tracker_new[chosen_promoter]['current_strength']
        genome_tracker_new[chosen_promoter]['current_strength'] = genome_tracker_new[chosen_promoter]['current_strength'] * np.random.normal(1, 0.1)
        while genome_tracker_new[chosen_promoter]['current_strength'] < min_promoter_strength or genome_tracker_new[chosen_promoter]['current_strength'] > max_promoter_strength:
            genome_tracker_new[chosen_promoter]['current_strength'] = genome_tracker_new[chosen_promoter]['previous_strength'] * np.random.normal(1, 0.1)

    if chosen_prom_modification == 'remove':
        promoter_possibilities = ['promoter2', 'promoter3']
        chosen_promoter = random.choice(promoter_possibilities)
        #Removing promoters
        genome_tracker_new[chosen_promoter]['start'] = 0
        genome_tracker_new[chosen_promoter]['stop'] = 0
        #Resetting strengths
        genome_tracker_new[chosen_promoter]['previous_strength'] = genome_tracker_new[chosen_promoter]['current_strength']
        genome_tracker_new[chosen_promoter]['current_strength'] = 0.0

    with open(output_dir, 'w') as prom_yaml:
        yaml.dump(genome_tracker_new, prom_yaml, default_flow_style=False)


def modify_rnase(genome_tracker_new, output_dir):
    """
    Rnases are either added with a starting polymerase strength, removed from the genome all together, or has its degredation rate modified.
    """

    rnase_modification = ['remove', 'add', 'modify']
    chosen_rnase_modification = random.choice(rnase_modification)
    rnase_possibilities = ["rnase1", "rnase2", "rnase3"]
    chosen_rnase = random.choice(rnase_possibilities)
    rnase_slots = ['A', 'B', 'C']

    if chosen_rnase_modification == 'add':
        if chosen_rnase == 'rnase2':
            regionA = 'region2a'
            regionB = 'region2b'
            regionC = 'region2c'
            promoter = 'promoter2'
            terminator = 'terminator1'
            region_endA = genome_tracker_new['region2a']['start']
            region_endB = genome_tracker_new['region2b']['start']
            region_endC = genome_tracker_new['region2c']['start']
        elif chosen_rnase == 'rnase3':
            regionA = 'region3a'
            regionB = 'region3b'
            regionC = 'region3c'
            promoter = 'promoter3'
            terminator = 'terminator2'
            region_endA = genome_tracker_new['region3a']['start']
            region_endB = genome_tracker_new['region3b']['start']
            region_endC = genome_tracker_new['region3c']['start']

        if chosen_rnase == "rnase1":
            #Adds rnase after first promoter
            rnase1_start = genome_tracker_new['rnase1']['start'] = random.randint(genome_tracker_new['region1']['start'], 15)
            rnase1_stop = genome_tracker_new['rnase1']['stop'] = rnase1_start + rnase_offset
        else:
            items = [promoter, terminator]
            for item in items:
            #Adds rnase after first gene
                if (genome_tracker_new[item]['start'] >= genome_tracker_new[regionA]['start']) and (genome_tracker_new[item]['start'] <= genome_tracker_new[regionA]['stop']):
                    rnase_slots.remove('A')
                if (genome_tracker_new[item]['start'] >= genome_tracker_new[regionB]['start']) and (genome_tracker_new[item]['start'] <= genome_tracker_new[regionB]['stop']):
                    rnase_slots.remove('B')
            if (genome_tracker_new[terminator]['start'] >= genome_tracker_new[regionC]['start']) and (genome_tracker_new[terminator]['start'] <= genome_tracker_new[regionC]['stop']):
                rnase_slots.remove('C')
            if rnase_slots == []:
                genome_tracker_new[chosen_rnase]['start'] = 0
                genome_tracker_new[chosen_rnase]['stop'] = 0
            else:
                available_slot = random.choice(rnase_slots)
                if available_slot == 'A':
                    region = regionA
                    region_end = region_endA
                if available_slot == 'B':
                    region = regionB
                    region_end = region_endB
                if available_slot == 'C':
                    region = regionC
                    region_end = region_endC
                rnase_start = genome_tracker_new[chosen_rnase]['start'] = random.randint(genome_tracker_new[region]['start'], region_end)
                rnase_stop = genome_tracker_new[chosen_rnase]['stop'] = rnase_start + rnase_offset

        genome_tracker_new[chosen_rnase]['previous_strength'] = genome_tracker_new[chosen_rnase]['current_strength']
        genome_tracker_new[chosen_rnase]['current_strength'] = min_rnase_strength

    if chosen_rnase_modification == 'remove':
        genome_tracker_new[chosen_rnase]['start'] = 0
        genome_tracker_new[chosen_rnase]['stop'] = 0
        genome_tracker_new[chosen_rnase]['previous_strength'] = genome_tracker_new[chosen_rnase]['current_strength']
        genome_tracker_new[chosen_rnase]['current_strength'] = 0

    if chosen_rnase_modification == 'modify':
        rnase_possibilities = []
        if genome_tracker_new['rnase1']['current_strength'] > 0:
            rnase_possibilities.append('rnase1')
        if genome_tracker_new['rnase2']['current_strength'] > 0:
            rnase_possibilities.append('rnase2')
        if genome_tracker_new['rnase3']['current_strength'] > 0:
            rnase_possibilities.append('rnase3')
        if rnase_possibilities != []:
            chosen_rnase = random.choice(rnase_possibilities)
            genome_tracker_new[chosen_rnase]['previous_strength'] = genome_tracker_new[chosen_rnase]['current_strength']
            genome_tracker_new[chosen_rnase]['current_strength'] = genome_tracker_new[chosen_rnase]['current_strength'] * np.random.normal(1, 0.1)
            while genome_tracker_new[chosen_rnase]['current_strength'] <= min_rnase_strength or genome_tracker_new[chosen_rnase]['current_strength'] >= max_rnase_strength:
                genome_tracker_new[chosen_rnase]['current_strength'] = genome_tracker_new[chosen_rnase]['previous_strength'] * np.random.normal(1, 0.1)

    with open(output_dir, 'w') as rnase_yaml:
        yaml.dump(genome_tracker_new, rnase_yaml, default_flow_style=False)


def modify_terminator(genome_tracker_new, output_dir):
    """
    Terminators are either added with a starting polymerase strength, removed from the genome all together, or has its terminator efficiency value modified.
    """

    terminator_modification = ['remove', 'add', 'modify']
    chosen_term_modification = random.choice(terminator_modification)
    terminator_possibilities = ["terminator1", "terminator2", "terminator3"]
    chosen_terminator = random.choice(terminator_possibilities)
    terminator_slots = ['A', 'B', 'C']

    if chosen_term_modification == 'add':
        if chosen_terminator == 'terminator1':
            regionA = 'region2a'
            regionB = 'region2b'
            regionC = 'region2c'
            rnase = 'rnase1'
            promoter = 'promoter2'
            region_endA = genome_tracker_new['region2a']['start']
            region_endB = genome_tracker_new['region2b']['start']
            region_endC = genome_tracker_new['region2c']['start']
        elif chosen_terminator == 'terminator2':
            regionA = 'region3a'
            regionB = 'region3b'
            regionC = 'region3c'
            rnase = 'rnase2'
            promoter = 'promoter3'
            region_endA = genome_tracker_new['region3a']['start']
            region_endB = genome_tracker_new['region3b']['start']
            region_endC = genome_tracker_new['region3c']['start']

        if chosen_terminator == "terminator3":
            #Adds terminator after third gene
            genome_tracker_new[chosen_terminator]['start'] = genome_tracker_new['geneZ']['stop']
            genome_tracker_new[chosen_terminator]['stop'] = genome_tracker_new['geneZ']['stop'] + terminator_offset
        else:
            items = [promoter, rnase]
            for item in items:
                #Adds terminator after first gene
                if (genome_tracker_new[item]['start'] >= genome_tracker_new[regionA]['start']) and (genome_tracker_new[item]['start'] <= genome_tracker_new[regionA]['stop']):
                    terminator_slots.remove('A')
                if (genome_tracker_new[item]['start'] >= genome_tracker_new[regionB]['start']) and (genome_tracker_new[item]['start'] <= genome_tracker_new[regionB]['stop']):
                    terminator_slots.remove('B')
            if (genome_tracker_new[rnase]['start'] >= genome_tracker_new[regionC]['start']) and (genome_tracker_new[rnase]['start'] <= genome_tracker_new[regionC]['stop']):
                terminator_slots.remove('C')
            if terminator_slots == []:
                genome_tracker_new[chosen_terminator]['start'] = 0
                genome_tracker_new[chosen_terminator]['stop'] = 0
            else:
                available_slot = random.choice(terminator_slots)
                if available_slot == 'A':
                    region = regionA
                    region_end = region_endA
                if available_slot == 'B':
                    region = regionB
                    region_end = region_endB
                if available_slot == 'C':
                    region = regionC
                    region_end = region_endC
                term_start = genome_tracker_new[chosen_terminator]['start'] = random.randint(genome_tracker_new[region]['start'], region_end)
                term_stop = genome_tracker_new[chosen_terminator]['stop'] = term_start + terminator_offset

        genome_tracker_new[chosen_terminator]['previous_strength'] = genome_tracker_new[chosen_terminator]['current_strength']
        genome_tracker_new[chosen_terminator]['current_strength'] = terminator_starting_strength

    #Altering the terminator efficiency rate
    if chosen_term_modification == 'modify':
        terminator_possibilities = []
        if genome_tracker_new['terminator1']['start'] > 0:
            terminator_possibilities.append('terminator1')
        if genome_tracker_new['terminator2']['start'] > 0:
            terminator_possibilities.append('terminator2')
        if genome_tracker_new['terminator3']['start'] > 0:
            terminator_possibilities.append('terminator3')
        if terminator_possibilities != []:
            chosen_terminator = random.choice(terminator_possibilities)
            genome_tracker_new[chosen_terminator]['previous_strength'] = genome_tracker_new[chosen_terminator]['current_strength']
            genome_tracker_new[chosen_terminator]['current_strength'] = genome_tracker_new[chosen_terminator]['current_strength'] * np.random.normal(1, 0.1)
            while genome_tracker_new[chosen_terminator]['current_strength'] <= min_terminator_strength or genome_tracker_new[chosen_terminator]['current_strength'] >= max_terminator_strength:
                genome_tracker_new[chosen_terminator]['current_strength'] = genome_tracker_new[chosen_terminator]['previous_strength'] * np.random.normal(1, 0.1)

    if chosen_term_modification == 'remove':
        genome_tracker_new[chosen_terminator]['start'] = 0
        genome_tracker_new[chosen_terminator]['stop'] = 0
        genome_tracker_new[chosen_terminator]['previous_strength'] = genome_tracker_new[chosen_terminator]['current_strength']
        genome_tracker_new[chosen_terminator]['current_strength'] = 0.0

    with open(output_dir, 'w') as term_yaml:
        yaml.dump(genome_tracker_new, term_yaml, default_flow_style=False)
