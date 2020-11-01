import copy
import file_setup
import mutation_analysis
import mutation_choices
import numpy as np
import os
import pandas as pd
import sys
import yaml

sys.path.insert(1, '../manuscript_figures_help/')

import count_successful_sims


def cleanup_genome(output_dir, genome_tracker_old, target_df, mutation_number):
    """
    Removes elements that have low strengths that are statistically insignificant to the overal gene expression pattern produced from the best found architecture.
    Input(s):
    output_dir is a string containing information of the path to the directory in which all the saved files are stored by the program.
    genome_tracker_old is the dataframe containing the most recently accepted genomic data.
    target_df is the user-inputted tsv dataframe containing transcript abundances for each gene.
    Output(s):
    rmse_best refers to the normalized RMSE value of the best genome architecture found.
    Saves the best found genomic architecture to the output directory.
    """
    
    remove_elements = []

    #Open genome architecture file
    genome_tracker_best = copy.deepcopy(genome_tracker_old)
    genome_tracker_saved = copy.deepcopy(genome_tracker_best)
    #Calculate RMSE for best genome and get a range of values
    nrmse_best = mutation_analysis.analyze_mutation(genome_tracker_best, output_dir, target_df, mutation_number, 1, True)
    print("nrmse_best:", nrmse_best[0])

    #Iterate through each element and remove them and compare RMSE values of new architecture to the best architecture previously found
    for gene in range(genome_tracker_old['num_genes']+1):
        promoter = 'promoter_{}'.format(gene)
        terminator = 'terminator_{}'.format(gene)
        rnase = 'rnase_{}'.format(gene)
        #If the region selected is not before the first gene or after the last gene, analyze the significance of promoters and terminators to the genome
        if gene != 0 and gene != genome_tracker_old['num_genes']:
        	#Determine importance of promoters to the genome's fitness
            if genome_tracker_best[promoter]['start'] > 0:
                genome_tracker_saved = mutation_choices.remove_element(genome_tracker_saved, output_dir, genome_tracker_old['num_genes'], promoter)
                #Get RMSE values range to compare to the best found architecture
                nrmse_comp = mutation_analysis.analyze_mutation(genome_tracker_saved, output_dir, target_df, mutation_number, 1, True)
                #If the RMSE values are insignificantly different from one another, remove element
                if (abs(nrmse_comp[0] - nrmse_best[0]) <= 0.01):
                    remove_elements.append((promoter, nrmse_comp[0]))
                elif nrmse_comp[0] < nrmse_best[0]:
                    remove_elements.append((promoter, nrmse_comp[0]))
                #Set the genome architecture file back to original state
                genome_tracker_saved = copy.deepcopy(genome_tracker_best)
            #Determine importance of terminators to the genome's fitness
            if genome_tracker_best[terminator]['start'] > 0:
                genome_tracker_saved = mutation_choices.remove_element(genome_tracker_saved, output_dir, genome_tracker_old['num_genes'], terminator)
                #Get RMSE values range to compare to the best found architecture
                nrmse_comp = mutation_analysis.analyze_mutation(genome_tracker_saved, output_dir, target_df, mutation_number, 1, True)
                #If the RMSE values are insignificantly different from one another, remove element
                if (abs(nrmse_comp[0] - nrmse_best[0]) <= 0.01):
                    remove_elements.append((terminator, nrmse_comp[0]))
                elif nrmse_comp[0] < nrmse_best[0]:
                    remove_elements.append((terminator, nrmse_comp[0]))
                #Set the genome architecture file back to original state
                genome_tracker_saved = copy.deepcopy(genome_tracker_best)
        #If the region selected is not after the last gene, analyze the significance of RNAses to the genome
        if gene != genome_tracker_old['num_genes']:
        	#Determine importance of RNAses to the genome's fitness
            if genome_tracker_best[rnase]['start'] > 0:
                genome_tracker_saved = mutation_choices.remove_element(genome_tracker_saved, output_dir, genome_tracker_old['num_genes'], rnase)
                #Get RMSE values range to compare to the best found architecture
                nrmse_comp = mutation_analysis.analyze_mutation(genome_tracker_saved, output_dir, target_df, mutation_number, 1, True)
                #If the RMSE values are insignificantly different from one another, remove element
                if (abs(nrmse_comp[0] - nrmse_best[0]) <= 0.01):
                    remove_elements.append((rnase, nrmse_comp[0]))
                elif nrmse_comp[0] < nrmse_best[0]:
                    remove_elements.append((rnase, nrmse_comp[0]))
                #Set the genome architecture file back to original state
                genome_tracker_saved = copy.deepcopy(genome_tracker_best)



    #Remove all elements that did not significantly alter the gene expression pattern produced
    remove_elements.sort(key=lambda x: x[1])
    for element in remove_elements:
        genome_tracker_saved = mutation_choices.remove_element(genome_tracker_saved, output_dir, genome_tracker_old['num_genes'], element[0])
        nrmse_comp = mutation_analysis.analyze_mutation(genome_tracker_saved, output_dir, target_df, mutation_number, 1, True)
        if (abs(nrmse_comp[0] - nrmse_best[0]) <= 0.01):
            print('removed')
            print("nrmse_comp", nrmse_comp[0])
            genome_tracker_best = copy.deepcopy(genome_tracker_saved)
        elif nrmse_comp[0] < nrmse_best[0]:
            print('removed2')
            print("nrmse_comp", nrmse_comp[0])
            genome_tracker_best = copy.deepcopy(genome_tracker_saved)
        else:
            break

    #Saves best genome architecture found
    with open(output_dir+'final/gene_clean.yml', 'w') as save_yaml:
        yaml.dump(genome_tracker_best, save_yaml)
    rmse_fin = mutation_analysis.analyze_mutation(genome_tracker_best, output_dir, target_df, mutation_number, 1)
    save_df = pd.read_csv(output_dir+"expression_pattern.tsv", header=0, sep='\t')
    save_df.to_csv(output_dir+"final/expression_pattern_clean.tsv", sep='\t', index=False)
    os.remove(output_dir+"/expression_pattern.tsv")

    return rmse_fin

def main():

    #User inputted data to determine successful simulations
    target_file_name = input("Please input the name of the target file with the extension: ")
    target_file = pd.read_csv("../../../data/targets/{}".format(target_file_name), header=0, sep='\t')
    target_df = file_setup.rearrange_file(target_file, target_file.iloc[-1]['time'], 3)
    num_folders = int(input("Please input number of directories to access: "))

    successes = count_successful_sims.calc_success(target_df, target_file_name, num_folders, True)

    for i in (successes):
        print(i)
        output_dir = '../../../results/{}/rep{}/'.format(target_file_name.strip('.tsv'), i)
        with open(output_dir+'final/gene_best.yml', 'r') as gene_elements:
            genome_tracker_old = yaml.safe_load(gene_elements)
        mutation_number = 10
        print("nrmse_final:", cleanup_genome(output_dir, genome_tracker_old, target_df, mutation_number))

if __name__ == '__main__':
    main()