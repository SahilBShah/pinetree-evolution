#common imports
import pandas as pd
from scipy.stats import entropy
import sim_success
import sys
import yaml

sys.path.insert(1, '../')

#lib imports
from lib.file_setup import rearrange_file

def get_architecture(genome_tracker):
	"""
	This function returns a string representation of the desired genome architecture.
	Input(s):
	genome_tracker is the inputted yaml dataframe containing genome infromation
	Output(s):
	String representation of genome architecture.
	"""

	genome_elements = []
	arch_string = ''

	#Iterate through each region on the genome
	for n_genes in range(genome_tracker['num_genes'] + 1):
		#Enumarate possible genome elements present
		promoter = 'promoter_{}'.format(n_genes)
		terminator = 'terminator_{}'.format(n_genes)
		rnase = 'rnase_{}'.format(n_genes)
		gene = 'gene_{}'.format(n_genes)

		#If not the last region
		if n_genes != genome_tracker['num_genes']:
			#If promoter is present on genome add to list
			if genome_tracker[promoter]['start'] > 0:
				genome_elements.append((genome_tracker[promoter]['start'], 'p'))
			#If rnase is present on genome add to list
			if genome_tracker[rnase]['start'] > 0:
				genome_elements.append((genome_tracker[rnase]['start'], 'r'))
		#If not the first region (region upstream gene X)
		if n_genes != 0:
			#If terminator is present on genome add to list
			if genome_tracker[terminator]['start'] > 0:
					genome_elements.append((genome_tracker[terminator]['start'], 't'))
			genome_elements.append((genome_tracker[gene]['start'], 'g'))
	#Sort elements based on which occur first on the genome
	genome_elements.sort()

	#Create string representing genome architecture
	for element in genome_elements:
		arch_string+=element[1]

	return arch_string

def calc_string_prob(string_archs):
	"""
	Calculate the probabilities that a genome architecture occurs within the list of architectures.
	Input(s):
	string_archs is the list of strings representing genome architectures.
	Output(s):
	arch_probs is a list of probabilities each architecture occurs.
	"""

	arch_probs = []
	used_archs = []
	string_count = []

	#Iterate through list of genetic architectures
	for n_archs in range(len(string_archs)):
		#If architecture has not already been analyzed
		if string_archs[n_archs] not in used_archs:
			#Get probability architecture occurs in the list of architectures
			arch_probs.append(string_archs.count(string_archs[n_archs]) / len(string_archs))
			#Once architecture is analyzed, add to list
			used_archs.append(string_archs[n_archs])
			string_count.append(string_archs.count(string_archs[n_archs]))

	print(string_count)
	print(used_archs)

	return(arch_probs)


def main():

	successful_archs = []

	#User inputted data to determine successful simulations
	target_file_name = input("Please input the name of the target file with the extension: ")
	target_file = pd.read_csv("../../../data/targets/{}".format(target_file_name), header=0, sep='\t')
	target_df = rearrange_file(target_file, target_file.iloc[-1]['time'], 3)
	num_folders = int(input("Please input number of directories to access: "))

	successes = sim_success.calc_success(target_df, target_file_name, num_folders, True)

	#Iterate through directories containing yaml files
	for i in (successes):
		output_dir = '../../../results/{}/rep{}/'.format(target_file_name.strip('.tsv'), i)
		with open(output_dir+'final/gene_clean.yml', 'r') as gene_elements:
			genome_tracker = yaml.safe_load(gene_elements)
		#Get string representation of architectures and add to list
		arch_string = get_architecture(genome_tracker)
		successful_archs.append(arch_string)
	#Calculate probabilites that each architecture occurs within the ones found
	probabilities = calc_string_prob(successful_archs)

	#Print out entropy value of the genome architectures
	print(entropy(probabilities, base=2))


if __name__ == '__main__':
	main()

