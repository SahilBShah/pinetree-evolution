#Common imports
import sys

sys.path.insert(1, '../lib/')

import pandas as pd
import copy
import yaml
from scipy import stats
import numpy as np
import os

#lib imports
import file_setup
import root_mean_square_error
import mutation_analysis
from mutation_choices import remove_element

def calc_success(target_df, target_name, n_folders, successes_out=False):

	success = 0
	#max_rmse = root_mean_square_error.calc_accepted_rmse_range(target_df, 3)
	max_rmse = 0.10
	success_list = []

	print("Successes:")
	for i in range(1, n_folders+1):
		if os.path.isdir("../../../results/{}/rep{}".format(target_name.strip(".tsv"), i)):
			#final_file = pd.read_csv("../../../results/2020_7_9/{}_rep{}_nmut10/final/expression_pattern_best.tsv".format(target_name.strip(".tsv"), i), header=0, sep='\t')
			#final_df = file_setup.rearrange_file(final_file, final_file.iloc[-1]['time'], 3)
			rmse_df = pd.read_csv("../../../results/{}/rep{}/final/rmse_data.tsv".format(target_name.strip(".tsv"), i), header=0, sep='\t')
			#df_rmse = root_mean_square_error.calc_nrmse(target_df, final_df)
			#Get index with lowest sum of squared error value
			rmse_df = rmse_df[rmse_df['Accepted'] == 'yes']
			min_rmse_df = rmse_df[rmse_df.NRMSE == rmse_df.NRMSE.min()]
			min_rmse = min_rmse_df.iloc[-1]['NRMSE']
			#print(min_rmse)
			if min_rmse <= max_rmse:
				success+=1
				success_list.append(i)

	if successes_out:
		return(success_list)
	else:
		print(success_list)
		return ((success / n_folders) * 100)


def main():

	target_file_name = input("Please input the name of the target file with the extension: ")
	target_file = pd.read_csv("../../../data/targets/{}".format(target_file_name), header=0, sep='\t')
	target_df = file_setup.rearrange_file(target_file, target_file.iloc[-1]['time'], 3)
	num_folders = int(input("Please input number of directories to access: "))

	print()
	print('Percent successful:', str(calc_success(target_df, target_file_name, num_folders)) + '%')

if __name__ == '__main__':
	main()