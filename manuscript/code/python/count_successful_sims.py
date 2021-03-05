#Common imports
import sys

sys.path.insert(1, '../')

import pandas as pd
import copy
import yaml
from scipy import stats
import numpy as np
import os

#lib imports
from lib.file_setup import rearrange_file
from lib.root_mean_square_error import calc_accepted_rmse_range
from lib.mutation_choices import remove_element

def calc_success(target_df, target_name, n_folders, successes_out=False):

	success = 0
	#max_rmse = calc_accepted_rmse_range(target_df, 3)
	max_rmse = 0.10
	success_list = []

	print("Successes:")
	for i in range(1, n_folders+1):
		if os.path.isdir("../../../results/{}/rep{}".format(target_name.strip(".tsv"), i)):
			rmse_df = pd.read_csv("../../../results/{}/rep{}/final/rmse_data.tsv".format(target_name.strip(".tsv"), i), header=0, sep='\t')
			#Get index with lowest sum of squared error value
			rmse_df = rmse_df[rmse_df['Accepted'] == 'yes']
			min_rmse_df = rmse_df[rmse_df.NRMSE == rmse_df.NRMSE.min()]
			min_rmse = min_rmse_df.iloc[-1]['NRMSE']
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
	target_df = rearrange_file(target_file, target_file.iloc[-1]['time'], 3)
	num_folders = int(input("Please input number of directories to access: "))

	print()
	print('Percent successful:', str(calc_success(target_df, target_file_name, num_folders)) + '%')

if __name__ == '__main__':
	main()