import pandas as pd
import numpy as np
import math

all_vals = []

for i in range(151, 161):

	deleterious = 0

	df = pd.read_csv('../../../results/nrmse_evaluation/paper_data1_rep{}_nmut10/final/rmse_data.tsv'.format(i), header=0, sep='\t')
	df = df[df['Accepted'] == 'yes']
	for i in range(1, len(df)):
		if df.iloc[i]['NRMSE'] > df.iloc[i-1]['NRMSE']:
			deleterious+=1
	all_vals.append(deleterious / len(df))
	#print(deleterious / len(df))
print(np.average(all_vals) / math.sqrt(len(all_vals)))