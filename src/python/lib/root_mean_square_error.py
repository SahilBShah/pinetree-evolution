import numpy as np
from statsmodels.tools.eval_measures import rmse

#Minimizes distance between observed values and regression line from given data
def calc_nrmse(df_target, df_new):
    """
    Calculates the normalized root mean square error of the target input dataframe and the simulated output dataframe to determine if the rmse decreased with a new mutation or not.
    Input(s):
    df_target is the user-inputted tsv file containing transcript abundances for each gene.
    df_new is the simulator-generated tsv file containng transcript abundances for each gene.
    Output(s):
    RMSE is a floating point number that refers to the root mean square error calculated.
    """

    #Confirms that the dataframes are the same shape
    assert df_target.shape == df_new.shape
    assert all(df_target.columns == df_new.columns)
    assert all(df_target.index == df_new.index)

    norm_errs = []

    #Performs a normalized RMSE to help determine the fitness of the new genome
    for column in df_target.columns:
        nrmse = rmse(df_target[column], df_new[column])/\
                        np.mean(df_target[column])
        norm_errs.append(nrmse)

    return np.mean(norm_errs)

def calc_accepted_rmse_range(df_target):
    """
    Calculates the highest RMSE value allowed before the program has found the a suitable architecture to reproduce the target data.
    Input(s):
    df_target is the user-inputted tsv file containing transcript abundances for each gene.
    Output(s):
    max_sse is a floating point number that refers to the highest sum of squared error value that is considered a success.
    """

    #Decreases the transcript abundances of the target file by 10%
    df_altered = df_target[['protein1', 'protein2', 'protein3']] * 0.95
    #Calculates the sum of squared error between the reduced values and the target values and returns it
    max_rmse = calc_nrmse(df_target, df_altered)
    return max_rmse
