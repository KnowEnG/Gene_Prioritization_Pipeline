"""
Created on Fri Sep 23 16:39:35 2016
@author: The Gene Prioritization dev team
"""

import numpy as np
import pandas as pd
from scipy.stats import pearsonr as pcc
from sklearn.preprocessing import normalize
import knpackage.toolbox as kn

def perform_pearson_correlation(spreadsheet, drug_response):
    """ Find pearson correlation coefficient(PCC) for each gene expression (spreadsheet row)
    with the drug response.
    Args:
        spreadsheet: genes x samples gene expression data
        drug_response: one drug response for each sample
    """
    pc_array = np.zeros(spreadsheet.shape[0])
    for row in range(0, spreadsheet.shape[0]):
        pcc_value = pcc(spreadsheet[row, :], drug_response)[0]
        pc_array[row] = pcc_value

    return pc_array

def run_gene_correlation(run_parameters):
    ''' wrapper: call sequence to perform gene prioritization
    Args:
        spreadsheet_df_full_path: spreadsheet_df path and file name.
        drug_response_full_path: drug_response path and file name.
    Returns:
        result_df: result dataframe of gene prioritization. Values are pearson
        correlation coefficient values in a descending order.
    '''

    #run_parameters['samples_file_name'] = run_parameters["spreadsheet_name_full_path"]
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])
    #run_parameters['samples_file_name'] = run_parameters["drug_response_full_path"]
    drug_response = kn.get_spreadsheet_df(run_parameters["drug_response_full_path"])
    pc_array = perform_pearson_correlation(spreadsheet_df.values, drug_response.values[0])
    result_df = pd.DataFrame(pc_array, index=spreadsheet_df.index.values,
                             columns=['PCC']).abs().sort_values("PCC", ascending=0)
    file_name = kn.create_timestamped_filename("pcc_result")
    file_name = file_name + '.df'
    result_df.to_csv(file_name, header=True, index=True, sep='\t')

    return result_df

def run_net_gene_correlation(run_parameters):
    print('run_net_gene_correlation called')

    return


def run_bootstrap_correlation(run_parameters):
    print('run_bootstrap_correlation called')

    return


def run_bootstrap_net_correlation(run_parameters):
    print('run_bootstrap_net_correlation called')

    return