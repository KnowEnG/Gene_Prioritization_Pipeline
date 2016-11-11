"""
@author: The KnowEnG dev team
"""
import os
import numpy as np
import pandas as pd

import knpackage.toolbox as kn

from gene_prioritization_toolbox import run_net_correlation

def correlation(run_parameters):
    """ gene prioritization """
    from gene_prioritization_toolbox import run_correlation
    run_correlation(run_parameters)

def net_correlation(run_parameters):
    """ gene prioritization """
    from gene_prioritization_toolbox import run_net_correlation
    run_net_correlation(run_parameters)

def bootstrap_correlation(run_parameters):
    """ gene prioritization """
    from gene_prioritization_toolbox import run_bootstrap_correlation
    run_bootstrap_correlation(run_parameters)

def bootstrap_net_correlation(run_parameters):
    """ gene prioritization """
    from gene_prioritization_toolbox import run_bootstrap_net_correlation
    run_bootstrap_net_correlation(run_parameters)

SELECT = {
    "correlation": correlation,
    "net_correlation": net_correlation,
    "bootstrap_correlation": bootstrap_correlation,
    "bootstrap_net_correlation": bootstrap_net_correlation}


def run_multi_drug_correlation(run_parameters):
    """  write correlation results for multi-drug (with NAs) input  """
    drug_response_df = kn.get_spreadsheet_df(run_parameters["drug_response_full_path"])
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])
    genes_list, drugs_list, consolodated_df = get_consolodated_dataframe(spreadsheet_df, drug_response_df)

    run_parameters['tmp_data'] = kn.create_dir(run_parameters['results_directory'], 'run_multi_drug_correlation_tmp', kn.get_timestamp())
    run_parameters['method'] = run_parameters['multi_drug_method']
    for drug_name in drugs_list:
        pass_parameters = write_spreadsheet_drug_dataframes(consolodated_df, genes_list, drug_name, run_parameters)
        SELECT[run_parameters["method"]](pass_parameters)

    kn.remove_dir(run_parameters['tmp_data'])
    return

def get_consolodated_dataframe(gene_samples_df, drug_samples_df):
    """  assemble the features X samples and data X samples dataframes into one  """
    genes_list = gene_samples_df.index.values.tolist()
    drugs_list = drug_samples_df.index.values.tolist()
    consolodated_df = pd.DataFrame(
        gene_samples_df.as_matrix(), index=genes_list, columns=gene_samples_df.columns.values)
    consolodated_df = consolodated_df.append(drug_samples_df)

    return genes_list, drugs_list, consolodated_df

def get_data_for_drug(consolodated_dataframe, genes_list, drug_name, run_parameters):
    """     """
    names_selector = genes_list.copy()
    names_selector.insert(0, drug_name)
    tmp_df = consolodated_dataframe.loc[names_selector]
    if run_parameters['drop_method'] is 'drop_NA':
        tmp_df = tmp_df.dropna(axis=1)
    else:
        tmp_df = tmp_df.fillna(0)

    drug_array = np.array([tmp_df.as_matrix()[0, :]])
    tmp_df = tmp_df.drop(drug_name)
    spreadsheet_matrix = tmp_df.as_matrix()
    sample_names = tmp_df.columns.values

    return drug_array, spreadsheet_matrix, sample_names

def write_spreadsheet_drug_dataframes(consolodated_df, genes_list, drug_name, run_parameters):
    """   """

    file_name_prefix = drug_name + '_' + run_parameters['drop_method']

    drug_array, spreadsheet_matrix, sample_names = get_data_for_drug(consolodated_df, genes_list, drug_name, run_parameters)

    spreadsheet_for_drug_df = pd.DataFrame(spreadsheet_matrix, index=genes_list, columns=sample_names)
    spreadsheet_file_name = file_name_prefix + '_gene_expression.tsv'
    spreadsheet_file_name = os.path.join(run_parameters['tmp_data'], spreadsheet_file_name)
    spreadsheet_for_drug_df.to_csv(spreadsheet_file_name, sep='\t')
    run_parameters['spreadsheet_name_full_path'] = spreadsheet_file_name

    drug_for_spreadsheet_df = pd.DataFrame(drug_array, index={drug_name}, columns=sample_names)
    drug_file_name = file_name_prefix + '_single_drug_response.tsv'
    drug_file_name = os.path.join(run_parameters['tmp_data'], drug_file_name)
    drug_for_spreadsheet_df.to_csv(drug_file_name, sep='\t')
    run_parameters["drug_response_full_path"] = drug_file_name


    return run_parameters