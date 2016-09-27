"""
Created on Fri Sep 23 16:39:35 2016
@author: The Gene Prioritization dev team
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr as pcc
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
        run_parameters: dict object with keys:
                run_parameters["spreadsheet_name_full_path"]
                run_parameters["drug_response_full_path"]
    Returns:
        result_df:  result - dataframe of gene prioritization. Values are pearson
                    correlation coefficient values in descending order.
    '''
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])
    drug_response = kn.get_spreadsheet_df(run_parameters["drug_response_full_path"])
    pc_array = perform_pearson_correlation(spreadsheet_df.values, drug_response.values[0])
    result_df = pd.DataFrame(pc_array, index=spreadsheet_df.index.values,
                             columns=['PCC']).abs().sort_values("PCC", ascending=0)

    target_file_base_name = "gene_drug_correlation"
    target_file_base_name = os.path.join(run_parameters["results_directory"], target_file_base_name)
    file_name = kn.create_timestamped_filename(target_file_base_name)
    file_name = file_name + '.df'
    result_df.to_csv(file_name, header=True, index=True, sep='\t')

    return result_df

def run_net_correlation(run_parameters):
    ''' wrapper: call sequence to perform gene prioritization

    Args:
        run_parameters: dict object with keys:
                run_parameters["spreadsheet_name_full_path"]    (samples x genes spreadsheet)
                run_parameters["drug_response_full_path"]       (one drug response spreadsheet)
                run_parameters['gg_network_name_full_path']     (gene, gene, weight,...   network)

    Returns:
        result_df: result dataframe of gene prioritization. Values are pearson
        correlation coefficient values in descending order.
    '''
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])
    drug_response = kn.get_spreadsheet_df(run_parameters["drug_response_full_path"])
    network_df = kn.get_network_df(run_parameters['gg_network_name_full_path'])
    node_1_names, node_2_names = kn.extract_network_node_names(network_df)
    unique_gene_names = kn.find_unique_node_names(node_1_names, node_2_names)
    genes_lookup_table = kn.create_node_names_dict(unique_gene_names)

    network_df = kn.map_node_names_to_index(network_df, genes_lookup_table, 'node_1')
    network_df = kn.map_node_names_to_index(network_df, genes_lookup_table, 'node_2')

    network_df = kn.symmetrize_df(network_df)
    network_mat = kn.convert_network_df_to_sparse(
        network_df, len(unique_gene_names), len(unique_gene_names))

    #network_mat = kn.normalize_sparse_mat_by_diagonal(network_mat)
    spreadsheet_df = kn.update_spreadsheet_df(spreadsheet_df, unique_gene_names)
    spreadsheet_mat = spreadsheet_df.as_matrix()
    sample_names = spreadsheet_df.columns

    sample_smooth, iterations = kn.smooth_matrix_with_rwr(
        spreadsheet_mat, network_mat, run_parameters)

    pc_array = perform_pearson_correlation(sample_smooth, drug_response.values[0])
    result_df = pd.DataFrame(pc_array, index=spreadsheet_df.index.values,
                             columns=['PCC']).abs().sort_values("PCC", ascending=0)

    target_file_base_name = "gene_drug_network_correlation"
    target_file_base_name = os.path.join(run_parameters["results_directory"], target_file_base_name)
    file_name = kn.create_timestamped_filename(target_file_base_name)
    file_name = file_name + '.df'
    result_df.to_csv(file_name, header=True, index=True, sep='\t')

    return result_df


def run_bootstrap_correlation(run_parameters):
    """ wrapper: call sequence to perform gene prioritization using bootstrap sampling
    Args:
        run_parameters: dict object with keys:
            run_parameters["spreadsheet_name_full_path"]    (samples x genes spreadsheet)
            run_parameters["drug_response_full_path"]       (one drug response spreadsheet)

    Returns:
        result_df: result dataframe of gene prioritization. Values are pearson
        correlation coefficient values in descending order.
    """
    print('run_bootstrap_correlation called')

    return


def run_bootstrap_net_correlation(run_parameters):
    """ wrapper: call sequence to perform gene prioritization using bootstrap sampling and network smoothing
    Args:
        run_parameters: dict object with keys:
            run_parameters["spreadsheet_name_full_path"]    (samples x genes spreadsheet)
            run_parameters["drug_response_full_path"]       (one drug response spreadsheet)
            run_parameters['gg_network_name_full_path']     (gene, gene, weight,...   network)

    Returns:
        result_df: result dataframe of gene prioritization. Values are pearson
        correlation coefficient values in descending order.
    """
    print('run_bootstrap_net_correlation called')

    return