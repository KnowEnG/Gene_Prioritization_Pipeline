"""
Created on Fri Sep 23 16:39:35 2016
@author: The Gene Prioritization dev team
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr as pcc
from sklearn.preprocessing import normalize
from sklearn.linear_model import LassoCV
import knpackage.toolbox as kn


def perform_lasso_cv_regression(spreadsheet, drug_response, normalize=True, max_iter=1e6):
    """ perform lasso with cross validation

    Args:
        spreadsheet:        features x samples matrix (np.array)
        drug_response:      1 x samples array (np.array)
        normalize_lasso:    (default=true) normalize inputs before
        max_iter:           (default 1000000) integer limit for iterations to converge

    Returns:
        coef_list:          features coefficient list
        lasso_predictor:    lasso_obj.predict (np.array)
    """
    lasso_cv_obj = LassoCV(normalize=normalize, max_iter=max_iter)
    lasso_cv_residual = lasso_cv_obj.fit(spreadsheet.T, drug_response[0])
    coef_array = lasso_cv_residual.coef_
    lasso_predictor = lasso_cv_residual.predict(spreadsheet.T)

    return coef_array, lasso_predictor

def perform_pearson_correlation(spreadsheet, drug_response):
    """ Find pearson correlation coefficient(PCC) for each gene expression (spreadsheet row)
            with the drug response.

    Args:
        spreadsheet: genes x samples gene expression data
        drug_response: one drug response for each sample

    Returns:
        pc_array: row-spreadsheet-features-ordered array of pearson correlation coefficients
    """
    pc_array = np.zeros(spreadsheet.shape[0])
    for row in range(0, spreadsheet.shape[0]):
        pcc_value = pcc(spreadsheet[row, :], drug_response)[0]
        pc_array[row] = pcc_value

    return pc_array


def run_gene_correlation(run_parameters):
    ''' pearson cc:  call sequence to perform gene prioritization

    Args:
        run_parameters: dict object with keys:
                run_parameters["spreadsheet_name_full_path"]
                run_parameters["drug_response_full_path"]

    Returns:
        result_df:  result - dataframe of gene prioritization. Values are pearson
                    correlation coefficient values in descending order.
    '''
    result_df = perform_gene_correlation(run_parameters)
    target_file_base_name = "gene_drug_correlation"
    target_file_base_name = os.path.join(run_parameters["results_directory"], target_file_base_name)
    file_name = kn.create_timestamped_filename(target_file_base_name)
    file_name = file_name + '.df'
    result_df.to_csv(file_name, header=True, index=True, sep='\t')

    return


def perform_gene_correlation(run_parameters):
    """ find the pearson correlation coefficient per gene

    Args:
        run_parameters: dict with keys:
                        spreadsheet_name_full_path
                        drug_response_full_path
    Returns:
        result_df: dataframe with genes and correlation columns
    """
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])
    drug_response = kn.get_spreadsheet_df(run_parameters["drug_response_full_path"])
    pc_array = perform_pearson_correlation(spreadsheet_df.values, drug_response.values[0])
    result_df = pd.DataFrame(pc_array, index=spreadsheet_df.index.values,
                             columns=['PCC']).abs().sort_values("PCC", ascending=0)
    return result_df


def run_net_correlation(run_parameters):
    ''' pearson cc:  call sequence to perform gene prioritization with network smoothing

    Args:
        run_parameters: dict object with keys:
                run_parameters["spreadsheet_name_full_path"]    (samples x genes spreadsheet)
                run_parameters["drug_response_full_path"]       (one drug response spreadsheet)
                run_parameters['gg_network_name_full_path']     (gene, gene, weight,...   network)

    Returns:
        result_df: result dataframe of gene prioritization. Values are pearson
        correlation coefficient values in descending order.
    '''
    result_df = perform_net_correlation(run_parameters)
    target_file_base_name = "gene_drug_network_correlation"
    target_file_base_name = os.path.join(run_parameters["results_directory"], target_file_base_name)
    file_name = kn.create_timestamped_filename(target_file_base_name)
    file_name = file_name + '.df'
    result_df.to_csv(file_name, header=True, index=True, sep='\t')

    return


def perform_net_correlation(run_parameters):
    ''' perform gene prioritization with network smoothing and Pearson correlation

    Args:
        run_parameters: dict object with keys:
                run_parameters["spreadsheet_name_full_path"]    (samples x genes spreadsheet)
                run_parameters["drug_response_full_path"]       (one drug response spreadsheet)
                run_parameters['gg_network_name_full_path']     (gene, gene, weight,...   network)

    Returns:
        result_df: result dataframe of gene prioritization. Values are pearson
        correlation coefficient values in descending order.
    '''
    spreadsheet_df, sample_smooth = get_smooth_spreadsheet_matrix(run_parameters)
    drug_response = kn.get_spreadsheet_df(run_parameters["drug_response_full_path"])

    pc_array = perform_pearson_correlation(sample_smooth, drug_response.values[0])
    result_df = pd.DataFrame(pc_array, index=spreadsheet_df.index.values,
                             columns=['PCC']).abs().sort_values("PCC", ascending=0)

    return result_df


def run_bootstrap_correlation(run_parameters):
    """ pearson cc: call sequence to perform gene prioritization using bootstrap sampling

    Args:
        run_parameters: dict object with keys:
            run_parameters["spreadsheet_name_full_path"]    (samples x genes spreadsheet)
            run_parameters["drug_response_full_path"]       (one drug response spreadsheet)

    Returns:
        result_df: result dataframe of gene prioritization. Values are pearson
        correlation coefficient values in descending order.
    """
    result_df = perform_bootstrap_correlation(run_parameters)
    target_file_base_name = "gene_drug_bootstrap_correlation"
    target_file_base_name = os.path.join(run_parameters["results_directory"], target_file_base_name)
    file_name = kn.create_timestamped_filename(target_file_base_name)
    file_name = file_name + '.df'
    result_df.to_csv(file_name, header=True, index=True, sep='\t')

    return


def perform_bootstrap_correlation(run_parameters):
    """ perform gene prioritization using bootstrap sampling with Pearson correlation

    Args:
        run_parameters: dict object with keys:
            run_parameters["spreadsheet_name_full_path"]    (samples x genes spreadsheet)
            run_parameters["drug_response_full_path"]       (one drug response spreadsheet)

    Returns:
        result_df: result dataframe of gene prioritization. Values are pearson
        correlation coefficient values in descending order.
    """
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])
    drug_response = kn.get_spreadsheet_df(run_parameters["drug_response_full_path"])
    sample_smooth = spreadsheet_df.values

    borda_count = np.int_(np.zeros(sample_smooth.shape[0]))
    for bootstrap_number in range(0, int(run_parameters["number_of_bootstraps"])):
        sample_random, sample_permutation = kn.sample_a_matrix(
            sample_smooth, np.float64(run_parameters["rows_sampling_fraction"]),
            np.float64(run_parameters["cols_sampling_fraction"]))

        D = drug_response.values[0, None]
        D = D[0, sample_permutation]
        pc_array = perform_pearson_correlation(sample_random, D)
        pc_array_idx = np.argsort(pc_array)[::-1]
        borda_count = sum_vote_perm_to_borda_count(borda_count, pc_array_idx)

    borda_count = borda_count / max(borda_count)
    result_df = pd.DataFrame(borda_count, index=spreadsheet_df.index.values,
                             columns=['PCC']).abs().sort_values("PCC", ascending=0)

    return result_df


def run_bootstrap_net_correlation(run_parameters):
    """ pearson cc: call sequence to perform gene prioritization using bootstrap sampling and network smoothing

    Args:
        run_parameters: dict object with keys:
            run_parameters["spreadsheet_name_full_path"]    (samples x genes spreadsheet)
            run_parameters["drug_response_full_path"]       (one drug response spreadsheet)
            run_parameters['gg_network_name_full_path']     (gene, gene, weight,...   network)

    Returns:
        result_df: result dataframe of gene prioritization. Values are pearson
        correlation coefficient values in descending order.
    """
    result_df = perform_bootstrap_net_correlation(run_parameters)
    target_file_base_name = "gene_drug_network_bootstrap_correlation"
    target_file_base_name = os.path.join(run_parameters["results_directory"], target_file_base_name)
    file_name = kn.create_timestamped_filename(target_file_base_name)
    file_name = file_name + '.df'
    result_df.to_csv(file_name, header=True, index=True, sep='\t')

    return


def perform_bootstrap_net_correlation(run_parameters):
    """ perform gene prioritization using bootstrap sampling and network smoothing before Pearson correlation
    Args:
        run_parameters: dict object with keys:
            run_parameters["spreadsheet_name_full_path"]    (samples x genes spreadsheet)
            run_parameters["drug_response_full_path"]       (one drug response spreadsheet)
            run_parameters['gg_network_name_full_path']     (gene, gene, weight,...   network)

    Returns:
        result_df: result dataframe of gene prioritization. Values are pearson
                    correlation coefficient values in descending order.
    """
    spreadsheet_df, sample_smooth = get_smooth_spreadsheet_matrix(run_parameters)
    drug_response = kn.get_spreadsheet_df(run_parameters["drug_response_full_path"])

    borda_count = np.int_(np.zeros(sample_smooth.shape[0]))
    for bootstrap_number in range(0, int(run_parameters["number_of_bootstraps"])):
        sample_random, sample_permutation = kn.sample_a_matrix(
            sample_smooth, np.float64(run_parameters["rows_sampling_fraction"]),
            np.float64(run_parameters["cols_sampling_fraction"]))

        D = drug_response.values[0, None]
        D = D[0, sample_permutation]
        pc_array = perform_pearson_correlation(sample_random, D)
        pc_array_idx = np.argsort(pc_array)[::-1]
        borda_count = sum_vote_perm_to_borda_count(borda_count, pc_array_idx)

    borda_count = borda_count / max(borda_count)
    result_df = pd.DataFrame(borda_count, index=spreadsheet_df.index.values,
                             columns=['PCC']).abs().sort_values("PCC", ascending=0)

    return result_df


def sum_vote_perm_to_borda_count(borda_count, vote_rank, vote_perm=None):
    """ incrementally update count by borda weighted vote in a subsample of the full borda size

    Args:
        vote_rank: (np.int_(vote_rank)) python rank array same size as borda_count (0 == first, 1 == second,...)
        borda_count: (np.int_(borda_count) ) the current running total of Borda weighted votes
        vote_perm: (None) or: the sample permutation of the vote ranking

    Returns:
        borda_count: input borda_count with borda weighted vote rankings added
    """
    rank_array = np.int_(sorted(np.arange(0, vote_rank.size) + 1, reverse=True))
    if vote_perm is None:
        borda_count += rank_array[np.int_(vote_rank)]
    else:
        borda_count[vote_perm] += rank_array[np.int_(vote_rank)]

    return borda_count


def get_smooth_spreadsheet_matrix(run_parameters):
    """ common task collection for reading files and performaing random walk with restart

    Args:
        run_parameters: dict with keys:
                        spreadsheet_name_full_path
                        gg_network_name_full_path
    Returns:
        spreadsheet_df: the spreadsheet dataframe
        sample_smooth: the spreadsheet matrix smoothed with the network
    """
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])
    network_df = kn.get_network_df(run_parameters['gg_network_name_full_path'])
    node_1_names, node_2_names = kn.extract_network_node_names(network_df)
    unique_gene_names = kn.find_unique_node_names(node_1_names, node_2_names)
    genes_lookup_table = kn.create_node_names_dict(unique_gene_names)

    network_df = kn.map_node_names_to_index(network_df, genes_lookup_table, 'node_1')
    network_df = kn.map_node_names_to_index(network_df, genes_lookup_table, 'node_2')

    network_df = kn.symmetrize_df(network_df)
    network_mat = kn.convert_network_df_to_sparse(
        network_df, len(unique_gene_names), len(unique_gene_names))

    network_mat = normalize(network_mat, norm="l1", axis=0)
    spreadsheet_df = kn.update_spreadsheet_df(spreadsheet_df, unique_gene_names)
    spreadsheet_mat = spreadsheet_df.as_matrix()

    sample_smooth, iterations = kn.smooth_matrix_with_rwr(spreadsheet_mat, network_mat, run_parameters)

    return spreadsheet_df, sample_smooth
