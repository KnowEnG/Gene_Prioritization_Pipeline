"""
@author: The KnowEnG dev team
"""
import os
import numpy as np
import pandas as pd

from scipy.stats import pearsonr as pcc
from scipy.stats.mstats import zscore
from sklearn.preprocessing import normalize
from sklearn.linear_model import LassoCV

import knpackage.toolbox as kn

def run_correlation(run_parameters):
    ''' perform gene prioritization

    Args:
        run_parameters: parameter set dictionary.
    '''
    drug_response_df = kn.get_spreadsheet_df(run_parameters["drug_response_full_path"])
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])

    pc_array = get_correlation(spreadsheet_df.as_matrix(), drug_response_df.values[0], run_parameters)

    result_df = pd.DataFrame(pc_array, index=spreadsheet_df.index.values,
                             columns=['PCC']).abs().sort_values("PCC", ascending=0)

    write_results_dataframe(result_df, run_parameters["results_directory"], "gene_drug_correlation")

    generate_correlation_output(pc_array, drug_response_df.index.values, spreadsheet_df.index, run_parameters)

    return

def generate_correlation_output(pc_array, drug_name, gene_name_list, run_parameters):
    """ Save final output of correlation
    
    Args:
        pc_array:
        drug_name:
        gene_name_list:
        run_parameters:
    """
    drug_name_list = np.repeat(drug_name, len(gene_name_list))
    output_val = np.column_stack(
        (drug_name_list, gene_name_list, np.abs(pc_array), np.abs(pc_array), pc_array))

    df_header = ['Response', 'Gene ENSEMBL ID', 'quantitative sorting score', 'visualization score', 'baseline score']
    result_df = pd.DataFrame(output_val, columns=df_header).sort_values("visualization score", ascending=0)
    result_df.index = range(result_df.shape[0])
    target_file_base_name = os.path.join(run_parameters["results_directory"], "correlation_final_result")
    file_name = kn.create_timestamped_filename(target_file_base_name) + '.df'
    result_df.to_csv(file_name, header=True, index=False, sep='\t')

    return 

def run_bootstrap_correlation(run_parameters):
    """ perform gene prioritization using bootstrap sampling

    Args:
        run_parameters: parameter set dictionary.
    """
    drug_response_df = kn.get_spreadsheet_df(run_parameters["drug_response_full_path"])
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])

    borda_count = np.int_(np.zeros(spreadsheet_df.shape[0]))
    for bootstrap_number in range(0, run_parameters["number_of_bootstraps"]):
        sample_random, sample_permutation = sample_a_matrix_pearson(
            spreadsheet_df.as_matrix(), run_parameters["rows_sampling_fraction"],
            run_parameters["cols_sampling_fraction"])

        drug_response = drug_response_df.values[0, None]
        drug_response = drug_response[0, sample_permutation]
        pc_array = get_correlation(sample_random, drug_response, run_parameters)

        borda_count = sum_vote_to_borda_count(borda_count, np.abs(pc_array))

    borda_count = borda_count / max(borda_count)
    result_df = pd.DataFrame(borda_count, index=spreadsheet_df.index.values,
                             columns=['PCC']).sort_values("PCC", ascending=0)
    write_results_dataframe(result_df, run_parameters["results_directory"], "gene_drug_bootstrap_correlation")
    return


def run_net_correlation(run_parameters):
    """ perform gene prioritization with network smoothing

    Args:
        run_parameters: parameter set dictionary.
    """
    drug_response_df = kn.get_spreadsheet_df(run_parameters["drug_response_full_path"])
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])
    spreadsheet_df = zscore_dataframe(spreadsheet_df)
    spreadsheet_genes_as_input = spreadsheet_df.index.values

    network_df = kn.get_network_df(run_parameters['gg_network_name_full_path'])
    node_1_names, node_2_names = kn.extract_network_node_names(network_df)
    unique_gene_names = kn.find_unique_node_names(node_1_names, node_2_names)

    unique_gene_names = sorted(unique_gene_names)

    genes_lookup_table = kn.create_node_names_dict(unique_gene_names)

    network_df = kn.map_node_names_to_index(network_df, genes_lookup_table, 'node_1')
    network_df = kn.map_node_names_to_index(network_df, genes_lookup_table, 'node_2')

    network_df = kn.symmetrize_df(network_df)
    network_mat_sparse = kn.convert_network_df_to_sparse(
        network_df, len(unique_gene_names), len(unique_gene_names))

    network_mat = normalize(network_mat_sparse, norm="l1", axis=0)
    spreadsheet_df = kn.update_spreadsheet_df(spreadsheet_df, unique_gene_names)

    sample_smooth, iterations = kn.smooth_matrix_with_rwr(spreadsheet_df.as_matrix(), network_mat.T, run_parameters)

    baseline_array = np.ones(network_mat.shape[0]) / network_mat.shape[0]
    baseline_array = kn.smooth_matrix_with_rwr(baseline_array, network_mat, run_parameters)[0]

    pc_array = get_correlation(sample_smooth,  drug_response_df.values[0], run_parameters)
    Pearson_array = pc_array.copy()
    pc_array[~np.in1d(spreadsheet_df.index, spreadsheet_genes_as_input)] = 0.0
    pc_array = np.abs(trim_to_top_beta(pc_array, run_parameters["top_beta_of_sort"]))
    pc_array = pc_array / sum(pc_array)
    pc_array = kn.smooth_matrix_with_rwr(pc_array, network_mat, run_parameters)[0]

    pc_array = pc_array - baseline_array
    result_df = pd.DataFrame(pc_array, index=spreadsheet_df.index.values,
                             columns=['run_net_correlation']).sort_values("run_net_correlation",
                                                                                    ascending=0)
    result_df = result_df.loc[result_df.index.isin(spreadsheet_genes_as_input)]

    write_results_dataframe(result_df, run_parameters["results_directory"], "gene_drug_network_correlation")

    generate_net_correlation_output(
        Pearson_array, pc_array, drug_response_df.index.values, spreadsheet_df.index, spreadsheet_genes_as_input, run_parameters)

    return

def generate_net_correlation_output(Pearson_array, pc_array, drug_name, gene_name_list, gene_orig_list, run_parameters):
    """ Save final output of correlation
    
    Args:
        Pearson_array:
        pc_array:
        drug_name:
        gene_name_list:
        gene_orig_list:
        run_parameters:
    """
    mask = np.in1d(gene_name_list, gene_orig_list)
    pc_array = pc_array[mask]
    Pearson_array = Pearson_array[mask]
    min_max_pc = (pc_array-min(pc_array))/(max(pc_array)-min(pc_array))

    gene_name_list = gene_name_list[mask]
    drug_name_list = np.repeat(drug_name, len(gene_name_list))
    output_val = np.column_stack(
        (drug_name_list, gene_name_list, pc_array, min_max_pc, Pearson_array, np.ones(len(gene_name_list), dtype=np.int)))

    df_header = ['Response', 'Gene ENSEMBL ID', 'quantitative sorting score', 'visualization score', \
    'baseline score', 'Percent of appearing in restart set']
    result_df = pd.DataFrame(output_val, columns=df_header).sort_values("visualization score", ascending=0)


    target_file_base_name = os.path.join(run_parameters["results_directory"], "net_correlation_final_result")
    file_name = kn.create_timestamped_filename(target_file_base_name) + '.df'
    result_df.to_csv(file_name, header=True, index=False, sep='\t')

    return 

def run_bootstrap_net_correlation(run_parameters):
    """ perform gene prioritization using bootstrap sampling and network smoothing

    Args:
        run_parameters: parameter set dictionary.
    """
    drug_response_df = kn.get_spreadsheet_df(run_parameters["drug_response_full_path"])
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])
    spreadsheet_df = zscore_dataframe(spreadsheet_df)
    spreadsheet_genes_as_input = spreadsheet_df.index.values

    network_df = kn.get_network_df(run_parameters['gg_network_name_full_path'])
    node_1_names, node_2_names = kn.extract_network_node_names(network_df)
    unique_gene_names = kn.find_unique_node_names(node_1_names, node_2_names)

    unique_gene_names = sorted(unique_gene_names)

    genes_lookup_table = kn.create_node_names_dict(unique_gene_names)

    network_df = kn.map_node_names_to_index(network_df, genes_lookup_table, 'node_1')
    network_df = kn.map_node_names_to_index(network_df, genes_lookup_table, 'node_2')

    network_df = kn.symmetrize_df(network_df)
    network_mat_sparse = kn.convert_network_df_to_sparse(
        network_df, len(unique_gene_names), len(unique_gene_names))

    network_mat = normalize(network_mat_sparse, norm="l1", axis=0)
    spreadsheet_df = kn.update_spreadsheet_df(spreadsheet_df, unique_gene_names)

    sample_smooth, iterations = kn.smooth_matrix_with_rwr(spreadsheet_df.as_matrix(), network_mat.T, run_parameters)

    baseline_array = np.ones(network_mat.shape[0]) / network_mat.shape[0]
    baseline_array = kn.smooth_matrix_with_rwr(baseline_array, network_mat, run_parameters)[0]

    borda_count = np.zeros(sample_smooth.shape[0])
    for bootstrap_number in range(0, run_parameters["number_of_bootstraps"]):
        sample_random, sample_permutation = sample_a_matrix_pearson(
            sample_smooth, run_parameters["rows_sampling_fraction"],
            run_parameters["cols_sampling_fraction"])

        drug_response = drug_response_df.values[0, None]
        pc_array = get_correlation(sample_random, drug_response[0, sample_permutation], run_parameters)
        pc_array[~np.in1d(spreadsheet_df.index, spreadsheet_genes_as_input)] = 0.0
        pc_array = trim_to_top_beta(pc_array, run_parameters["top_beta_of_sort"])
        pc_array = kn.smooth_matrix_with_rwr(pc_array, network_mat, run_parameters)[0]
        pc_array = pc_array - baseline_array
        borda_count = sum_vote_to_borda_count(borda_count, pc_array)

    borda_count = borda_count / max(borda_count)
    result_df = pd.DataFrame(borda_count, index=spreadsheet_df.index.values,
                             columns=['run_bootstrap_net_correlation']).sort_values("run_bootstrap_net_correlation",
                                                                                    ascending=0)
    result_df = result_df.loc[result_df.index.isin(spreadsheet_genes_as_input)]

    write_results_dataframe(result_df, run_parameters["results_directory"], "gene_drug_network_bootstrap_correlation")
    return


def get_correlation(spreadsheet, drug_response, run_parameters, normalize=True, max_iter=1e6):
    """
    Args:
        spreadsheet:
        drug_response:
        run_parameters:
        normalize:
        max_iter:
    Returns:
        correlation_array
    """
    correlation_array = np.zeros(spreadsheet.shape[0])
    if 'correlation_method' in run_parameters:
        if run_parameters['correlation_method'] == 'pearson':
            spreadsheet = zscore(spreadsheet, axis=1, ddof=0)
            for row in range(0, spreadsheet.shape[0]):
                correlation_array[row] = pcc(spreadsheet[row, :], drug_response)[0]
            correlation_array[~(np.isfinite(correlation_array))] = 0
            return correlation_array

        if run_parameters['correlation_method'] == 'lasso':
            drug_response = np.array([drug_response])
            drug_response[ ~(np.isfinite(drug_response)) ] = 0
            lasso_cv_obj = LassoCV(normalize=normalize, max_iter=max_iter)
            lasso_cv_residual = lasso_cv_obj.fit(spreadsheet.T, drug_response[0])
            correlation_array = lasso_cv_residual.coef_

    return correlation_array


def sum_vote_to_borda_count(borda_count, corr_array):
    """ incrementally update borda count by borda weighted vote from correlation array

    Args:
        borda_count: (np.int_(borda_count) ) the current running total of Borda weighted votes
        corr_array:  correlation ranking largest value gets largest borda vote score

    Returns:
        borda_count: input borda_count with borda weighted vote rankings added
    """
    borda_count[np.argsort(corr_array)] += np.int_(sorted(np.arange(0, corr_array.size) + 1))
    return borda_count


def write_results_dataframe(result_df, run_dir, write_file_name):
    """ Writes a dataframe with a header and row names to a tab separated text file

    Args:
        result_df:          a dataframe with row and column names
        run_dir:            the directory to write in
        write_file_name:    the file name to join with the directory
    """
    target_file_base_name = os.path.join(run_dir, write_file_name)
    file_name = kn.create_timestamped_filename(target_file_base_name) + '.txt'
    result_df.to_csv(file_name, header=True, index=True, sep='\t')

    return

def sample_a_matrix_pearson(spreadsheet_mat, rows_fraction, cols_fraction):
    """ percent_sample x percent_sample random sample, from spreadsheet_mat.

    Args:
        spreadsheet_mat: gene x sample spread sheet as matrix.
        percent_sample: decimal fraction (slang-percent) - [0 : 1].

    Returns:
        sample_random: A specified precentage sample of the spread sheet.
        sample_permutation: the array that correponds to columns sample.
    """
    features_size = int(np.round(spreadsheet_mat.shape[0] * (1-rows_fraction)))
    features_permutation = np.random.permutation(spreadsheet_mat.shape[0])
    features_permutation = features_permutation[0:features_size].T

    patients_size = int(np.round(spreadsheet_mat.shape[1] * cols_fraction))
    sample_permutation = np.random.permutation(spreadsheet_mat.shape[1])
    sample_permutation = sample_permutation[0:patients_size]

    sample_random = spreadsheet_mat[:, sample_permutation]
    sample_random[features_permutation[:, None], :] = 0

    return sample_random, sample_permutation


def trim_to_top_beta(corr_arr, Beta):
    """ Preserve corr_arr order: set the top Beta members of the correlation array to one and set the rest to zero

    Args:
        corr_arr: an array of sortable values
        Beta:     a fraction between 0 and 1 to designate the top percentage of corr_arr to select

    Returns:
        corr_arr: the correlation array as binary with ones int the top Beta percent
    """
    Beta = min(corr_arr.size, np.round(corr_arr.size * Beta)) - 1
    Beta = int(max(Beta, 0))
    corr_arr[np.abs(corr_arr) < sorted(np.abs(corr_arr))[::-1][Beta]] = 0

    return corr_arr

def zscore_dataframe(gxs_df):
    """ zscore by rows for genes x samples dataframe
    Args:
        spreads_df
    Returns:
        spreadsheet_df: rows add up to zero, normalized to the mean and std deveiation
    """
    zscore_df = (gxs_df.sub(gxs_df.mean(axis=1), axis=0)).truediv(np.maximum(gxs_df.std(axis=1), 1e-12), axis=0)
    return zscore_df
