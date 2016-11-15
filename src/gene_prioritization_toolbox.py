"""
@author: The KnowEnG dev team
"""
import os
import numpy as np
import pandas as pd

from scipy.stats import pearsonr as pcc
from scipy.stats.mstats import zscore
from scipy.stats.mstats import gmean
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

    generate_correlation_output(pc_array, drug_response_df.index.values[0], spreadsheet_df.index, run_parameters)

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
        (drug_name_list, gene_name_list, abs(pc_array), abs(pc_array), pc_array))

    df_header = ['Response', 'Gene ENSEMBL ID', 'quantitative sorting score', 'visualization score', 'baseline score']
    result_df = pd.DataFrame(output_val, columns=df_header).sort_values("visualization score", ascending=0)
    result_df.index = range(result_df.shape[0])
    target_file_base_name = os.path.join(run_parameters["results_directory"], drug_name + '_' + "correlation_final_result")
    file_name = kn.create_timestamped_filename(target_file_base_name) + '.tsv'
    result_df.to_csv(file_name, header=True, index=False, sep='\t')

    return 

def run_bootstrap_correlation(run_parameters):
    """ perform gene prioritization using bootstrap sampling

    Args:
        run_parameters: parameter set dictionary.
    """
    if run_parameters["number_of_bootstraps"] <= 1:
        run_correlation(run_parameters)
        return

    drug_response_df = kn.get_spreadsheet_df(run_parameters["drug_response_full_path"])
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])

    pearson_array = get_correlation(spreadsheet_df.as_matrix(), drug_response_df.values[0], run_parameters)
    pc_array = np.int_(np.zeros(spreadsheet_df.shape[0]))
    for bootstrap_number in range(0, run_parameters["number_of_bootstraps"]):
        sample_random, sample_permutation = sample_a_matrix_pearson(
            spreadsheet_df.as_matrix(), run_parameters["rows_sampling_fraction"],
            run_parameters["cols_sampling_fraction"])

        drug_response = drug_response_df.values[0, None]
        drug_response = drug_response[0, sample_permutation]
        pc_array = get_correlation(sample_random, drug_response, run_parameters)
        save_a_sample_correlation(pc_array, run_parameters)

    pcc_gm_array, borda_count= get_bootstrap_correlation_score(run_parameters, pc_array.size)
    kn.remove_dir(run_parameters["pc_array_tmp_dir"])
    run_parameters['out_filename'] = 'bootstrap_correlation'
    generate_bootstrap_correlation_output(borda_count, pcc_gm_array, pearson_array, drug_response_df.index.values[0],
                                          spreadsheet_df.index, run_parameters)
    return

def generate_bootstrap_correlation_output(borda_count, pcc_gm_array, pc_array, drug_name, gene_name_list, run_parameters):
    """ Save final output of correlation

    Args:
        pc_array:
        drug_name:
        gene_name_list:
        run_parameters:
    """
    drug_name_list = np.repeat(drug_name, len(gene_name_list))
    output_val = np.column_stack(
        (drug_name_list, gene_name_list, borda_count, pcc_gm_array, pc_array))

    df_header = ['Response', 'Gene ENSEMBL ID', 'quantitative sorting score', 'visualization score', 'baseline score']
    result_df = pd.DataFrame(output_val, columns=df_header).sort_values("quantitative sorting score", ascending=0)
    result_df.index = range(result_df.shape[0])
    target_file_base_name = os.path.join(run_parameters["results_directory"], drug_name + '_' + run_parameters['out_filename'])
    file_name = kn.create_timestamped_filename(target_file_base_name) + '.tsv'
    result_df.to_csv(file_name, header=True, index=False, sep='\t')

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
    pearson_array = pc_array.copy()
    pc_array[~np.in1d(spreadsheet_df.index, spreadsheet_genes_as_input)] = 0.0
    pc_array = np.abs(trim_to_top_beta(pc_array, run_parameters["top_beta_of_sort"]))
    restart_accumulator = pc_array.copy()
    restart_accumulator[restart_accumulator != 0] = 1

    pc_array = pc_array / sum(pc_array)
    pc_array = kn.smooth_matrix_with_rwr(pc_array, network_mat, run_parameters)[0]

    pc_array = pc_array - baseline_array
    min_max_pc = (pc_array - min(pc_array)) / (max(pc_array) - min(pc_array))

    generate_net_correlation_output(
        pearson_array, pc_array, min_max_pc, restart_accumulator, drug_response_df.index.values[0],
        spreadsheet_df.index, spreadsheet_genes_as_input, run_parameters)

    return

def generate_net_correlation_output(pearson_array, pc_array, min_max_pc, restart_accumulator,
                                    drug_name, gene_name_list, gene_orig_list, run_parameters):
    """ Save final output of correlation
    
    Args:
        pearson_array:
        pc_array:
        drug_name:
        gene_name_list:
        gene_orig_list:
        run_parameters:
    """
    mask = np.in1d(gene_name_list, gene_orig_list)
    pc_array = pc_array[mask]
    pearson_array = pearson_array[mask]
    min_max_pc = min_max_pc[mask]
    gene_name_list = gene_name_list[mask]
    restart_accumulator = restart_accumulator[mask]
    drug_name_list = np.repeat(drug_name, len(gene_name_list))

    output_val = np.column_stack(
        (drug_name_list, gene_name_list, pc_array, min_max_pc, pearson_array, restart_accumulator))

    df_header = ['Response', 'Gene ENSEMBL ID', 'quantitative sorting score', 'visualization score',
                 'baseline score', 'Percent appearing in restart set']
    result_df = pd.DataFrame(output_val, columns=df_header).sort_values("visualization score", ascending=0)

    target_file_base_name = os.path.join(run_parameters["results_directory"], drug_name + '_' + "net_correlation_final_result")
    file_name = kn.create_timestamped_filename(target_file_base_name) + '.tsv'
    result_df.to_csv(file_name, header=True, index=False, sep='\t')

    return 

def run_bootstrap_net_correlation(run_parameters):
    """ perform gene prioritization using bootstrap sampling and network smoothing

    Args:
        run_parameters: parameter set dictionary.
    """
    if run_parameters["number_of_bootstraps"] <= 1:
        run_net_correlation(run_parameters)
        return

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

    restart_accumulator = np.zeros(network_mat.shape[0])
    baseline_array = np.ones(network_mat.shape[0]) / network_mat.shape[0]
    baseline_array = kn.smooth_matrix_with_rwr(baseline_array, network_mat, run_parameters)[0]

    pearson_array = get_correlation(sample_smooth,  drug_response_df.values[0], run_parameters)
    for bootstrap_number in range(0, run_parameters["number_of_bootstraps"]):
        sample_random, sample_permutation = sample_a_matrix_pearson(
            sample_smooth, run_parameters["rows_sampling_fraction"],
            run_parameters["cols_sampling_fraction"])

        drug_response = drug_response_df.values[0, None]
        drug_response = drug_response[0, sample_permutation]
        pc_array = get_correlation(sample_random, drug_response, run_parameters)
        pc_array[~np.in1d(spreadsheet_df.index, spreadsheet_genes_as_input)] = 0.0
        pc_array = np.abs(trim_to_top_beta(pc_array, run_parameters["top_beta_of_sort"]))
        restart_accumulator[pc_array != 0] += 1

        pc_array = pc_array / sum(pc_array)
        pc_array = kn.smooth_matrix_with_rwr(pc_array, network_mat, run_parameters)[0]
        pc_array = pc_array - baseline_array
        save_a_sample_correlation(pc_array, run_parameters)

    restart_accumulator = restart_accumulator / run_parameters["number_of_bootstraps"]
    run_parameters['out_filename'] = 'bootstrap_net_correlation'
    pcc_gm_array, borda_count = get_bootstrap_net_correlation_score(run_parameters, pearson_array.size)
    kn.remove_dir(run_parameters["pc_array_tmp_dir"])

    generate_net_correlation_output(
        pearson_array, borda_count, pcc_gm_array, restart_accumulator, drug_response_df.index.values[0],
        spreadsheet_df.index, spreadsheet_genes_as_input, run_parameters)
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

    uniq_vals, val_counts = np.unique(corr_array, return_counts=True)
    if uniq_vals.size < corr_array.size:
        for k in range(0, uniq_vals.size):
            if val_counts[k] > 1:
                borda_count[corr_array == uniq_vals[k]] = sum(borda_count[corr_array == uniq_vals[k]]) / val_counts[k]

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

def save_a_sample_correlation(pc_array, run_parameters):
    """ Save a correlation array to the pc_array_tmp_dir
    Args:
        pc_array:
        run_paramters: with key 'pc_array_tmp_dir'
    """
    tmp_dir = run_parameters['pc_array_tmp_dir']
    os.makedirs(tmp_dir, mode=0o755, exist_ok=True)
    tmp_filename = kn.create_timestamped_filename('sampled_pc_array', name_extension=None, precision=1e12)

    pc_array_filename = os.path.join(tmp_dir, tmp_filename)
    with open(pc_array_filename , 'wb') as tmp_file_handle:
        pc_array.dump(tmp_file_handle)

    return


def get_bootstrap_correlation_score(run_parameters, n_rows):
    """

    """
    tmp_dir = run_parameters['pc_array_tmp_dir']
    dir_list = os.listdir(tmp_dir)
    n_cols = len(dir_list)

    pc_array_vectors = np.zeros((n_rows, n_cols))
    borda_count = np.zeros(n_rows)
    current_column = 0
    for tmp_f in dir_list:
        if tmp_f[0:16] == 'sampled_pc_array':
            pc_name = os.path.join(tmp_dir, tmp_f)
            corr_array = np.load(pc_name)
            sum_vote_to_borda_count(borda_count, np.abs(corr_array))
            pc_array_vectors[:,current_column] = corr_array
            current_column += 1

    borda_count = borda_count / max(borda_count)

    pc_array_vectors = np.abs(pc_array_vectors[:,0:current_column - 1])
    pc_array_vectors_gm_max = max(gmean(pc_array_vectors, axis=1))
    pcc_gm_array = gmean(pc_array_vectors / pc_array_vectors_gm_max, axis=1)

    return pcc_gm_array, borda_count

def get_bootstrap_net_correlation_score(run_parameters, n_rows):
    """
    pcc_gm_array, borda_count = get_bootstrap_net_correlation_score(run_parameters, n_rows)
    """
    tmp_dir = run_parameters['pc_array_tmp_dir']
    dir_list = os.listdir(tmp_dir)
    n_cols = len(dir_list)
    pc_array_vectors = np.zeros((n_rows, n_cols))
    borda_count = np.zeros(n_rows)
    current_column = 0
    for tmp_f in dir_list:
        if tmp_f[0:16] == 'sampled_pc_array':
            pc_name = os.path.join(tmp_dir, tmp_f)
            corr_array = np.load(pc_name)
            sum_vote_to_borda_count(borda_count, corr_array)
            corr_array = corr_array - min(corr_array)
            pc_array_vectors[:,current_column] = corr_array / max(corr_array)
            current_column += 1

    borda_count = borda_count / max(borda_count)

    pc_array_vectors = pc_array_vectors[:,0:current_column - 1]
    pcc_gm_array = gmean(pc_array_vectors, axis=1)

    return pcc_gm_array, borda_count
