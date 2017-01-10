"""
@author: The KnowEnG dev team
"""
import gc
import os
import numpy as np
import pandas as pd

from scipy.stats import ttest_ind
from scipy.stats import pearsonr as pcc
from scipy.stats.mstats import zscore
from sklearn.preprocessing import normalize

import knpackage.toolbox as kn
import knpackage.distributed_computing_utils as dstutil
import knpackage.data_cleanup_toolbox as datacln

EPSILON_0 = 1e-7


def get_consolodated_dataframe(gene_samples_df, drug_samples_df):
    """  assemble the features X samples and data X samples dataframes into one

    Args:
        gene_samples_df:  spreadsheet with the same column headers as the drug_samples_df
        drug_samples_df:  spreadsheet of drug response rows with same column header

    Returns:
        genes_list:       the list of genes that are the rows of the spreadsheet in the consolodated dataframe
        drugs_list:       the list of drugs that for row names in the consolodated dataframe
        consolodated_df:  input spreadsheet dataframe with the drug dataframe appended
    """
    genes_list = gene_samples_df.index.values.tolist()
    drugs_list = drug_samples_df.index.values.tolist()
    consolodated_df = pd.DataFrame(
        gene_samples_df.as_matrix(), index=genes_list, columns=gene_samples_df.columns.values)
    consolodated_df = consolodated_df.append(drug_samples_df)

    return genes_list, drugs_list, consolodated_df


def get_data_for_drug(consolodated_dataframe, genes_list, drug_name, run_parameters):
    """ get spreadsheet and drug dataframes with the NA columns removed

    Args:
        consolodated_dataframe: the spreadsheet dataframe appended with the drug dataframe
        genes_list:             the row names that constitute the spreadsheet
        drug_name:              the individual drug to select
        run_parameters:         parameters dict

    Returns:
        drug_df:        single drug dataframe with NA columns removed
        spreadsheet_df: spreadsheet dataframe with same columns as drug_df
    """
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

    drug_df = pd.DataFrame(drug_array, index={drug_name}, columns=sample_names)
    spreadsheet_df = pd.DataFrame(spreadsheet_matrix, index=genes_list, columns=sample_names)

    return drug_df, spreadsheet_df


def run_correlation(run_parameters):
    """ perform gene prioritization

    Args:
        run_parameters: parameter set dictionary.
    """
    phenotype_df = kn.get_spreadsheet_df(run_parameters["drug_response_full_path"])
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])

    number_of_jobs = len(phenotype_df.index)
    list_of_job_id = range(0, number_of_jobs)
    zipped_arguments = dstutil.zip_parameters(run_parameters, spreadsheet_df, phenotype_df, list_of_job_id)
    dstutil.parallelize_processes_locally(worker_for_run_correlation, zipped_arguments, number_of_jobs)


def worker_for_run_correlation(run_parameters, spreadsheet_df, phenotype_df, i):
    """ core function for parallel run_correlation

    Args:
        run_parameters:  dict of parameters
        spreadsheet_df: spreadsheet - drug consolodated data frame
        genes_list:      ordered list of genes in consolodated_df
        drugs_list:      ordered list of drugs in consolodated_df
        i:               paralell iteration number
    """
    # selects the ith row in phenotype_df
    phenotype_df = phenotype_df.iloc[[i], :]

    spreadsheet_df_trimmed, phenotype_df_trimmed, err_msg = datacln.check_input_value_for_gene_prioritazion(
        spreadsheet_df, phenotype_df)

    pc_array = get_correlation(spreadsheet_df_trimmed.as_matrix(), phenotype_df_trimmed.values[0], run_parameters)

    generate_correlation_output(pc_array, phenotype_df.index.values[0], spreadsheet_df_trimmed.index, run_parameters)


def generate_correlation_output(pc_array, drug_name, gene_name_list, run_parameters):
    """ Save final output of correlation
    
    Args:
        pc_array: pearson correlation coefficient array
        drug_name: name of the drug
        gene_name_list: list of the genes correlated (size of pc_array
        run_parameters: dictionary of run parameters with key 'results_directory'
    """
    drug_name_list = np.repeat(drug_name, len(gene_name_list))
    output_val = np.column_stack(
        (drug_name_list, gene_name_list, abs(pc_array), abs(pc_array), pc_array))

    df_header = ['Response', 'Gene ENSEMBL ID', 'quantitative sorting score', 'visualization score', 'baseline score']
    result_df = pd.DataFrame(output_val, columns=df_header).sort_values("visualization score", ascending=0)
    result_df.index = range(result_df.shape[0])
    target_file_base_name = os.path.join(run_parameters["results_directory"],
                                         drug_name + '_' + run_parameters['out_filename'])
    file_name = kn.create_timestamped_filename(target_file_base_name) + '.tsv'
    result_df.to_csv(file_name, header=True, index=False, sep='\t')

    download_result_df = pd.DataFrame(data=None, index=None, columns=['Response', 'Gene ENSEMBL ID'])
    download_result_df['Response'] = result_df['Response']
    download_result_df['Gene ENSEMBL ID'] = result_df['Gene ENSEMBL ID']

    download_target_file_base_name = os.path.join(run_parameters["results_directory"],
                                                  drug_name + '_' + run_parameters['out_filename'])

    download_file_name = kn.create_timestamped_filename(download_target_file_base_name) + '_download' + '.tsv'
    download_result_df.to_csv(download_file_name, header=True, index=False, sep='\t')


def run_bootstrap_correlation(run_parameters):
    """ perform gene prioritization using bootstrap sampling

    Args:
        run_parameters: parameter set dictionary.
    """
    drug_response_df = kn.get_spreadsheet_df(run_parameters["drug_response_full_path"])
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])

    n_bootstraps = run_parameters["number_of_bootstraps"]

    number_of_jobs = len(drug_response_df.index)
    list_of_job_id = range(0, number_of_jobs)
    zipped_arguments = dstutil.zip_parameters(run_parameters, spreadsheet_df, drug_response_df, n_bootstraps, list_of_job_id)
    dstutil.parallelize_processes_locally(worker_for_run_bootstrap_correlation, zipped_arguments, number_of_jobs)


def worker_for_run_bootstrap_correlation(run_parameters, spreadsheet_df, phenotype_df, n_bootstraps, i):
    """  core function for parallel run_bootstrap_correlation

    Args:
        run_parameters:  dict of parameters
        consolodated_df: spreadsheet - drug consolodated data frame
        genes_list:      ordered list of genes in consolodated_df
        n_bootstraps:    number of bootstrap samples to use
        drugs_list:      ordered list of drugs in consolodated_df
        i:               paralell iteration number
    """
    phenotype_df = phenotype_df.iloc[[i], :]
    spreadsheet_df_trimmed, phenotype_df_trimmed, ret_msg = datacln.check_input_value_for_gene_prioritazion(
        spreadsheet_df, phenotype_df)

    pearson_array = get_correlation(spreadsheet_df_trimmed.as_matrix(), phenotype_df_trimmed.values[0], run_parameters)
    borda_count = np.zeros(spreadsheet_df.shape[0])
    gm_accumulator = np.ones(spreadsheet_df.shape[0])
    for bootstrap_number in range(0, n_bootstraps):
        sample_random, sample_permutation = sample_a_matrix_pearson(
            spreadsheet_df_trimmed.as_matrix(), run_parameters["rows_sampling_fraction"],
            run_parameters["cols_sampling_fraction"])
        drug_response = phenotype_df_trimmed.values[0, None]
        drug_response = drug_response[0, sample_permutation]
        pc_array = get_correlation(sample_random, drug_response, run_parameters)
        borda_count = sum_array_ranking_to_borda_count(borda_count, np.abs(pc_array))
        gm_accumulator = (np.abs(pc_array) + EPSILON_0) * gm_accumulator
    pcc_gm_array = gm_accumulator ** (1 / n_bootstraps)
    borda_count = borda_count / n_bootstraps
    generate_bootstrap_correlation_output(borda_count, pcc_gm_array, pearson_array,
                                          phenotype_df_trimmed.index.values[0],
                                          spreadsheet_df_trimmed.index, run_parameters)


def generate_bootstrap_correlation_output(borda_count, pcc_gm_array, pc_array, drug_name, gene_name_list,
                                          run_parameters):
    """ Save final output of correlation

    Args:
        pc_array: pearson correlation coefficient array
        drug_name: name of the drug
        gene_name_list: list of the genes correlated (size of pc_array
        run_parameters: dictionary of run parameters with key 'results_directory'
    """
    drug_name_list = np.repeat(drug_name, len(gene_name_list))
    output_val = np.column_stack(
        (drug_name_list, gene_name_list, borda_count, pcc_gm_array, pc_array))

    df_header = ['Response', 'Gene ENSEMBL ID', 'quantitative sorting score', 'visualization score', 'baseline score']
    result_df = pd.DataFrame(output_val, columns=df_header).sort_values("quantitative sorting score", ascending=0)
    result_df.index = range(result_df.shape[0])
    target_file_base_name = os.path.join(run_parameters["results_directory"],
                                         drug_name + '_' + run_parameters['out_filename'])
    file_name = kn.create_timestamped_filename(target_file_base_name) + '.tsv'
    result_df.to_csv(file_name, header=True, index=False, sep='\t')

    download_result_df = pd.DataFrame(data=None, index=None, columns=['Response', 'Gene ENSEMBL ID'])
    download_result_df['Response'] = result_df['Response']
    download_result_df['Gene ENSEMBL ID'] = result_df['Gene ENSEMBL ID']

    download_target_file_base_name = os.path.join(run_parameters["results_directory"],
                                                  drug_name + '_' + run_parameters['out_filename'])

    download_file_name = kn.create_timestamped_filename(download_target_file_base_name) + '_download' + '.tsv'
    download_result_df.to_csv(download_file_name, header=True, index=False, sep='\t')


def run_net_correlation(run_parameters):
    """ perform gene prioritization with network smoothing

    Args:
        run_parameters: parameter set dictionary.
    """
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

    del network_df
    del network_mat_sparse
    del node_1_names
    del node_2_names
    del genes_lookup_table
    gc.collect()

    drug_response_df = kn.get_spreadsheet_df(run_parameters["drug_response_full_path"])
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])
    spreadsheet_genes_as_input = spreadsheet_df.index.values

    spreadsheet_df = kn.update_spreadsheet_df(spreadsheet_df, unique_gene_names)
    spreadsheet_df = zscore_dataframe(spreadsheet_df)

    sample_smooth, iterations = kn.smooth_matrix_with_rwr(spreadsheet_df.as_matrix(), network_mat.T, run_parameters)
    spreadsheet_df = pd.DataFrame(sample_smooth, index=spreadsheet_df.index, columns=spreadsheet_df.columns)

    baseline_array = np.ones(network_mat.shape[0]) / network_mat.shape[0]
    baseline_array = kn.smooth_matrix_with_rwr(baseline_array, network_mat, run_parameters)[0]

    genes_list, drugs_list, consolodated_df = get_consolodated_dataframe(spreadsheet_df, drug_response_df)

    del spreadsheet_df
    del drug_response_df
    gc.collect()

    number_of_drugs = len(drugs_list)
    zipped_arguments = dstutil.zip_parameters(run_parameters, consolodated_df, genes_list, network_mat,
                                              spreadsheet_genes_as_input,
                                              baseline_array, drugs_list, range(0, number_of_drugs))
    dstutil.parallelize_processes_locally(worker_for_run_net_correlation, zipped_arguments, number_of_drugs)


def worker_for_run_net_correlation(run_parameters, consolodated_df, genes_list,
                                   network_mat, spreadsheet_genes_as_input, baseline_array, drugs_list, i):
    drug_response_df, spreadsheet_df = get_data_for_drug(consolodated_df, genes_list, drugs_list[i], run_parameters)
    sample_smooth = spreadsheet_df.as_matrix()

    pc_array = get_correlation(sample_smooth, drug_response_df.values[0], run_parameters)
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


def run_bootstrap_net_correlation(run_parameters):
    """ perform gene prioritization using bootstrap sampling and network smoothing

    Args:
        run_parameters: parameter set dictionary.
    """
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

    del network_df
    del network_mat_sparse
    del node_1_names
    del node_2_names
    del genes_lookup_table
    gc.collect()

    drug_response_df = kn.get_spreadsheet_df(run_parameters["drug_response_full_path"])
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])
    spreadsheet_genes_as_input = spreadsheet_df.index.values

    spreadsheet_df = kn.update_spreadsheet_df(spreadsheet_df, unique_gene_names)
    spreadsheet_df = zscore_dataframe(spreadsheet_df)
    sample_smooth, iterations = kn.smooth_matrix_with_rwr(spreadsheet_df.as_matrix(), network_mat.T, run_parameters)
    spreadsheet_df = pd.DataFrame(sample_smooth, index=spreadsheet_df.index, columns=spreadsheet_df.columns)

    baseline_array = np.ones(network_mat.shape[0]) / network_mat.shape[0]
    baseline_array = kn.smooth_matrix_with_rwr(baseline_array, network_mat, run_parameters)[0]

    genes_list, drugs_list, consolodated_df = get_consolodated_dataframe(spreadsheet_df, drug_response_df)

    del spreadsheet_df
    del drug_response_df
    gc.collect()

    number_of_drugs = len(drugs_list)
    zipped_arguments = dstutil.zip_parameters(run_parameters, consolodated_df, genes_list, network_mat,
                                              spreadsheet_genes_as_input,
                                              baseline_array, drugs_list, range(0, number_of_drugs))
    dstutil.parallelize_processes_locally(worker_for_run_bootstrap_net_correlation, zipped_arguments, number_of_drugs)


def worker_for_run_bootstrap_net_correlation(run_parameters, consolodated_df, genes_list, network_mat,
                                         spreadsheet_genes_as_input, baseline_array, drugs_list, i):
    """ worker for drug level parallelization.

    Args:
        run_parameters:
        consolodated_df:
        genes_list:
        unique_gene_names:
        network_mat:
        baseline_array:
        drugs_list:
        i:

    Returns:

    """
    restart_accumulator = np.zeros(network_mat.shape[0])
    gm_accumulator = np.ones(network_mat.shape[0])
    borda_count = np.zeros(network_mat.shape[0])

    drug_response_df, spreadsheet_df = get_data_for_drug(consolodated_df, genes_list, drugs_list[i], run_parameters)

    sample_smooth = spreadsheet_df.as_matrix()

    del consolodated_df
    gc.collect()

    pearson_array = get_correlation(sample_smooth, drug_response_df.values[0], run_parameters)
    n_bootstraps = run_parameters["number_of_bootstraps"]
    for bootstrap_number in range(0, n_bootstraps):
        sample_random, sample_permutation = sample_a_matrix_pearson(
            sample_smooth, run_parameters["rows_sampling_fraction"],
            run_parameters["cols_sampling_fraction"])

        drug_response = drug_response_df.values[0, None]
        drug_response = drug_response[0, sample_permutation]
        pc_array = get_correlation(sample_random, drug_response, run_parameters)
        pc_array[~np.in1d(spreadsheet_df.index, spreadsheet_genes_as_input)] = 0.0
        pc_array = np.abs(trim_to_top_beta(pc_array, run_parameters["top_beta_of_sort"]))
        restart_accumulator[pc_array != 0] += 1.0

        pc_array = pc_array / sum(pc_array)
        pc_array = kn.smooth_matrix_with_rwr(pc_array, network_mat, run_parameters)[0]
        pc_array = pc_array - baseline_array

        borda_count = sum_array_ranking_to_borda_count(borda_count, pc_array)
        gm_accumulator = (np.abs(pc_array) + EPSILON_0) * gm_accumulator

    restart_accumulator = restart_accumulator / n_bootstraps
    borda_count = borda_count / n_bootstraps
    pcc_gm_array = gm_accumulator ** (1 / n_bootstraps)

    generate_net_correlation_output(pearson_array, borda_count, pcc_gm_array, restart_accumulator,
                                    drug_response_df.index.values[0],
                                    spreadsheet_df.index, spreadsheet_genes_as_input, run_parameters)


def generate_net_correlation_output(pearson_array, pc_array, min_max_pc, restart_accumulator,
                                    drug_name, gene_name_list, gene_orig_list, run_parameters):
    """ Save final output of correlation with network

    Args:
        pearson_array: pearson correlation coefficient array
        pc_array: correlation score used to sort genes
        drug_name: name of the drug correlated
        gene_name_list: list of genes correlated (pc_array & all will be trimmed to these)
        gene_orig_list: original list of genes - size of pc_array
        run_parameters: with key 'results_directory'
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
    result_df = pd.DataFrame(output_val, columns=df_header).sort_values('quantitative sorting score', ascending=0)

    target_file_base_name = os.path.join(run_parameters["results_directory"],
                                         drug_name + '_' + run_parameters['out_filename'])

    file_name = kn.create_timestamped_filename(target_file_base_name) + '.tsv'
    result_df.to_csv(file_name, header=True, index=False, sep='\t')

    download_result_df = pd.DataFrame(data=None, index=None, columns=['Response', 'Gene ENSEMBL ID'])
    download_result_df['Response'] = result_df['Response']
    download_result_df['Gene ENSEMBL ID'] = result_df['Gene ENSEMBL ID']

    download_target_file_base_name = os.path.join(run_parameters["results_directory"],
                                                  drug_name + '_' + run_parameters['out_filename'])

    download_file_name = kn.create_timestamped_filename(download_target_file_base_name) + '_download' + '.tsv'
    download_result_df.to_csv(download_file_name, header=True, index=False, sep='\t')


def get_correlation(spreadsheet, drug_response, run_parameters):
    """ correlation function definition for all run methods

    Args:
        spreadsheet: genes x samples
        drug_response: one x samples
        run_parameters: with key 'correlation_measure'
        normalize: for lasso only
        max_iter: for lasso only

    Returns:
        correlation_array: genes x one
    """
    correlation_array = np.zeros(spreadsheet.shape[0])
    if 'correlation_measure' in run_parameters:
        if run_parameters['correlation_measure'] == 'pearson':
            spreadsheet = zscore(spreadsheet, axis=1, ddof=0)
            for row in range(0, spreadsheet.shape[0]):
                correlation_array[row] = pcc(spreadsheet[row, :], drug_response)[0]
            correlation_array[~(np.isfinite(correlation_array))] = 0
            return correlation_array

        if run_parameters['correlation_measure'] == 't_test':
            for row in range(0, spreadsheet.shape[0]):
                d = np.int_(drug_response)
                a = spreadsheet[row, d != 0]
                b = spreadsheet[row, d == 0]
                correlation_array[row] = np.abs(ttest_ind(a, b, axis=None, equal_var=False)[0])

            return correlation_array

    return correlation_array


def sum_array_ranking_to_borda_count(borda_count, corr_array):
    """ sum to borda count with a contigous array added to borda count
    Args:
        borda_count: the current borda count - same size as correlation array
        corr_array:  the correlation array to rank and add to the count
    Returns:
        borda_count: the ranking of the correlation array added to the input borda count
    """
    num_elem = borda_count.size

    # either assign (no duplicate case) or enumerate the correlation array
    if num_elem == (np.unique(corr_array)).size:
        borda_count[np.argsort(corr_array)] += np.int_(sorted(np.arange(0, corr_array.size) + 1))
        return borda_count

    # enumerate the borda vote
    borda_add = np.zeros(num_elem)
    enum_value = 1
    sort_order = np.argsort(corr_array)
    current_value = corr_array[sort_order[0]]
    for k in range(0, num_elem):
        if corr_array[sort_order[k]] != current_value:
            enum_value += 1
            current_value = corr_array[sort_order[k]]
        borda_add[sort_order[k]] = enum_value

    # scale to the number of elements in the array -- philosopical choice here --
    borda_add = borda_add + (num_elem - enum_value)

    return borda_count + borda_add


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


def sample_a_matrix_pearson(spreadsheet_mat, rows_fraction, cols_fraction):
    """ percent_sample x percent_sample random sample, from spreadsheet_mat.

    Args:
        spreadsheet_mat: gene x sample spread sheet as matrix.
        percent_sample: decimal fraction (slang-percent) - [0 : 1].

    Returns:
        sample_random: A specified precentage sample of the spread sheet.
        sample_permutation: the array that correponds to columns sample.
    """
    features_size = int(np.round(spreadsheet_mat.shape[0] * (1 - rows_fraction)))
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
    Beta = max(min(corr_arr.size, Beta) - 1, 0)
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


def write_phenotype_data_all(run_parameters, top_n=100):
    """ Post Processing: writes rows as genes, cols as drugs, data is gene in top n for the drug T or F.

    Args:
        run_parameters: with field: 'results_directory'
        top_n:          number of genes to rank (default=100)

    Returns: (writes consolodation file)
    """
    if 'top_beta_of_sort' in run_parameters:   top_n = run_parameters['top_beta_of_sort']
    dirList = sorted(os.listdir(run_parameters["results_directory"]))
    download_list = []
    for fileName in dirList:
        if (fileName[0:4] != 'all_') & (fileName[-12:] == 'download.tsv'):
            download_list.append(fileName)

    if len(download_list) == 0:
        return

    StartFileName = os.path.join(run_parameters["results_directory"], download_list[0])
    src_df = pd.read_csv(StartFileName, sep='\t', header=0, index_col=None)
    index_list = src_df['Gene ENSEMBL ID'].values

    all_phenotypes_df = pd.DataFrame(data=None, index=index_list)

    for fileName in download_list:
        tFileName = os.path.join(run_parameters["results_directory"], fileName)
        src_df = pd.read_csv(tFileName, sep='\t', header=0, index_col=None)
        all_phenotypes_df.insert(all_phenotypes_df.shape[1], src_df['Response'][1], [0] * all_phenotypes_df.shape[0],
                                 allow_duplicates=True)
        top_n_list = src_df['Gene ENSEMBL ID'][0:top_n]
        all_phenotypes_df[src_df['Response'][1]].loc[top_n_list] = 1

    download_target_file_base_name = os.path.join(run_parameters["results_directory"],
                                                  'all_phenotypes' + '_' + run_parameters['out_filename'])
    write_file_name = kn.create_timestamped_filename(download_target_file_base_name) + '_download.tsv'
    all_phenotypes_df.to_csv(write_file_name, header=True, index=True, sep='\t')
