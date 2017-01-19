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


def run_correlation(run_parameters):
    """ perform gene prioritization

    Args:
        run_parameters: parameter set dictionary.
    """
    phenotype_df = kn.get_spreadsheet_df(run_parameters["drug_response_full_path"])
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])

    number_of_jobs = len(phenotype_df.index)
    jobs_id = range(0, number_of_jobs)
    zipped_arguments = dstutil.zip_parameters(run_parameters, spreadsheet_df, phenotype_df, jobs_id)
    dstutil.parallelize_processes_locally(worker_for_run_correlation, zipped_arguments, number_of_jobs)


def worker_for_run_correlation(run_parameters, spreadsheet_df, phenotype_df, job_id):
    """ core function for parallel run_correlation

    Args:
        run_parameters:  dict of parameters
        spreadsheet_df:  spreadsheet data frame
        phenotype_df:    drug data frame
        job_id:          parallel iteration number
    """
    # selects the ith row in phenotype_df
    phenotype_df = phenotype_df.iloc[[job_id], :]

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

    df_header = ['Response', 'Gene_ENSEMBL_ID', 'quantitative_sorting_score', 'visualization_score', 'baseline_score']
    result_df = pd.DataFrame(output_val, columns=df_header).sort_values("visualization score", ascending=0)
    result_df.index = range(result_df.shape[0])

    result_df.to_csv(get_output_file_name(run_parameters, drug_name), header=True, index=False, sep='\t')

    download_result_df = pd.DataFrame(data=None, index=None, columns=['Gene_ENSEMBL_ID'])
    # download_result_df['Response'] = result_df['Response']
    download_result_df['Gene_ENSEMBL_ID'] = result_df['Gene_ENSEMBL_ID']

    download_result_df.to_csv(get_output_file_name(run_parameters, drug_name, 'download'),
                              header=True, index=False, sep='\t')


def run_bootstrap_correlation(run_parameters):
    """ perform gene prioritization using bootstrap sampling

    Args:
        run_parameters: parameter set dictionary.
    """
    drug_response_df = kn.get_spreadsheet_df(run_parameters["drug_response_full_path"])
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])

    n_bootstraps = run_parameters["number_of_bootstraps"]

    number_of_jobs = len(drug_response_df.index)
    jobs_id = range(0, number_of_jobs)
    zipped_arguments = dstutil.zip_parameters(run_parameters, spreadsheet_df, drug_response_df, n_bootstraps, jobs_id)
    dstutil.parallelize_processes_locally(worker_for_run_bootstrap_correlation, zipped_arguments, number_of_jobs)


def worker_for_run_bootstrap_correlation(run_parameters, spreadsheet_df, phenotype_df, n_bootstraps, job_id):
    """  core function for parallel run_bootstrap_correlation

    Args:
        run_parameters:  dict of parameters
        spreadsheet_df:  spreadsheet data frame
        phenotype_df:    phenotype data frame
        n_bootstraps:    number of bootstrap samples to use
        job_id:          parallel iteration number
    """
    phenotype_df = phenotype_df.iloc[[job_id], :]
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

    df_header = ['Response', 'Gene_ENSEMBL_ID', 'quantitative_sorting_score', 'visualization_score', 'baseline_score']
    result_df = pd.DataFrame(output_val, columns=df_header).sort_values("quantitative sorting score", ascending=0)
    result_df.index = range(result_df.shape[0])

    result_df.to_csv(get_output_file_name(run_parameters, drug_name), header=True, index=False, sep='\t')

    download_result_df = pd.DataFrame(data=None, index=None, columns=['Gene_ENSEMBL_ID'])
    # download_result_df['Response'] = result_df['Response']
    download_result_df['Gene_ENSEMBL_ID'] = result_df['Gene_ENSEMBL_ID']

    download_result_df.to_csv(get_output_file_name(run_parameters, drug_name, 'download'), header=True, index=False, sep='\t')


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

    number_of_jobs = len(drug_response_df.index)
    jobs_id = range(0, number_of_jobs)
    zipped_arguments = dstutil.zip_parameters(run_parameters, spreadsheet_df, drug_response_df, network_mat,
                                              spreadsheet_genes_as_input, baseline_array, jobs_id)
    dstutil.parallelize_processes_locally(worker_for_run_net_correlation, zipped_arguments, number_of_jobs)


def worker_for_run_net_correlation(run_parameters, spreadsheet_df, phenotype_df, network_mat,
                                   spreadsheet_genes_as_input, baseline_array, job_id):

    phenotype_df = phenotype_df.iloc[[job_id], :]
    spreadsheet_df_trimmed, phenotype_df_trimmed, ret_msg = datacln.check_input_value_for_gene_prioritazion(
        spreadsheet_df, phenotype_df)

    sample_smooth = spreadsheet_df_trimmed.as_matrix()

    pc_array = get_correlation(sample_smooth, phenotype_df_trimmed.values[0], run_parameters)
    pearson_array = pc_array.copy()
    pc_array[~np.in1d(spreadsheet_df_trimmed.index, spreadsheet_genes_as_input)] = 0.0
    pc_array = np.abs(trim_to_top_beta(pc_array, run_parameters["top_beta_of_sort"]))
    restart_accumulator = pc_array.copy()
    restart_accumulator[restart_accumulator != 0] = 1

    pc_array = pc_array / sum(pc_array)
    pc_array = kn.smooth_matrix_with_rwr(pc_array, network_mat, run_parameters)[0]

    pc_array = pc_array - baseline_array
    min_max_pc = (pc_array - min(pc_array)) / (max(pc_array) - min(pc_array))

    generate_net_correlation_output(
        pearson_array, pc_array, min_max_pc, restart_accumulator, phenotype_df_trimmed.index.values[0],
        spreadsheet_df_trimmed.index, spreadsheet_genes_as_input, run_parameters)


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

    number_of_jobs = len(drug_response_df.index)
    jobs_id = range(0, number_of_jobs)
    zipped_arguments = dstutil.zip_parameters(run_parameters, spreadsheet_df, drug_response_df, network_mat,
                                              spreadsheet_genes_as_input, baseline_array, jobs_id)
    dstutil.parallelize_processes_locally(worker_for_run_bootstrap_net_correlation, zipped_arguments, number_of_jobs)


def worker_for_run_bootstrap_net_correlation(run_parameters, spreadsheet_df, phenotype_df, network_mat,
                                             spreadsheet_genes_as_input, baseline_array, job_id):
    """ worker for drug level parallelization.

    Args:
        run_parameters:
        spreadsheet_df:
        phenotype_df:
        network_mat:
        baseline_array:
        job_id:

    Returns:

    """
    restart_accumulator = np.zeros(network_mat.shape[0])
    gm_accumulator = np.ones(network_mat.shape[0])
    borda_count = np.zeros(network_mat.shape[0])

    phenotype_df = phenotype_df.iloc[[job_id], :]
    spreadsheet_df_trimmed, phenotype_df_trimmed, ret_msg = datacln.check_input_value_for_gene_prioritazion(
        spreadsheet_df, phenotype_df)

    sample_smooth = spreadsheet_df_trimmed.as_matrix()

    pearson_array = get_correlation(sample_smooth, phenotype_df_trimmed.values[0], run_parameters)
    n_bootstraps = run_parameters["number_of_bootstraps"]
    for bootstrap_number in range(0, n_bootstraps):
        sample_random, sample_permutation = sample_a_matrix_pearson(
            sample_smooth, run_parameters["rows_sampling_fraction"],
            run_parameters["cols_sampling_fraction"])

        drug_response = phenotype_df_trimmed.values[0, None]
        drug_response = drug_response[0, sample_permutation]
        pc_array = get_correlation(sample_random, drug_response, run_parameters)
        pc_array[~np.in1d(spreadsheet_df_trimmed.index, spreadsheet_genes_as_input)] = 0.0
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
                                    phenotype_df_trimmed.index.values[0],
                                    spreadsheet_df_trimmed.index, spreadsheet_genes_as_input, run_parameters)


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

    df_header = ['Response', 'Gene_ENSEMBL_ID', 'quantitative_sorting_score', 'visualization_score',
                 'baseline_score', 'Percent_appearing_in_restart_set']
    result_df = pd.DataFrame(output_val, columns=df_header).sort_values('quantitative_sorting_score', ascending=0)

    result_df.to_csv(get_output_file_name(run_parameters, drug_name), header=True, index=False, sep='\t')

    curr_dir = run_parameters["results_directory"]
    new_dir = os.path.join(run_parameters["results_directory"], 'tmp')
    if not os.path.isdir(new_dir):
        os.mkdir(new_dir)
    run_parameters["results_directory"] = new_dir

    download_result_df = pd.DataFrame(data=None, index=None, columns=[drug_name])
    download_result_df[drug_name] = result_df['Gene_ENSEMBL_ID']
    download_result_df.to_csv(get_output_file_name(run_parameters, drug_name, 'download'), header=True, index=False, sep='\t')
    
    top_genes = download_result_df.values[: run_parameters['top_beta_of_sort']]
    update_orig_result_df = pd.DataFrame(np.in1d(gene_orig_list, top_genes).astype(int), index=gene_orig_list, columns=[drug_name])
    update_orig_result_df.to_csv(get_output_file_name(run_parameters, drug_name, 'original'), header=True, index=True, sep='\t')

    run_parameters["results_directory"] = curr_dir

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
            #linear relationship between two (normally distributed) datasets. same slope = 1, cross = 0, opposite = -1
            spreadsheet = zscore(spreadsheet, axis=1, ddof=0)
            for row in range(0, spreadsheet.shape[0]):
                correlation_array[row] = pcc(spreadsheet[row, :], drug_response)[0]
            correlation_array[~(np.isfinite(correlation_array))] = 0
            return correlation_array

        if run_parameters['correlation_measure'] == 't_test':
            #ttest_ind: two sided test for null hypothesis that 2 independent samples have identical (expected) vlaues
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
    tmp_dir = os.path.join(run_parameters["results_directory"], 'tmp')
    dirList = sorted(os.listdir(tmp_dir))
    download_list = []
    original_list = []
    for fileName in dirList:
        if (fileName[-12:] == 'download.tsv'):
            download_list.append(fileName)
        if (fileName[-12:] == 'original.tsv'):
            original_list.append(fileName)

    if (len(download_list) == 0 or len(original_list) == 0):
        return

    StartFileName = os.path.join(tmp_dir, original_list[0])
    src_df = pd.read_csv(StartFileName, sep='\t', header=0, index_col=0)
    index_list = src_df.index.values

    all_phenotypes_download_df = pd.DataFrame(data=None, index=None)
    all_phenotypes_original_df = pd.DataFrame(data=None, index=index_list)


    for fileName in download_list:
        tFileName = os.path.join(tmp_dir, fileName)
        src_df = pd.read_csv(tFileName, sep='\t', header=0, index_col=None)
        drug_name = src_df.columns.values[0]
        all_phenotypes_download_df[drug_name] = src_df[drug_name]

    for fileName in original_list:
        tFileName = os.path.join(tmp_dir, fileName)
        src_df = pd.read_csv(tFileName, sep='\t', header=0, index_col=0)          
        drug_name = src_df.columns.values[0]
        all_phenotypes_original_df[drug_name] = src_df[drug_name]


    all_phenotypes_download_df.to_csv(get_output_file_name(run_parameters, 'all_phenotypes', 'download'), header=True, index=True, sep='\t')
    all_phenotypes_original_df.to_csv(get_output_file_name(run_parameters, 'all_phenotypes', 'original'), header=True, index=True, sep='\t')
    kn.remove_dir(tmp_dir)

# def write_phenotype_data_all(run_parameters, top_n=100):
#     """ Post Processing: writes rows as genes, cols as drugs, data is gene in top n for the drug T or F.

#     Args:
#         run_parameters: with field: 'results_directory'
#         top_n:          number of genes to rank (default=100)

#     Returns: (writes consolodation file)
#     """
#     if 'top_beta_of_sort' in run_parameters:   top_n = run_parameters['top_beta_of_sort']
#     dirList = sorted(os.listdir(run_parameters["results_directory"]))
#     download_list = []
#     for fileName in dirList:
#         if (fileName[0:4] != 'all_') & (fileName[-12:] == 'download.tsv'):
#             download_list.append(fileName)

#     if len(download_list) == 0:
#         return

#     StartFileName = os.path.join(run_parameters["results_directory"], download_list[0])
#     src_df = pd.read_csv(StartFileName, sep='\t', header=0, index_col=None)
#     index_list = src_df['Gene_ENSEMBL_ID'].values

#     all_phenotypes_df = pd.DataFrame(data=None, index=index_list)

#     for fileName in download_list:
#         tFileName = os.path.join(run_parameters["results_directory"], fileName)
#         src_df = pd.read_csv(tFileName, sep='\t', header=0, index_col=None)
#         all_phenotypes_df.insert(all_phenotypes_df.shape[1], src_df['Response'][1], [0] * all_phenotypes_df.shape[0],
#                                  allow_duplicates=True)
#         top_n_list = src_df['Gene_ENSEMBL_ID'][0:top_n]
#         all_phenotypes_df[src_df['Response'][1]].loc[top_n_list] = 1

#     all_phenotypes_df.to_csv(get_output_file_name(run_parameters, 'all_phenotypes', 'download'), header=True, index=True, sep='\t')


def get_output_file_name(run_parameters, prefix_string, suffix_string='', type_suffix='tsv'):
    """ get the full directory / filename for writing
    Args:
        run_parameters: dictionary with keys: "results_directory", "method" and "correlation_measure"
        prefix_string:  the first letters of the ouput file name
        suffix_string:  the last letters of the output file name before '.tsv'

    Returns:
        output_file_name:   full file and directory name suitable for file writing
    """
    output_file_name = os.path.join(run_parameters["results_directory"], prefix_string + '_' +
                                    run_parameters['method'] + '_' + run_parameters["correlation_measure"])

    output_file_name = kn.create_timestamped_filename(output_file_name) + '_' + suffix_string + '.' + type_suffix
    return output_file_name
