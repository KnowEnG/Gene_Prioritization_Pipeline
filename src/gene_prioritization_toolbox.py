"""
@author: The KnowEnG dev team
"""
import os
import numpy as np
import pandas as pd

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

    run_parameters["results_tmp_directory"] = kn.create_dir(run_parameters["results_directory"], 'tmp')

    results_tmp_directory      = run_parameters["results_tmp_directory"     ]
    phenotype_name_full_path   = run_parameters["phenotype_name_full_path"  ]
    spreadsheet_name_full_path = run_parameters["spreadsheet_name_full_path"]

    spreadsheet_df             = kn.get_spreadsheet_df(spreadsheet_name_full_path)
    phenotype_df               = kn.get_spreadsheet_df(phenotype_name_full_path  )
    phenotype_df               = phenotype_df.T

    number_of_jobs             = len(phenotype_df.index)
    jobs_id                    = range(0, number_of_jobs)
    zipped_arguments           = dstutil.zip_parameters( run_parameters
                                                       , spreadsheet_df
                                                       , phenotype_df
                                                       , jobs_id
                                                       )

    dstutil.parallelize_processes_locally( run_correlation_worker
                                         , zipped_arguments
                                         , number_of_jobs
                                         )

    write_phenotype_data_all(run_parameters       )
    kn.remove_dir           (results_tmp_directory)


def run_correlation_worker(run_parameters, spreadsheet_df, phenotype_df, job_id):
    """ core function for parallel run_correlation

    Args:
        run_parameters:  dict of parameters
        spreadsheet_df:  spreadsheet data frame
        phenotype_df:    phenotype data frame
        job_id:          parallel iteration number
    """

    np.random.seed(job_id) # selects a row from the phenotype_df spreadsheet

    phenotype_df           = phenotype_df.iloc[[job_id], :]
    spreadsheet_df,phenotype_df,msg= datacln.check_input_value_for_gene_prioritazion(spreadsheet_df, phenotype_df)
    pc_array               = get_correlation(spreadsheet_df.as_matrix(), phenotype_df.values[0], run_parameters)
    gene_name_list         = spreadsheet_df.index
    phenotype_name         = phenotype_df.index.values[0]

    generate_correlation_output(pc_array, phenotype_name, gene_name_list, run_parameters)


def generate_correlation_output(pc_array, phenotype_name, gene_name_list, run_parameters):
    """ Save final output of correlation

    Args:
        pc_array: pearson correlation coefficient array
        phenotype_name: name of the phenotype
        gene_name_list: list of the genes correlated (size of pc_array
        run_parameters: dictionary of run parameters with key 'results_directory'
    """

    phenotype_name_list = np.repeat(phenotype_name, len(gene_name_list))

    baseline_score      =     pc_array
    pc_array            = abs(pc_array)
    viz_score           =    (pc_array - min(pc_array)) / (max(pc_array) - min(pc_array))

    baseline_score      = np.round( baseline_score, 8)
    pc_array            = np.round( pc_array,       8)
    viz_score           = np.round( viz_score,      8)

    output_val          = np.column_stack( (phenotype_name_list, gene_name_list, pc_array, viz_score, baseline_score))
    df_header           = ['Response', 'Gene_ENSEMBL_ID', 'quantitative_sorting_score', 'visualization_score', 'baseline_score']
    result_df           = pd.DataFrame( output_val , columns=df_header)

    result_df           =       result_df.sort_values("visualization_score", ascending=0)
    result_df.index     = range(result_df.shape[0])

    write_one_phenotype(result_df, phenotype_name, gene_name_list, run_parameters)


def run_bootstrap_correlation(run_parameters):
    """ perform gene prioritization using bootstrap sampling

    Args:
        run_parameters: parameter set dictionary.
    """
    run_parameters["results_tmp_directory"] = kn.create_dir(run_parameters["results_directory"], 'tmp')

    phenotype_df = kn.get_spreadsheet_df(run_parameters["phenotype_name_full_path"])
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])
    phenotype_df = phenotype_df.T
    n_bootstraps = run_parameters["number_of_bootstraps"]

    number_of_jobs = len(phenotype_df.index)
    jobs_id = range(0, number_of_jobs)
    zipped_arguments = dstutil.zip_parameters(run_parameters, spreadsheet_df, phenotype_df, n_bootstraps, jobs_id)
    dstutil.parallelize_processes_locally(run_bootstrap_correlation_worker, zipped_arguments, number_of_jobs)

    write_phenotype_data_all(run_parameters)
    kn.remove_dir(run_parameters["results_tmp_directory"])


def run_bootstrap_correlation_worker(run_parameters, spreadsheet_df, phenotype_df, n_bootstraps, job_id):
    """  core function for parallel run_bootstrap_correlation

    Args:
        run_parameters:  dict of parameters
        spreadsheet_df:  spreadsheet data frame
        phenotype_df:    phenotype data frame
        n_bootstraps:    number of bootstrap samples to use
        job_id:          parallel iteration number
    """

    np.random.seed(job_id)

    phenotype_df = phenotype_df.iloc[[job_id], :]
    spreadsheet_df, phenotype_df, msg = datacln.check_input_value_for_gene_prioritazion(spreadsheet_df, phenotype_df)

    pearson_array = get_correlation(spreadsheet_df.as_matrix(), phenotype_df.values[0], run_parameters)
    borda_count = np.zeros(spreadsheet_df.shape[0])
    gm_accumulator = np.ones(spreadsheet_df.shape[0])
    for bootstrap_number in range(0, n_bootstraps):
        sample_random, sample_permutation = sample_a_matrix_pearson(
            spreadsheet_df.as_matrix(), 1.0, run_parameters["cols_sampling_fraction"])
        phenotype_response = phenotype_df.values[0, None]
        phenotype_response = phenotype_response[0, sample_permutation]
        pc_array = get_correlation(sample_random, phenotype_response, run_parameters)
        borda_count = sum_array_ranking_to_borda_count(borda_count, np.abs(pc_array))
        gm_accumulator = (np.abs(pc_array) + EPSILON_0) * gm_accumulator
    pcc_gm_array = gm_accumulator ** (1 / n_bootstraps)
    borda_count = borda_count / n_bootstraps

    phenotype_name = phenotype_df.index.values[0]
    gene_name_list = spreadsheet_df.index
    viz_score = (borda_count - min(borda_count)) / (max(borda_count) - min(borda_count))

    generate_bootstrap_correlation_output(borda_count, viz_score, pearson_array,
                                          phenotype_name, gene_name_list, run_parameters)


def generate_bootstrap_correlation_output(borda_count, viz_score, pearson_array, 
                                          phenotype_name, gene_name_list, run_parameters):
    """ Save final output of correlation

    Args:
        pearson_array: pearson correlation coefficient array
        phenotype_name: name of the phenotype
        gene_name_list: list of the genes correlated (size of pearson_array
        run_parameters: dictionary of run parameters with key 'results_directory'
    """
    phenotype_name_list = np.repeat(phenotype_name, len(gene_name_list))
    viz_score = np.round(viz_score, 8)
    borda_count = np.round(borda_count, 8)
    pearson_array = np.round(pearson_array, 8)

    output_val = np.column_stack(
        (phenotype_name_list, gene_name_list, borda_count, viz_score, pearson_array))

    df_header = ['Response', 'Gene_ENSEMBL_ID', 'quantitative_sorting_score', 'visualization_score', 'baseline_score']
    result_df = pd.DataFrame(output_val, columns=df_header).sort_values("visualization_score", ascending=0)
    result_df.index = range(result_df.shape[0])

    write_one_phenotype(result_df, phenotype_name, gene_name_list, run_parameters)


def run_net_correlation(run_parameters):
    """ perform gene prioritization with network smoothing

    Args:
        run_parameters: parameter set dictionary.
    """
    run_parameters["results_tmp_directory"] = kn.create_dir(run_parameters["results_directory"], 'tmp')
    gg_network_name_full_path = run_parameters['gg_network_name_full_path']
    network_mat, unique_gene_names = kn.get_sparse_network_matrix(gg_network_name_full_path)

    network_mat = normalize(network_mat, norm="l1", axis=0)

    phenotype_df = kn.get_spreadsheet_df(run_parameters["phenotype_name_full_path"])
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])
    spreadsheet_genes_as_input = spreadsheet_df.index.values
    phenotype_df = phenotype_df.T

    spreadsheet_df = kn.update_spreadsheet_df(spreadsheet_df, unique_gene_names)
    spreadsheet_df = zscore_dataframe(spreadsheet_df)

    sample_smooth, iterations = kn.smooth_matrix_with_rwr(spreadsheet_df.as_matrix(), network_mat.T, run_parameters)
    spreadsheet_df = pd.DataFrame(sample_smooth, index=spreadsheet_df.index, columns=spreadsheet_df.columns)

    baseline_array = np.ones(network_mat.shape[0]) / network_mat.shape[0]
    baseline_array = kn.smooth_matrix_with_rwr(baseline_array, network_mat, run_parameters)[0]

    number_of_jobs = len(phenotype_df.index)
    jobs_id = range(0, number_of_jobs)
    zipped_arguments = dstutil.zip_parameters(run_parameters, spreadsheet_df, phenotype_df, network_mat,
                                              spreadsheet_genes_as_input, baseline_array, jobs_id)
    dstutil.parallelize_processes_locally(run_net_correlation_worker, zipped_arguments, number_of_jobs)

    write_phenotype_data_all(run_parameters)
    kn.remove_dir(run_parameters["results_tmp_directory"])


def run_net_correlation_worker(run_parameters, spreadsheet_df, phenotype_df, network_mat,
                                   spreadsheet_genes_as_input, baseline_array, job_id):
    """  core function for parallel run_net_correlation

    Args:
        run_parameters:  dict of parameters
        spreadsheet_df:  spreadsheet data frame
        phenotype_df:    phenotype data frame
        network_mat:     adjacency matrix
        spreadsheet_genes_as_input: list of genes 
        baseline_array:  network smooted baseline array
        job_id:          parallel iteration number
    """

    np.random.seed(job_id)

    phenotype_df = phenotype_df.iloc[[job_id], :]
    spreadsheet_df, phenotype_df, msg = datacln.check_input_value_for_gene_prioritazion(spreadsheet_df, phenotype_df)

    sample_smooth = spreadsheet_df.as_matrix()

    pc_array = get_correlation(sample_smooth, phenotype_df.values[0], run_parameters)
    pearson_array = pc_array.copy()
    pc_array[~np.in1d(spreadsheet_df.index, spreadsheet_genes_as_input)] = 0.0
    pc_array = np.abs(trim_to_top_beta(pc_array, run_parameters["top_beta_of_sort"]))
    restart_accumulator = pc_array.copy()
    restart_accumulator[restart_accumulator != 0] = 1

    pc_array = pc_array / max(sum(pc_array), EPSILON_0)
    pc_array = kn.smooth_matrix_with_rwr(pc_array, network_mat, run_parameters)[0]

    pc_array = pc_array - baseline_array
    quantitative_score = pc_array
    viz_score = (pc_array - min(pc_array)) / (max(pc_array) - min(pc_array))

    phenotype_name = phenotype_df.index.values[0]
    gene_name_list = spreadsheet_df.index
    gene_orig_list = spreadsheet_genes_as_input

    generate_net_correlation_output(pearson_array, quantitative_score, viz_score, restart_accumulator,
                                    phenotype_name, gene_name_list, gene_orig_list, run_parameters)


def run_bootstrap_net_correlation(run_parameters):
    """ perform gene prioritization using bootstrap sampling and network smoothing

    Args:
        run_parameters: parameter set dictionary.
    """
    run_parameters["results_tmp_directory"] = kn.create_dir(run_parameters["results_directory"], 'tmp')
    gg_network_name_full_path = run_parameters['gg_network_name_full_path']
    network_mat, unique_gene_names = kn.get_sparse_network_matrix(gg_network_name_full_path)

    network_mat = normalize(network_mat, norm="l1", axis=0)

    phenotype_df = kn.get_spreadsheet_df(run_parameters["phenotype_name_full_path"])
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])
    spreadsheet_genes_as_input = spreadsheet_df.index.values
    phenotype_df = phenotype_df.T

    spreadsheet_df = kn.update_spreadsheet_df(spreadsheet_df, unique_gene_names)
    spreadsheet_df = zscore_dataframe(spreadsheet_df)
    sample_smooth, iterations = kn.smooth_matrix_with_rwr(spreadsheet_df.as_matrix(), network_mat.T, run_parameters)
    spreadsheet_df = pd.DataFrame(sample_smooth, index=spreadsheet_df.index, columns=spreadsheet_df.columns)

    baseline_array = np.ones(network_mat.shape[0]) / network_mat.shape[0]
    baseline_array = kn.smooth_matrix_with_rwr(baseline_array, network_mat, run_parameters)[0]

    number_of_jobs = len(phenotype_df.index)
    jobs_id = range(0, number_of_jobs)
    zipped_arguments = dstutil.zip_parameters(run_parameters, spreadsheet_df, phenotype_df, network_mat,
                                              spreadsheet_genes_as_input, baseline_array, jobs_id)
    dstutil.parallelize_processes_locally(run_bootstrap_net_correlation_worker, zipped_arguments, number_of_jobs)

    write_phenotype_data_all(run_parameters)
    kn.remove_dir(run_parameters["results_tmp_directory"])


def run_bootstrap_net_correlation_worker(run_parameters, spreadsheet_df, phenotype_df, network_mat,
                                             spreadsheet_genes_as_input, baseline_array, job_id):
    """ worker for bootstrap network parallelization.

    Args:
        run_parameters:  dict of parameters
        spreadsheet_df:  spreadsheet data frame
        phenotype_df:    phenotype data frame
        network_mat:     adjacency matrix
        spreadsheet_genes_as_input: list of genes 
        baseline_array:  network smooted baseline array
        job_id:          parallel iteration number
    """

    np.random.seed(job_id)

    restart_accumulator = np.zeros(network_mat.shape[0])
    gm_accumulator = np.ones(network_mat.shape[0])
    borda_count = np.zeros(network_mat.shape[0])

    phenotype_df = phenotype_df.iloc[[job_id], :]
    spreadsheet_df, phenotype_df, msg = datacln.check_input_value_for_gene_prioritazion(spreadsheet_df, phenotype_df)

    sample_smooth = spreadsheet_df.as_matrix()

    pearson_array = get_correlation(sample_smooth, phenotype_df.values[0], run_parameters)
    n_bootstraps = run_parameters["number_of_bootstraps"]
    for bootstrap_number in range(0, n_bootstraps):
        sample_random, sample_permutation = sample_a_matrix_pearson(
            sample_smooth, 1.0, run_parameters["cols_sampling_fraction"])

        phenotype_response = phenotype_df.values[0, None]
        phenotype_response = phenotype_response[0, sample_permutation]
        pc_array = get_correlation(sample_random, phenotype_response, run_parameters)

        pc_array[~np.in1d(spreadsheet_df.index, spreadsheet_genes_as_input)] = 0.0
        pc_array = np.abs(trim_to_top_beta(pc_array, run_parameters["top_beta_of_sort"]))
        restart_accumulator[pc_array != 0] += 1.0

        pc_array = pc_array / max(sum(pc_array), EPSILON_0)
        pc_array = kn.smooth_matrix_with_rwr(pc_array, network_mat, run_parameters)[0]
        pc_array = pc_array - baseline_array

        borda_count = sum_array_ranking_to_borda_count(borda_count, pc_array)
        gm_accumulator = (np.abs(pc_array) + EPSILON_0) * gm_accumulator

    restart_accumulator = restart_accumulator / n_bootstraps
    borda_count = borda_count / n_bootstraps
    # pcc_gm_array = gm_accumulator ** (1 / n_bootstraps)
    viz_score = (borda_count - min(borda_count)) / (max(borda_count) - min(borda_count))

    phenotype_name = phenotype_df.index.values[0]
    gene_name_list = spreadsheet_df.index
    gene_orig_list = spreadsheet_genes_as_input
    quantitative_score = borda_count
    generate_net_correlation_output(pearson_array, quantitative_score, viz_score, restart_accumulator,
                                    phenotype_name, gene_name_list, gene_orig_list, run_parameters)


def generate_net_correlation_output(pearson_array, quantitative_score, viz_score, restart_accumulator,
                                    phenotype_name, gene_name_list, gene_orig_list, run_parameters):
    """ Save final output of correlation with network

    Args:
        pearson_array: pearson correlation coefficient array
        quantitative_score: correlation score used to sort genes
        phenotype_name: name of the phenotype correlated
        gene_name_list: list of genes correlated (quantitative_score & all will be trimmed to these)
        gene_orig_list: original list of genes - size of quantitative_score
        run_parameters: with key 'results_directory'
    """
    mask = np.in1d(gene_name_list, gene_orig_list)
    quantitative_score = quantitative_score[mask]
    pearson_array = pearson_array[mask]
    viz_score = viz_score[mask]
    gene_name_list = gene_name_list[mask]
    restart_accumulator = restart_accumulator[mask]
    phenotype_name_list = np.repeat(phenotype_name, len(gene_name_list))

    quantitative_score = np.round(quantitative_score, 8)
    viz_score = np.round(viz_score, 8)
    pearson_array = np.round(pearson_array, 8)

    output_val = np.column_stack(
        (phenotype_name_list, gene_name_list, quantitative_score, viz_score, pearson_array, restart_accumulator))
    
    df_header = ['Response', 'Gene_ENSEMBL_ID', 'quantitative_sorting_score', 'visualization_score',
                 'baseline_score', 'Percent_appearing_in_restart_set']
    result_df = pd.DataFrame(output_val, columns=df_header).sort_values('visualization_score', ascending=0)

    write_one_phenotype(result_df, phenotype_name, gene_orig_list, run_parameters)



def get_correlation(spreadsheet_mat, phenotype_response, run_parameters):
    """ correlation function definition for all run methods

    Args:
        spreadsheet_mat: genes x samples
        phenotype_response: one x samples
        run_parameters: with key 'correlation_measure'

    Returns:
        correlation_array: genes x one
    """
    correlation_array = np.zeros(spreadsheet_mat.shape[0])
    if 'correlation_measure' in run_parameters:
        if run_parameters['correlation_measure'] == 'pearson':

            spreadsheet_mat = spreadsheet_mat - spreadsheet_mat.mean(axis=1).reshape((-1, 1))
            phenotype_response = phenotype_response - phenotype_response.mean()
            spreadsheet_mat_var = np.std(spreadsheet_mat, axis=1)
            phenotype_response_var = np.std(phenotype_response)
            numerator = spreadsheet_mat.dot(phenotype_response)
            denominator = spreadsheet_mat_var * phenotype_response_var * spreadsheet_mat.shape[1]
            with np.errstate(divide='ignore', invalid='ignore'):
                correlation_array = np.true_divide(numerator, denominator)
                correlation_array[denominator==0] = 0

            return correlation_array

        if run_parameters['correlation_measure'] == 't_test':
        
            a = spreadsheet_mat[:, phenotype_response!=0]
            b = spreadsheet_mat[:, phenotype_response==0]
            d = np.mean(a, axis=1) - np.mean(b, axis=1)
            denom = np.sqrt(np.var(a, axis=1, ddof=1)/a.shape[1] + np.var(b, axis=1, ddof=1)/b.shape[1])
            with np.errstate(divide='ignore', invalid='ignore'):
                correlation_array = np.divide(d, denom)
                correlation_array[np.isnan(denom)] = 0

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
    abs_corr_arr = np.abs(corr_arr)
    abs_corr_arr_cutoff_value = sorted(abs_corr_arr)[::-1][Beta]
    corr_arr[abs_corr_arr < abs_corr_arr_cutoff_value] = 0
    return corr_arr


def zscore_dataframe(genes_by_sample_df):
    """ zscore by rows for genes x samples dataframe

    Args:
        genes_by_sample_df: zscore along rows for genes x phenotypes dataframe

    Returns:
        spreadsheet_df: rows add up to zero, normalized to the mean and std deveiation
    """
    zscore_df = (genes_by_sample_df.sub(genes_by_sample_df.mean(axis=1), axis=0)).truediv(
                    np.maximum(genes_by_sample_df.std(axis=1), 1e-12), axis=0)
    return zscore_df


def write_one_phenotype(result_df, phenotype_name, gene_name_list, run_parameters):
    """ write the phenotype output file to the results directory and the temporary directory files

    Args:
        result_df:
        phenotype_name:
        gene_name_list:
        run_parameters:
        
    Output:
        {phenotype}_{method}_{correlation_measure}_{timestamp}_viz.tsv
    """
    result_df.to_csv(get_output_file_name(run_parameters, 'results_directory', phenotype_name, 'viz'), header=True, index=False, sep='\t')

    download_result_df = pd.DataFrame(data=None, index=None, columns=[phenotype_name])
    download_result_df[phenotype_name] = result_df['Gene_ENSEMBL_ID']
    download_result_df.to_csv(
        get_output_file_name(run_parameters, 'results_tmp_directory', phenotype_name, 'download'), header=True, index=False, sep='\t')

    top_genes = download_result_df.values[: run_parameters['top_beta_of_sort']]
    update_orig_result_df = pd.DataFrame(np.in1d(gene_name_list, top_genes).astype(int), index=gene_name_list, columns=[phenotype_name])
    update_orig_result_df.to_csv(
        get_output_file_name(run_parameters, 'results_tmp_directory', phenotype_name, 'original'), header=True, index=True, sep='\t')


def write_phenotype_data_all(run_parameters):
    """ Post Processing: writes rows as genes, cols as phenotypes, data is gene in top n for the phenotype T or F.

    Args:
        run_parameters: with field: 'results_directory'
        
    Output:
        ranked_genes_per_phenotype_{method}_{correlation_measure}_{timestamp}_download.tsv
        top_genes_per_phenotype_{method}_{correlation_measure}_{timestamp}_download.tsv
    """  
    tmp_dir = run_parameters["results_tmp_directory"]
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
        phenotype_name = src_df.columns.values[0]
        all_phenotypes_download_df[phenotype_name] = src_df[phenotype_name]

    for fileName in original_list:
        tFileName = os.path.join(tmp_dir, fileName)
        src_df = pd.read_csv(tFileName, sep='\t', header=0, index_col=0)          
        phenotype_name = src_df.columns.values[0]
        all_phenotypes_original_df[phenotype_name] = src_df[phenotype_name]

    all_phenotypes_download_df.index = range(1, all_phenotypes_download_df.shape[0]+1)
    all_phenotypes_download_df.to_csv(
        get_output_file_name(run_parameters, 'results_directory', 'ranked_genes_per_phenotype', 'download'), header=True, index=True, sep='\t')
    all_phenotypes_original_df.to_csv(
        get_output_file_name(run_parameters, 'results_directory', 'top_genes_per_phenotype', 'download'), header=True, index=True, sep='\t')


def get_output_file_name(run_parameters, dir_name_key, prefix_string, suffix_string='', type_suffix='tsv'):
    """ get the full directory / filename for writing
    Args:
        run_parameters: dictionary with keys: dir_name_key, "method" and "correlation_measure"
        dir_name_key:   run_parameters dictionary key for the output directory
        prefix_string:  the first letters of the ouput file name
        suffix_string:  the last letters of the output file name before type_suffix
        type_suffix:    the file type extenstion (default 'tsv') without period character

    Returns:
        output_file_name:   full file and directory name suitable for file writing
    """
    output_file_name = os.path.join(run_parameters[dir_name_key], prefix_string + '_' +
                                    run_parameters['method'] + '_' + run_parameters["correlation_measure"])

    output_file_name = kn.create_timestamped_filename(output_file_name) + '_' + suffix_string + '.' + type_suffix
    return output_file_name
