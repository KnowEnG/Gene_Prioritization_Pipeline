import time
import pandas as pd
import numpy as np

#import knpackage.toolbox as kn
import gene_prioritization_toolbox as gptbx


def benchmark_test_run_gene_prioritization(run_parameters):
    """

    Args:
         run_parameters: dict of run parameters with keys:
            benchmark_results_full_path, ...

    Returns:
        result_df: dataframe of correlation coefficients for each gene with each drug
    """

    start_time = 'hack_time_' + time.strftime("%a_%d_%b_%Y_%H_%M_%S", time.localtime())
    print('{}\n   using result test: {}\n'.format(start_time, run_parameters["benchmark_results_full_path"]))

    result_df = gptbx.run_gene_correlation(run_parameters)

    return result_df


def trim_spreadsheet_and_drug_df(spreadsheet_df, drug_df):
    """ remove the columns of spreadsheet_df where drug_df has NA's and remove corresponding rows of the drug_df
    """
    good_cols_list = get_no_NAs_cols(drug_df)
    spreadsheet_df = spreadsheet_df[good_cols_list]
    drug_df = drug_df.loc[good_cols_list]

    return spreadsheet_df, drug_df


def get_no_NAs_cols(one_col_df):
    """
    Args:
        one_col_df: a one column dataframe with rows indices corresponding to spreadsheet_df columns
    Returns:
        good_cols_list: list of good columns suitable for indexing
    """
    col = one_col_df.columns[0]
    good_cols = one_col_df[col][one_col_df[col] != 'NA'].index
    good_cols_list = good_cols.values

    return good_cols_list