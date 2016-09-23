"""
Created on Fri Sep 23 16:39:35 2016
@author: The Gene Prioritization dev team
"""

import time
def gene_prioritization_pcc(spreadsheet_df_full_path, drug_response_full_path):
    """Pearson correlation coefficient(PCC) gene prioritization

    Args:
        spreadsheet_df_full_path: spreadsheet_df path and file name.
        drug_response_full_path: drug_response path and file name.
    """
    import gene_prioritization_toolbox as tl
    t0 = time.time()
    tl.run_gene_prioritization(spreadsheet_df_full_path, drug_response_full_path)
    return time.time() - t0


def main():
    """
    This the main function to perform gene prioritization.
    """
    f_time = gene_prioritization_pcc("../data/gene_expression.csv", "../data/drug_response.csv")

if __name__ == "__main__":
    main()