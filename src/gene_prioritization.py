"""
Created on Fri Sep 23 16:39:35 2016
@author: The Gene Prioritization dev team
"""

import time

def correlation(run_parameters):
    """ gene prioritization """
    from gene_prioritization_toolbox import run_gene_prioritization
    run_gene_prioritization(run_parameters)

def net_correlate(run_parameters):
    """ gene prioritization """
    from gene_prioritization_toolbox import run_net_gene_prioritization
    run_net_gene_prioritization(run_parameters)

def bootstrap_correlate(run_parameters):
    """ gene prioritization """
    from gene_prioritization_toolbox import run_bootstrap_correlate
    run_bootstrap_correlate(run_parameters)

def bootstrap_net_correlate(run_parameters):
    """ gene prioritization """
    from gene_prioritization_toolbox import run_bootstrap_net_correlate
    run_bootstrap_net_correlate(run_parameters)

SELECT = {
    "correlate": correlation,
    "net_correlate": net_correlate,
    "bootstrap_correlate": bootstrap_correlate,
    "bootstrap_net_correlate": bootstrap_net_correlate}

def main():
    """
    This the main function to perform gene prioritization.
    """
    #f_time = gene_prioritization_pcc("../data/gene_expression.csv", "../data/drug_response.csv")
    import sys
    from knpackage.toolbox import get_run_directory_and_file
    from knpackage.toolbox import get_run_parameters
    run_directory, run_file = get_run_directory_and_file(sys.argv)
    run_parameters = get_run_parameters(run_directory, run_file)
    SELECT[run_parameters["method"]](run_parameters)

if __name__ == "__main__":
    main()

