"""
Created on Fri Sep 23 16:39:35 2016
@author: The KnowEnG dev team
"""

def correlation(run_parameters):
    """ gene prioritization """
    from gene_prioritization_toolbox import run_correlation
    run_correlation(run_parameters)

def net_correlation(run_parameters):
    """ gene prioritization """
    from gene_prioritization_toolbox import run_net_correlation
    run_net_correlation(run_parameters)

def bootstrap_correlation(run_parameters):
    """ gene prioritization """
    from gene_prioritization_toolbox import run_bootstrap_correlation
    run_bootstrap_correlation(run_parameters)

def bootstrap_net_correlation(run_parameters):
    """ gene prioritization """
    from gene_prioritization_toolbox import run_bootstrap_net_correlation
    run_bootstrap_net_correlation(run_parameters)

SELECT = {
    "correlation": correlation,
    "net_correlation": net_correlation,
    "bootstrap_correlation": bootstrap_correlation,
    "bootstrap_net_correlation": bootstrap_net_correlation}

def main():
    """
    This is the main function to perform gene prioritization.
    """
    import sys
    from knpackage.toolbox import get_run_directory_and_file
    from knpackage.toolbox import get_run_parameters

    run_directory, run_file = get_run_directory_and_file(sys.argv)
    run_parameters = get_run_parameters(run_directory, run_file)
    SELECT[run_parameters["method"]](run_parameters)

    #from gene_prioritization_toolbox import write_phenotype_data_all
    #write_phenotype_data_all(run_parameters)

if __name__ == "__main__":
    main()
