"""
lanier4@illinois.edu

"""
import os
import filecmp
import time

GP_options_dict = {
                'run_pearson'                 : 'BENCHMARK_1_GP_pearson',
                'run_net_pearson'             : 'BENCHMARK_3_GP_net_pearson',
                'run_bootstrap_pearson'       : 'BENCHMARK_2_GP_bootstrap_pearson',
                'run_bootstrap_net_pearson'   : 'BENCHMARK_4_GP_bootstrap_net_pearson',
                'run_t_test'                  : 'BENCHMARK_5_GP_t_test',
                'run_net_t_test'              : 'BENCHMARK_7_GP_net_t_test',
                'run_bootstrap_t_test'        : 'BENCHMARK_6_GP_bootstrap_t_test',
                'run_bootstrap_net_t_test'    : 'BENCHMARK_8_GP_bootstrap_net_t_test'}


verify_root_dir = '../data/verification'
results_dir = './run_dir/results'


def verify_results_directory(v_dir, results_dir):
    """ Compare the verification files in v_dir with the result files in results_dir """
    number_of_equals = 0
    v_files_list = os.listdir(v_dir)
    result_files_list = os.listdir(results_dir)
    for file_name in v_files_list:
        file_prefix = file_name[0:-4]
        for res_file_name in result_files_list:
            if res_file_name[0:len(file_prefix)] == file_prefix:
                res_file_full_name = os.path.join(results_dir, res_file_name)
                veri_file_name = os.path.join(v_dir, file_name)
                if filecmp.cmp(res_file_full_name, veri_file_name) == False:
                    print('NOT EQUAL:\n', file_name, '\n', res_file_name, '\n')
                else:
                    number_of_equals += 1

    return number_of_equals, len(v_files_list)


def run_all_BENCHMARKs_and_TESTs():
    """ run the make file targes for all yaml files and compre the results with their verification files """
    t0 = time.time()
    directory_methods_dict = {v: k for k, v in (GP_options_dict).items()}
    verification_directory_list = sorted(directory_methods_dict.keys(), reverse=True)

    for test_directory in verification_directory_list:
        verification_directory = os.path.join(verify_root_dir, test_directory)
        verification_method = 'make' + ' ' + directory_methods_dict[test_directory]
        print('\n\n\nRun Method:\t\t\t', verification_method, '\n', verification_directory)
        os.system(verification_method)

        n_eq, n_tot = verify_results_directory(verification_directory, results_dir)
        tt = '%0.2f'%(time.time() - t0)
        if n_eq == n_tot:
            print('\n\t\t', tt, 'sec\t', test_directory, n_eq, 'of', n_tot, 'equal\t', 'PASS')
        else:
            print('\n\t\t', tt, 'sec\t', n_eq, 'of', n_tot, 'equal\t', test_directory, '<-- FAIL')
        print('clear results and run_files')
        for tmp_file_name in os.listdir(results_dir):
            if os.path.isfile(os.path.join(results_dir, tmp_file_name)):
                os.remove(os.path.join(results_dir, tmp_file_name))


def main():
    try:
        os.system('make env_setup')
    except:
        pass

    run_all_BENCHMARKs_and_TESTs()
    print('\n')


if __name__ == "__main__":
    main()
