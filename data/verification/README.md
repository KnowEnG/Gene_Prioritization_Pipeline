# Verification directory files are benchmark result files
The Makefile in the test directory contains the targes, needed to build the **Gene Prioritization Pipeline** benchmarks.


* Follow the instructions on the **Gene Prioritization Pipeline** landing page to set up the environmet:
```
    cd Gene_Prioritization_Pipeline/test
    make decompress_input_data
    make env_setup
```
### 1. Run the small data test
```
    make run_small_data_test
```

### 2. Compare the results with verification file: TEST_1_bootstrap_net_correlstion_RESULTS.txt
```
    open the verification file in exel or a text editor, compare the result file with a name like this:
    run_dir/results/TEST_1_bootstrap_net_correlation_XXXXX..... .txt
```

