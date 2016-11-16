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

* Compare the results with verification file: **TEST_1_bootstrap_net_correlstion_RESULTS.txt**

* The result file will be in **run_dir/results/**  TEST_1_bootstrap_net_correlation_XXX_a_timestamp_XXX.txt

* Results files are all based on the yaml parameter files preserved in **verification_yaml_set_used.zip**

* If the results are not equal, unzip the file in **data/verification** and replace the .yml files in **test/run_dir**

### 2. Run the other targets in the **test/Makefile**

```
    make run_correlation
        or
    make run_bootstrap_correlation
```
* Results will match **data/verification/17-AAG_correlation_final_result.txt**


```
    make run_bootstrap_net_correlation
        or
    make run_net_correlation
```
* Results will match **data/verification/17-AAG_net_correlation_final_result.txt**


```
    make run_multi_drug_small
```
* Results will match **data/verification/TEST_1_multi_drug_debug_RESUTS.zip** (when unzipped)


```
    make run_multi_drug_bootstrap_correlation
```
* The first five files will match the five unzipped from **data/verification/CCLE_firs_5_bootstrap_net_correlation.zip.txt**
