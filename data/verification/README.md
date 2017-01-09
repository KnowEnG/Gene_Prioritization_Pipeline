# Verification directory files are benchmark result files
The Makefile in the test directory contains the targets, needed to build the **Gene Prioritization Pipeline** benchmarks.
For estimating base memory requirments use the above heuristic _gene_prioritization_p_memory_estimator.xlsx_


* Follow the instructions on the **Gene Prioritization Pipeline** landing page to set up the environmet:
```
    cd Gene_Prioritization_Pipeline/test
    make env_setup
```
### 1. Run the single drug small data test
```
    make run_small_data_test
```

* Compare the results with verification file: **TEST_1_bootstrap_net_correlstion_RESULTS.txt**

* The result file will be in **run_dir/results/**  TEST_1_bootstrap_net_correlation_XXX_a_timestamp_XXX.txt

* Results files are all based on the yaml parameter files preserved in **verification_yaml_set_used.zip**

* If the results are not equal, unzip the file in **data/verification** and replace the .yml files in **test/run_dir**

### 2. Run the multi-drug samll data test **test/Makefile**

```
    make run_small_data_multidrug_test
```
* Results will match **data/verification/TEST_1_multi_drug_debug_RESUTS.zip** (when unzipped)
