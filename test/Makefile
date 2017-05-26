# ----------------------------------------------------------------
# Variables setup Section                                        -
# ----------------------------------------------------------------
SCRIPT =             ../src/gene_prioritization.py
RUN_DIR =            ./run_dir
RUN_FILES =          ../data/run_files
RESULTS_DIR =        $(RUN_DIR)/results
INPUT_DATA_DIR =     ./data
UNIT_TEST_TEMP_DIR = ./unit/tmp
DATA_DIR =           ../data/networks ../data/spreadsheets
INPUT_DATA =         $(foreach dir,$(DATA_DIR),$(wildcard $(dir)/*.gz))
INPUT_DATA_SUFFIX =  $(suffix $(INPUT_DATA))

RMR =                rm -rf
GZIP =               gzip
GZIP_D =             gzip -d
MKDIR_P =            mkdir -p

# ----------------------------------------------------------------
# Benchmarks setup Section                                       -
# ----------------------------------------------------------------
env_setup: create_run_dir copy_run_files decompress_input_data

decompress_input_data:
	@echo Decompress input data...
	$(GZIP_D) $(INPUT_DATA)

create_run_dir:
	$(MKDIR_P) $(RESULTS_DIR) 

copy_run_files:
	cp $(RUN_FILES)/*.yml $(RUN_DIR) 

run_all_methods: run_all_small run_all_pearson run_all_t_test
# ----------------------------------------------------------------
# Large Benchmarks RUN Section                                   -
# ----------------------------------------------------------------
run_all_pearson: run_pearson run_bootstrap_pearson run_net_pearson run_bootstrap_net_pearson

run_pearson:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_1_GP_pearson.yml

run_net_pearson:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_3_GP_net_pearson.yml

run_bootstrap_pearson:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_2_GP_bootstrap_pearson.yml

run_bootstrap_net_pearson:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_4_GP_bootstrap_net_pearson.yml

# ----------------------------------------------------------------
run_all_t_test: run_t_test run_bootstrap_t_test run_net_t_test run_bootstrap_net_t_test

run_t_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_5_GP_t_test.yml

run_net_t_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_7_GP_net_t_test.yml

run_bootstrap_t_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_6_GP_bootstrap_t_test.yml

run_bootstrap_net_t_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_8_GP_bootstrap_net_t_test.yml

# ----------------------------------------------------------------
# Small Benchmarks RUN Section                                   -
# ----------------------------------------------------------------
run_all_small: run_single_drug_pearson run_multidrug_pearson run_single_drug_t_test run_multidrug_t_test

run_single_drug_pearson:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_1_GP_single_drug_pearson.yml

run_multidrug_pearson:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_2_GP_many_drugs_pearson.yml

run_single_drug_t_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_3_GP_single_drug_t_test.yml

run_multidrug_t_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_4_GP_many_drugs_t_test.yml

# ----------------------------------------------------------------
# - VERIFICATION TESTS RUN SECTION                               -
# ----------------------------------------------------------------
verification_tests:
	python3 ./integration/verify_benchmarks.py

# ----------------------------------------------------------------
# - UNIT TESTS RUN SECTION                                       -
# ----------------------------------------------------------------
unit_test_setup:
	$(MKDIR_P) $(UNIT_TEST_TEMP_DIR)

unit_tests:
	cd unit; make all_unit_tests

# ----------------------------------------------------------------
# - CLEAN UP SECTION                                             -
# ----------------------------------------------------------------	
final_clean: clean_dir_recursively 

clean_dir_recursively:
	$(RMR) $(RUN_DIR)	

unit_test_cleanup:
	$(RMR) $(UNIT_TEST_TEMP_DIR)
