# ----------------------------------------------------------------
# Variables setup Section                                        -
# ----------------------------------------------------------------
SCRIPT = 		../src/gene_prioritization.py
RUN_DIR = 		./run_dir
RUN_FILES = 		../data/run_files
RESULTS_DIR = 		$(RUN_DIR)/results
INPUT_DATA_DIR = 	../data
UNIT_TEST_TEMP_DIR = 	./unit/tmp

RMR = 			rm -rf
GZIP = 			gzip
GZIP_D = 		gzip -d
MKDIR_P = 		mkdir -p

# ----------------------------------------------------------------
# Benchmarks setup Section                                       -
# ----------------------------------------------------------------
env_setup: create_run_dir copy_run_files decompress_input_data

decompress_input_data:
	$(GZIP_D) $(INPUT_DATA_DIR)/networks/* $(INPUT_DATA_DIR)/spreadsheets/*
 
compress_input_data:
	$(GZIP) $(INPUT_DATA_DIR)/networks/* $(INPUT_DATA_DIR)/spreadsheets/*

create_run_dir:
	$(MKDIR_P) $(RESULTS_DIR) 

copy_run_files:
	cp $(RUN_FILES)/*.yml $(RUN_DIR) 

# ----------------------------------------------------------------
# Large Benchmarks RUN Section                                   -
# ----------------------------------------------------------------
# pearson

run_correlation_pearson:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_1_correlation_pearson.yml

run_net_correlation_pearson:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_1_net_correlation_pearson.yml

run_bootstrap_correlation_pearson:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_1_bootstrap_correlation_pearson.yml

run_bootstrap_net_correlation_pearson:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_1_bootstrap_net_correlation_pearson.yml

# ----------------------------------------------------------------
# t_test

run_correlation_t_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_1_correlation_t_test.yml

run_net_correlation_t_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_1_net_correlation_t_test.yml

run_bootstrap_correlation_t_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_1_bootstrap_correlation_t_test.yml

run_bootstrap_net_correlation_t_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_1_bootstrap_net_correlation_t_test.yml

# ----------------------------------------------------------------
# Small Benchmarks RUN Section                                   -
# ----------------------------------------------------------------
run_small_data_multidrug_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_1_multidrug_small_data_pearson.yml

run_small_data_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_1_single_drug_small_data_pearson.yml

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
final_clean: clean_dir_recursively compress_input_data

clean_dir_recursively:
	$(RMR) $(RUN_DIR)	

unit_test_cleanup:
	$(RMR) $(UNIT_TEST_TEMP_DIR)
