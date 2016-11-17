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
run_correlation:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file correlation.yml 

run_net_correlation:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file net_correlation.yml 

run_bootstrap_correlation:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file bootstrap_correlation.yml 

run_bootstrap_net_correlation:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file bootstrap_net_correlation.yml 

run_all_correlations: run_correlation run_net_correlation run_bootstrap_correlation run_bootstrap_net_correlation

# ----------------------------------------------------------------
# Small Benchmarks RUN Section                                   -
# ----------------------------------------------------------------
run_small_data_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_1_bootstrap_net_correlation.yml 

run_debug_data_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_1_debug.yml 

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
