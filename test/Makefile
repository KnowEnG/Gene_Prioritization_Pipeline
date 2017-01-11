# ----------------------------------------------------------------
# Variables setup Section                                        -
# ----------------------------------------------------------------
SCRIPT = 		../src/gene_prioritization.py
RUN_DIR = 		./run_dir
RUN_FILES = 		../data/run_files
RESULTS_DIR = 		$(RUN_DIR)/results
INPUT_DATA_DIR = 	../data
UNIT_TEST_TEMP_DIR = 	./unit/tmp
DATA_DIR = ../data/networks ../data/spreadsheets
INPUT_DATA = $(foreach dir,$(DATA_DIR),$(wildcard $(dir)/*))
INPUT_DATA_SUFFIX = $(suffix $(INPUT_DATA))

RMR = 			rm -rf
GZIP = 			gzip
GZIP_D = 		gzip -d
MKDIR_P = 		mkdir -p

# ----------------------------------------------------------------
# Benchmarks setup Section                                       -
# ----------------------------------------------------------------
env_setup: create_run_dir copy_run_files decompress_input_data

decompress_input_data:
	@echo Decompress input data...
ifeq ($(filter-out .gz, $(INPUT_DATA_SUFFIX)),)
	$(GZIP_D) $(INPUT_DATA)
else
	@echo No gzip file needs to be decompressed.
endif	
 
compress_input_data:
	@echo Compress input data...
ifeq ($(filter .gz, $(INPUT_DATA_SUFFIX)),)
	$(GZIP) $(INPUT_DATA)
else
	@echo Target file in already in gzip format. Aborting compressing input data...
endif

create_run_dir:
	$(MKDIR_P) $(RESULTS_DIR) 

copy_run_files:
	cp $(RUN_FILES)/*.yml $(RUN_DIR) 

# ----------------------------------------------------------------
# Large Benchmarks RUN Section                                   -
# ----------------------------------------------------------------
# pearson

run_correlation_pearson:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_1_GP_pearson.yml

run_net_correlation_pearson:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_3_GP_net_pearson.yml

run_bootstrap_correlation_pearson:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_2_GP_bootstrap_pearson.yml

run_bootstrap_net_correlation_pearson:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_4_GP_bootstrap_net_pearson.yml

# ----------------------------------------------------------------
# t_test

run_correlation_t_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_5_GP_t_test.yml

run_net_correlation_t_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_7_GP_net_t_test.yml

run_bootstrap_correlation_t_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_6_GP_bootstrap_t_test.yml

run_bootstrap_net_correlation_t_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file BENCHMARK_8_GP_bootstrap_net_t_test.yml

# ----------------------------------------------------------------
# Small Benchmarks RUN Section                                   -
# ----------------------------------------------------------------
run_small_data_multidrug_pearson:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_2_GP_many_drugs_pearson.yml

run_small_data_pearson:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_1_GP_single_drug_pearson.yml

run_small_data_multidrug_t_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_4_GP_many_drugs_t_test.yml

run_small_data_t_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file TEST_3_GP_single_drug_t_test.yml

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
