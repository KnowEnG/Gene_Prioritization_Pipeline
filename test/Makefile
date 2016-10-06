MKDIR_P = mkdir -p
RMR = rm -rf
GZIP = gzip
GZIP_D = gzip -d
INPUT_DATA_DIR = ../data
RUN_DIR = ./run_dir
UNIT_TEST_TEMP_DIR = ./unit/tmp
RESULTS_DIR = $(RUN_DIR)/results
RUN_FILES = ../data/run_files
SCRIPT = ../src/gene_prioritization.py


env_setup: create_run_dir copy_run_files

clean_dir_recursively:
	$(RMR) $(RUN_DIR)

run_all_correlations: run_gene_correlation run_net_correlation run_bootstrap_correlation run_bootstrap_net_correlation

run_gene_correlation:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file correlation_run_file.yml 

run_net_correlation:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file net_correlation_run_file.yml 

run_bootstrap_correlation:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file bootstrap_correlation_run_file.yml 

run_bootstrap_net_correlation:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file bootstrap_net_correlation_run_file.yml 

run_small_data_test:
	python3 $(SCRIPT) -run_directory $(RUN_DIR) -run_file run_parameters_small_test.yml 

decompress_input_data:
	$(GZIP_D) $(INPUT_DATA_DIR)/*
 
compress_input_data:
	$(GZIP) $(INPUT_DATA_DIR)/*

create_input_data:
	$(MKDIR_P) $(INPUT_DATA_DIR)	

create_run_dir:
	$(MKDIR_P) $(RESULTS_DIR) 

copy_run_files:
	cp $(RUN_FILES)/*.yml $(RUN_DIR) 

unit_test_setup:
	$(MKDIR_P) $(UNIT_TEST_TEMP_DIR)

unit_tests:
	cd unit; make all_unit_tests
	
unit_test_cleanup:
	$(RMR) $(UNIT_TEST_TEMP_DIR)

benchmarks_run: benchmarks_setup benchmarks_copy_files

benchmarks_setup:
	$(MKDIR_P) $(RESULTS_DIR)

benchmarks_copy_files:
	cp $(RUN_FILES)/*.yml $(RUN_DIR)

run_gene_correlation_benchmark:
	PYTHONPATH='../src/' python3 ./benchmarks/gene_prioritization_benchmarks.py -run_directory $(RUN_DIR) -run_file correlation_run_benchmarks.yml

run_small_data_benchmark:
	PYTHONPATH='../src/' python3 ./benchmarks/gene_prioritization_benchmarks.py -run_directory $(RUN_DIR) -run_file run_parameters_small_test.yml 
