### The compressed data set ./results was generated using "bootstrap_corrleation.yml"
Assuming you have cloned the Gene_Prioritization_Pipeline to your computer as described in the root directory; you may run this yaml file from a terminal, by changing to the "Gene_Prioritization_Pipeline/test/" directory and typing the commands listed below.
* make env_setup
* cp ../data/verification/bootstrap_correlation.yml ./run_dir/bootstrap_correlation.yml
* python3 gene_prioritization.py -run_directory ./run_dir -run_file bootstrap_correlation.yml

It may take more than 60 seconds depending on your computer hardware.

### The output files in "results" are compressed and may be uncompressed by changing to "./Dec_22_2016_24drug_output/results" and entering the unzip command.
* gzip -d *.gz

The unzipped *.tsv files may be viewed as spreadsheet files to compare with your *.tsv results.

