# KnowEnG's Gene Prioritization Pipeline
This is the Knowledge Engine for Genomics (KnowEnG), an NIH, BD2K Center of Excellence, Gene Prioritization Pipeline.

This pipeline **ranks** the rows of a given spreadsheet, where spreadsheet's rows correspond to gene-labels and columns correspond to sample-labels. The ranking is based on correlating gene expression data (network smoothed) against pheno-type data.

There are four prioritization methods that one can choose from:


| **Options**                                        | **Method**                           | **Parameters**            |
| -------------------------------------------------- | -------------------------------------| ------------------------- |
| Correlation                                        | correlation                          | correlation               |
| Bootstrap Correlation                              | bootstrap sampling correlation       | net_correlation           |
| Correlation with network regularization            | network-based correlation            | bootstrap_correlation     |
| Bootstrap Correlation with network regularization  | bootstrapping w network correlation  | bootstrap_net_correlation |


Note: all of the correlation methods mentioned above use the Pearson correlation coefficient.

* * * 
## How to run this pipeline with Our data
* * * 
###1. Get Access to KnowEnG-Research Repository:
Email omarsobh@illinois.edu infrastructure team (IST) lead to:

* __Access__ KnowEnG-Research github repo

###2. Clone the Gene_Prioritization_Pipeline Repo
```
 git clone https://github.com/KnowEnG-Research/Gene_Prioritization_Pipeline.git
```
 
###3. Install the following (Ubuntu or Linux)
  ```
 apt-get install -y python3-pip
 apt-get install -y libblas-dev liblapack-dev libatlas-base-dev gfortran
 pip3 install numpy==1.11.1
 pip3 install pandas==0.18.1
 pip3 install scipy==0.18.0
 pip3 install scikit-learn==0.17.1
 apt-get install -y libfreetype6-dev libxft-dev
 pip3 install matplotlib==1.4.2
 pip3 install pyyaml
 pip3 install knpackage
```

###4. Change directory to Gene_Prioritization_Pipeline

```
cd Gene_Prioritization_Pipeline
```

###5. Change directory to test

```
cd test
```
 
###6. Create a local directory "run_dir" and place all the run files in it
```
make env_setup
```

###7. Use one of the following "make" commands to select and run a clustering option:


| **Command**                        | **Option**                                        | 
|:---------------------------------- |:------------------------------------------------- | 
| make run_correlation          | correlation                                       |
| make run_net_correlation           | bootstrap sampling correlation                    |
| make run_bootstrap_correlation     | correlation with network regularization           |
| make run_bootstrap_net_correlation | bootstrap correlation with network regularization |

 
* * * 
## How to run this pipeline with Your data
* * * 

__***Follow steps 1-4 above then do the following:***__

### * Create your run directory

 ```
 mkdir run_directory
 ```

### * Change directory to the run_directory

 ```
 cd run_directory
 ```

### * Create your results directory

 ```
 mkdir results_directory
 ```
 
### * Create run_paramters file  (YAML Format)
 ``` 
Look for examples of run_parameters in ./Gene_Prioritization_Pipeline/data/run_files/run_parameters_template.yml
 ```
### * Modify run_paramters file  (YAML Format)
```
set the spreadsheet, network and drug_response (phenotype data) file names to point to your data
```

### * Run the Samples Clustering Pipeline:

  * Update PYTHONPATH enviroment variable
   ``` 
   export PYTHONPATH='../src':$PYTHONPATH    
   ```
   
  * Run
   ```
  python3 ../src/gene_prioritization.py -run_directory ./ -run_file run_parameters_template.yml
   ```

* * * 
## Description of "run_parameters" file
* * * 

| **Key**                   | **Value** | **Comments** |
| ------------------------- | --------- | ------------ |
| method                    | correlation or net_correlation or bootstrap_correlation or bootstrap_net_correlation | Choose gene prioritization method |
| correlation_measure       | pearson or t_test | Choose correlation measure method |
| gg_network_name_full_path | directory+gg_network_name |Path and file name of the 4 col network file|
| spreadsheet_name_full_path | directory+spreadsheet_name|  Path and file name of user supplied gene sets |
| drug_response_full_path | directory+drug_response_full_path| Path and file name of user supplied drug response file |
| results_directory | directory | Directory to save the output files |
| number_of_bootstraps | 5 | Number of bootstraps |
| cols_sampling_fraction | 0.9 | Select 90% of spreadsheet columns |
| rwr_max_iterations | 100| Maximum number of iterations without convergence in random walk with restart |
| rwr_convergence_tolerence | 1.0e-2 | Frobenius norm tolerence of spreadsheet vector in random walk|
| rwr_restart_probability | 0.5 | alpha in `V_(n+1) = alpha * N * Vn + (1-alpha) * Vo` |
| top_beta_of_sort| 100| Number of top genes selected 
gg_network_name_full_path = STRING_experimental_gene_gene.edge</br>
spreadsheet_name = CCLE_Expression_ensembl.df</br>
drug_response_full_path = CCLE_drug_ec50_cleaned_NAremoved.txt

* * * 
## Description of Output files saved in results directory
* * * 

* Output files of all four methods save separate files per phenotype with name {phenotype}\_{method}\_{correlation_measure}\_{timestamp}.tsv. Genes are sorted in descending order based on `quantitative_sorting_score`. </br>  

 | **Response** | **Gene_ENSEMBL_ID** | **quantitative_sorting_score** | **visualization_score** | **baseline_score** |
 |:-------------:|:------------:|:---------:|:--------------:|:--------------:|
 |   phenotype 1      |   gene 1     |    float    |    float         |   float          | 
 |    ...      |   ...     |    ...    |    ...         |   ...          | 
 |   phenotype 1      |   gene n     |    float    |    float         |   float          | 


* Output files of all four methods save sorted genes for each phenotype with name all_phenotypes\_{method}\_{correlation_measure}\_{timestamp}_download.tsv.

 | **phenotype 1** |**phenotype 2**|**...**|**phenotype n**|
 | :--------------------: |:--------------------:|---|:--------------------:|
 | gene </br> (most significant) |gene </br> (most significant)|...|gene </br> (most significant)|
 |...| ... |...|...|...|
 |gene </br> (least significant) |gene </br> (least significant)|...|gene </br> (least significant)|
 
 
 
* Output files of all four methods save spreadsheet with top ranked genes per phenotype with name  all_phenotypes\_{method}\_{correlation_measure}\_{timestamp}\_original.tsv.

 |**Genes**| **phenotype 1**|**...**|**phenotype n**|
 | :--------------------: |:--------------------:|---|:--------------------:|
 | gene 1 |1/0 |...|1/0|
 | ... |...|...|...|
 | gene n | 1/0|...|1/0|
