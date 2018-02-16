# KnowEnG's Gene Prioritization Pipeline
This is the Knowledge Engine for Genomics (KnowEnG), an NIH, BD2K Center of Excellence, Gene Prioritization Pipeline.

This pipeline **ranks** the rows of a given spreadsheet, where spreadsheet's rows correspond to gene-labels and columns correspond to sample-labels. The ranking is based on correlating gene expression data (network smoothed) against pheno-type data.

There are four prioritization methods, using either pearson or t-test as the measure of correlation:


| **Options**                                        | **Method**                           | **Parameters**            |
| -------------------------------------------------- | -------------------------------------| ------------------------- |
| Simple Correlation                                 | simple correlation                          | correlation               |
| Bootstrap Correlation                              | bootstrap sampling correlation       | bootstrap_correlation     |
| Correlation with network regularization            | network-based correlation            | net_correlation     |
| Bootstrap Correlation with network regularization  | bootstrapping w network correlation  | bootstrap_net_correlation |


Note: all of the correlation methods mentioned above use the Pearson or t-test correlation measure method.

* * * 
## How to run this pipeline with Our data
* * * 

### 1. Clone the Gene_Prioritization_Pipeline Repo
```
 git clone https://github.com/KnowEnG/Gene_Prioritization_Pipeline.git
```

### 2. Install the following (Ubuntu or Linux)
  ```
 apt-get install -y python3-pip
 apt-get install -y libblas-dev liblapack-dev libatlas-base-dev gfortran
 pip3 install numpy==1.11.1
 pip3 install pandas==0.18.1
 pip3 install scipy==0.19.1
 pip3 install scikit-learn==0.17.1
 apt-get install -y libfreetype6-dev libxft-dev
 pip3 install matplotlib==1.4.2
 pip3 install pyyaml
 pip3 install knpackage
```

### 3. Change directory to Gene_Prioritization_Pipeline

```
cd Gene_Prioritization_Pipeline
```

### 4. Change directory to test

```
cd test
```
 
### 5. Create a local directory "run_dir" and place all the run files in it
```
make env_setup
```

### 6. Use one of the following "make" commands to select and run a clustering option:


| **Command**                        | **Option**                                        | 
|:---------------------------------- |:------------------------------------------------- | 
| make run_pearson          | pearson correlation                                       |
| make run_bootstrap_pearson | bootstrap sampling with pearson correlation                    |
| make run_net_pearson     | pearson correlation with network regularization           |
| make run_bootstrap_net_pearson | bootstrap pearson correlation with network regularization |
| make run_t_test          | t-test correlation                                       |
| make run_bootstrap_t_test | bootstrap sampling with t-test correlation                    |
| make run_net_t_test     | t-test correlation with network regularization           |
| make run_bootstrap_net_t_test | bootstrap t-test correlation with network regularization |

 
* * * 
## How to run this pipeline with Your data
* * * 

__***Follow steps 1-3 above then do the following:***__

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
Look for examples of run_parameters in ./Gene_Prioritization_Pipeline/data/run_files/zTEMPLATE_GP_BENCHMARKS.yml
 ```
### * Modify run_paramters file  (YAML Format)
```
set the spreadsheet, network and phenotype data file names to point to your data
```

### * Run the Gene Prioritization Pipeline:

  * Update PYTHONPATH enviroment variable
   ``` 
   export PYTHONPATH='../src':$PYTHONPATH    
   ```
   
  * Run (in test directory with env_setup as described above)
   ```
  python3 ../src/gene_prioritization.py -run_directory ./run_dir -run_file zTEMPLATE_GP_BENCHMARKS.yml
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
| phenotype_name_full_path | directory+phenotype_response| Path and file name of user supplied phenotype response file |
| results_directory | directory | Directory to save the output files |
| number_of_bootstraps | 5 | Number of random samplings |
| cols_sampling_fraction | 0.9 | Select 90% of spreadsheet columns |
| rwr_max_iterations | 100| Maximum number of iterations without convergence in random walk with restart |
| rwr_convergence_tolerence | 1.0e-2 | Frobenius norm tolerence of spreadsheet vector in random walk|
| rwr_restart_probability | 0.5 | alpha in `V_(n+1) = alpha * N * Vn + (1-alpha) * Vo` |
| top_beta_of_sort| 100| Number of top genes selected |
| top_gamma_of_sort| 50| Number of top genes reported |

gg_network_name = STRING_experimental_gene_gene.edge</br>
spreadsheet_name = CCLE_Expression_ensembl.df</br>
phenotype_name = CCLE_drug_ec50_cleaned_NAremoved_pearson.txt

* * * 
## Description of Output files saved in results directory
* * * 

* Any method saves separate files per phenotype with name {phenotype}\_{method}\_{correlation_measure}\_{timestamp}\_viz.tsv. Genes are sorted in descending order based on `visualization_score`. </br>  

 | **Response** | **Gene_ENSEMBL_ID** | **quantitative_sorting_score** | **visualization_score** | **baseline_score** |
 |:-------------:|:------------:|:---------:|:--------------:|:--------------:|
 |   phenotype 1      |   gene 1     |    float    |    float         |   float          | 
 |    ...      |   ...     |    ...    |    ...         |   ...          | 
 |   phenotype 1      |   gene n     |    float    |    float         |   float          | 


* Any method saves sorted genes for each phenotype with name ranked_genes_per_phenotype\_{method}\_{correlation_measure}\_{timestamp}\_download.tsv.

 |**Ranking**| **phenotype 1** |**phenotype 2**|**...**|**phenotype n**|
 |:----:| :--------------------: |:--------------------:|---|:--------------------:|
 |1| gene </br> (most significant) |gene </br> (most significant)|...|gene </br> (most significant)|
 |...|...| ... |...|...|...|
 |n|gene </br> (least significant) |gene </br> (least significant)|...|gene </br> (least significant)|
 
 
 
* Any method saves spreadsheet with top ranked genes per phenotype with name  top_genes_per_phenotype\_{method}\_{correlation_measure}\_{timestamp}\_download.tsv.

 |**Genes**| **phenotype 1**|**...**|**phenotype n**|
 | :--------------------: |:--------------------:|---|:--------------------:|
 | gene 1 |1/0 |...|1/0|
 | ... |...|...|...|
 | gene n | 1/0|...|1/0|
