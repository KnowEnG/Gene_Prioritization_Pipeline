# KnowEnG's Gene Prioritization Pipeline
This is the Knowledge Engine for Genomics (KnowEnG), an NIH, BD2K Center of Excellence, Gene Prioritization Pipeline.

This pipeline **ranks** the columns of a given spreadsheet, where spreadsheet's columns correspond to sample-labels and rows correspond to gene-labels.

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
| make run_gene_correlation          | correlation                                       |
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
