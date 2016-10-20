# Building The Gene Prioritization Pipeline Docker Image
The Dockefile in this directory contains all the commands, in order, needed to build the **Gene Prioritization Pipeline** docker image.


* Run the "make" command to build the **Gene Prioritization Pipeline** docker image (output: docker image called "gene_prioritization_pipeline" and a tag with today's date and time):
```
    make build_docker_image
```

* Login to docker hub. When prompted, enter your password and press enter:
```
    make login_to_dockerhub username=*enter your docker login here* email=*enter your email here*
```

* Upload your image to docker hub:
```
    make push_to_dockerhub
```

* * * 
## How to run this docker image
* * * 

1 Change directory to the directory  where you want to run.

2 docker run -v \`pwd\`:/home/test/run_dir -it knowengdev/gene_prioritization_pipeline:10_19_2016

3 cd test

3 make env_setup

4 Edit changes in the (correlation) .yml file in the run_dir created by env_setup.

5 make correlation

* The make options in Gene_Prioritization_Pipeline/README.md apply.

* Check on docker.hub to get the latest image. 

* Your data will be in the directory where you ran the command in line 2 when you exit docker.
