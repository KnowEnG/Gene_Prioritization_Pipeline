# Building The Gene Prioritization Pipeline Docker Image
The Dockefile in this directory contains all the commands, in order, needed to build the **Gene Prioritization Pipeline** docker image.


* Run the "make" command to build the **Gene Prioritization Pipeline** docker image (output: docker image called "gene_prioritization_pipeline" and a tag with today's date and time):
```
    make build_docker_image
```

* Login to docker hub. When prompted, enter your password and press enter:
```
    make login_to_dockerhub username=<u>enter your docker login here</u> email=<u>enter your email here</u>
```

* Upload your image to docker hub:
```
    make push_to_dockerhub
```

* * * 
## How to run this docker image
* * * 

### 1. Run the following command with the specified docker image:
```
docker run -v \`pwd\`:/home/test/run_dir/ -it knowengdev/gene_prioritization_pipeline:10_19_2016 
```

### 2. Change directory to the "test" directory
```
cd test
```

### 3. Create the local directory "run_dir" and place all the run files in it
```
make env_setup
```

### 4. Run the Gene Prioritization Pipeline
```
make correlation
```
