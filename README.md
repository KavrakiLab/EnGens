# EnGeNs
Private repository for the development of Ensemble Generation Notebooks (EnGeNs)

## Instructions for running

### Step 1 - pull the docker image 

`docker pull ac121/engens:latest`


### Step 2 - run the following command from the working directory with this code

`docker run -it --rm -v ${pwd}:/home/engen/ -p 8888:8888 ac121/engens:latest jupyter-lab --allow-root --ip=0.0.0.0 --port=8888`


___

P.S. Folder **./baudry-data-scripts** contains various scripts used to analyse this data: generating trajectory from pdbqt files, extracting the binding site.

