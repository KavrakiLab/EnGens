# EnGeNs
Private repository for the development of Ensemble Generation Notebooks (EnGeNs)

## Instructions for running

### Step 1 - build/pull the docker image 

option 1 - build docker image by running:

`docker build image .`

option 2 - pull the image

`docker pull ac121/engens`


### Step 2 - run the following command from the working directory with this code

`docker run -it --rm -v ${pwd}:/home/engen/ -p 8888:8888 ac121/engens:latest jupyter-lab --ip=0.0.0.0 --port=8888`


