# EnGens 

![Alt text](./preprint/logo.svg)


Repository for the computational framework for generation and analysis of representative protein conformational ensembles.
___

## Demo 

Try runnning our notebooks on Binder:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/KavrakiLab/EnGens/binder?labpath=Workflow1-FeatureExtraction.ipynb)

Try running a short demo on Google Colab: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1rVeWH8CdUtbvmVCTZkxleRCRTe8dW5LN?usp=sharing)


## Installation instructions

#### Docker image pull
The prefered and easiest is by pulling the docker image made available publicly.

**prerequisites:** [docker](https://docs.docker.com/get-docker/)

Just pull the image:

```
docker pull ac121/engens:latest
```
You're all set!

_Note: this step should not take longer than 15min. On Windows the PowerShell sometimes gets stuck - do a right click in the terminal to check the progress after 10-15min._

For other installation options check the section [Advanced Installation](#advanced-installation) bellow.

## Running EnGens

### Linux: run the following command from the working directory with this code

`docker run -it --rm -v $(pwd):/home/engen/ -p 8888:8888 ac121/engens:latest jupyter notebook --ip=0.0.0.0 --port=8888`

### Windows: run the following command from the working directory with this code

`docker run -it --rm -v ${pwd}:/home/engen/ -p 8888:8888 ac121/engens:latest jupyter notebook --ip=0.0.0.0 --port=8888`


Here we provide a pipeline for generating ensembles of conformations from molecular dynamics trajectories as a first step towards ensemble docking.
The pipeline consists of the following steps:

  1. **Featurization of the trajectory** - as a first step you need to computationally represent your trajectory. In this step we heavilty rely on the <a href=http://www.emma-project.org/latest/> PyEMMA </a> package. You can use any of the featurizations they provide - or implement your own in order to represent the trajectories. As part of this step you can also choose a subset of the structure for which you want to do the analysis (e.g. binding site). 
  2. **Dimensionality reduction of the derived trajectory features** In order to proceed with the clustering step you need to reduce the dimensionality of the derived features. We propose and implement a couple of approaches <a href=http://www.emma-project.org/latest/api/generated/pyemma.coordinates.tica.html#pyemma.coordinates.tica`>TICA</a>, <a=https://github.com/hsidky/srv>HDE</a> and PCA. Note that TICA and HDE use the variational approach to conformational dynamics (VAC) and are thus better suited for the analysis of the MD data.
  3. **Clustering the trajectory and extracting cluster representatives** Finally we propose two methods (<a href=https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html>KMeans</a> and <a href=https://scikit-learn.org/stable/modules/mixture.html>Gaussian Mixture Models</a>) to preform the clustering of the data with the reduced dimensionality. The cluster centers are then defined as the representative structures of the MD and are extracted to form the ensemble that can later be used for ensemble docking.
  

## Workflows

You can find the three workflows inside the `./notebooks/` directory. To perform the full analysis you need to run the workflows in sequential order (Workflow1, Workflow2, Workflow3). Here we provide the overview of the notebooks and point to the steps where you need to make important choicec.

* *Workflow1-FeatureExtraction.ipynb*: in this notebook you will be featurizing your trajectory. 
 
    * First, you need to load the trajectory file and the reference PDB file of your MD. Then you need to decide if you want to apply your features to the full PDB structure or only an important subset (for example the binding site). If you know the binding site this is a good place to indicate it because it will reduce the computational complexity of steps that follow. You will use the <a href=https://mdtraj.org/1.9.4/atom_selection.html>mdtraj atom selection syntax </a> to specify the region of the structure. 

    * Next, you will have to choose a good featurization for the trajectory. We provide all options available in  <a href=http://www.emma-project.org/latest/api/generated/pyemma.coordinates.featurizer.html>PyEMMA featurizers</a>. Alternativelly you can use our default initialization of the featurizers and get three featurizers: one with amino-acid pairwose distances, one using torsion angles and one that includes both the distances and torsion angles. 
    
    * If you list more then one potential featurizer (like we do three in out default setting) you will have to choose between them. For this we provide the analysis using the VAMP2 score. You can see more details on how to choose features based on the VAMP2 score in this <a href=http://www.emma-project.org/latest/tutorials/notebooks/00-pentapeptide-showcase.html>PyEMMA tutorial</a>. You want to choose the featurization with the highest VAMP2 score in this analysis but you also want to vary the two important VAMP parameters (dimensionality and lag time) to see the robustness of the score and pick the featurization that has a consistently high score.

    * Finally based on the above analysis you can choose a featurization and save your result for the Workflow2.

* *Workflow2-DimensionalityReduction.ipynb*: in this notebook you can use different techniques to reduce the dimensionality of the features you generated in Workflow1.
    
    *   First you will load the data from Wrokflow1.
    *   Next, you need to choose the dimensionality reduction tecnhique (we provide HDE, TICA, PCA)
    *   If you choose a VAC approach (HDE or TICA) it will be important to choose the lag time correctly. You want to choose the lowest lag at which the timecale has converged. 
    *   If you choose TICA or PCA it will be important to correctly choose the number of dimensions to proceed with by looking at the amount of variance explained (or kinetic variance in case of TICA).
    *   Finally, after choosing appropriate lag time and the number of dimensions you can save this data  and proceed to Workflow3
    
* *Workflow3-Clustering.ipynb*: in this notebook you will use a clustering method tu cluster the trajectory and pick the representative samples

    * First you load the data from the Workflow2.
    * Next, you choose a clustering method (we provide Kmeans and GMMs but you can easily implement your own method).
    * Then you will have to find the best clustering that fits your data. You might want to change different parameters of your clustering. For example for Kmeans you can change the number of clusters, while for GMMs you can vary the number of components. To evaluate the clustering and choose the best number of clusters you can use the elbow method or the silhouette analysis.
    * Once you choose the best clustering you can also choose a subset of clusters that you want to explore. Perhaps you want your ensemble to contain representatives from clusters that have large number of components while ignoring the small clusters. 
    * Finally, after choosing the clusters you can find the cluster centers and structures closest to the cluster centers will be extracted to form your ensemble.

## Code

All the code and classes used in the notebooks are found in the directory `./EnGeNs/engens_code/engens/core/`

## Advanced Installation

#### 1. Docker image build
You can clone this repo and build the docker image yourself.

**prerequisites:**  [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git), [docker](https://docs.docker.com/get-docker/)


1. Clone the github repo:

```
git clone https://github.com/KavrakiLab/EnGens.git
```

2. Build the image:

```
cd EnGens
docker build -t test_engens:latest .
```

You're all set!

#### 2. Conda environment build
If you don't want to use docker, you can clone this repo and install using conda (or mamba which will be faster).

**prerequisites:** [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html#) or [mamba](https://mamba.readthedocs.io/en/latest/installation.html)

1. Clone the github repo:

```
git clone https://github.com/KavrakiLab/EnGens.git
```

2. Install with conda (or mamba)

```
cd EnGens
conda env create -f ./environment.yaml
#mamba create -f ./environment.yml

conda activate engens
#mamba activate engens

./linux_setup.sh
#or ./windows_setup.sh
```

___

P.S. Folder **./baudry-data-scripts** contains various scripts used to analyse this data: generating trajectory from pdbqt files, extracting the binding site.

