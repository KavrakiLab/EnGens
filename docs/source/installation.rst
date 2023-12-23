Installation instructions
===================================

.. _installation:

Installation
------------

There are a couple of options for installing EnGens locally.

1. Pulling the Docker image (all platforms)
2. Building the Docker image (all platforms)
3. Conda installation (LinuxOS)

For running EnGens demos online, checkout our Google Colab and Binder notebooks.

|Binder Link| 

.. |Binder Link| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/KavrakiLab/EnGens/binder?labpath=Workflow1-FeatureExtraction.ipynb

|Colab Link| 

.. |Colab Link| image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/drive/1rVeWH8CdUtbvmVCTZkxleRCRTe8dW5LN?usp=sharing


Docker pull
`````````````
Make sure you have the `Docker`_ installed.

.. _Docker: https://docs.docker.com/engine/install/

Open the console and simply run:

.. code-block:: console

   $ docker pull ac121/engens:latest
   
   

Docker build
`````````````
You can clone this repo and build the docker image yourself.

**prerequisites:**  `git`_ , `docker`_ 

.. _git: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git


1. Clone the github repo:

.. code-block:: console
   
   $ git clone https://github.com/KavrakiLab/EnGens.git


2. Build the image:

.. code-block:: console

   $ cd EnGens
   $ docker build -t test_engens:latest .

You're all set!


Conda build
`````````````

If you don't want to use docker, you can clone this repo and install using conda (or mamba which will be faster).

**prerequisites:** `conda`_ or `mamba`_

.. _conda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html#

.. _mamba: https://mamba.readthedocs.io/en/latest/installation.html

1. Clone the github repo:

.. code-block:: console

   $ git clone https://github.com/KavrakiLab/EnGens.git


2. Install with conda (or mamba)

.. code-block:: console

   $ cd EnGens
   $ conda env create -f ./environment.yaml
   $ #mamba env create -f ./environment.yml

   $ conda activate engens
   $ #mamba activate engens

   $ ./linux_setup.sh


If the command ./linux_setup.sh fails due to not having pypatch - do pip install pypatch.    


   
   
   
   
