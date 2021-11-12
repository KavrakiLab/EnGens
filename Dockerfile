FROM python:3.6.7

LABEL MAINTAINER="Anja Conev <ac121@rice.edu>"

WORKDIR /var/www/

ENV PACKAGES="\
	wget \ 
	vim \
"

ENV PYTHON_PACKAGES="\
    numpy \
    matplotlib \
    scipy \
    scikit-learn \
    pandas \
	notebook \
	jupyterhub \
" 

ENV BIO_PACKAGES="\
    mdtraj \
    pyemma \
    nglview \
    biopython \
" 

RUN pip install --upgrade pip 
	
RUN pip install --no-cache-dir $PYTHON_PACKAGES 

RUN pip install --no-cache-dir $BIO_PACKAGES 

RUN pip install --no-cache-dir jupyterlab 
RUN pip install --no-cache-dir wget 
RUN pip install --no-cache-dir pytraj 
	
ARG NB_USER=engen
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
	
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
WORKDIR ${HOME}