FROM mambaorg/micromamba:1.4.2
LABEL MAINTAINER="Anja Conev <ac121@rice.edu>"

# Create the environment:
COPY --chown=$MAMBA_USER:$MAMBA_USER ./environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python will not be found)
RUN echo $(python --version)

# Additional installations
RUN pip install -i https://test.pypi.org/simple/ engens
RUN git clone https://github.com/hsidky/hde.git
COPY hde.patch .
RUN cp hde.patch ./hde/hde.patch
RUN cd hde && echo $(ls) && git apply hde.patch
RUN pip install ./hde

COPY --chown=$MAMBA_USER:$MAMBA_USER ./dependencies/mTM-align.tar.bz2 /tmp/mTM-align.tar.bz2

RUN tar -xvf /tmp/mTM-align.tar.bz2 -C /tmp/
RUN cp /tmp/mTM-align/src/mTM-align ${CONDA_PREFIX}/bin/mTM-align \
&& rm /tmp/mTM-align.tar.bz2

RUN pypatch apply ./pdbfixer.patch pdbfixer

ADD --chown=$MAMBA_USER:$MAMBA_USER ./notebooks/ ${HOME}