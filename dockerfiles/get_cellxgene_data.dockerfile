# ############################################################
# # Dockerfile for BROAD GRO share-seq-pipeline
# # Based on Debian slim
# ############################################################

FROM python:3.8.16-slim@sha256:1083b6d3ebe9a9cf9fd7ad311643ae1a1d51cb8f13b6fd465de2fc5c5b1990ce

LABEL maintainer = "Zhijian Li"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="Download data from cellxgene"

## Create new user 
ENV USER=shareseq
WORKDIR /home/$USER

RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

# Install other libraries
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    python3-pip && \
    rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install cellxgene-census

ENV PYTHONPATH="/usr/local/python:$PYTHONPATH"

COPY --chown=$USER:$USER src/python/get_cellxgene_data.py /usr/local/bin/
COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin

USER ${USER}