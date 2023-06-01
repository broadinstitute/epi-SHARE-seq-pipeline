# ############################################################
# # Dockerfile for BROAD GRO share-seq-pipeline
# # Based on Debian slim
# ############################################################

FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f

LABEL maintainer = "Zhijian Li"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="Download data from cellxgene"

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive

## Create new user 
ENV USER=shareseq
WORKDIR /home/$USER
RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

# Install other libraries
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    gcc \
    git \
    python3 \
    python3-dev \
    python3-pip && \
    rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install tiledb tiledbsoma
# cellxgene-census

# COPY --chown=$USER:$USER src/python/get_cellxgene_data.py /usr/local/bin
# COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin

# USER ${USER}