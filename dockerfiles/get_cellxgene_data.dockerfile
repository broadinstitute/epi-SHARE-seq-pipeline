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
RUN groupadd -r $USER && \
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

# Install python 3
RUN apt-get update
RUN apt-get install -y --no-install-recommends python3 python3-pip
RUN rm -rf /var/lib/apt/lists/*

# Install cellxgene package
RUN python3 -m pip install cellxgene-census

# Copy script to /use/local/bin
COPY src/python/get_cellxgene_data.py /usr/local/bin
COPY src/bash/monitor_script.sh /usr/local/bin
