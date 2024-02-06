############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on Debian slim
############################################################

#FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f
FROM ubuntu:latest

LABEL maintainer = "Eugenio Mattei"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="qc-atac"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    gcc \
    git \
    pigz \
    python3 \
    python3-dev \
    python3-pip \
    r-base \
    tabix \
    wget \
    zlib1g-dev &&\
    rm -rf /var/lib/apt/lists/*

# Install packages for python3 scripts (pysam, SAMstats)
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install --no-cache-dir --ignore-installed numpy matplotlib pandas plotnine pysam xopen snapatac2

# Install R packages
RUN R --no-echo --no-restore --no-save -e "install.packages('ggplot2')"

# Create and setup new user
ENV USER=shareseq
WORKDIR /home/$USER

RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

RUN wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/bedClip && chmod a+x bedClip && mv bedClip /usr/local/bin

# Copy the compiled software from the builder
COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin
COPY --chown=$USER:$USER src/python/merge_atac_barcode_metadata.py /usr/local/bin
COPY --chown=$USER:$USER src/python/plot_insert_size_hist.py /usr/local/bin
COPY --chown=$USER:$USER src/python/snapatac2-tss-enrichment.py /usr/local/bin
COPY --chown=$USER:$USER src/python/compute_tss_enrichment.py /usr/local/bin
COPY --chown=$USER:$USER src/python/compute_tss_enrichment_bulk.py /usr/local/bin
COPY --chown=$USER:$USER src/R/atac_qc_plots.R /usr/local/bin
COPY --chown=$USER:$USER src/R/barcode_rank_functions.R /usr/local/bin

USER ${USER}

