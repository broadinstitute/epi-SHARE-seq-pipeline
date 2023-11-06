############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on Ubuntu 18.04.3
############################################################

# Set the base image to Ubuntu 18.04.3
#FROM ubuntu:focal
FROM ubuntu@sha256:d1d454df0f579c6be4d8161d227462d69e163a8ff9d20a847533989cf0c94d90

MAINTAINER Neva Durand

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \ 
    libncurses5-dev libcurl4-openssl-dev zlib1g-dev liblzma-dev libbz2-dev \
    python3 python3-setuptools python3-pip \
    git wget xmlstarlet \
    openjdk-8-jre \
    && rm -rf /var/lib/apt/lists/*

# Make directory for all softwares
RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"

# Install samtools 1.9
RUN git clone --branch 1.9 --single-branch https://github.com/samtools/samtools.git && \
    git clone --branch 1.9 --single-branch https://github.com/samtools/htslib.git && \
    cd samtools && make && make install && cd ../ && rm -rf samtools* htslib*

# Install system/math python packages (python3)
RUN pip3 install --no-cache-dir python-Levenshtein==0.12.2 pysam requests oauth2client

# Install Picard 2.26.11
RUN wget https://github.com/broadinstitute/picard/releases/download/2.26.11/picard.jar && chmod +x picard.jar

# Install gsutil
# Downloading gcloud package
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-453.0.0-linux-x86_64.tar.gz

# Installing the package
RUN tar -xf google-cloud-cli-453.0.0-linux-x86_64.tar.gz && \
    ./google-cloud-sdk/install.sh && \
    rm google-cloud-cli-453.0.0-linux-x86_64.tar.gz

# Adding the package path to local
ENV PATH $PATH:/software/google-cloud-sdk/bin

# Copy the external scripts inside
COPY src/python/bam_to_raw_fastq.py /software
COPY src/python/flexible_import_entities_standard.py /software
COPY src/python/write_terra_tables.py /software
COPY src/bash/monitor_script.sh /software