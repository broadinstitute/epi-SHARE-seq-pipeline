############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on Ubuntu 18.04.3
############################################################

# Set the base image to Ubuntu 18.04.3
#FROM ubuntu:focal
FROM ubuntu@sha256:d1d454df0f579c6be4d8161d227462d69e163a8ff9d20a847533989cf0c94d90 as builder

LABEL maintainer = "Neva Durand"

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive
ENV SAMTOOLS_VERSION 1.9

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    libncurses5-dev \
    libcurl4-openssl-dev \
    zlib1g-dev \
    liblzma-dev \
    libbz2-dev \
    python3 \
    python3-setuptools \
    python3-pip \
    git \
    wget \
    xmlstarlet \
    openjdk-8-jre \
    && rm -rf /var/lib/apt/lists/*

# Install samtools 1.9
RUN git clone --branch ${SAMTOOLS_VERSION} --single-branch https://github.com/samtools/samtools.git && \
    git clone --branch ${SAMTOOLS_VERSION} --single-branch https://github.com/samtools/htslib.git && \
    cd samtools && make && make install && cd ../ && rm -rf samtools* htslib*

# Install system/math python packages (python3)
RUN pip3 install --ignore-installed -t /usr/local/python --no-cache-dir python-Levenshtein==0.12.2 pysam

# Install Picard 2.26.11
RUN wget https://github.com/broadinstitute/picard/releases/download/2.26.11/picard.jar && \
    chmod +x picard.jar && \
    mv picard.jar /usr/local/bin


#############################################################
FROM google/cloud-sdk@sha256:29dfab2e106fbd68ce826db03437c5f9ee8b154bb1c3f8c5eaf104c71c2883f1

LABEL maintainer = "Eugenio Mattei"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="preprocess"

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install --no-install-recommends -y \
    python3 \
    openjdk-11-jre &&\
    rm -rf /var/lib/apt/lists/*

# Create and setup new user
ENV USER=shareseq
WORKDIR /home/$USER
RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

WORKDIR /software
ENV PATH="/software:${PATH}"
ENV PYTHONPATH="/usr/local/python:$PYTHONPATH"

# Copy the compiled software from the builder
COPY --from=builder --chown=$USER:$USER /usr/local/bin/* /usr/local/bin/
COPY --from=builder --chown=$USER:$USER /usr/local/python/ /usr/local/python
COPY --chown=$USER:$USER src/python/bam_fastq.py /software
COPY --chown=$USER:$USER src/python/write_terra_tables.py /software
COPY --chown=$USER:$USER src/python/flexible_import_entities_standard.py /software

USER ${USER}
