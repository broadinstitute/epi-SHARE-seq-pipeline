############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
############################################################

FROM ubuntu@sha256:d1d454df0f579c6be4d8161d227462d69e163a8ff9d20a847533989cf0c94d90

LABEL maintainer="Mei Knudson"

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive
ENV SAMTOOLS_VERSION 1.9

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    autoconf \
    binutils \
    build-essential \
    gcc \
    git \
    libcurl4-openssl-dev \
    libjpeg-dev \
    liblz4-dev \
    liblzma-dev \
    libncurses5-dev \
    libbz2-dev \
    openjdk-11-jre \
    python2.7 \
    python3 \ 
    python-pip \
    python3-pip \
    r-base \
    wget \
    zlib1g-dev &&\
    rm -rf /var/lib/apt/lists/*

# Install packages for python2 scripts
RUN pip2 install --ignore-installed --no-cache-dir matplotlib pysam numpy==1.16.0

# Install packages for python3 scripts
RUN python3 -m pip install --ignore-installed matplotlib numpy==1.19.5 pandas plotnine pysam

# Install Picard 2.20.7
RUN wget https://github.com/broadinstitute/picard/releases/download/2.20.7/picard.jar && chmod +x picard.jar && mv picard.jar /usr/local/bin

# Install samtools 1.9
RUN git clone --branch ${SAMTOOLS_VERSION} --single-branch https://github.com/samtools/samtools.git && \
    git clone --branch ${SAMTOOLS_VERSION} --single-branch https://github.com/samtools/htslib.git && \
    cd samtools && make && make install && cd ../ && rm -rf samtools* htslib*

# Create and setup new user
ENV USER=shareseq
WORKDIR /home/$USER
RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

ENV PYTHONPATH="/usr/local/python:$PYTHONPATH"

COPY --chown=$USER:$USER src/python/make-tss-pileup-jbd.py /usr/local/bin
COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin
COPY --chown=$USER:$USER src/python/plot_insert_size_hist.py /usr/local/bin

USER ${USER}
