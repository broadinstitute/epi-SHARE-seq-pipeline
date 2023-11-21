############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on Ubuntu latest
############################################################

FROM ubuntu:latest

LABEL maintainer = "Eugenio Mattei"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="make-track"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    bc \
    build-essential \
    libkrb5-3 \
    libcurl4 \
    libbz2-dev \
    liblz4-dev \
    liblzma-dev \
    libncurses5-dev \
    python3 \
    pigz \
    tabix \
    wget \
    zlib1g-dev &&\
    rm -rf /var/lib/apt/lists/*

RUN wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && chmod +x bedGraphToBigWig && mv bedGraphToBigWig /usr/local/bin
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && tar -zxvf bedtools-2.31.1.tar.gz && cd bedtools2 && make && make install
#RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static && mv bedtools.static bedtools && chmod a+x bedtools && mv bedtools /usr/local/bin


# Create and setup new user
ENV USER=shareseq
WORKDIR /home/$USER

RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin
COPY --chown=$USER:$USER src/bash/make_track.sh /usr/local/bin


USER ${USER}





