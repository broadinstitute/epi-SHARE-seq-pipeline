############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on Debian slim
############################################################

FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f

LABEL maintainer = "Eugenio Mattei"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="merge_fragments"

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    autoconf \
    automake \
    build-essential \
    libcurl4-openssl-dev \
    liblz4-dev \
    liblzma-dev \
    libncurses5-dev \
    libbz2-dev \
    zlib1g-dev &&\
    rm -rf /var/lib/apt/lists/*

# Create and setup new user
ENV USER=shareseq
WORKDIR /home/$USER
RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

# Copy scripts
COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin

USER ${USER}
