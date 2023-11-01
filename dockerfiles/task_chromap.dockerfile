############################################################
# Dockerfile for chromap
# Based on Debian slim
############################################################

FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f as builder

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    liblz4-dev \
    liblzma-dev \
    libncurses5-dev \
    libbz2-dev \
    unzip \
    wget \
    zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

# Make directory for all softwares
RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"

RUN wget https://github.com/haowenz/chromap/archive/refs/heads/master.zip && unzip master.zip && \ 
    cd chromap-master && make STATIC_BUILD=1

FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f

LABEL maintainer = "Eugenio Mattei"
LABEL software = "combinomics Single Cell pipeline"
LABEL software.version="2.0.0"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="Yes"
LABEL software.task="chromap"

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    tabix \
    gcc \
    pigz && \
    rm -rf /var/lib/apt/lists/*

# Create and setup new user
ENV USER=shareseq
WORKDIR /home/$USER

RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

# Add folder with software to the path
ENV PATH="/software/:${PATH}"

# Copy the compiled software from the builder
COPY --from=builder --chown=$USER:$USER /software/chromap-master /software/
COPY --from=builder --chown=$USER:$USER /usr/include/* /usr/include/
COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin

USER $USER
