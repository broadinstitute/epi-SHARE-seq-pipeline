############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on Debian slim
############################################################

FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f as builder

ENV SAMBAMBA_VERSION 0.6.6
ENV PICARD_VERSION 2.27.5

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    autoconf \
    build-essential \
    git \
    libcurl4-openssl-dev \
    liblz4-dev \
    liblzma-dev \
    libncurses5-dev \
    libbz2-dev \
    python \
    unzip \
    wget \
    zlib1g-dev &&\
    rm -rf /var/lib/apt/lists/*


# Make directory for all softwares
RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"

# Install sambamba 0.6.6
RUN wget https://github.com/lomereiter/sambamba/releases/download/v${SAMBAMBA_VERSION}}/sambamba_v0.6.6_linux.tar.bz2 && \
    tar -xvjf sambamba_v${SAMBAMBA_VERSION}}_linux.tar.bz2 && \
    mv sambamba_v${SAMBAMBA_VERSION}} /usr/local/bin/sambamba && \
    rm -rf sambamba_*

# Install Picard 2.20.7
RUN wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar && chmod +x picard.jar && mv picard.jar /usr/local/bin

FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f

LABEL maintainer = "Eugenio Mattei"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="merge"

RUN apt-get update && apt-get install -y \
    openjdk-11-jre &&\
    rm -rf /var/lib/apt/lists/*

# Create and setup new user
ENV USER=shareseq
WORKDIR /home/$USER

RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

# Add folder with software to the path
ENV PATH="/software:${PATH}"

# Copy the compiled software from the builder
COPY --from=builder --chown=$USER:$USER /usr/local/bin/* /usr/local/bin/
COPY --from=builder --chown=$USER:$USER /lib/x86_64-linux-gnu/* /lib/x86_64-linux-gnu/
COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin


USER ${USER}


