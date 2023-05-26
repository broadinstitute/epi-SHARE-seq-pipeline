############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on Debian slim
############################################################

FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f as builder

ENV BEDTOOLS_VERSION v2.29.0
ENV PICARD_VERSION 2.27.5
ENV SAMTOOLS_VERSION 1.16
ENV SAMBAMBA_VERSION 0.6.6

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    autoconf \
    automake \
    build-essential \
    git \
    libcurl4-openssl-dev \
    liblz4-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    perl \
    python \
    unzip \
    xz-utils \
    wget \
    zlib1g-dev &&\
    rm -rf /var/lib/apt/lists/*

# Make directory for all softwares
RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"

# Install bedtools 2.29.0
RUN git clone --branch ${BEDTOOLS_VERSION} --single-branch https://github.com/arq5x/bedtools2.git && \
    cd bedtools2 && make && make install && cd ../ && rm -rf bedtools2*

# Install sambamba 1.32.2
#RUN wget https://github.com/ldc-developers/ldc/releases/download/v${SAMBAMBA_VERSION}/ldc2-${SAMBAMBA_VERSION}-linux-x86_64.tar.xz && \
#    tar -xvf ldc2-${SAMBAMBA_VERSION}-linux-x86_64.tar.xz && \
#    echo "export LIBRARY_PATH=${PWD}/ldc2-1.7.0-linux-x86_64/lib" >> ${HOME}/.bashrc && \
#    echo "export LIBRARY_PATH=${PWD}/ldc2-1.7.0-linux-x86_64/lib" >> ${HOME}/.bashrc && \
#    source ${HOME}/.bashrc \
#    git clone --branch v${SAMBAMBA_VERSION} --recursive https://github.com/biod/sambamba.git && \
#    cd sambamba && make && make install && cd .. &&  rm -rf sambamba_* && rm -rf ldc2-${SAMBAMBA_VERSION}-linux-x86_64

# Install sambamba 0.6.6
RUN wget https://github.com/lomereiter/sambamba/releases/download/v${SAMBAMBA_VERSION}/sambamba_v${SAMBAMBA_VERSION}_linux.tar.bz2 && \
    tar -xvjf sambamba_v${SAMBAMBA_VERSION}_linux.tar.bz2 && \
    mv sambamba_v${SAMBAMBA_VERSION} /usr/local/bin/sambamba && \
    rm -rf sambamba_*

# Install samtools 1.16
RUN git clone --branch ${SAMTOOLS_VERSION} --single-branch https://github.com/samtools/htslib.git && \
    cd htslib && git submodule update --init --recursive && autoreconf -i && make && make install && cd ../ && \
    git clone --branch ${SAMTOOLS_VERSION} --single-branch https://github.com/samtools/samtools.git && \
    cd samtools && make && make install && cd ../ && rm -rf samtools* && rm -rf htslib*
    

# Install Picard 2.20.7
RUN wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar && chmod +x picard.jar && mv picard.jar /usr/local/bin



FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f

LABEL maintainer = "Eugenio Mattei"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="filter"

RUN apt-get update && apt-get install -y \
    gcc \
    libcurl4-openssl-dev \
    python3 \
    python3-dev \
    python3-pip \
    openjdk-11-jre \
    zlib1g-dev &&\
    rm -rf /var/lib/apt/lists/*

# Install packages for python3 scripts
RUN python3 -m pip install --no-cache-dir --ignore-installed pysam

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
COPY --chown=$USER:$USER src/python/filter_mito_reads.py /usr/local/bin
COPY --chown=$USER:$USER src/python/bam_to_fragments.py /usr/local/bin
COPY --chown=$USER:$USER src/python/assign_multimappers.py /usr/local/bin





USER ${USER}


