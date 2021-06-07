############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on Ubuntu 18.04.3
############################################################

# Set the base image to Ubuntu 20.04.3
#FROM ubuntu:focal
FROM ubuntu@sha256:86ac87f73641c920fb42cc9612d4fb57b5626b56ea2a19b894d0673fd5b4f2e9

MAINTAINER Eugenio Mattei

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    libncurses5-dev libcurl4-openssl-dev libfreetype6-dev zlib1g-dev \
    python python-setuptools python3 python3-setuptools python3-pip \
    git wget unzip ghostscript pkg-config libboost-dev \
    openjdk-8-jre apt-transport-https tabix \
    r-base pigz \ 
    && rm -rf /var/lib/apt/lists/*

# Make directory for all softwares
RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"

# Install system/math python packages (python3)
RUN pip3 install --no-cache-dir common python-Levenshtein==0.12.2 biopython==1.79

# Copy the preprocess python file inside
COPY src/python /software
COPY src/R /software

# Install R packages including spp
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.r-project.org'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e "install.packages('Matrix')" && \
    Rscript -e "install.packages('matrixStats')" && \
    Rscript -e "install.packages('reshape2')" && \
    Rscript -e "install.packages('data.table')" && \
    Rscript -e "install.packages('tidyr')"

# Install Picard 2.20.7
RUN wget https://github.com/broadinstitute/picard/releases/download/2.20.7/picard.jar && chmod +x picard.jar

# Install samtools 1.9
RUN git clone --branch 1.9 --single-branch https://github.com/samtools/samtools.git && \
    git clone --branch 1.9 --single-branch https://github.com/samtools/htslib.git && \
    cd samtools && make && make install && cd ../ && rm -rf samtools* htslib*

# Install Bowtie2 2.3.4.3
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.4/bowtie2-2.4.4-source.zip && \
    unzip bowtie2-2.4.4-source.zip && cd bowtie2-2.4.4 && make static-libs && make STATIC_BUILD=1 && mv bowtie2* .. && \
    cd .. && rm -rf bowtie2-2.4.4*
