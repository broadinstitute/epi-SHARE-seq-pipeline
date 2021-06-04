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
COPY ../src/python/fastq.process.py3.v0.6.py /software
