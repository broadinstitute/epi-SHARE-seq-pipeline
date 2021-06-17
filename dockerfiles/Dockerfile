############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on Ubuntu 18.04.3
############################################################

# Set the base image to Ubuntu 18.04.3
#FROM ubuntu:focal
FROM ubuntu@sha256:d1d454df0f579c6be4d8161d227462d69e163a8ff9d20a847533989cf0c94d90

MAINTAINER Eugenio Mattei

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    libncurses5-dev libcurl4-openssl-dev libfreetype6-dev zlib1g-dev \
    python python-setuptools python-pip python3 python3-setuptools python3-pip \
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

# Install packages for python2 scripts
RUN pip2 install --no-cache-dir numpy matplotlib pysam

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
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.3/bowtie2-2.4.3-source.zip && \
    unzip bowtie2-2.4.3-source.zip && cd bowtie2-2.4.3 && make static-libs && make STATIC_BUILD=1 && \
    cp bowtie2* .. && \
    cd .. && rm -rf bowtie2-2.4.3*

# Install bedtools 2.29.0
RUN git clone --branch v2.29.0 --single-branch https://github.com/arq5x/bedtools2.git && \
    cd bedtools2 && make && make install && cd ../ && rm -rf bedtools2*

# Install STAR 2.5.1b
RUN wget https://github.com/alexdobin/STAR/archive/2.5.1b.tar.gz && tar -xzf 2.5.1b.tar.gz
RUN cd STAR-2.5.1b && make STAR && rm ../2.5.1b.tar.gz
ENV PATH="/software/STAR-2.5.1b/bin/Linux_x86_64:${PATH}"

# Install featureCounts
RUN wget https://sourceforge.net/projects/subread/files/subread-2.0.2/subread-2.0.2-source.tar.gz && tar -zxf subread-2.0.2-source.tar.gz
RUN cd subread-2.0.2-source/src && make -f Makefile.Linux && mv ../bin/* /software/ && cd ../.. && rm -r subread-2.0.2-source

# Copy the external scripts inside
COPY src/python /software
COPY src/R /software
