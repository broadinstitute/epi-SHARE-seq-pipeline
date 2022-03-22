############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on Debian slim
############################################################

FROM r-base@sha256:fff003a52d076e963396876b83cfa88c4f40a8bc27e341339cd3cc0236c1db79 as builder

ENV R_VERSION=4.1.2

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive
ENV SAMTOOLS_VERSION 1.9

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
    libbz2-dev \
    wget \
    zlib1g-dev &&\
    rm -rf /var/lib/apt/lists/*

# Install htslib 1.9
RUN git clone --branch ${SAMTOOLS_VERSION} --single-branch https://github.com/samtools/htslib.git && \
    cd htslib && autoreconf -i && make && make install && cd ../ && rm -rf htslib*

ENV R_LIBS_USER=/usr/local/lib/R

# Install Tidyverse
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.r-project.org'; r['Ncpus'] = 4; options(repos = r); Sys.setenv('MAKE'='make -j4')" > ~/.Rprofile && Rscript -e "install.packages(c('tidyr','tibble'))"

#############################################################
FROM r-base@sha256:fff003a52d076e963396876b83cfa88c4f40a8bc27e341339cd3cc0236c1db79

LABEL maintainer = "Eugenio Mattei"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="count_atac"

# Create and setup new user
ENV USER=shareseq
WORKDIR /home/$USER
RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

ENV R_LIBS_USER=/usr/local/lib/R


# Copy the compiled software from the builder
COPY --chown=$USER:$USER src/R/sum_reads.R /usr/local/bin
COPY --from=builder --chown=$USER:$USER /usr/local/bin/* /usr/local/bin/
COPY --from=builder --chown=$USER:$USER ${R_LIBS_USER} ${R_LIBS_USER}

USER ${USER}
