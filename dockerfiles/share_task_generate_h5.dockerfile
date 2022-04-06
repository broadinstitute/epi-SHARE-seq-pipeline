############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on Debian slim
############################################################

FROM r-base@sha256:fff003a52d076e963396876b83cfa88c4f40a8bc27e341339cd3cc0236c1db79 as builder

ENV R_VERSION=4.1.2

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
    wget \
    zlib1g-dev &&\
    rm -rf /var/lib/apt/lists/*

ENV R_LIBS_USER=/usr/local/lib/R

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.r-project.org'; r['Ncpus'] = 4; options(repos = r); Sys.setenv('MAKE'='make -j4')" > ~/.Rprofile && Rscript -e "install.packages(c('matrixStats','reshape2','Matrix','dplyr','data.table','BiocManager'));BiocManager::install(c('SummarizedExperiment','rhdf5','DropletUtils'), update=T, ask=F)"

#############################################################
FROM r-base@sha256:fff003a52d076e963396876b83cfa88c4f40a8bc27e341339cd3cc0236c1db79

LABEL maintainer = "Eugenio Mattei"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="generate_h5"

# Create and setup new user
ENV USER=shareseq
WORKDIR /home/$USER
RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

ENV R_LIBS_USER=/usr/local/lib/R

# Copy the compiled software from the builder
COPY --chown=$USER:$USER src/R/UMI_gene_perCell_plot_v2.R /usr/local/bin/
COPY --from=builder --chown=$USER:$USER ${R_LIBS_USER} ${R_LIBS_USER}

USER ${USER}
