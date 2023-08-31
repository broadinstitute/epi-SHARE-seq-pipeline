############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on Debian slim
############################################################

FROM r-base@sha256:fff003a52d076e963396876b83cfa88c4f40a8bc27e341339cd3cc0236c1db79 as builder

LABEL maintainer = "Siddarth Wekhande"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="archr"

RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site

ENV R_LIBS_USER=/usr/local/lib/R
ENV RETICULATE_MINICONDA_ENABLED=FALSE

RUN apt-get update -qq && \
    apt-get install -y -qq --no-install-recommends\
    binutils \
    gtk-doc-tools \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libgsl-dev \
    libharfbuzz-dev \
    libhdf5-dev \
    libjpeg-dev \
    libmpfr-dev \
    libpng-dev \
    libssl-dev \
    libtiff5-dev \
    libxml2-dev \
    libxt-dev \
    libmagick++-dev \
    libgeos-dev \
    meson \
    python3 \
    python3-pip && \
    rm -rf /var/lib/apt/lists/*

ENV USER=shareseq
WORKDIR /home/$USER

RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

RUN R --no-echo --no-restore --no-save -e "install.packages(c('devtools','hdf5r','IRkernel','BiocManager','Cairo','magick'))"

RUN R --no-echo --no-restore --no-save -e "BiocManager::install(c('GenomeInfoDbData','GenomicRanges','Rsamtools'), update=F, ask=F)"

RUN R --no-echo --no-restore --no-save -e "devtools::install_github('GreenleafLab/ArchR@v1.0.1', repos = BiocManager::repositories());ArchR::installExtraPackages()"

RUN R --no-echo --no-restore --no-save -e "devtools::install_github('immunogenomics/presto')"

RUN R --no-echo --no-restore --no-save -e "remotes::install_version('Seurat', version = '4.3.0')"

RUN R --no-echo --no-restore --no-save -e "install.packages(c('logr','hexbin', 'ggpointdensity'))"

RUN python3 -m pip install --break-system-packages jupyter papermill

COPY src/jupyter_nb/archr_notebook.ipynb /usr/local/bin/

COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin


RUN R -e "IRkernel::installspec()"
