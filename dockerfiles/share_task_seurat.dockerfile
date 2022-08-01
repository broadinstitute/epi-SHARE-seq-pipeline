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
LABEL software.task="seurat"

RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site

ENV R_LIBS_USER=/usr/local/lib/R
ENV RETICULATE_MINICONDA_ENABLED=FALSE

RUN apt-get update -qq && \
	apt-get install -y -qq --no-install-recommends\
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
    libgeos-dev \
    meson \
    pkg-config \
    python3 \
    python3-pip && \
    rm -rf /var/lib/apt/lists/*
        
RUN R --no-echo --no-restore --no-save -e "install.packages(c('hdf5r','remotes','IRkernel','logr', 'BiocManager'))"

RUN R --no-echo --no-restore --no-save -e "remotes::install_version('Seurat', version = '4.1.1')"

RUN R --no-echo --no-restore --no-save -e "BiocManager::install(c('glmGamPoi'), update=F, ask=F)"

RUN python3 -m pip install jupyter papermill

COPY src/jupyter_nb/seurat_notebook.ipynb /usr/local/bin/

RUN R -e "IRkernel::installspec()"