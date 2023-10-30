############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on Debian slim
############################################################

FROM ubuntu@sha256:2fdb1cf4995abb74c035e5f520c0f3a46f12b3377a59e86ecca66d8606ad64f9

LABEL maintainer = "Zhijian Li"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="cell-annotation"

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive
ENV RETICULATE_MINICONDA_ENABLED=FALSE

## Create new user 
ENV USER=shareseq
WORKDIR /home/$USER
RUN groupadd -r $USER && \
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

# Install libraries
RUN apt-get update
RUN apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    gfortran \
    patch \
    build-essential \
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
    libgeos-dev \
    meson \
    libblas-dev \
    liblapack-dev \
    libbz2-dev
    
# Install python and R
RUN apt-get install -y --no-install-recommends \
    python3 python3-pip python3-dev python3-venv r-base

RUN rm -rf /var/lib/apt/lists/*

RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site
ENV R_LIBS_USER=/usr/local/lib/R

RUN R --no-echo --no-restore --no-save -e "install.packages('hdf5r')"
RUN R --no-echo --no-restore --no-save -e "install.packages('remotes')"
RUN R --no-echo --no-restore --no-save -e "install.packages('IRkernel')"
RUN R --no-echo --no-restore --no-save -e "install.packages('logr')"
RUN R --no-echo --no-restore --no-save -e "install.packages('BiocManager')"
RUN R --no-echo --no-restore --no-save -e "install.packages('glue')"
RUN R --no-echo --no-restore --no-save -e "install.packages('Matrix')"
#RUN R --no-echo --no-restore --no-save -e "install.packages('SeuratObject')"
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('bnprks/BPCells')"
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('satijalab/seurat', 'seurat5')"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install('rhdf5', update=F, ask=F)"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install('EnsDb.Mmusculus.v79', update=F, ask=F)"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install('EnsDb.Hsapiens.v86', update=F, ask=F)"
RUN R --no-echo --no-restore --no-save -e "install.packages('optparse')"
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('mojaveazure/seurat-object', 'seurat5')"
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('stuart-lab/signac', 'seurat5')"
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('satijalab/azimuth', 'seurat5')"
RUN R --no-echo --no-restore --no-save -e "install.packages('sctransform')"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install('glmGamPoi', update=F, ask=F)"

RUN python3 -m pip install anndata cellxgene-census

COPY src/bash/monitor_script.sh /usr/local/bin
COPY src/python/get_cellxgene_data.py /usr/local/bin
COPY src/R/cell_annotation.R /usr/local/bin/
COPY src/R/h5ad_to_seurat.R /usr/local/bin/

