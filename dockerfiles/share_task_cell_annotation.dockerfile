############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on Debian slim
############################################################

FROM r-base@sha256:fff003a52d076e963396876b83cfa88c4f40a8bc27e341339cd3cc0236c1db79 as builder

LABEL maintainer = "Zhijian Li"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="cell-annotation"

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
    libgeos-dev \
    meson \
    pkg-config \
    python3 \
    python3-pip && \
    rm -rf /var/lib/apt/lists/*

ENV USER=shareseq
WORKDIR /home/$USER

RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

RUN R --no-echo --no-restore --no-save -e "install.packages('hdf5r')"
RUN R --no-echo --no-restore --no-save -e "install.packages('remotes')"
RUN R --no-echo --no-restore --no-save -e "install.packages('IRkernel')"
RUN R --no-echo --no-restore --no-save -e "install.packages('logr')"
RUN R --no-echo --no-restore --no-save -e "install.packages('BiocManager')"
RUN R --no-echo --no-restore --no-save -e "install.packages('glue')"
RUN R --no-echo --no-restore --no-save -e "install.packages('Matrix')"
RUN R --no-echo --no-restore --no-save -e "install.packages('SeuratObject')"
RUN R --no-echo --no-restore --no-save -e "remotes::install_version('Seurat', version = '4.1.1')"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install('rhdf5', update=F, ask=F)"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install('EnsDb.Mmusculus.v79', update=F, ask=F)"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install('EnsDb.Hsapiens.v86', update=F, ask=F)"
RUN R --no-echo --no-restore --no-save -e "install.packages('devtools')"
RUN R --no-echo --no-restore --no-save -e "devtools::install_github('chanzuckerberg/cellxgene-census/api/r/cellxgene.census')"


COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin

RUN python3 -m pip install --break-system-packages jupyter papermill anndata

COPY src/jupyter_nb/cell_annotation_notebook.ipynb /usr/local/bin/

RUN R -e "IRkernel::installspec()"
