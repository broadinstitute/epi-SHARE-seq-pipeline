# ############################################################
# # Dockerfile for BROAD GRO share-seq-pipeline
# # Based on Debian slim
# ############################################################

FROM debian@sha256:8d498c9133965638a6c161f541829352e4a9907969e6b0fd3f9efa1b3acae80b

LABEL maintainer = "Zhijian Li"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="cell-annotation"

## install R 4.3.0
ENV R_BASE_VERSION 4.3.0

## Use Debian unstable via pinning -- new style via APT::Default-Release
RUN echo "deb http://http.debian.net/debian sid main" > /etc/apt/sources.list.d/debian-unstable.list

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libopenblas0-pthread \
    r-base=${R_BASE_VERSION}-* \
    r-base-dev=${R_BASE_VERSION}-* \
    r-base-core=${R_BASE_VERSION}-* \
    r-recommended=${R_BASE_VERSION}-* \
	&& chown root:staff "/usr/local/lib/R/site-library" \
	&& chmod g+ws "/usr/local/lib/R/site-library"

## Create new user 
RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site
ENV R_LIBS_USER=/usr/local/lib/R
ENV USER=shareseq
WORKDIR /home/$USER

RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

## Install other libraries
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    cmake \
    git \
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

RUN R --no-echo --no-restore --no-save -e "install.packages('hdf5r')"
RUN R --no-echo --no-restore --no-save -e "install.packages('tiledb')"
RUN R --no-echo --no-restore --no-save -e "install.packages('arrow')"
RUN R --no-echo --no-restore --no-save -e "install.packages('fs')"
RUN R --no-echo --no-restore --no-save -e "install.packages('urltools')"
RUN R --no-echo --no-restore --no-save -e "install.packages('dplyr')"
RUN R --no-echo --no-restore --no-save -e "install.packages('data.table')"
RUN R --no-echo --no-restore --no-save -e "install.packages('remotes')"
RUN R --no-echo --no-restore --no-save -e "install.packages('knitr')"

## Install TileDB-SOMA
RUN git clone https://github.com/single-cell-data/TileDB-SOMA.git
RUN cd TileDB-SOMA
RUN scripts/bld --prefix=/usr/local
#RUN cd apis/r
#RUN ./configure && R CMD INSTALL .
#RUN cd /home/$USER

#RUN R --no-echo --no-restore --no-save -e "devtools::install_github('single-cell-data/TileDB-SOMA/apis/r')"

# # RUN R --no-echo --no-restore --no-save -e "install.packages('IRkernel')"
# # RUN R --no-echo --no-restore --no-save -e "install.packages('logr')"
# # RUN R --no-echo --no-restore --no-save -e "install.packages('BiocManager')"
# # RUN R --no-echo --no-restore --no-save -e "install.packages('glue')"
# # RUN R --no-echo --no-restore --no-save -e "install.packages('Matrix')"
# # RUN R --no-echo --no-restore --no-save -e "install.packages('SeuratObject')"
# # RUN R --no-echo --no-restore --no-save -e "devtools::install_version('Seurat', version = '4.3.0')"
# # RUN R --no-echo --no-restore --no-save -e "BiocManager::install('rhdf5', update=F, ask=F)"
# # RUN R --no-echo --no-restore --no-save -e "BiocManager::install('EnsDb.Mmusculus.v79', update=F, ask=F)"
# # RUN R --no-echo --no-restore --no-save -e "BiocManager::install('EnsDb.Hsapiens.v86', update=F, ask=F)"
# # RUN R --no-echo --no-restore --no-save -e "install.packages('anndata')"

# # COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin

# # RUN python3 -m pip install --break-system-packages jupyter papermill anndata

# # COPY src/jupyter_nb/cell_annotation_notebook.ipynb /usr/local/bin/
# # COPY src/R/cell_annotation_helper_functions.R /usr/local/bin/

# # RUN R -e "IRkernel::installspec()"
