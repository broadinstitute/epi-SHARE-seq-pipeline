############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
############################################################

FROM r-base@sha256:fff003a52d076e963396876b83cfa88c4f40a8bc27e341339cd3cc0236c1db79

LABEL maintainer = "Eugenio Mattei"
LABEL software = "combinomics pipeline"
LABEL software.version="2.0.0"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="joint-qc"

RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site

ENV R_LIBS_USER=/usr/local/lib/R

RUN apt-get update -qq && \
    apt-get install -y -qq --no-install-recommends\
    binutils \
    gtk-doc-tools \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libfreetype-dev \
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

# Install R packages
RUN R --no-echo --no-restore --no-save -e "install.packages(c('ggplot2','remotes'))"
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('LKremer/ggpointdensity')"
# Install python packages
RUN python3 -m pip install --break-system-packages matplotlib numpy pandas plotnine

COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin
COPY --chown=$USER:$USER src/python/joint_cell_plotting.py /usr/local/bin
COPY --chown=$USER:$USER src/R/joint_cell_plotting_density.R /usr/local/bin
COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin

USER ${USER}