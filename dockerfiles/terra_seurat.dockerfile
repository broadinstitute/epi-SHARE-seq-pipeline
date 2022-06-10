############################################################
# Dockerfile for Terra to support Seurat
# Based on Debian slim
############################################################

FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-r:2.1.3

LABEL maintainer = "Siddarth Wekhande"
LABEL software = "Seurat on Terra"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="seurat"

USER root

RUN R --no-echo --no-restore --no-save -e "install.packages(c('hdf5r','remotes'))"

RUN R --no-echo --no-restore --no-save -e "remotes::install_version('Seurat', version = '4.1.1')"

ENV USER jupyter
USER $USER 

ENTRYPOINT ["/bin/bash"]