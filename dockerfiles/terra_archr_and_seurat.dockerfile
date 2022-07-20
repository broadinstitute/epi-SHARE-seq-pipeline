############################################################
# Dockerfile for Terra to support ArchR
# Based on Debian slim
############################################################

FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-r:2.1.3

LABEL maintainer = "Siddarth Wekhande"
LABEL software = "ArchR on Terra"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="archr"

USER root

RUN R --no-echo --no-restore --no-save -e "install.packages(c('hdf5r','remotes'))"

RUN R --no-echo --no-restore --no-save -e "remotes::install_version('Seurat', version = '4.1.1')"

RUN R --no-echo --no-restore --no-save -e "remotes::install_github('GreenleafLab/ArchR@v1.0.1', repos = BiocManager::repositories());ArchR::installExtraPackages()"

RUN R --no-echo --no-restore --no-save -e "remotes::install_github('immunogenomics/presto')"

ENV USER jupyter
USER $USER 

ENTRYPOINT ["/bin/bash"]