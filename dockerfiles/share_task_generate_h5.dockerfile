############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on Debian slim
############################################################

FROM python@sha256:7ad180fdf785219c4a23124e53745fbd683bd6e23d0885e3554aff59eddbc377

LABEL maintainer = "Eugenio Mattei"
LABEL software = "Share-seq pipeline"
LABEL software.version="1.0.0"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="Yes"
LABEL software.task="generate_h5"

# Install python packages
RUN pip install --no-cache-dir h5py scipy

# Create and setup new user
ENV USER=shareseq
WORKDIR /home/$USER
RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

# Copy scripts
COPY --chown=$USER:$USER src/python/generate_h5_rna.py /usr/local/bin/
COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin

USER ${USER}
