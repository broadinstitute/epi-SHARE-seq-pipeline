############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on python slim
############################################################

FROM python@sha256:7ad180fdf785219c4a23124e53745fbd683bd6e23d0885e3554aff59eddbc377

LABEL maintainer = "Eugenio Mattei"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="generate_h5"

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive

# Install python packages
RUN python3 -m pip install --no-cache-dir h5py scipy==1.10.0

# Create and setup new user
ENV USER=shareseq
WORKDIR /home/$USER
RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

# Copy scripts
COPY --chown=$USER:$USER src/python/generate_h5_rna.py /usr/local/bin/
COPY --chown=$USER:$USER src/python/merge_rna_counts.py /usr/local/bin/
COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin

USER ${USER}
