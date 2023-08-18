############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on Python slim
############################################################

FROM python@sha256:7ad180fdf785219c4a23124e53745fbd683bd6e23d0885e3554aff59eddbc377

LABEL maintainer = "Eugenio Mattei"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="correct_fastq"

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    pigz && \
    rm -rf /var/lib/apt/lists/*

# Install python packages
RUN python3 -m pip install --no-cache-dir --break-system-packages xopen

# Create and setup new user
ENV USER=shareseq
WORKDIR /home/$USER
RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

# Copy scripts
COPY --chown=$USER:$USER src/python/correct_fastq.py /usr/local/bin/
COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin

USER ${USER}
