############################################################
# Dockerfile for BROAD GRO combinomics pipeline
# Based on Python
############################################################

FROM python@sha256:fd0fa50d997eb56ce560c6e5ca6a1f5cf8fdff87572a16ac07fb1f5ca01eb608

LABEL maintainer="Eugenio Mattei"
LABEL software="Combinomics pipeline"
LABEL software.version="2.0.0.rc"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="qc-atac"
LABEL software.description="Quality control for ATAC-seq data"

# Install the required packages
RUN pip install --upgrade pip
# RUN pip install click matplotlib numpy pandas snapatac2[recommend]

RUN mkdir /software
COPY src/python/qc_atac /software
RUN cd /software && pip install --editable .

# Create and setup new user
ENV USER=combinomics
WORKDIR /home/$USER

RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

# Copy the compiled software from the builder
COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin
#COPY --chown=$USER:$USER src/python/snapatac2_qc_generation.py /usr/local/bin

#RUN chmod +x /usr/local/bin/snapatac2_qc_generation.py

USER ${USER}

