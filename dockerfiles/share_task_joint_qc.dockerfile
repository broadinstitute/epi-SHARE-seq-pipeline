############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
############################################################

#FROM ubuntu@sha256:d1d454df0f579c6be4d8161d227462d69e163a8ff9d20a847533989cf0c94d90
FROM python:3.8-buster@sha256:7e7f4c5508b85268a93b573566c8eb321a6fdb466e3b60c663a42300c73a7400

LABEL maintainer="Mei Knudson"

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive
ENV SAMTOOLS_VERSION 1.9

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    r-base  &&\
    rm -rf /var/lib/apt/lists/*

# Install packages for python3 scripts
RUN python3 -m pip install matplotlib numpy pandas plotnine

# Install packages for R scripts
RUN R -e "install.packages(c('ggplot2', 'remotes'))"
RUN R -e "remotes::install_github('LKremer/ggpointdensity')"

# Create and setup new user
ENV USER=shareseq
WORKDIR /home/$USER
RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

ENV PYTHONPATH="/usr/local/python:$PYTHONPATH"
ENV R_LIBS_USER=/usr/local/lib/R

COPY --chown=$USER:$USER src/python/joint_cell_plotting.py /usr/local/bin
COPY --chown=$USER:$USER src/R/joint_cell_plotting_density.R /usr/local/bin
COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin

USER ${USER}
