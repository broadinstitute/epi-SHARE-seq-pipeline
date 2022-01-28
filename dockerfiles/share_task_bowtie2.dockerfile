############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on Debian slim
############################################################

FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f as builder

ENV BOWTIE2_VERSION 2.4.3

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    build-essential \
    cpanminus \
    liblz4-dev \
    liblzma-dev \
    unzip \
    wget \
    zlib1g-dev &&\
    rm -rf /var/lib/apt/lists/*


# Make directory for all softwares
RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"

RUN cpanm Sys::Hostname


# Install Bowtie2 2.3.4.3
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-source.zip && \
    unzip bowtie2-${BOWTIE2_VERSION}-source.zip && cd bowtie2-${BOWTIE2_VERSION} && make static-libs && make STATIC_BUILD=1 && \
    cp bowtie2* .. && \
    cd .. && rm -rf bowtie2-${BOWTIE2_VERSION}*

FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f

LABEL maintainer = "Eugenio Mattei"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="Align ATAC using Bowtie2"

# Create and setup new user
ENV USER=shareseq
WORKDIR /home/$USER

RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

# Add folder with software to the path
ENV PATH="/software:${PATH}"

# Copy the compiled software from the builder
COPY --from=builder --chown=$USER:$USER /software/bowtie2* /software/
COPY --from=builder /usr/lib/x86_64-linux-gnu/perl/5.32 /usr/lib/x86_64-linux-gnu/perl/5.32/



USER $USER
