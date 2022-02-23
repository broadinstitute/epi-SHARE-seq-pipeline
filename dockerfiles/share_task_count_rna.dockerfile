############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on Debian slim
############################################################

FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f as builder

ENV SAMTOOLS_VERSION 1.9
ENV FEATURECOUNTS_VERSION 2.0.2

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    liblz4-dev \
    liblzma-dev \
    libncurses5-dev \
    libbz2-dev \
    wget \
    zlib1g-dev &&\
    rm -rf /var/lib/apt/lists/*


# Make directory for all softwares
RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"


# Install samtools 1.9
RUN git clone --branch ${SAMTOOLS_VERSION} --single-branch https://github.com/samtools/samtools.git && \
    git clone --branch ${SAMTOOLS_VERSION} --single-branch https://github.com/samtools/htslib.git && \
    cd samtools && make && make install && cd ../ && rm -rf samtools* htslib*

# Install featureCounts
RUN wget https://sourceforge.net/projects/subread/files/subread-${FEATURECOUNTS_VERSION}/subread-${FEATURECOUNTS_VERSION}-source.tar.gz && tar -zxf subread-${FEATURECOUNTS_VERSION}-source.tar.gz
RUN cd subread-${FEATURECOUNTS_VERSION}-source/src && make -f Makefile.Linux && mv ../bin/* /usr/local/bin/ && cd ../.. && rm -r subread-${FEATURECOUNTS_VERSION}-source

FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f

LABEL maintainer = "Eugenio Mattei"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="count_rna"

# Create and setup new user
ENV USER=shareseq
WORKDIR /home/$USER

RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

# Add folder with software to the path
ENV PATH="/software:${PATH}"

# Copy the compiled software from the builder
COPY --from=builder --chown=$USER:$USER /usr/local/bin/* /usr/local/bin/
COPY --from=builder --chown=$USER:$USER /lib/x86_64-linux-gnu/libncurses.so.6 /lib/x86_64-linux-gnu/


USER $USER
