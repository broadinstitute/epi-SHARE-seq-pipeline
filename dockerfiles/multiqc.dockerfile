FROM ubuntu:jammy

ARG FASTQC_VER="0.12.1"

# install dependencies; cleanup apt garbage
RUN apt-get update && apt-get install -y --no-install-recommends \
    unzip \
    python3-pip \
    wget \
    perl \
    default-jre \
    ca-certificates \
    procps && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

# Install FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VER}.zip && \
    unzip fastqc_v${FASTQC_VER}.zip && \
    rm fastqc_v${FASTQC_VER}.zip && \
    chmod +x FastQC/fastqc && \
    mkdir /common


# Install MultiQC
RUN pip install multiqc

# set PATH and working directory
ENV USER=combinomics
WORKDIR /home/$USER

RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

ENV PATH="${PATH}:/FastQC/"

# Copy SHARE-seq specific contaminants TSV file.
COPY --chown=$USER:$USER common/share_contaminants.tsv /common/share_contaminants.tsv

# Default command
CMD ["bash"]


