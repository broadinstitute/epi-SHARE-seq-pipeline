FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f

LABEL maintainer = "Eugenio Mattei"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="Trim ATAC fastqs"

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    autoconf \
    automake \
    binutils \
    build-essential \
    libcurl4-openssl-dev \
    liblz4-dev \
    liblzma-dev \
    libncurses5-dev \
    libbz2-dev \
    python3-dev \
    python3-pip \ 
    wget \
    zlib1g-dev &&\
    rm -rf /var/lib/apt/lists/*

# Install python packages
RUN pip install --no-cache-dir dnaio Levenshtein
# Install fastp
RUN wget http://opengene.org/fastp/fastp.0.20.1 && mv fastp.0.20.1 fastp && chmod a+x ./fastp && mv ./fastp /usr/local/bin

# Create and setup new user
ENV USER=shareseq
WORKDIR /home/$USER

RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

# Add folder with software to the path
ENV PATH="/software:${PATH}"

# Copy the compiled software from the builder
COPY --chown=$USER:$USER src/python/trim_fastq.py /usr/local/bin
COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin

USER ${USER}
