FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f

LABEL maintainer = "Eugenio Mattei"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="10x preprocess"

RUN apt-get update && apt-get install -y \
    gcc \
    git \
    python3 \
    python3-dev \
    python3-pip \
    zlib1g-dev \
    wget &&\
    rm -rf /var/lib/apt/lists/*

# Install packages for python3 scripts (pysam, SAMstats)
RUN python3 -m pip install --no-cache-dir --ignore-installed --break-system-packages numpy pandas pybind11 --editable=git+https://github.com/GreenleafLab/matcha.git#egg=matcha

# Create and setup new user
ENV USER=shareseq
WORKDIR /home/$USER

RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

# Add folder with software to the path
ENV PATH="/software:${PATH}"

# Copy the compiled software from the builder
COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin
COPY --chown=$USER:$USER src/python/barcode_revcomp_detect.py /usr/local/bin
COPY --chown=$USER:$USER src/python/match_barcodes.py /usr/local/bin

USER ${USER}
