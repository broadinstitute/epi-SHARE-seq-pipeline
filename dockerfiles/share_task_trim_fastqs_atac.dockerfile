FROM debian@sha256:946dcd5fcb0c6d072a8b12c25abdba542a40126c6e11c432128e8fe5dce511a7

LABEL maintainer = "Eugenio Mattei"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="Trim ATAC fastqs"

RUN apt-get update && apt-get install -y \
    wget &&\
    rm -rf /var/lib/apt/lists/*

# Install packages for python3 scripts (pysam, SAMstats)
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
COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin

USER ${USER}
