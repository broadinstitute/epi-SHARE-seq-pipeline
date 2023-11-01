FROM python@sha256:7ad180fdf785219c4a23124e53745fbd683bd6e23d0885e3554aff59eddbc377

LABEL maintainer = "Eugenio Mattei"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="Trim ATAC fastqs"

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    pigz \
    wget &&\
    rm -rf /var/lib/apt/lists/*

# Install python packages
RUN pip install --break-system-packages --no-cache-dir dnaio Levenshtein
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
