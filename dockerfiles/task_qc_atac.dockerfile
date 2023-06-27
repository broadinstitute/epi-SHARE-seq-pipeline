############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on Debian slim
############################################################

FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f

LABEL maintainer = "Eugenio Mattei"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="qc-atac"

RUN apt-get update && apt-get install -y \
    gcc \
    git \
    python3 \
    python3-dev \
    python3-pip \
    r-base \
    tabix \
    zlib1g-dev &&\
    rm -rf /var/lib/apt/lists/*

# Install packages for python3 scripts (pysam, SAMstats)
RUN python3 -m pip install --break-system-packages --no-cache-dir --ignore-installed numpy matplotlib pandas plotnine pysam xopen

# Add folder with software to the path
ENV PATH="/software:${PATH}"

# Copy the compiled software from the builder
COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin
COPY --chown=$USER:$USER src/python/pbc_stats.py /usr/local/bin
COPY --chown=$USER:$USER src/python/qc_atac_compute_tss_enrichment.py /usr/local/bin
COPY --chown=$USER:$USER src/python/qc_atac_count_duplicates_per_barcode.py /usr/local/bin
COPY --chown=$USER:$USER src/python/qc_atac_compute_reads_in_peaks.py /usr/local/bin
COPY --chown=$USER:$USER src/python/plot_insert_size_hist.py /usr/local/bin
COPY --chown=$USER:$USER src/R/barcode_rank_functions.R /usr/local/bin
COPY --chown=$USER:$USER src/R/atac_qc_plots.R /usr/local/bin
COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin


USER ${USER}