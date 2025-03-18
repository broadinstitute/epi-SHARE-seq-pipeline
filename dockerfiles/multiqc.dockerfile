# Use an official Python runtime as a parent image
FROM python:3.8-slim

# Set the working directory
WORKDIR /usr/src/app

# Install dependencies
RUN apt-get update && apt-get install -y \
    default-jre \
    wget \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Install FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip \
    && unzip fastqc_v0.12.1.zip \
    && rm fastqc_v0.12.1.zip \
    && chmod +x FastQC/fastqc \
    && ln -s /usr/src/app/FastQC/fastqc /usr/local/bin/fastqc

# Install MultiQC
RUN pip install multiqc

# Add FastQC to PATH
ENV PATH="/usr/src/app/FastQC:${PATH}"

# Default command
CMD ["bash"]


