FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f

LABEL maintainer = "Eugenio Mattei"
LABEL software.organization="Broad Institute of MIT and Harvard"

RUN apt-get update && apt-get install -y \
    pigz &&\
    rm -rf /var/lib/apt/lists/*

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
