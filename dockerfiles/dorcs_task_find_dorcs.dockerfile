############################################################
# Dockerfile for BROAD GRO share-seq-pipeline
# Based on Debian slim
############################################################

FROM r-base@sha256:fff003a52d076e963396876b83cfa88c4f40a8bc27e341339cd3cc0236c1db79 as builder

RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site

ENV R_LIBS_USER=/usr/local/lib/R
ENV RETICULATE_MINICONDA_ENABLED=FALSE

RUN apt-get update -qq && \
	apt-get install -y -qq --no-install-recommends\
		gtk-doc-tools \
        libssl-dev \
		libcairo2-dev \
		libcurl4-openssl-dev \
		libfreetype6-dev \
		libfribidi-dev \
		libgsl-dev \
		libharfbuzz-dev \
		libhdf5-dev \
		libjpeg-dev \
		libmpfr-dev \
		libpng-dev \
		libssl-dev \
		libtiff5-dev \
		libxml2-dev \
		libxt-dev \
        libgeos-dev \
		meson \
		pkg-config \
        python3 \
        python3-pip 
        
RUN R --no-echo --no-restore --no-save -e "install.packages(c('dplyr','patchwork','ggplot2','ggrepel','reshape2','circlize','networkD3','GGally','igraph','network','foreach','iterators','hdf5r','ggrastr','BiocManager','remotes','pbmcapply','doSNOW','Rmpfr', 'glue','magrittr','pillar','RcppArmadillo','reticulate','rlang','yaml','rpart','IRkernel','data.table', 'tidyft','qlcMatrix','logr'))"

RUN R --no-echo --no-restore --no-save -e "remotes::install_version('Seurat', version = '4.1.1')"

RUN R --no-echo --no-restore --no-save -e "BiocManager::install(c('Biostrings','rtracklayer','GenomicRanges','motifmatchr','ComplexHeatmap','chromVAR'), update=T, ask=F)"

RUN R --no-echo --no-restore --no-save -e  "remotes::install_github('caleblareau/BuenColors')"

ENV USER=shareseq
WORKDIR /home/$USER
RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER
        
RUN python3 -m pip install jupyter papermill

RUN chown $USER:$USER /usr/local/lib/R

COPY --chown=$USER:$USER src/jupyter_nb/dorcs_jplot_notebook.ipynb /usr/local/bin/

RUN mkdir -p /home/R/

COPY --chown=$USER:$USER src/R/DORCS_helper_functions_optimized.R src/R/TSSRanges.RData /home/R/

USER ${USER}

RUN R -e "IRkernel::installspec()"