# syntax=docker/dockerfile:1

From r-base:4.1.2

RUN apt-get update && apt-get install -y python3 python3-pip libssl-dev pkg-config libxml2-dev libcurl4-openssl-dev libxt-dev libcairo2-dev libgsl-dev libfribidi-dev meson gtk-doc-tools libharfbuzz-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libhdf5-dev libmpfr-dev && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install jupyter papermill

COPY src/R/install_dorcs_packages.R src/R/DORCS_helper_functions.R src/R/TSSRanges.RData R/

RUN Rscript R/install_dorcs_packages.R

COPY src/jupyter_nb/dorcs_jplot_notebook.ipynb dorcs_jplot_notebook.ipynb

