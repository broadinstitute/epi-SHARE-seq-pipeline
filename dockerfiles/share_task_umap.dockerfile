FROM bioconductor/bioconductor_docker:RELEASE_3_14

RUN apt-get update -qq && \
	apt-get install -y -qq --no-install-recommends\
		gtk-doc-tools \
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
		meson \
		pkg-config \
		python3 \
		python3-pip

RUN python3 -m pip install jupyter papermill

COPY install.R .

RUN Rscript install.R