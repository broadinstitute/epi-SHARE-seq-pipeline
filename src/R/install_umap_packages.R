#!/usr/bin/env Rscript

# Use all available CPU cores for the build
library(parallel)
options(Ncpus = detectCores())

# helper to install and load packages one by one
# to verify they were installed successfully;
# otherwise, R treats installation errors as warnings
setup <- function(installer, packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      installer(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}

install.packages("magrittr")

setup(install.packages, c(
  "IRkernel",
  "hdf5r",
#   "ggplot2",
#   "magrittr",
#   "plyr",
  
  # Cistopics
  "densityClust",
  "scatterplot3d",
  "umap",

  # for Bioconductor packages
  "BiocManager",
  
  # for development packages
  "devtools"
))

setup(BiocManager::install, c(
  "AnVIL"
))

setup(AnVIL::install, c(
  "biovizBase",
  "BSgenome.Hsapiens.UCSC.hg38",
  "EnsDb.Hsapiens.v86",
  "GenomeInfoDbData",
  "GenomicRanges",
  "Rsamtools",
  "Seurat"
))

# Jupyter kernel
IRkernel::installspec()

# Signac installation
if (!require("Signac", character.only = TRUE)){
    setRepositories(ind=1:2)
    install.packages("Signac")
}

# ArchR installation
if (!require("ArchR", character.only = TRUE)) {
  devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
  devtools::reload(pkg = pkgload::inst("pillar"), quiet = FALSE)
  devtools::reload(pkg = pkgload::inst("magrittr"), quiet = FALSE)
  library(ArchR)
  ArchR::installExtraPackages()
}

# Cistopic installation
if (!require("cisTopic", character.only = TRUE)) {
  Sys.setenv(LIBARROW_MINIMAL = "false")
  devtools::install_github("aertslab/cisTopic")
  library(cisTopic)
}
