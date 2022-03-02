#!/usr/bin/env Rscript

#Possibly not needed anymore

# Use all available CPU cores for the build
library(parallel)
options(Ncpus = 4)

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
  
  # Cistopics
  "densityClust",
  "scatterplot3d",
  "umap",

  # for Bioconductor packages
  "BiocManager",
  
  # for development packages
  "devtools"
))

#setup(BiocManager::install, c(
#  "AnVIL"
#))

setup(BiocManager::install, c(
  "biovizBase",
  #"BSgenome.Hsapiens.UCSC.hg38",
  "EnsDb.Hsapiens.v86",
  "GenomeInfoDbData",
  "GenomicRanges",
  "Rsamtools",
  "Seurat"
))

# Jupyter kernel
#Run this on the production image
#IRkernel::installspec()

# Signac installation

# if (!require("Signac", character.only = TRUE)){
#     setRepositories(ind=1:2)
#     install.packages("Signac")
# }

# ArchR installation
if (!require("ArchR", character.only = TRUE)) {
  devtools::install_github("GreenleafLab/ArchR@v1.0.1", repos = BiocManager::repositories())
  devtools::reload(pkg = pkgload::inst("pillar"), quiet = FALSE)
  devtools::reload(pkg = pkgload::inst("magrittr"), quiet = FALSE)
  library(ArchR)
  ArchR::installExtraPackages()
}

# Cistopic installation

# if (!require("cisTopic", character.only = TRUE)) {
#   Sys.setenv(LIBARROW_MINIMAL = "false")
#   devtools::install_github("aertslab/cisTopic")
#   library(cisTopic)
# }
